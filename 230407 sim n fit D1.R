# Conditions: Change line 9 if necessary for other samples sizes (and at the end)
#
#

rm(list=ls()); gc()
setwd("C:/Users/nuria/OneDrive/Escritorio/Universidad/TFM-Paper/Sim v4")

# 100 175 250 325 400 475 550
n <- 550
sup_data <- readRDS(paste0("sup_data_",n,".rds"))                                      #!#
par <- readRDS("par.rds")
n_sim <- length(sup_data)

# 1. Prepare data for D1

# 1.1 Add ID cohorts

P <- 2 / c(rep(15,2), rep(15*2,11))
breaks <- cumsum(round(P*n))[1:12]

for (j in 1:n_sim){
  p1 <- sup_data[[j]]
  for (i in seq_along(p1)) {
    # Get the index of the break that the current data frame belongs to
    index <- findInterval(i, breaks)
    # Add a new column with the corresponding number
    p1[[i]][["coh"]] <- index
  }
  sup_data[[j]] <- p1
}

# 1.2 Select rows for each cohort. D1
sup_data_coh <- list()

for (j in 1:n_sim){
p2 <- sup_data[[j]]

for (i in seq_along(p2)) {
  col_value <- unique(p2[[i]]$coh)+1
  row_indices <- c(col_value, col_value + 2) # D1
  p2[[i]] <- p2[[i]][row_indices, ]
}

sup_data_coh[[j]] <- p2
}

# 2. Fit and save parameters

## Specify State-Space Model ----
library(OpenMx)
opmxL <- list()

# Autoregressive dynamics A(t)
opmxL$amat <- mxMatrix("Full", nrow=4, ncol=4, 
                       free = c(T,T,F,F,
                                T,T,F,F,
                                F,F,F,F,
                                F,F,F,F),
                       values = c(par$A_ct[1,1], par$A_ct[2,1], 0, 0,
                                  par$A_ct[1,2], par$A_ct[2,2], 0, 0,
                                  1, 0, 0, 0,
                                  0, 1, 0, 0),
                       labels = c("a11", "a21", NA, NA,
                                  "a12", "a22", NA,NA,
                                  NA,NA,NA,NA,
                                  NA,NA,NA,NA),
                       lbound = c(-1, -1, NA, NA,
                                  -1, -1, NA, NA,
                                  NA, NA, NA, NA,
                                  NA, NA, NA, NA),
                       dimnames = list(c("x0", "y0", "xA", "yA"),
                                       c("x0", "y0", "xA", "yA")), 
                       name="A")

# B·u[t] # Input effects on the latent variables
opmxL$bmat <- mxMatrix("Zero", 4, 1, name = "B")

# C·x[t] # Factor loadings in the measurement model
opmxL$cmat <- mxMatrix("Full", 2, 4, free = FALSE, 
                       values = c(1,0,0,1,
                                  0,0,0,0), 
                       dimnames = list(c("x", "y"),
                                       c("x0", "y0", "xA", "yA")),
                       labels = c(NA, NA, NA, NA,
                                  NA, NA, NA, NA),
                       name = "C")

# D·u[t] # Input effects on the observed variables
opmxL$dmat <- mxMatrix("Zero", 2, 1, name = "D")

# Dynamic error (0) q[t]
opmxL$qmat <- mxMatrix("Zero", 4, 4, name = "Q")

# Measurement error r[t]
opmxL$rmat <- mxMatrix("Full", 2, 2, free = TRUE,
                       values = c(par$xmer, par$mercv,
                                  par$mercv, par$ymer),
                       labels = c("xmer", "mercv",
                                  "mercv", "ymer"),
                       name = "R")

# x(t) # Mean vector of latent variables at t=0
opmxL$xmat <- mxMatrix("Full", 4, 1, free = TRUE, values = par$Mu, 
                       labels= c("x0mn", "y0mn", "xAmn", "yAmn"),
                       name = "x0")

# P0 # Covariance matrix of latent variables at t=0
opmxL$pmat <- mxMatrix("Symm", 4, 4, values= as.vector(round(par$Sigma, 6)), # 0.001 ,
                       labels = c(  "x0v",     "x0y0cv",  "x0xAcv",  "x0yAcv", # x0
                                    "x0y0cv",  "y0v",     "y0xAcv",  "y0yAcv", # y0
                                    "x0xAcv",  "y0xAcv",  "xAv",     "xAyAcv", # xA
                                    "x0yAcv",  "y0yAcv",  "xAyAcv",  "yAv"),   # yA
                       free=TRUE,
                       lbound=c(0, NA, NA, NA,
                                NA, 0, NA, NA,
                                NA, NA, 0, NA,
                                NA, NA, NA, 0), 
                       name="P0")

# u[t] # Covariates
opmxL$umat <- mxMatrix("Zero", 1, 1, name = "u")

# Specification of the time index [t]
opmxL$tmat <- mxMatrix('Full', 1, 1, name='time', labels='data.time')

# Store matrices and algebra
opmxL$modataL_ct <- with(opmxL, list(opmxL$amat, opmxL$bmat, 
                                     opmxL$cmat, opmxL$dmat,
                                     opmxL$qmat, opmxL$rmat, 
                                     opmxL$umat, opmxL$tmat,
                                     opmxL$xmat, opmxL$pmat))
# Select expectation
opmxL$expSSCT <- mxExpectationStateSpaceContinuousTime(A = "A", B = "B",
                                                       C = "C", D = "D",
                                                       Q = "Q", R = "R",
                                                       x0 = "x0", P0 = "P0",
                                                       u = "u", t = "time")

################################################################################

## Create SSM-CT multisubject models
modNames <- paste0("indiv", 1:n)

# GET FITS

# Set up a parallel backend
library(parallel)
library(snow)
cl <- makeCluster(7)

# Get functions
ex <- Filter(function(x) {is.function(get(x, .GlobalEnv))},
             c(ls(.GlobalEnv, all.names = TRUE), ls(package:OpenMx)))

# Export objects that are needed to workers
clusterExport(cl, c(ex, "opmxL", "n", "modNames"))

# Parallelization

BLCS.fit <- parLapply(cl, sup_data_coh, function(data) {
  indivmodels <- list()
  for (k in 1:n) {
    DataSetForSubjectK <- data[[k]]
    indivmodels[[k]] <- mxModel(name = modNames[k],
                                opmxL$modataL_ct,
                                opmxL$expSSCT,
                                mxFitFunctionML(),
                                mxData(DataSetForSubjectK, type = 'raw'))
  }
  BLCS.model <- mxModel(name = "BLCS_CT", indivmodels,
                        mxFitFunctionMultigroup(modNames))
  mxRun(BLCS.model)
}) ; beepr::beep(5)

# Stop the parallel backend
stopCluster(cl)

setwd("C:/Users/nuria/OneDrive/Escritorio/Universidad/TFM-Paper/Sim v4/Results D1")

saveRDS(BLCS.fit, paste0("BLCS.fit_D1_",n,".rds"))                                     #!#

# GET SUMMARIES

# Split the input into four parts
chunks <- split(BLCS.fit, rep(1:4, each = 50, length.out = length(BLCS.fit)))

# Set up a new parallel backend for summaries
cl <- makeCluster(7)

# Process the first chunk of 50 elements
summaries1 <- parLapply(cl, chunks[[1]], summary); beepr::beep(1)

# Process the second chunk of 50 elements
summaries2 <- parLapply(cl, chunks[[2]], summary); beepr::beep(2)

# Process the third chunk of 50 elements
summaries3 <- parLapply(cl, chunks[[3]], summary); beepr::beep(11)

# Process the fourth chunk of 50 elements
summaries4 <- parLapply(cl, chunks[[4]], summary); beepr::beep(4)

# Stop the parallel backend
stopCluster(cl)

summaries <- unlist(list(summaries1, summaries2, summaries3, summaries4), recursive = FALSE)

setwd("C:/Users/nuria/OneDrive/Escritorio/Universidad/TFM-Paper/Sim v4/Results D1")
saveRDS(summaries, paste0("summaries_D1_",n,".rds"))                                   #!#


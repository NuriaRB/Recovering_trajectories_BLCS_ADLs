# 1. Get pars in CT
# 2. Simulate latent values
# 3. Select and add measurement error. (Switch T/F to select random or fixed weeks)
# 4. Plot
rm(list=ls()); gc()
setwd("C:/Users/nuria/OneDrive/Escritorio/Universidad/TFM-Paper/Sim v4")
source("functionsv2.R")

# 1. Get original pars and simulation conditions ###############################
par <- readRDS("parCT.rds")

# par$Mu <- with(par, c(x0mn, y0mn, xAmn, yAmn))
# par$Sigma <- with(par, matrix(c(
#   # x0     y0       xA       yA
#   x0v,     x0y0cv,  x0xAcv,  x0yAcv,         # x0
#   x0y0cv,  y0v,     y0xAcv,  y0yAcv,         # y0
#   x0xAcv,  y0xAcv,  xAv,     xAyAcv,         # xA
#   x0yAcv,  y0yAcv,  xAyAcv,  yAv), nrow=4) ) # yA
par$Mu
par$Sigma 

par$A_ct <- with(par, matrix(c(a11, a21, a12, a22), nrow=2)) 

parsdt <- CT2DT(par$Mu, par$Sigma, par$A_ct, delta=1/52)

# 2. Simulate latent values ####################################################

n_sim <- 200

parsdt$n = 550   # Subjects                                                  #!#
parsdt$years = 15
parsdt$delta = 1/52 # 1/delta = number of measurements per year

latent <- list()
for (i in 1:n_sim){
  latent[[i]] <- BLCSsim(parsdt$Mu, parsdt$Sigma, parsdt$A_dt, parsdt)
}

# 3. Sample observed values ####################################################

parsdt$xmer <- par$xmer
parsdt$ymer <- par$ymer
parsdt$mercv <- par$mercv

manifest <- list()
for (i in 1:n_sim){
  manifest[[i]] <- OBSsample(parsdt, latent[[i]], random=T)
}


# Pair time points and manifest values for each subject

sup_data <- list()
col <- ncol(manifest[[1]][["values"]])
for (j in 1:n_sim){
  data <- list()
  for (i in 1:parsdt$n){
    data[[i]] <- data.frame(time = manifest[[j]][[1]][[i]], 
                            x = manifest[[j]][[2]][i, 1:(col/2)],
                            y= manifest[[j]][[2]][i, (col/2+1):col])
  }
  sup_data[[j]] <- data
}


# 4. Prepare data for ALD D0

# sup_data # No changes needed
saveRDS(sup_data, "sup_data_550.rds")                                       #!#
saveRDS(par, "par.rds")

# Fit observed measurements  #######################################

n <- parsdt$n

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
                       values = c(parsdt$xmer, parsdt$mercv,
                                  parsdt$mercv, parsdt$ymer),
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


# Set up a parallel backend
library(parallel)
library(snow)
cl <- makeCluster(7)

# Get functions
ex <- Filter(function(x) {is.function(get(x, .GlobalEnv))},
             c(ls(.GlobalEnv, all.names = TRUE), ls(package:OpenMx)))

# Export objects that are needed to workers
clusterExport(cl, c(ex, "opmxL", "modNames", "n"))

# Parallelization

BLCS.fit <- parLapply(cl, sup_data, function(data) {
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
}) ; beepr::beep(2)

# Stop the parallel backend
stopCluster(cl)

setwd("C:/Users/nuria/OneDrive/Escritorio/Universidad/TFM-Paper/Sim v4/Results D0")
saveRDS(BLCS.fit, "BLCS.fit_D0_550.rds")                                     #!#

# Set up a new parallel backend for summaries
cl <- makeCluster(7)
# Export objects that are needed to workers
clusterExport(cl, c("BLCS.fit"))

summaries <- parLapply(cl, BLCS.fit, summary);beepr::beep(9)

# Stop the parallel backend
stopCluster(cl)

setwd("C:/Users/nuria/OneDrive/Escritorio/Universidad/TFM-Paper/Sim v4/Results D0")
saveRDS(summaries, "summaries_D0_550.rds")                                   #!#


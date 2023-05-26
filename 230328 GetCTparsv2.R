# 1. Take Pablo's parameters in DT
# 2. Simulate various datasets
# 3. Fit each dataset and obtain CT parameters
# 4. Average CT parameters

rm(list=ls());gc()
setwd("C:/Users/nuria/OneDrive/Escritorio/Universidad/TFM-Paper/Sim v4")
source("functionsv2.R")

# 1. Take Pablo's parameters in DT ---------------------------------------------
# Values from Cáncer & Estrada (2023) Table 1

pablo <- list(# Mean vector
              x0mn = -.8, y0mn = -.2, xAmn = .8, yAmn = 1,
              # Dynamics
              a11 = -.35, a22 = -.25, a12 = .1, a21 = .2,
              # Var-covar
              x0v = .3, y0v = .6, xAv = .02, yAv = .08,
              x0xAcv = .046476, x0y0cv = .212132, x0yAcv = .030984,
              y0xAcv = .032863, y0yAcv = .131453, xAyAcv = .016,
              # Measurement error
              xmer = .1, ymer = .2, mercv = .02,
              # Simulation conditions
              seed = 35,
              n = 250,   # Subjects
              years = 15,
              delta = 1  # Delta = 1 (comes from DT)
              )

pablo$Mu <- with(pablo, c(x0mn, y0mn, xAmn, yAmn))
pablo$Sigma <- with(pablo, matrix(c(
  # x0     y0       xA       yA
  x0v,     x0y0cv,  x0xAcv,  x0yAcv,         # x0
  x0y0cv,  y0v,     y0xAcv,  y0yAcv,         # y0
  x0xAcv,  y0xAcv,  xAv,     xAyAcv,         # xA
  x0yAcv,  y0yAcv,  xAyAcv,  yAv), nrow=4) ) # yA

pablo$A_dt <- with(pablo, matrix(c(a11, a21, a12, a22), nrow=2)) 


# 2. Simulate various datasets -------------------------------------------------
n_sim <- 20
n <- pablo$n

I <- diag(1, 15)
m.err <- with(pablo, matrix(c(xmer, mercv, mercv, ymer), nrow=2))
psi <- kronecker(m.err, I)

latent <- list()
manifest <- list()


for (i in 1:n_sim){
  cat("\r", i/n_sim*100, "%")
  latent[[i]] <- BLCSsim(pablo$Mu, pablo$Sigma, pablo$A_dt, pablo)

  m.errors <- MASS::mvrnorm(n = n,
                            mu = rep(0, 30),
                            Sigma = psi)

  manifest[[i]] <-  latent[[i]]+m.errors
}


# 3. Fit each dataset and obtain CT parameters ---------------------------------

# Prepare the data for various subjects estimation model
# standarize such data on the way 

Tmax <- pablo$years/pablo$delta
time <- seq(0, by=pablo$delta, to=pablo$years-1) 
sup_data <- list()

for (j in 1:n_sim){
  data <- list()
  for (i in 1:n){
    data[[i]] <- data.frame(time = time, 
                            x = (manifest[[j]][i, 
                                1:(ncol(manifest[[j]])/2)] - pablo$x0mn)/sqrt(pablo$x0v),
                            y = (manifest[[j]][i, 
                                (ncol(manifest[[j]])/2+1):ncol(manifest[[j]])] - pablo$y0mn)/sqrt(pablo$y0v))
  }
    sup_data[[j]] <- data
}


# Plot random simulation to check
library(data.table);library(ggplot2)
combined_df <- rbindlist(sup_data[[sample(1:n_sim,1)]], idcol = "subject")

ggplot(combined_df, aes(x = time, y = y, group = subject)) +
  geom_line(size=.1, colour = "blue") +
  labs(x = "Time points", y = "Y")
ggplot(combined_df, aes(x = time, y = x, group = subject)) +
  geom_line(size=.1, colour = "red") +
  labs(x = "Time points", y = "X")


## Specify State-Space Model ---------------------------------------------------
A = pablo$A_dt
Mu = pablo$Mu
Sigma = pablo$Sigma


library(OpenMx)
opmxL <- list()

# Autoregressive dynamics A(t)
opmxL$amat <- mxMatrix("Full", nrow=4, ncol=4, 
                       free = c(T,T,F,F,
                                T,T,F,F,
                                F,F,F,F,
                                F,F,F,F),
                       values = c(A[1,1], A[2,1], 0, 0,
                                  A[1,2], A[2,2], 0, 0,
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
                       values = c(pablo$xmer, pablo$mercv,
                                  pablo$mercv, pablo$ymer),
                       labels = c("xmer", "mercv",
                                  "mercv", "ymer"),
                       name = "R")

# x(t) # Mean vector of latent variables at t=0
opmxL$xmat <- mxMatrix("Full", 4, 1, free = TRUE, values = Mu, 
                       labels= c("x0mn", "y0mn", "xAmn", "yAmn"),
                       name = "x0")

# P0 # Covariance matrix of latent variables at t=0
opmxL$pmat <- mxMatrix("Symm", 4, 4, values = as.vector(round(Sigma,6)), # 0.001 ,
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
## End of SSM specification ----------------------------------------------------
modNames <- paste0("indiv", 1:n)

# Set up a parallel backend
library(parallel)
library(snow)
cl <- makeCluster(7)

# Get functions
ex <- Filter(function(x) {is.function(get(x, .GlobalEnv))},
             c(ls(.GlobalEnv, all.names = TRUE), ls(package:OpenMx)))

# Export objects that are needed to workers
clusterExport(cl, c("opmxL", "modNames","n",ls(package:OpenMx)))

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
}) ; beepr::beep(8)

# Stop the parallel backend
stopCluster(cl)

# Set up a new parallel backend for summaries
cl <- makeCluster(8)
# Export objects that are needed to workers
clusterExport(cl, c("BLCS.fit"))

summaries <- parLapply(cl, BLCS.fit, summary);beepr::beep(8)

# Stop the parallel backend
stopCluster(cl)


# 4. Average CT parameters -----------------------------------------------------

# Save each estimation from each model (in total: n_sim) in a list

saveRDS(BLCS.fit, "BLCS.fit_10_std.rds")
saveRDS(summaries, "summaries_10_std.rds")
 
# BLCS.fit <- readRDS("BLCS.fit.rds")
# summaries <- readRDS("summaries_10_std.rds")
n_sim <- length(summaries)
# Get estimates
df_list <- list()
for(j in 1:n_sim){
  df_list[[j]] <- summaries[[j]][["parameters"]][c(1,5)]
}

# Merge in one dataframe
library(tidyverse)
merged_df <- reduce(df_list, full_join, by = "name")

# Get mean value of all simulations
mean_estimate <- rowMeans(merged_df[-1])

# Save in one list for future uses
names <- merged_df[[1]]
parCT <- list()
for(i in 1:length(names)){
  parCT[[names[i]]] <- mean_estimate[i]
}

# Correct covariances so correlations remain

pablo$corM <- round(cov2cor(pablo$Sigma), 3)

parCT$Sigma <- psych::cor2cov(rho = pablo$corM,
               sigma = sapply(parCT[c("x0v", "y0v", "xAv", "yAv")],
                              sqrt) )
dimnames(parCT$Sigma) <- list(c("x0v", "y0v", "xAv", "yAv"),
                         c("x0v", "y0v", "xAv", "yAv") )

# parCT$x0v    <- 
# parCT$y0v    <- 
# parCT$x0sd   <- 
# parCT$y0sd   <- 
# parCT$xAv    <- 
# parCT$yAv    <- 
# parCT$xAsd   <- 
# parCT$yAsd   <- 
# parCT$x0y0cv <- 
# parCT$x0xAcv <- 
# parCT$x0yAcv <- 
# parCT$y0xAcv <- 
# parCT$y0yAcv <- 
# parCT$xAyAcv <- 



parCT$Mu <- with(parCT, c(x0mn, y0mn, xAmn, yAmn))

par$A_ct <- with(parCT, matrix(c(a11, a21, a12, a22), nrow=2)) 

parCT 

saveRDS(parCT, "parCT.rds")



# TESTS 
################################################################################

pablo <- list(# Mean vector
  x0mn = -.8, y0mn = -.2, xAmn = .8, yAmn = 1,
  # Dynamics
  a11 = -.35, a22 = -.25, a12 = .1, a21 = .2,
  # Var-covar
  x0v = .3, y0v = .6, xAv = .02, yAv = .08,
  x0xAcv = .046476, x0y0cv = .212132, x0yAcv = .030984,
  y0xAcv = .032863, y0yAcv = .131453, xAyAcv = .016,
  # Measurement error
  xmer = .1, ymer = .2, mercv = .02,
  # Simulation conditions
  seed = 35,
  n = 200,   # Subjects
  years = 15,
  delta = 1  # Delta = 1 (comes from DT)
)

setwd("C:/Users/nuria/OneDrive/Escritorio/Universidad/TFM-Paper/Sim v4")
parCT <- readRDS("parCT.rds")

pablo$x0y0r <-with(pablo, x0y0cv/sqrt(x0v*y0v))
pablo$x0xAr <-with(pablo, x0xAcv/sqrt(x0v*xAv))
pablo$x0yAr <-with(pablo, x0yAcv/sqrt(x0v*yAv))
pablo$y0xAr <-with(pablo, y0xAcv/sqrt(y0v*xAv))
pablo$y0yAr <-with(pablo, y0yAcv/sqrt(y0v*yAv))
pablo$xAyAr <-with(pablo, xAyAcv/sqrt(xAv*yAv))


parCT$x0y0r <-with(parCT, x0y0cv/sqrt(x0v*y0v))
parCT$x0xAr <-with(parCT, x0xAcv/sqrt(x0v*xAv))
parCT$x0yAr <-with(parCT, x0yAcv/sqrt(x0v*yAv))
parCT$y0xAr <-with(parCT, y0xAcv/sqrt(y0v*xAv))
parCT$y0yAr <-with(parCT, y0yAcv/sqrt(y0v*yAv))
parCT$xAyAr <-with(parCT, xAyAcv/sqrt(xAv*yAv))


pablo$x0y0r;parCT$x0y0r
pablo$x0xAr;parCT$x0xAr
pablo$x0yAr;parCT$x0yAr
pablo$y0xAr;parCT$y0xAr
pablo$y0yAr;parCT$y0yAr
pablo$xAyAr;parCT$xAyAr

pablo$Sigma <- with(pablo, matrix(c(
  # x0     y0       xA       yA
  x0v,     x0y0cv,  x0xAcv,  x0yAcv,         # x0
  x0y0cv,  y0v,     y0xAcv,  y0yAcv,         # y0
  x0xAcv,  y0xAcv,  xAv,     xAyAcv,         # xA
  x0yAcv,  y0yAcv,  xAyAcv,  yAv), nrow=4) ) # yA

parCT$Sigma <- with(parCT, matrix(c(
  # x0     y0       xA       yA
  x0v,     x0y0cv,  x0xAcv,  x0yAcv,         # x0
  x0y0cv,  y0v,     y0xAcv,  y0yAcv,         # y0
  x0xAcv,  y0xAcv,  xAv,     xAyAcv,         # xA
  x0yAcv,  y0yAcv,  xAyAcv,  yAv), nrow=4) ) # yA


cov2cor(pablo$Sigma)
cov2cor(parCT$Sigma)


# 1. Convert CT to DT

CT2DT <- function(Mu, Sigma, A_ct, delta){
  
  # Initial level means
  mn0 <- Mu[c(1,2)]
  # Initial level var-covar
  v0 <- Sigma[1:2,1:2]
  
  # Add. comp. means
  b <- Mu[c(3,4)]
  # Add. comp, var-covar
  B <- Sigma[3:4,3:4]
  
  # Cross covariances
  v <- t(Sigma[3:4,1:2])
  
  ## DT Conversion ## Voelkle & Oud 2015 page 370##
  library(MASS)
  library(Matrix)
  I <- diag(1, 2)
  
  # Beta and gamma
  A_dt <- expm(A_ct*delta)-I 
  conversion <- solve(A_ct) %*% (expm(A_ct*delta)-I)
  
  # Intercept vector
  b_dt <- conversion %*% b
  
  # Covariance matrix between traits
  B_dt <- conversion %*% B %*% t(conversion)
  
  # Covariates between latent traits and initial time points
  v_dt <- v %*% t(conversion) 
  
  # Build V_dt matrix
  V_dt <- matrix(c(Sigma[1,1], Sigma[1,2], v_dt[1], v_dt[3],
                   Sigma[2,1], Sigma[2,2], v_dt[2], v_dt[4],
                   v_dt[1],    v_dt[2],    B_dt[1], B_dt[3],
                   v_dt[3],    v_dt[4],    B_dt[2], B_dt[4]), ncol=4, byrow=T)
  
  res <- list(Mu = c(mn0, as.vector(b_dt)),
              Sigma = V_dt,
              A_dt = A_dt,
              delta = delta)
}

# 2. Simulate BLCS data

BLCSsim <- function(Mu_dt, Sigma_dt, A_dt, con){
  # con is a list that contains all parameters and simulation conditions
  
  # Total repeated measurements
  con$Tmax = con$years/con$delta
  
  con$L <- MASS::mvrnorm(n = con$n,
                         mu = Mu_dt,
                         Sigma = Sigma_dt)
  
  # Create the matrices to store the repeated measures
  xmat <- cbind(con$L[,1], matrix(NA, ncol=con$Tmax-1, nrow=con$n))
  ymat <- cbind(con$L[,2], matrix(NA, ncol=con$Tmax-1, nrow=con$n))
  
  # Calculate latent variables and fill the matrices with the results:
  for(i in 1:(con$Tmax-1)){
    
    # x_i = x_i-1 + c.a.x_i + b_x*x_i-1 + g_x*y_i-1
    xmat[,(i+1)] <- xmat[,i] + con$L[,3] + A_dt[1,1]*xmat[,i] + A_dt[1,2]*ymat[,i]
    
    # y_i = y_i-1 + c.a.y_i + b_y*y_i-1 + g_y*x_i-1
    ymat[,(i+1)] <- ymat[,i] + con$L[,4] + A_dt[2,2]*ymat[,i] + A_dt[2,1]*xmat[,i]
  }
  latent <- cbind(xmat, ymat)
  colnames(latent) <- c(paste0("X", 1:con$Tmax), paste0("Y", 1:con$Tmax))
  
  return(latent)
}

# 3. Sample weeks from latent data and add measurement error

OBSsample <- function(con, latent, random=T) {
  # con is a list that contains all parameters and simulation conditions
  
  # Total repeated measurements
  con$Tmax = con$years/con$delta
  
  if (random==F){
    # First week of every year
    selected <- seq(from=1, to=con$Tmax, by=52)
    subj_week <- list()
    for (j in 1:con$n){
      subj_week[[j]] <- selected
    }
  }else{
    # Random week of each year
    subj_week <- list()
    for (j in 1:con$n){
      selected <- numeric(con$years)
      for (i in 1:con$years) {
        selected[i] <- sample(52, 1) + (i - 1) * 52
      }
      subj_week[[j]] <- selected
    }
  }
  
  # Matrix with selected latent values
  sel_lat <- matrix(NA, nrow = con$n, ncol=length(selected)*2)
  for (i in 1:con$n) {
    sel_lat[i, ] <- latent[i, c(subj_week[[i]], subj_week[[i]]+con$Tmax)]
  }
  
  # Number of observations
  obs <- length(selected)
  
  # Add measurement error to those selected
  I <- diag(1, obs)
  m.err <- with(con, matrix(c(xmer, mercv, mercv, ymer), nrow=2))
  psi <- kronecker(m.err, I)
  m.errors <- MASS::mvrnorm(n = con$n,
                            mu = rep(0, obs*2),
                            Sigma = psi)
  
  manifest <-  sel_lat+m.errors
  
  time <- seq(from=con$delta, by=con$delta, to=con$years)
  time_list <- list()
  for (i in 1:con$n){
    time_list[[i]] <- time[subj_week[[i]]]
  }
  
  return(list(time=time_list, values=manifest))
}

################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
#
# This file contains R functions that realize the main algorithms presented in the paper
#
#
# HBEF works here in its simplest form: without the L_o term.
#
# ``Hierarchical Bayes Ensemble Kalman Filtering''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
#
#
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 14 July 2015
################################################################################
# functions.R
# File with definitions of functions
#
# ==================================
# General comments.
#
# 1. By a "universe" we mean a particular realization of the structural sequences {F_k} and {sigma_k}.
# 2. By a "world"    we mean a particular realization of the "truth" {x_k} (given {F_k} and {sigma_k}) and observations.
# ==================================
#
# list of functions:
#      -- create_parameters_universe_world
#      -- update_parameters
#      -- generate_universe
#      -- parameters_kf
#      -- kf
#      -- parameters_var
#      -- var
#      -- parameters_enkf
#      -- enkf
#      -- parameters_henkf
#      -- henkf
#      -- parameters_hbef
#      -- hbef
#      -- rms

#====================================================================
# This function specifies parameters that govern the model of "truth"

create_parameters_universe_world <- function(){
  
  # External parameters that govern the model of ``truth'' and observation-error statistics
  
  tau_x     <- 12                        # e-folding time for x if F_0 = const
  #                                        tau_x defines F_0 by F_0 = exp(-1/tau_x) 
  #
  tau_F     <- 18                        # e-folding time for F
  #                                        tau_F defines mu by mu = exp(-1/tau_A)
  #
  tau_Sigma <- 18                        # e-folding time for Sigma
  #                                        tau_Sigma defines kappa by kappa = exp(-1/tau_Sigma)
  #
  pi        <- 0.05                      # Pr(F > 1). pi and F_0 define std_F and var_epsilon_F
  #
  Sigma_0   <- 0                         # mean of Sigma (this is rather a constant than a parameter, don't change it)
  #
  std_Sigma <- 0.5                       # standart deviation of Sigma
  #                           
  std_eta   <- 9                         # standart deviation of eta (observation-error st.dev.)
  ## -------------------------------------------------------------------------------
  
  # Internal (derived) parameters to be used in the model equations
  
  # Model of {F_k}: 
  # F_k=mu*(F_{k-1} - F_0) + sqrt(var_epsilon_F) * epsilon_F_k
  
  F_0               <- exp(-1/tau_x)               # mean of F
  #
  mu                <- exp(-1/tau_F)               # mu  
  #
  std_F             <- (1 - F_0) / qnorm(1 - pi)   # standart deviation of F
  #
  var_epsilon_F     <- std_F^2*(1-mu^2)            # var of epsilon_F
  
  # Model of {Sigma_k}: 
  # Sigma_k=kappa*(Sigma_{k-1} - F_0) + \sqrt(var_epsilon_Sigma) * epsilon_Sigma_k
  
  kappa                <- exp(-1/tau_Sigma)           # kappa
  #
  var_epsilon_Sigma <- std_Sigma^2*(1-kappa^2)        # var of epsilon_Sigma
  ## ===============================================================================
  
  # Other parameters that define the setup of a numerical experiment
  
  time              <- 200000                     # number of time moments in the time series (n_{time})
  N                 <- 5                          # ensemble size    
  L                 <- 1                          # number of independent runs (worlds) under the same {F_k} and {Sigma_k} (universe)
  distort_Q         <- 1                          # coefficient of Q distortion
  # ------------------------------
  seed_for_universe1     <- 42433425              # random seed for initiation the {F_k} time series
  seed_for_universe2     <- 91435564              # random seed for initiation the {SIgma_k} time series
  seed_for_world         <- 24674353              # random seed for initiation of world
  seed_for_filters       <- 56345352              # random seed for initiation of other pseudo-random sources
  
  # list of parameters
  list <- data.frame(seed_for_universe1, seed_for_universe2, seed_for_filters, seed_for_world, distort_Q, 
                     time, tau_x, tau_F, std_F, mu, var_epsilon_F, pi, tau_Sigma, 
                     std_Sigma, kappa, var_epsilon_Sigma, F_0, Sigma_0, std_eta, N, L)
  print(list)
  return(list)
}


update_parameters <- function(parameters){
  # Internal (derived) parameters to be used in the model equations
  
  # Model of {F_k}: 
  # F_k=mu*(F_{k-1} - F_0) + sqrt(var_epsilon_F) * epsilon_F_k
  
  parameters$F_0               <- exp(-1/parameters$tau_x)               # mean of F
  #
  parameters$mu                <- exp(-1/parameters$tau_F)               # mu  
  #
  parameters$std_F             <- (1 - parameters$F_0) / qnorm(1 - parameters$pi)   # standart deviation of F
  #
  parameters$var_epsilon_F     <- parameters$std_F^2*(1-parameters$mu^2)            # var of epsilon_F
  
  # Model of {Sigma_k}: 
  # Sigma_k=kappa*(Sigma_{k-1} - F_0) + \sqrt(var_epsilon_Sigma) * epsilon_Sigma_k
  
  parameters$kappa                <- exp(-1/parameters$tau_Sigma)           # kappa
  #
  parameters$var_epsilon_Sigma <- parameters$std_Sigma^2*(1-parameters$kappa^2)        # var of epsilon_Sigma
  ## ===============================================================================
  return(parameters)
}



#====================================================================
# This function simulates a realization of the sequences {F_k}, {Sigma_k}, and $Q_k$.

generate_universe <- function(parameters){
  set.seed(parameters$seed_for_universe1)
  F_0           <- rnorm(1, parameters$F_0, sd = parameters$std_F)
  Sigma_0       <- rnorm(1, parameters$Sigma_0, sd = parameters$std_Sigma)
  epsilon_F     <- rnorm(parameters$time, 0, sd = sqrt(parameters$var_epsilon_F))
  
  set.seed(parameters$seed_for_universe2)
  epsilon_Sigma <- rnorm(parameters$time, 0, sd = sqrt(parameters$var_epsilon_Sigma))
  
  F             <- c(1:parameters$time)
  Sigma         <- c(1:parameters$time)
  
  F[1]          <- parameters$F_0 + parameters$mu * (F_0 - parameters$F_0) + epsilon_F[1]
  Sigma[1]      <- parameters$Sigma_0 + parameters$kappa * (Sigma_0 - parameters$Sigma_0) + epsilon_Sigma[1]
  for(i in (2:parameters$time)){
    F[i]        <- parameters$F_0 + parameters$mu * (F[i-1] - parameters$F_0) + epsilon_F[i]
    Sigma[i]    <- parameters$Sigma_0 + parameters$kappa * (Sigma[i-1] - parameters$Sigma_0) + epsilon_Sigma[i]
  }
  Q <- (exp(Sigma))^2*parameters$distort_Q
  return(data.frame(F, Sigma, Q, F_0))
}



#====================================================================
# This function simulates a realization of the "truth" {x_k} and observations {y_k}.

generate_world <- function(universe, parameters){
  ## ============================================================
  ## ==================== Generating of the world ===============
  ## ============================================================
  set.seed(parameters$seed_for_world)
  
  X_0           <- rnorm(1,0,sd=sqrt(1/(1-(parameters$std_F)^2)))
  X             <- c(1:parameters$time)
  X_obs         <- c(1:parameters$time)
  epsilon       <- rnorm(parameters$time, 0, 1)
  eta           <- rnorm(parameters$time, 0, sd = parameters$std_eta)
  
  X[1]          <- universe$F[1] *X_0 + exp(universe$Sigma[1]) * epsilon[1]
  X_obs[1]      <- X[1] + eta[1]
  
  for(i in (2:parameters$time)){
    X[i] <- universe$F[i]*X[i-1] + exp(universe$Sigma[i]) * epsilon[i]
    X_obs[i] <- X[i] + eta[i]
  }
  return(data.frame(X, X_obs))
}



#====================================================================
# # Here, we specify additional parameters for the KF.

parameters_kf <- function(){
  B_f_0 <- 6.97                  # The assumed B_f at the start of the filtering 
  #  (close to the steady-state background-error variance).
  return(data.frame(B_f_0))
}



#====================================================================
# This function runs the KF.

filter_kf <- function(world, universe, parameters, parameters_kf){
  set.seed(parameters$seed_for_filters)
  # Kalman filter
  
  B_f <- c(1:parameters$time)
  A   <- c(1:parameters$time)
  K   <- c(1:parameters$time)
  X_a <- c(1:parameters$time)
  X_f <- c(1:parameters$time)
  
  
  A_0    <- parameters$std_eta^2 * parameters_kf$B_f_0 / (parameters_kf$B_f_0 + parameters$std_eta^2)
  K_0    <- parameters_kf$B_f_0 / (parameters_kf$B_f_0 + parameters$std_eta^2)
  
  X_f[1] <- universe$F_0[1] * rnorm(1,0,sd=sqrt(1/(1-(parameters$std_F)^2)))
  B_f[1] <- universe$F[1] * A_0 * universe$F[1] + universe$Q[1]
  A[1]   <- parameters$std_eta^2 * B_f[1] / (B_f[1] + parameters$std_eta^2)
  K[1]   <- B_f[1] / (B_f[1] + parameters$std_eta^2)
  X_a[1] <- X_f[1] + K[1]*(world$X_obs[1] - X_f[1])  
  
  for(i in (2:parameters$time)){
    X_f[i] <- universe$F[i] * X_a[i-1]
    B_f[i] <- universe$F[i] * A[i-1] * universe$F[i] + universe$Q[i]
    A[i] <- parameters$std_eta^2 * B_f[i] / (B_f[i] + parameters$std_eta^2)
    K[i]   <- B_f[i] / (B_f[i] + parameters$std_eta^2)
    X_a[i] <- X_f[i] + K[i]*(world$X_obs[i] - X_f[i]) 
  }
  B_a <- B_f
  return(data.frame(X_a, X_f, A, B_f, B_a))
}



#====================================================================
# Here, we specify additional parameters for Var.

parameters_var <- function(){
  mean_B <- 6.97        # The assumed (constant) B at the start of the filtering
  #  (close to the steady-state KF's background-error variance).
  return(data.frame(mean_B))
}



#====================================================================
# This function runs the Var filter.

filter_var <- function(world, universe, parameters, parameters_var){
  
  set.seed(parameters$seed_for_filters)
  
  # Variational assimilation filter
  
  B_var <- parameters_var$mean_B
  B_f <- rep(B_var, parameters$time)
  B_a <- B_f
  X_a <- c(1:parameters$time)
  X_f <- c(1:parameters$time)
  K_var   <- B_var / (B_var + parameters$std_eta^2)
  X_f[1] <- universe$F_0[1] * rnorm(1,0,sd=sqrt(1/(1-(parameters$std_F)^2)))
  for(i in (1:parameters$time)){
    if(i>1)  X_f[i] <- universe$F[i] * X_a[i-1]
    X_a[i] <- X_f[i] + K_var*(world$X_obs[i] - X_f[i])
  }
  return(data.frame(X_a, X_f, B_f, B_a))
}



#====================================================================
# Here, we specify additional parameters for the EnKF.

parameters_enkf <-function(){
  inflation             <- 1.02               # parameter of ensemble pertubations' inflation
  return(data.frame(inflation))
}



#====================================================================
# This function runs the EnKF.

filter_enkf <- function(world, universe, parameters, parameters_enkf){
  set.seed(parameters$seed_for_filters)
  
  # Ensemble Kalman Filter
  
  epsilon_en    <- matrix(rnorm(parameters$time * parameters$N, 0, 1), ncol = parameters$time, nrow = parameters$N)
  eta_en        <- matrix(rnorm(parameters$time * parameters$N, 0, sd = parameters$std_eta), ncol = parameters$time, nrow = parameters$N)
  
  X_en   <- matrix(ncol = parameters$time, nrow = parameters$N)
  X_en_f <- matrix(ncol = parameters$time, nrow = parameters$N)
  X_en_a <- matrix(ncol = parameters$time, nrow = parameters$N)
  
  X_en_0    <- rnorm(parameters$N,0,sd=sqrt(1/(1-(parameters$std_F)^2)))
  
  
  B_f <- c(1:parameters$time)
  B_a <- c(1:parameters$time)
  K   <- c(1:parameters$time)
  X_a <- c(1:parameters$time)
  X_f <- c(1:parameters$time)
  
  X_en_f[,1]  <- universe$F[1] * X_en_0 + sqrt(universe$Q[1]) * epsilon_en[,1]
  X_f[1]  <- universe$F[1] * rnorm(1,0,sd=sqrt(1/(1-(parameters$std_F)^2))) 
  X_en_f[,1]  <- rep(X_f[1], parameters$N) + (X_en_f[,1] - rep(X_f[1], parameters$N))*parameters_enkf$inflation
  B_f[1]  <- var(X_en_f[,1])
  K[1]   <- B_f[1] / (B_f[1] + parameters$std_eta^2)
  X_a[1] <- X_f[1] + K[1]*(world$X_obs[1] - X_f[1])  
  X_en_a[,1] <- X_en_f[,1] + K[1]*(rep(world$X_obs[1],parameters$N)+eta_en[,1]-X_en_f[,1])
  
  for(i in (2:parameters$time)){
    X_en_f[,i]  <- universe$F[i] * X_en_a[,i-1] + sqrt(universe$Q[i]) * epsilon_en[,i]
    X_f[i] <- universe$F[i] * X_a[i-1]
    X_en_f[,i]  <- rep(X_f[i], parameters$N) + (X_en_f[,i] - rep(X_f[i], parameters$N))*parameters_enkf$inflation
    B_f[i]  <- var(X_en_f[,i])
    K[i]   <- B_f[i] / (B_f[i] + parameters$std_eta^2)
    X_a[i] <- X_f[i] + K[i]*(world$X_obs[i] - X_f[i])  
    X_en_a[,i] <- X_en_f[,i] + K[i]*(rep(world$X_obs[i],parameters$N) + eta_en[,i]-X_en_f[,i]) 
  }
  B_a <- B_f
  return(data.frame(X_f, X_a, B_f, B_a))
}



#====================================================================
# Here, we specify additional parameters for the HEnKF.

parameters_henkf <- function(){
  mean_B        <- as.matrix(6.97)       # climatological  (time-mean) B
  theta         <- 10                    # dispersion parameter for the prior Inverse Gamma distribution of B
  return(data.frame(mean_B, theta))
}



#====================================================================
# This function runs the HEnKF.

filter_henkf <- function(world, universe, parameters, parameters_henkf){
  set.seed(parameters$seed_for_filters)
  
  # Hierarchical Ensemble Kalman Filter
  
  x_en_f  <- matrix(ncol = parameters$time, nrow = parameters$N)   # ensemble forecast
  x_en_a  <- matrix(ncol = parameters$time, nrow = parameters$N)   # ensemble forecast
  X_a     <- c(1:parameters$time)                                  # analysis with approximation
  B_f     <- c(1:parameters$time)                                  # background B
  B_a     <- c(1:parameters$time)                                  # analysis (posterior) B
  X_f     <- c(1:parameters$time)
  ## === generation of ensemble ====
  B_f[1] <- parameters_henkf$mean_B
  B_a[1] <- parameters_henkf$mean_B
  epsilon_en    <- matrix(rnorm(parameters$time * parameters$N, 0, 1), ncol = parameters$time, nrow = parameters$N)  
  x_en_f[,1] <- rnorm(parameters$N, 0, sqrt(B_f[1]))
  ## ===============================
  for(i in (1:parameters$time)){
    if(i == 1){
      B_f[i] <- parameters_henkf$mean_B
    }else{
      B_f[i] <- (1/B_a[i-1]+1/parameters$std_eta^2)^{-1}
    }  
    
    if(i==1){
      m_f <- 0
    }else{
      m_f <- universe$F[i] * X_a[i-1]  
    }
    X_f[i] <- m_f
    
    m_C <- function(C_arg){
      D <- C_arg + 1/(parameters$std_eta^2)
      result <- m_f + (world$X_obs[i] - m_f) / (D * parameters$std_eta^2) 
      return(result)
    }
    
    theta_hat   <- parameters_henkf$theta + parameters$N
    S           <- sum((x_en_f[,i]-rep(m_f, parameters$N))^2) / parameters$N
    B_a[i]      <- (parameters_henkf$theta*B_f[i] + parameters$N * S)/(parameters_henkf$theta + parameters$N) 
    X_a[i]      <- m_C(B_a[i]^(-1))
    
    if(i < parameters$time){
      C_sample <- rWISHART(parameters$N, df = theta_hat + 2,  S = (B_a[i] * theta_hat)^{-1})
      for(j in (1:parameters$N)){
        x_en_a[j,i] <- rnorm(1, m = m_C(C_sample[j]), sd = (C_sample[j]+1/(parameters$std_eta^2))^(-1/2))
      }
    }  
    if(i < parameters$time){
      x_en_f[,i+1] <- universe$F[i+1] * x_en_a[,i] + sqrt(universe$Q[i+1]) * epsilon_en[,i+1]
    }
  }
  return(data.frame(X_a, X_f, B_a, B_f))
}



#====================================================================
# Here, we specify additional parameters for the HBEF.

parameters_hbef <- function(){
  theta         <- 4                # dispersion parameter for the Inverse Gamma distribution P|Pi
  size_for_MC   <- 500              # size of the Monte-Carlo sample used to estimate the posterir mean m^a
  phi           <- 20               # dispersion parameter for the Inverse Gamma distribution Pi|Pi^f
  chi           <- 9                # dispersion parameter for the Inverse Gamma distribution Q|Q^f
  mean_A        <- 6.15             # starting value for the analysis-error variance A
  mean_Q        <- 1.64             # starting value for the model-error variance Q
  approximation <- FALSE            # use the approximated posterior (TRUE) or not (FALSE) ?
  use_L_o       <- FALSE             
  return(data.frame(theta, size_for_MC, phi, chi, mean_A, mean_Q, approximation, use_L_o))
}



#====================================================================
# This function runs the HBEF.

filter_hbef <- function(world, universe, parameters, parameters_hbef){
  ## ==========================================================
  ## ====================== HBEF ver 2 ========================
  ## ==========================================================
  # Hierarchical Bayes Ensemble Filter
  
  # Initialization (starting conditions)
  
  set.seed(parameters$seed_for_filters)
  
  x_en_f  <- matrix(ncol = parameters$time, nrow = parameters$N)   # ensemble forecast
  X_a     <- c(1:parameters$time)                                  # analysis with the approximated posterior
  X_f     <- c(1:parameters$time)
  B_f     <- c(1:parameters$time)                                  # B_f (B_tilde)
  B_a     <- c(1:parameters$time)                                  # B_a
  Q_f     <- c(1:parameters$time)
  Q_tilde <- c(1:parameters$time)
  Q_a     <- c(1:parameters$time)               
  Pi_f         <- c(1:parameters$time)
  Pi_tilde     <- c(1:parameters$time)
  P_a         <- c(1:parameters$time)
  D            <- c(1:parameters$time)
  S_pe         <- c(1:parameters$time)
  S_me         <- c(1:parameters$time)
  A            <- c(1:parameters$time)
  Pi           <- c(1:parameters$time)
  ## ================== Parameters of the Filter ===============
  ## 
  
  
  ## === generation of the ensembles ====
  epsilon_en    <- matrix(rnorm(parameters$time * parameters$N, 0, 1), ncol = parameters$time, nrow = parameters$N)  # X^{me}
  x_en_f[,1]    <- universe$F[1] * rnorm(parameters$N, 0, sqrt(parameters_hbef$mean_A)) ## for PI                    # X^{pe}
  ## ===============================
  
  A[1] <- parameters_hbef$mean_A
  Pi[1] <- universe$F[1]^2*A[1]
  
  for(i in (1:parameters$time)){
    
    if(i==1){
      m_f <- 0
      Pi_f[i]  <- (mean(universe$F))^2*parameters_hbef$mean_A
    }else{
      m_f <- universe$F[i] * X_a[i-1]
      Pi_f[i] <- P_a[i-1]
    }
    
    X_f[i] <- m_f
    
    m_C <- function(C_arg){
      D_m <- C_arg + 1/(parameters$std_eta^2)
      result <- m_f + (world$X_obs[i] - m_f) / (D_m * parameters$std_eta^2) 
      return(result)
    }
    
    S_me[i]     <- universe$Q[i] * mean(epsilon_en[,i]^2)
    S_pe[i]       <- mean(x_en_f[,i]^2)
    v          <- world$X_obs[i] - m_f
    
    if((parameters_hbef$approximation == FALSE) && (parameters_hbef$use_L_o == TRUE)){
      Pi_tilde[i] <- (parameters_hbef$phi * Pi_f[i] + parameters$N * S_pe[i])/(parameters_hbef$phi + parameters$N)
      if(i==1){
        Q_a[i] <- (parameters_hbef$chi*parameters_hbef$mean_Q + parameters$N*S_me[i])/(parameters_hbef$chi+parameters$N)
        Q_tilde[i] <- parameters_hbef$mean_Q
        Q_f[i] <- parameters_hbef$mean_Q
      }else{
        Q_f[i] <- Q_a[i-1]
        Q_tilde[i] <- (parameters_hbef$chi*Q_a[i-1] + parameters$N*S_me[i])/(parameters_hbef$chi+parameters$N)
        Q_samp <- rinvgamma(parameters_hbef$size_for_MC, shape = (parameters_hbef$chi+parameters$N)/2 + 1,  scale = Q_tilde[i] * (parameters_hbef$chi+parameters$N)/2)
        L_o <- function(Q){
          T <- Pi_tilde[i] + Q + parameters$std_eta^2
          result <- 1/sqrt(T) * exp(-v^2/(2*T))
        }
        L_o_samp <- L_o(Q_samp)
        Q_a[i] <- L_o_samp%*%Q_samp / sum(L_o_samp)
      }                                            
      P_samp <- rinvgamma(parameters_hbef$size_for_MC, shape = parameters_hbef$theta/2 + 1,  scale = Pi_tilde[i] * parameters_hbef$theta/2)
      L_o <- function(P){
        T <- P + Q_tilde[i] + parameters$std_eta^2
        result <- 1/sqrt(T) * exp(-v^2/(2*T))
      }
      L_o_samp <- L_o(P_samp)
      P_a[i] <- L_o_samp%*%P_samp / sum(L_o_samp)
    }
    
    if(parameters_hbef$approximation == TRUE && parameters_hbef$use_L_o == TRUE){
      Pi_tilde[i] <- (parameters_hbef$phi * Pi_f[i] + parameters$N * S_pe[i])/(parameters_hbef$phi + parameters$N)
      if(i==1){
        Q_a[i] <- (parameters_hbef$chi*parameters_hbef$mean_Q + parameters$N*S_me[i])/(parameters_hbef$chi+parameters$N)
        Q_tilde[i] <- parameters_hbef$mean_Q
        Q_f[i] <- parameters_hbef$mean_Q
      }else{
        Q_f[i] <- Q_a[i-1]
        Q_tilde[i] <- (parameters_hbef$chi*Q_f[i] + parameters$N*S_me[i])/(parameters_hbef$chi+parameters$N)
        Q_a[i] <- Q_tilde[i] + 1/ (parameters_hbef$chi + parameters$N)*(Q_tilde[i]/(Pi_tilde[i]+Q_tilde[i]+parameters$std_eta^2))^2*(v^2-(Pi_tilde[i]+Q_tilde[i]+parameters$std_eta^2))
      }
      P_a[i] <- Pi_tilde[i] + 1/parameters_hbef$theta*(Pi_tilde[i]/(Pi_tilde[i]+Q_tilde[i]+parameters$std_eta^2))^2*(v^2-(Pi_tilde[i]+Q_tilde[i]+parameters$std_eta^2))
    }
    
    if(parameters_hbef$use_L_o == FALSE){
      Pi_tilde[i] <- (parameters_hbef$phi * Pi_f[i] + parameters$N * S_pe[i])/(parameters_hbef$phi + parameters$N)
      if(i==1){
        Q_a[i] <- (parameters_hbef$chi*parameters_hbef$mean_Q + parameters$N*S_me[i])/(parameters_hbef$chi+parameters$N)
        Q_tilde[i] <- parameters_hbef$mean_Q
        Q_f[i] <- parameters_hbef$mean_Q
      }else{
        Q_f[i] <- Q_a[i-1]
        Q_tilde[i] <- (parameters_hbef$chi*Q_f[i] + parameters$N*S_me[i])/(parameters_hbef$chi+parameters$N)
        Q_a[i] <- Q_tilde[i]
      }
      P_a[i] <- Pi_tilde[i]
    }
    
    B_a[i]     <- P_a[i] + Q_a[i]
    B_f[i]     <- Pi_f[i] + Q_f[i]
    D[i]       <- (B_a[i])^(-1) + 1/(parameters$std_eta^2)
    X_a[i]     <- m_C(B_a[i]^(-1))
    
    if(i < parameters$time){
      x_en_f[,i+1] <- universe$F[i+1] * rnorm(parameters$N, 0, (D[i])^{-1/2})
      Pi[i+1] <- universe$F[i+1]^2*(D[i])^{-1}
    }
  }
  return(data.frame(X_f, X_a, Pi_f, Q_f, P_a, Q_a, B_f, B_a, Pi, Pi_tilde, Q_tilde, S_pe, S_me))
}

rms <- function(x) sqrt(sum(x^2)/length(x))  # calc RMS
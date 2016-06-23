################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
#
# This file contains the R functions that realize the main algorithms presented in the paper
#
# ``A Hierarchical Bayes Ensemble Kalman Filter''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
#
#
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 21 Dec 2015  Tsy
## 29 May 2016  Tsy
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
  pi        <- 0.05                      # Pr(F > 1). pi and F_0 define std_F and var_epsilon_F (nrm =0.05)
  #
  Sigma_0   <- 0                         # mean of Sigma (this is rather a constant than a parameter, don't change it)
  #
  std_Sigma <- 0.5                       # standart deviation of Sigma (nrm =0.5)
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
  
  time              <- 200000          # number of time moments in the time series (n_{time})
  N                 <- 5               # ensemble size    
  L                 <- 1               # number of independent runs (worlds) under the same {F_k} and {Sigma_k} (universe)
  distort_Q         <- 1               # coefficient of Q distortion (Q_flt <- Q_true * distort_Q)
  # ------------------------------
  seed_for_universe1 <- 42433425 # random seed for initiation the {F_k} time series
  seed_for_universe2 <- 91435564 # random seed for initiation the {SIgma_k} time series
  seed_for_world     <- 24674353 # random seed for initiation of world
  seed_for_filters   <- 56345352 # random seed for initiation of other pseudo-random sources
  
  # list of parameters
  list <- data.frame(seed_for_universe1, seed_for_universe2, seed_for_filters, seed_for_world, distort_Q, 
                     time, tau_x, tau_F, std_F, mu, var_epsilon_F, pi, tau_Sigma, 
                     std_Sigma, kappa, var_epsilon_Sigma, F_0, Sigma_0, std_eta, N, L)
  print(list)
  return(list)
}


#====================================================================
# This function updates internal parameters from the external ones.

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
  # Sigma_k=kappa*Sigma_{k-1} + \sqrt(var_epsilon_Sigma) * epsilon_Sigma_k
  
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
  Q <- (exp(Sigma))^2
  return(data.frame(F, Sigma, Q, F_0))
}



#====================================================================
# This function simulates a realization of the "truth" {x_k} and observations {y_k}.

generate_world <- function(universe, parameters){
  ## =====================================================================================
  ## ==================== Generating of the world (and two ensemble feeds) ===============
  ## =====================================================================================
  set.seed(parameters$seed_for_world)
  
  X_0           <- rnorm(1,0,sd=sqrt(1/(1-(parameters$std_F)^2)))
  X             <- c(1:parameters$time)
  X_obs         <- c(1:parameters$time)
  epsilon       <- rnorm(parameters$time, 0, 1)
  eta           <- rnorm(parameters$time, 0, sd = parameters$std_eta)
  
  X[1]          <- universe$F[1] *X_0 + exp(universe$Sigma[1]) * epsilon[1]
  X_obs[1]      <- X[1] + eta[1]
  
  # truth and obs
  
  for(i in (2:parameters$time)){
    X[i] <- universe$F[i]*X[i-1] + exp(universe$Sigma[i]) * epsilon[i]
    X_obs[i] <- X[i] + eta[i]
  }
  
  # ensembles of simulated model errors (epsilon_en) and simulated obs err (eta_en)

  epsilon_en <- matrix(rnorm(parameters$time * parameters$N, 0, 1), ncol = parameters$time, nrow = parameters$N)
  for(i in (1:parameters$time)){
    epsilon_en[,i] <- epsilon_en[,i] * sqrt(universe$Q[i]) * sqrt(parameters$distort_Q)
  }
   
  eta_en     <- matrix(rnorm(parameters$time * parameters$N, 0, sd = parameters$std_eta), 
                          ncol = parameters$time, nrow = parameters$N)
  
  return(list(X=X, X_obs=X_obs, epsilon_en=epsilon_en, eta_en=eta_en))
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
  
  X_f[1] <- 0 # universe$F_0[1] * rnorm(1,0,sd=sqrt(1/(1-(parameters$std_F)^2)))
  B_f[1] <- universe$F[1] * A_0 * universe$F[1] + universe$Q[1]
  A[1]   <- parameters$std_eta^2 * B_f[1] / (B_f[1] + parameters$std_eta^2)
  K[1]   <- B_f[1] / (B_f[1] + parameters$std_eta^2)
  X_a[1] <- X_f[1] + K[1]*(world$X_obs[1] - X_f[1])  
  
  for(i in (2:parameters$time)){
    X_f[i] <- universe$F[i] * X_a[i-1]
    B_f[i] <- universe$F[i] * A[i-1] * universe$F[i] + universe$Q[i]*parameters$distort_Q
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
  
  # Before running this filter, run the KF and replace parameters_var$mean_B by
  # the mean_B from the KF.
  
  B_var <- parameters_var$mean_B
  B_f <- rep(B_var, parameters$time)
  B_a <- B_f
  X_a <- c(1:parameters$time)
  X_f <- c(1:parameters$time)
  K_var   <- B_var / (B_var + parameters$std_eta^2)
  X_f[1] <- 0 # universe$F_0[1] * rnorm(1,0,sd=sqrt(1/(1-(parameters$std_F)^2)))
  for(i in (1:parameters$time)){
    if(i>1)  X_f[i] <- universe$F[i] * X_a[i-1]
    X_a[i] <- X_f[i] + K_var*(world$X_obs[i] - X_f[i])
  }
  return(data.frame(X_a, X_f, B_f, B_a))
}



#====================================================================
# Here, we specify additional parameters for the EnKF.

parameters_enkf <-function(){
  inflation     <- 1.005    # parameter of ensemble pertubations' inflation (applied to perturbations, not variance)
  mean_A        <- 6.0      # default starting value for the analysis-error variance A
  return(list(inflation=inflation, mean_A=mean_A))
}



#====================================================================
# This function runs the EnKF.

filter_enkf <- function(world, universe, parameters, parameters_enkf){
  set.seed(parameters$seed_for_filters)
  
  # Ensemble Kalman Filter
  
  # Before running this filter, run the KF and replace parameters_enkf$mean_A by
  # the mean_B from the KF.
  
  R <- parameters$std_eta^2
  
  X_en_f <- matrix(ncol = parameters$time, nrow = parameters$N)
  X_en_a <- matrix(ncol = parameters$time, nrow = parameters$N)
  
  X_me <- matrix(ncol = parameters$time, nrow = parameters$N)
  X_pe <- matrix(ncol = parameters$time, nrow = parameters$N)
  
  B_f <- c(1:parameters$time)
  B_a <- c(1:parameters$time)
  K   <- c(1:parameters$time)
  X_a <- c(1:parameters$time)
  X_f <- c(1:parameters$time)
  
  X_pe[,1]   <- 0 
  X_me[,1]   <- world$epsilon_en[,1]
  X_en_f[,1] <- X_pe[,1] + X_me[,1]
  x_fe_mean  <- mean(X_en_f[,1])
  X_f[1]     <- 0         # det fc
  #X_f[1]  <- x_fe_mean # ensm mean
  #X_en_f[,1]  <- rep(X_f[1], parameters$N) + (X_en_f[,1] - rep(X_f[1], parameters$N))*parameters_enkf$inflation
  X_en_f[,1]  <- rep(x_fe_mean, parameters$N) + (X_en_f[,1] - rep(x_fe_mean, parameters$N))*parameters_enkf$inflation
  B_f[1]  <- var(X_en_f[,1])
  K[1]   <- B_f[1] / (B_f[1] + R)
  
  X_en_a[,1] <- X_en_f[,1] + K[1]*(rep(world$X_obs[1],parameters$N) + world$eta_en[,1]-X_en_f[,1])
  X_a[1] <- X_f[1] + K[1]*(world$X_obs[1] - X_f[1]) 
  #X_a[1] <- mean(X_en_a[,1]) 
  
  for(i in (2:parameters$time)){
    X_pe[,i]   <- universe$F[i] * X_en_a[,i-1] 
    X_me[,i]   <- world$epsilon_en[,i]
    X_en_f[,i] <- X_pe[,i] + X_me[,i]  # bef infl
    x_fe_mean  <- mean(X_en_f[,i])
    
    X_f[i] <- universe$F[i] * X_a[i-1]  # det fc
    #X_f[i] <- x_fe_mean                 # ensm mean
    X_en_f[,i]  <- rep(x_fe_mean, parameters$N) + (X_en_f[,i] - rep(x_fe_mean, parameters$N))*parameters_enkf$inflation
    
    B_f[i]  <- var(X_en_f[,i])  # of the inflated ensm
    K[i]   <- B_f[i] / (B_f[i] + R)
    
    X_a[i] <- X_f[i] + K[i]*(world$X_obs[i] - X_f[i])  # det anls
    #X_a[i] <- mean(X_en_a[,i])                         # ensm mean
    
    # Stoch EnKF
    
    #X_en_a[,i] <- X_en_f[,i] + K[i]*(rep(world$X_obs[i],parameters$N) + world$eta_en[,i]-X_en_f[,i])
    
    # ETKF
    # anls perturbations f_j=c*e_j, where e_j are the fcst perturbtions
    # A=c^2*B_f   =>  c=sqrt(A/B_f)
    
    A <- K[i] * R    # anls err var
    c <- sqrt(A/B_f[i])
    X_en_a[,i] <- rep(X_a[i],parameters$N) + c*(X_en_f[,i] - rep(X_f[i],parameters$N)) 
    
  }
  B_a <- B_f
  return(list(X_f=X_f, X_a=X_a, B_f=B_f, B_a=B_a, X_pe=X_pe, X_me=X_me, X_en_f=X_en_f, 
              X_en_a=X_en_a, K=K ))
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
  
  # Before running this filter, run the KF and replace parameters_henkf$mean_B by
  # the mean_B from the KF.
  
  X_en_f  <- matrix(ncol = parameters$time, nrow = parameters$N)   # ensemble forecast
  X_en_a  <- matrix(ncol = parameters$time, nrow = parameters$N)   # ensemble anls
  X_a     <- c(1:parameters$time)                                  # det.analysis 
  B_f     <- c(1:parameters$time)                                  # prior mean B
  B_a     <- c(1:parameters$time)                                  # analysis (posterior) mean B
  X_f     <- c(1:parameters$time)                                  # det.fcst 
  
  R <- parameters$std_eta^2

  B_f[1] <- parameters_henkf$mean_B
  B_a[1] <- parameters_henkf$mean_B
  
  X_en_f[,1] <- rnorm(parameters$N, 0, sqrt(B_f[1])) 
  ## ===============================
  for(i in (1:parameters$time)){
    if(i == 1){
      B_f[i] <- parameters_henkf$mean_B
    }else{
      B_f[i] <- (1/B_a[i-1]+1/R)^{-1}  # =A_hat[i-1], following Myrseth & Omre.
    }  
    
    if(i==1){
      m_f <- 0
    }else{
      m_f <- universe$F[i] * X_a[i-1]  
    }
    x_fe_mean  <- mean(X_en_f[,i])
    
    X_f[i] <- m_f
    #X_f[i] <- x_fe_mean
    
    m_C <- function(C_arg){
      D <- C_arg + 1/R
      result <- X_f[i] + (world$X_obs[i] - X_f[i]) / (D * R) 
      return(result)
    }
    
    theta_hat   <- parameters_henkf$theta + parameters$N
    #S           <- sum((X_en_f[,i]-rep(m_f, parameters$N))^2) / parameters$N
    S           <- var(X_en_f[,i])
    B_a[i]      <- (parameters_henkf$theta*B_f[i] + parameters$N * S)/(parameters_henkf$theta + parameters$N) 
    X_a[i]      <- m_C(B_a[i]^(-1))
    
    # preparations for the next step
    # anls ensm
    
    if(i < parameters$time){
      C_sample <- rWISHART(parameters$N, df = theta_hat + 2,  S = (B_a[i] * theta_hat)^{-1})
      K_sample <- (1 + C_sample[] *R)^(-1)
      for(j in (1:parameters$N)){
        X_en_a[j,i] <- (1 - K_sample[j])*X_en_f[j,i] + K_sample[j]*(world$X_obs[i] + world$eta_en[j,i])
      }
    }  
    
    # next-step fc ensm
    
    if(i < parameters$time){
      X_en_f[,i+1] <- universe$F[i+1] * X_en_a[,i] + world$epsilon_en[,i+1]
    }
  }
  return(list(X_a=X_a, X_f=X_f, B_a=B_a, B_f=B_f, X_en_f=X_en_f))
}



#====================================================================
# Here, we specify additional parameters for the HBEF.

parameters_hbef <- function(){
  theta         <- 2                # dispersion parameter for the Inverse Gamma distribution P|Pi
  size_for_MC   <- 100              # size of the Monte-Carlo sample used to estimate the posterir mean m^a
  phi           <- 30               # dispersion parameter for the Inverse Gamma distribution Pi|Pi^f
  chi           <- 5                # dispersion parameter for the Inverse Gamma distribution Q|Q^f
  mean_A        <- 6                # starting value for the analysis-error variance A
  mean_Q        <- 2                # starting value for the model-error variance Q
  approximation <- FALSE            # use the approximated posterior (TRUE) or not (FALSE) ?
  use_L_o       <- FALSE            # multiplication by L_o(B) in the posterior: TRUE if yes.
  return(data.frame(theta, size_for_MC, phi, chi, mean_A, mean_Q, approximation, use_L_o))
}



#====================================================================
# This function runs the HBEF.

filter_hbef <- function(world, universe, parameters, parameters_hbef){
  ## ==========================================================
  ## ====================== HBEF ver 2 ========================
  ## ==========================================================
  # Hierarchical Bayes Ensemble Filter
  
  # Before running this filter, run the KF and replace 
  #  parameters_hbef$mean_A and parameters_hbef$mean_Q
  # by the respective values from the KF.

  # Initialization (starting conditions)
  
  set.seed(parameters$seed_for_filters)
  
  X_me <- matrix(ncol = parameters$time, nrow = parameters$N)
  X_pe <- matrix(ncol = parameters$time, nrow = parameters$N)
  X_ae <- matrix(ncol = parameters$time, nrow = parameters$N)
  X_a     <- c(1:parameters$time)                                  # analysis with the approximated posterior
  X_f     <- c(1:parameters$time)
  B_f     <- c(1:parameters$time)                                  # B_f (B_tilde)
  B_a     <- c(1:parameters$time)                                  # B_a
  Q_f     <- c(1:parameters$time)
  Q_tilde <- c(1:parameters$time)
  Q_a     <- c(1:parameters$time)   
  K       <- c(1:parameters$time)
  P_f         <- c(1:parameters$time)
  P_tilde     <- c(1:parameters$time)
  P_a         <- c(1:parameters$time)
  D            <- c(1:parameters$time)
  S_pe         <- c(1:parameters$time)
  S_me         <- c(1:parameters$time)
  A            <- c(1:parameters$time)
  Pi           <- c(1:parameters$time)
  ## ================== Parameters of the Filter ===============
  
  R <- parameters$std_eta^2
  ## ===============================
  X_pe[,1]    <- 0 # universe$F[1] * rnorm(parameters$N, 0, sqrt(parameters_hbef$mean_A)) ## for PI   #  ??!! no sd=                 # X^{pe} for t=1
  
  A[1] <- parameters_hbef$mean_A
  Pi[1] <- universe$F[1]^2*A[1]
  
  for(i in (1:parameters$time)){
    
    if(i==1){
      m_f <- 0
      Q_f[i] <- parameters_hbef$mean_Q
      P_f[i]  <- (mean(universe$F))^2*parameters_hbef$mean_A
    }else{
      m_f <- universe$F[i] * X_a[i-1]
      Q_f[i] <- Q_a[i-1]
      P_f[i] <- P_a[i-1]
    }
    
    B_f[i] <- P_f[i] + Q_f[i]
    X_f[i] <- m_f
    v      <- world$X_obs[i] - m_f
    
    m_B <- function(B_arg){
      K_B <- B_arg / (B_arg + R)
      result <- m_f + K_B * v  
      return(result)
    }
    
    X_me[,i] <- world$epsilon_en[,i]
    # X_pe[,i] is computed at the end of the previous cycle (or set to be =0 at i=1)
    
    S_me[i] <- mean(X_me[,i]^2) # var(X_me[,i]) 
    S_pe[i] <- var(X_pe[,i])

    # MCarlo
    
    if((parameters_hbef$approximation == FALSE) && (parameters_hbef$use_L_o == TRUE)){
      P_tilde[i] <- (parameters_hbef$phi * P_f[i] + parameters$N * S_pe[i])/(parameters_hbef$phi + parameters$N)
      if(i==1){
        Q_a[i] <- (parameters_hbef$chi*parameters_hbef$mean_Q + parameters$N*S_me[i])/(parameters_hbef$chi+parameters$N)
        Q_tilde[i] <- parameters_hbef$mean_Q
        P_samp <- rinvgamma(parameters_hbef$size_for_MC, shape = parameters_hbef$theta/2 + 1,  scale = P_tilde[i] * parameters_hbef$theta/2)
        L_o <- function(P){
          T <- P + Q_tilde[i] + parameters$std_eta^2
          result <- 1/sqrt(T) * exp(-v^2/(2*T))
        }
        L_o_samp <- L_o(P_samp)
        c <- 1 / sum(L_o_samp)
        P_a[i] <- c * t(L_o_samp)%*%P_samp 
        
        X_a_samp <- m_B(P_samp[] + Q_tilde[i])
        X_a[i] <- c * t(L_o_samp)%*%X_a_samp 
        
      }else{  # i>1  
        Q_tilde[i] <- (parameters_hbef$chi*Q_a[i-1] + parameters$N*S_me[i])/(parameters_hbef$chi+parameters$N)
        Q_samp <- rinvgamma(parameters_hbef$size_for_MC, shape = (parameters_hbef$chi+parameters$N)/2 + 1,  scale = Q_tilde[i] * (parameters_hbef$chi+parameters$N)/2)
        P_samp <- rinvgamma(parameters_hbef$size_for_MC, shape =  parameters_hbef$theta/2 + 1,              scale = P_tilde[i] *  parameters_hbef$theta/2)
        L_o <- function(X){
          T <- X + parameters$std_eta^2
          result <- 1/sqrt(T) * exp(-v^2/(2*T))
        }
        L_o_samp <- L_o(Q_samp+P_samp)
        c <- 1 / sum(L_o_samp)
        
        X_a_samp <- m_B(P_samp[] + Q_samp[])
        
        Q_a[i] <- c * t(L_o_samp)%*%Q_samp 
        P_a[i] <- c * t(L_o_samp)%*%P_samp 
        X_a[i] <- c * t(L_o_samp)%*%X_a_samp
      }                                            
      
      B_a[i]     <- P_a[i] + Q_a[i]
      B <- B_a[i]
      K[i] <- B/ (B+R)
      
      #X_a[i] <- X_f[i] + K[i] * v  # an alternative MC formulation (only P,Q are sampled)
      
    }
    
    # APPROX
    
    if(parameters_hbef$approximation == TRUE && parameters_hbef$use_L_o == TRUE){
      P_tilde[i] <- (parameters_hbef$phi * P_f[i] + parameters$N * S_pe[i])/(parameters_hbef$phi + parameters$N)
      if(i==1){
        Q_a[i] <- (parameters_hbef$chi*parameters_hbef$mean_Q + parameters$N*S_me[i])/(parameters_hbef$chi+parameters$N)
        Q_tilde[i] <- parameters_hbef$mean_Q
            
  }else{
        Q_tilde[i] <- (parameters_hbef$chi*Q_f[i] + parameters$N*S_me[i])/(parameters_hbef$chi+parameters$N)
        Q_a[i] <- Q_tilde[i] + 1/ (parameters_hbef$chi + parameters$N)*(Q_tilde[i]/(P_tilde[i]+Q_tilde[i]+parameters$std_eta^2))^2*(v^2-(P_tilde[i]+Q_tilde[i]+parameters$std_eta^2))
      }
      P_a[i] <- P_tilde[i] + 1/parameters_hbef$theta*(P_tilde[i]/(P_tilde[i]+Q_tilde[i]+parameters$std_eta^2))^2*(v^2-(P_tilde[i]+Q_tilde[i]+parameters$std_eta^2))
      
      B_a[i] <- P_a[i] + Q_a[i]
      B <- B_a[i]
      K[i] <- B/ (B+R)
      X_a[i] <- X_f[i] + K[i] * v
    }
    
    # SIMPLEST 
    
    if(parameters_hbef$use_L_o == FALSE){
      P_tilde[i] <- (parameters_hbef$phi * P_f[i] + parameters$N * S_pe[i])/(parameters_hbef$phi + parameters$N)
      P_a[i] <- P_tilde[i]
      
      Q_tilde[i] <- (parameters_hbef$chi * Q_f[i] + parameters$N * S_me[i])/(parameters_hbef$chi + parameters$N)
      Q_a[i] <- Q_tilde[i]
      
      B_a[i] <- P_a[i] + Q_a[i]
      B <- B_a[i]
      K[i] <- B/ (B+R)
      X_a[i] <- X_f[i] + K[i] * v
    }
    
    Q_a[i]
    P_a[i]
    B_a[i]
    D[i]       <- (B_a[i])^(-1) + 1/R
    D[i]
    A[i] <- 1/D[i]  # anls err var
    
    # Anls ensm
    
    if(i < parameters$time){
      
      #------------------- ANLS ENSM -----------------
      # Stoch-HBEF: Anls ensm pertbns: an_perbn(:) = R/(B+R) * fc_prtbn(:) + B/(B+R) * sim_obs_err
      
      X_ae[,i] <- R/(B+R) * (X_pe[,i]+ X_me[,i])  + B/(B+R) * world$eta_en[,i] # anls ensm prtbns
      #X_ae[,i] <- X_ae[,i] - rep(mean(X_ae[,i]), times=parameters$N)  # centering
      X_ae[,i] <- X_ae[,i] + X_a[i]            # add the deterministic anls
      
      #-------
      # Transform-HBEF: 
      # anls ensm perturbations f_j=c*e_j, where e_j are the fcst perturbtions
      # A=c^2*B_a   =>  c=sqrt(A/B_a)
      
      c <- sqrt(A[i] / B_a[i])
      X_ae[,i] <- rep(X_a[i],parameters$N) + c*(X_pe[,i]+ X_me[,i]) 
      #-------------------
      
      # Next-step predictability ensm
      
      X_pe[,i+1] <- universe$F[i+1] * X_ae[,i] # predictability ensm 
      
      Pi[i+1] <- universe$F[i+1]^2*A[i] # just for diagnostics
      
    }
  }
  return(list(X_f=X_f, X_a=X_a, P_f=P_f, Q_f=Q_f, P_a=P_a, Q_a=Q_a, B_f=B_f, B_a=B_a, K=K, 
              Pi=Pi, P_tilde=P_tilde, Q_tilde=Q_tilde, S_pe=S_pe, S_me=S_me, X_me=X_me, X_pe=X_pe))
}


#====================================================================


rms <- function(x) sqrt(sum(x^2)/length(x))  # calc RMS


#====================================================================
# bootstrap -- function for the confidence interval estimation of the  
#              variance of a sample
# vec -- the input sample
# nboot -- the number of the bootstrap samples desired
# conf -- level of confidence, P(X is inside the confidence interval).
# result -> two values -- the confidence interval (low, upp)

bootstrap <- function(vec, nboot, conf){
  v_samp <- 1:nboot
  m_samp <- 1:nboot
  
  p_1 <- (1 - conf) /2
  p_2 <- 1 - p_1
  for(i in 1:nboot){
    ind <- sample(1:length(vec), replace = TRUE)
    v_samp[i] <- var (vec[ind])
    m_samp[i] <- mean(vec[ind])
  }
  ci_mean=c(as.numeric(quantile(m_samp,p_1)), as.numeric(quantile(m_samp,p_2)))
  ci_var= c(as.numeric(quantile(v_samp,p_1)) ,as.numeric(quantile(v_samp,p_2)))
  return(list(ci_mean=ci_mean, ci_var= ci_var))
}


#====================================================================


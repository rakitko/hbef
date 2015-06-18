# File with definitions of functions
#
# ==================================
#
# list of functions:
#      -- create_parameters
#      -- 
#      --
#      -- 
#      --
#      -- 
#      --
#      -- 

create_parameters_universe_world <- function(){
  # initial parameters (from user)
  #
  # tau_A = 1/(ln(1/mu))                    -- e-folding time
  #
  # var_A = (var(epsilon_A)/(1-mu^2))       -- variation of A
  #
  # 
  ## ========= Parameters of the World =============================================
  tau_x     <- 12
  # e-folding time for x if A = F_0 = const
  # tau_x defines F_0 by F_0 = exp(-1/tau_x) 
  #
  tau_A     <- 18                        # e-folding time for F
  # tau_A defines mu by mu = exp(-1/tau_A)
  #
  tau_sigma <- 18                        # e-folding time for sigma
  # tau_A defines nu by nu = exp(-1/tau_sigma)
  #
  pi        <- 0.05                      # P(A > 1). pi and F_0 define std_A and var_epsilon_A
  #
  sigma_0   <- 0                         # mean of sigma
  #
  gamma     <- 0.5                       # proportional parameter for std_sigma
  #
  std_sigma <- 0.5                       # standart deviation of sigma 
  #                           
  std_eta   <- 9                      # sd of eta
  ## -------------------------------------------------------------------------------
  
  F_0               <- exp(-1/tau_x)               # mean of A
  #
  mu                <- exp(-1/tau_A)               # mu for A
  #
  std_A             <- (1 - F_0) / qnorm(1 - pi)   # standart deviation of A
  #
  var_epsilon_A     <- std_A^2*(1-mu^2)            # var of epsilon_A
  #
  nu                <- exp(-1/tau_sigma)           # mu for A
  #
  var_epsilon_sigma <- std_sigma^2*(1-nu^2)        # var of epsilon_A 
  ## ===============================================================================
  time <- 5000                     # number of realization
  prob_of_obs <- 1               # proportion of sucsessobservations
  kappa     <- 1.02
  blow      <- 1
  N         <- 5                        # ensemble size
  B_f_0     <- 10                     
  L         <- 1                    # number of iterations
  distort_Q <- 1
  # ------------------------------
  seed_for_universe1     <- 42433425
  seed_for_universe2     <- 91435564
  seed_for_filters   <- 24674353
  
  
  # list of parameters
  list <- data.frame(seed_for_universe1, seed_for_universe2, seed_for_filters, distort_Q, prob_of_obs, time, tau_x, tau_A, std_A, mu, var_epsilon_A, blow,pi,
                     tau_sigma, std_sigma, nu, var_epsilon_sigma, F_0, sigma_0, std_eta, N, B_f_0, kappa, L)
  print(list)
  return(list)
}

generate_universe <- function(parameters){
  # Recalculations
  ## -------------------------------------------------------------------------------
  
  parameters$F_0               <- exp(-1/parameters$tau_x)               # mean of A
  #
  parameters$mu                <- exp(-1/parameters$tau_A)               # mu for A
  #
  parameters$std_A             <- (1 - parameters$F_0) / qnorm(1 - parameters$pi)   # standart deviation of A
  #
  parameters$var_epsilon_A     <- parameters$std_A^2*(1-parameters$mu^2)            # var of epsilon_A
  #
  parameters$nu                <- exp(-1/parameters$tau_sigma)           # mu for A
  #
  parameters$var_epsilon_sigma <- parameters$std_sigma^2*(1-parameters$nu^2)        # var of epsilon_A 
  ## ===============================================================================
  set.seed(parameters$seed_for_universe1)
  F_0 <- rnorm(1, parameters$F_0, sd = parameters$std_A)
  Sigma_0 <- rnorm(1, parameters$sigma_0, sd = parameters$std_sigma)
  epsilon_A     <- rnorm(parameters$time, 0, sd = sqrt(parameters$var_epsilon_A))
  
  set.seed(parameters$seed_for_universe2)
  epsilon_sigma <- rnorm(parameters$time, 0, sd = sqrt(parameters$var_epsilon_sigma))
  
  F      <- c(1:parameters$time)
  Sigma  <- c(1:parameters$time)
  
  F[1]     <- parameters$F_0 + parameters$mu * (F_0 - parameters$F_0) + epsilon_A[1]
  Sigma[1] <- parameters$sigma_0 + parameters$nu * (Sigma_0 - parameters$sigma_0) + epsilon_sigma[1]
  for(i in (2:parameters$time)){
    F[i]     <- parameters$F_0 + parameters$mu * (F[i-1] - parameters$F_0) + epsilon_A[i]
    Sigma[i] <- parameters$sigma_0 + parameters$nu * (Sigma[i-1] - parameters$sigma_0) + epsilon_sigma[i]
  }
  Q <- (exp(Sigma))^2*parameters$distort_Q
  return(data.frame(F, Sigma, Q, F_0))
}

generate_world <- function(universe, parameters){
  # Recalculations
  ## -------------------------------------------------------------------------------

  parameters$F_0               <- exp(-1/parameters$tau_x)               # mean of A
  #
  parameters$mu                <- exp(-1/parameters$tau_A)               # mu for A
  #
  parameters$std_A             <- (1 - parameters$F_0) / qnorm(1 - parameters$pi)   # standart deviation of A
  #
  parameters$var_epsilon_A     <- parameters$std_A^2*(1-parameters$mu^2)            # var of epsilon_A
  #
  parameters$nu                <- exp(-1/parameters$tau_sigma)           # mu for A
  #
  parameters$var_epsilon_sigma <- parameters$std_sigma^2*(1-parameters$nu^2)        # var of epsilon_A 
  ## ===============================================================================


  ## ============================================================
  ## ==================== Generating of the world ===============
  ## ============================================================
  
  X_0 <- rnorm(1,0,sd=sqrt(1/(1-(parameters$std_A)^2)))

  X      <- c(1:parameters$time)
  X_obs  <- c(1:parameters$time)

  epsilon       <- rnorm(parameters$time, 0, 1)
  eta           <- rnorm(parameters$time, 0, sd = parameters$std_eta)
  ind_obs       <- c(1,1,1,rbinom(parameters$time - 3, 1, parameters$prob_of_obs))
  
  X[1]     <- universe$F[1] *X_0 + exp(universe$Sigma[1]) * epsilon[1]
  X_obs[1] <- X[1] + eta[1]
  for(i in (2:parameters$time)){
    X[i]     <- universe$F[i]*X[i-1] + exp(universe$Sigma[i]) * epsilon[i]
      if(ind_obs[i] == 1){
        X_obs[i] <- X[i] + eta[i]
      }else{
        X_obs[i] = NA
      }
  }
  return(data.frame(X, X_obs))
}

kalman_filter <- function(world, universe, parameters){
  # Kalman filtration
  
  B_f <- c(1:parameters$time)
  B_a <- c(1:parameters$time)
  K   <- c(1:parameters$time)
  X_a <- c(1:parameters$time)
  X_f <- c(1:parameters$time)
  
  
  B_a_0  <- parameters$std_eta^2 * parameters$B_f_0 / (parameters$B_f_0 + parameters$std_eta^2)
  K_0    <- parameters$B_f_0 / (parameters$B_f_0 + parameters$std_eta^2)
  
  X_f[1] <- universe$F_0[1] * rnorm(1,0,sd=sqrt(1/(1-(parameters$std_A)^2)))
  B_f[1] <- universe$F[1] * B_a_0 * universe$F[1] + universe$Q[1]
  B_a[1] <- parameters$std_eta^2 * B_f[1] / (B_f[1] + parameters$std_eta^2)
  K[1]   <- B_f[1] / (B_f[1] + parameters$std_eta^2)
  X_a[1] <- X_f[1] + K[1]*(world$X_obs[1] - X_f[1])  
  
  for(i in (2:parameters$time)){
    X_f[i] <- universe$F[i] * X_a[i-1]
    B_f[i] <- universe$F[i] * B_a[i-1] * universe$F[i] + universe$Q[i]
    if(is.na(world$X_obs[i])){
      B_a[i] <- B_f[i]
      K[i]   <- 0
      X_a[i] <- X_f[i]
    }else{
      B_a[i] <- parameters$std_eta^2 * B_f[i] / (B_f[i] + parameters$std_eta^2)
      K[i]   <- B_f[i] / (B_f[i] + parameters$std_eta^2)
      X_a[i] <- X_f[i] + K[i]*(world$X_obs[i] - X_f[i]) 
    }
  }
  return(data.frame(X_a, X_f, B_a, B_f))
}

var_filter <- function(world, universe, parameters, mean_B){
  B_var <- mean_B
  B_f <- rep(B_var, parameters$time)
  B_a <- B_f
  X_a <- c(1:parameters$time)
  X_f <- c(1:parameters$time)
  K_var   <- B_var / (B_var + parameters$std_eta^2)
  X_f[1] <- universe$F_0[1] * rnorm(1,0,sd=sqrt(1/(1-(parameters$std_A)^2)))
  for(i in (1:parameters$time)){
    if(i>1)  X_f[i] <- universe$F[i] * X_a[i-1]
    X_a[i] <- X_f[i] + K_var*(world$X_obs[i] - X_f[i])
  }
  return(data.frame(X_a, X_f, B_f, B_a))
}

ensemble_kalman_filter <- function(world, universe, parameters){
  # Ensemble Kalman Filter
  epsilon_en    <- matrix(rnorm(parameters$time * parameters$N, 0, 1), ncol = parameters$time, nrow = parameters$N)
  eta_en        <- matrix(rnorm(parameters$time * parameters$N, 0, sd = parameters$std_eta), ncol = parameters$time, nrow = parameters$N)
  
  X_en   <- matrix(ncol = parameters$time, nrow = parameters$N)
  X_en_f <- matrix(ncol = parameters$time, nrow = parameters$N)
  X_en_a <- matrix(ncol = parameters$time, nrow = parameters$N)
  
  X_en_0    <- rnorm(parameters$N,0,sd=sqrt(1/(1-(parameters$std_A)^2)))
  
  
  B_f <- c(1:parameters$time)
  B_a <- c(1:parameters$time)
  K   <- c(1:parameters$time)
  X_a <- c(1:parameters$time)
  X_f <- c(1:parameters$time)
  
  X_en_f[,1]  <- universe$F[1] * X_en_0 + sqrt(universe$Q[1]) * epsilon_en[,1]
  X_f[1]  <- universe$F[1] * rnorm(1,0,sd=sqrt(1/(1-(parameters$std_A)^2))) 
  X_en_f[,1]  <- rep(X_f[1], parameters$N) + (X_en_f[,1] - rep(X_f[1], parameters$N))*parameters$kappa
  B_f[1]  <- var(X_en_f[,1]) * parameters$blow
  K[1]   <- B_f[1] / (B_f[1] + parameters$std_eta^2)
  X_a[1] <- X_f[1] + K[1]*(world$X_obs[1] - X_f[1])  
  X_en_a[,1] <- X_en_f[,1] + K[1]*(rep(world$X_obs[1],parameters$N)+eta_en[,1]-X_en_f[,1])
  
  for(i in (2:parameters$time)){
    X_en_f[,i]  <- universe$F[i] * X_en_a[,i-1] + sqrt(universe$Q[i]) * epsilon_en[,i]
    X_f[i] <- universe$F[i] * X_a[i-1]
    X_en_f[,i]  <- rep(X_f[i], parameters$N) + (X_en_f[,i] - rep(X_f[i], parameters$N))*parameters$kappa
    B_f[i]  <- var(X_en_f[,i])* parameters$blow
    
    if(is.na(world$X_obs[i])){
      K[i]   <- 0
      X_a[i] <- X_f[i]
      X_en_a[,i]<- X_en_f[,i]
    }else{
      K[i]   <- B_f[i] / (B_f[i] + parameters$std_eta^2)
      X_a[i] <- X_f[i] + K[i]*(world$X_obs[i] - X_f[i])  
      X_en_a[,i] <- X_en_f[,i] + K[i]*(rep(world$X_obs[i],parameters$N) + eta_en[,i]-X_en_f[,i])
    }
  }
  B_a <- B_f
  return(data.frame(X_f, X_a, B_f, B_a))
}

parameters_henkf <- function(){
  B_clim        <- as.matrix(9)       # scaling matrix(mean over time)
  theta         <- 10                 # scale parameter for Wishart distribution
  size_for_MC   <- 500
  psi           <- 0                  # weight for averaging B_f
  w_mix         <- 1                  # weight of B_a[k-1] for mixing with climate
  return(data.frame(B_clim, theta, size_for_MC, psi, w_mix))
}

henkf <- function(world, universe, parameters, parameters_henkf){
  # Hierarchical Bayes Ensemble Filter
  
  C_f     <- c(1:parameters$time)                                  # covariance of forecasts
  x_en_f  <- matrix(ncol = parameters$time, nrow = parameters$N)   # ensemble forecast
  x_en_a  <- matrix(ncol = parameters$time, nrow = parameters$N)   # ensemble forecast
  X_a     <- c(1:parameters$time)                                  # analysis with approximation
  B_f     <- c(1:parameters$time)                                  # B_f
  B_a     <- c(1:parameters$time)                                  # B_a
  X_f     <- c(1:parameters$time)
  S_ar       <- c(1:parameters$time)
  B_tilde <- c(1:parameters$time)
  B_tilde_arr  <- c(1:parameters$time)
  ## === generation of ensemble ====
  C_f[1] <- rWISHART(1, df = (parameters_henkf$theta+2),  S = (parameters_henkf$B_clim * parameters_henkf$theta)^{-1})  ## sample form aprior for ensemble
  epsilon_en    <- matrix(rnorm(parameters$time * parameters$N, 0, 1), ncol = parameters$time, nrow = parameters$N)  
  x_en_f[,1] <- rnorm(parameters$N, 0, sqrt(1/C_f[1]))
  ## ===============================
  C_hat <- 1/parameters_henkf$B_clim
  
  B_a[1] <- parameters_henkf$B_clim
  
  for(i in (1:parameters$time)){
    if(i == 1){
      B_f[i] <- parameters_henkf$B_clim
    }else{
      B_f[i] <- (1/B_a[i-1]+1/parameters$std_eta^2)^{-1}
    }  
    
    if(i==1){
      m_f <- 0
    }else{
      m_f <- universe$F[i] * X_a[i-1]  
    }
    X_f[i] <- m_f
    
    if(is.na(world$X_obs[i])){
      X_a[i] <- m_f
      C_mode <- 1/B_f[i]      
      x_en_a[,i] <- x_en_f[,i-1]
    }else{
      m_C <- function(C_arg){
        D <- C_arg + 1/(parameters$std_eta^2)
        result <- m_f + (world$X_obs[i] - m_f) / (D * parameters$std_eta^2) 
        return(result)
      }
      
      theta_tilde <- parameters_henkf$theta + parameters$N
      theta_hat   <- theta_tilde
      S           <- sum((x_en_f[,i]-rep(m_f, parameters$N))^2) / parameters$N
      S_ar[i]     <- S
      
      B_hat     <- (parameters_henkf$theta*B_f[i] + parameters$N * S)/(parameters_henkf$theta + parameters$N)
      B_tilde   <- B_hat
      C_mode      <- B_hat^(-1)
      
      X_a[i] <- m_C(C_mode)
      
      
      if(i < parameters$time){
        
        C_sample <- rWISHART(parameters$N, df = theta_hat + 2,  S = (B_hat * theta_hat)^{-1})
        
        for(j in (1:parameters$N)){
          x_en_a[j,i] <- rnorm(1, m = m_C(C_sample[j]), sd = (C_sample[j]+1/(parameters$std_eta^2))^(-1/2))
        }
      }  
      
    }
    if(i < parameters$time){
      x_en_f[,i+1] <- universe$F[i+1] * x_en_a[,i] + sqrt(universe$Q[i+1]) * epsilon_en[,i+1]
    }
    
    C_hat <- C_mode
    B_a[i] <- 1/C_mode    
  }
  return(data.frame(X_a, X_f, B_a, B_f))
}

parameters_hbef <- function(){
  theta         <- 4                # scale parameter for Wishart distribution
  size_for_MC   <- 500
  phi           <- 20
  chi           <- 9
  A_clim        <- 6
  Q_clim        <- 3
  
  return(data.frame(theta, size_for_MC, phi, chi, A_clim, Q_clim))
}

hbef <- function(world, universe, parameters, parameters_hbef){
  ## ==========================================================
  ## ====================== HBEF ver 2 ========================
  ## ==========================================================
  # Hierarchical Bayes Ensemble Filter
  
  # Initialization of zero conditions
  
  
  C_f     <- c(1:parameters$time)                                  # covariance of forecasts
  x_en_f  <- matrix(ncol = parameters$time, nrow = parameters$N)   # ensemble forecast
  x_en_a  <- matrix(ncol = parameters$time, nrow = parameters$N)   # ensemble forecast
  X_a     <- c(1:parameters$time)                                  # analysis with approximation
  B_f     <- c(1:parameters$time)                                  # B_f (B_tilde)
  B_a     <- c(1:parameters$time)                                  # B_a
  B_tilde <- c(1:parameters$time)
  X_a     <- c(1:parameters$time)                                  # analysis
  X_f     <- c(1:parameters$time)
  Q_hat        <- c(1:parameters$time)
  Pi_hat       <- c(1:parameters$time)
  Pi_f         <- c(1:parameters$time)
  Pi           <- c(1:parameters$time)
  S       <- c(1:parameters$time)
  S_Q     <- c(1:parameters$time)
  A       <- c(1:parameters$time)
  D            <- c(1:parameters$time)
  mean_A       <- c(1:parameters$time)
  temp         <- c(1:parameters$time)
  Pi_tilde     <- c(1:parameters$time)
  Q_tilde      <- c(1:parameters$time)
  Q_f          <- c(1:parameters$time)
  
  ## ================== Parameters of the Filter ===============
  ## 
  
  
  Pi[1]         <- universe$F[1]^2*parameters_hbef$A_clim
  
  
  ## === generation of ensemble ====
  epsilon_en    <- matrix(rnorm(parameters$time * parameters$N, 0, 1), ncol = parameters$time, nrow = parameters$N)  
  x_en_f[,1] <- universe$F[1] * rnorm(parameters$N, 0, sqrt(parameters_hbef$A_clim)) ## for PI
  ## ===============================
  
  A[1] <- parameters_hbef$A_clim
  delta_B <-0
  
  for(i in (1:parameters$time)){
  
    if(i==1){
      m_f <- 0
      Pi_f[i]  <- (mean(universe$F))^2*parameters_hbef$A_clim
    }else{
      m_f <- universe$F[i] * X_a[i-1]
      Pi_f[i] <- Pi_hat[i-1]
    }
    
    X_f[i] <- m_f
    
    if(is.na(world$X_obs[i])){
      print("ALarm")
    }else{
      
      m_C <- function(C_arg){
        D_m <- C_arg + 1/(parameters$std_eta^2)
        result <- m_f + (world$X_obs[i] - m_f) / (D_m * parameters$std_eta^2) 
        return(result)
      }
      
      
      S_Q[i]     <- universe$Q[i] * mean(epsilon_en[,i]^2)
      S[i]       <- mean(x_en_f[,i]^2)
      
      v           <- world$X_obs[i] - m_f
      
      Pi_tilde[i] <- (parameters_hbef$phi * Pi_f[i] + parameters$N * S[i])/(parameters_hbef$phi + parameters$N)
      
      
      if(i==1){
        Q_hat[i] <- (parameters_hbef$chi*parameters_hbef$Q_clim + parameters$N*S_Q[i])/(parameters_hbef$chi+parameters$N)
        Q_tilde[i] <- parameters_hbef$Q_clim
        Q_f[i] <- parameters_hbef$Q_clim
      }else{
        Q_f[i] <- Q_hat[i-1]
        Q_tilde[i] <- (parameters_hbef$chi*Q_hat[i-1] + parameters$N*S_Q[i])/(parameters_hbef$chi+parameters$N)
        
        Q_samp <- rinvgamma(500, shape = (parameters_hbef$chi+parameters$N)/2 + 1,  scale = Q_tilde[i] * (parameters_hbef$chi+parameters$N)/2)
        L_o <- function(Q){
          T <- Pi_tilde[i] + Q + parameters$std_eta^2
          result <- 1/sqrt(T) * exp(-v^2/(2*T))
        }
        L_o_samp <- L_o(Q_samp)
        Q_hat[i] <- L_o_samp%*%Q_samp / sum(L_o_samp)
        
       
      }                                            
      
      B_tilde[i] <- Pi_tilde[i] + Q_tilde[i]
      
      P_samp <- rinvgamma(500, shape = parameters_hbef$theta/2 + 1,  scale = Pi_tilde[i] * parameters_hbef$theta/2)
      L_o <- function(P){
        T <- P + Q_tilde[i] + parameters$std_eta^2
        result <- 1/sqrt(T) * exp(-v^2/(2*T))
      }
      L_o_samp <- L_o(P_samp)
      Pi_hat[i] <- L_o_samp%*%P_samp / sum(L_o_samp)
      
      
      B_f[i]     <- Pi_hat[i] + Q_hat[i]
      B_a[i] <- B_f[i]
      
      C_mode      <- 1/ B_a[i]
      D[i]           <- (B_a[i])^(-1) + 1/(parameters$std_eta^2)
      
      
      
      
      #C_sample <- rWISHART(100, df = parameters_hbef$theta + 2,  S = (B_a[i] * parameters_hbef$theta)^{-1})
      #mean_A[i] <- mean((C_sample+1/(parameters$std_eta^2))^{-1})
      
      #temp[i] <- (mean(((C_sample+1/rep(parameters$std_eta^2,100))^{-1}-rep(D[i]^{-1}, 100))^2)*(X_f[i] - world$X_obs[i])^2/(parameters$std_eta^4) )/(D[i]^{-1})
      
      
      X_a[i] <- m_C(C_mode)
      
     
      
    }
    if(i < parameters$time){
      x_en_f[,i+1] <- universe$F[i+1] * rnorm(parameters$N, 0, (D[i])^{-1/2})
      Pi[i+1] <- universe$F[i+1]^2*(D[i])^{-1}
    }
  }
  return(data.frame(X_f, X_a, Pi_hat, Q_hat))
}

rms <- function(x) sqrt(sum(x^2)/length(x))

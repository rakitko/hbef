################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## RMSE.R
#
# This file contains the R script that computes the RMSEs for the state x
# for the 5 filters (KF, Var, EnKF, HEnKF, and HBEF) in their basic configurations, 
# as in the paper
#
# ``Hierarchical Bayes Ensemble Kalman Filtering''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
#
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 27 Sep 2015
## 21 Dec 2015
## 26 May 2016
## 29 May 2016  Tsy
################################################################################
library(mixAK)
library(MCMCpack)
source('functions.R')


parameters          <- create_parameters_universe_world()
universe            <- generate_universe(parameters)
world               <- generate_world(universe, parameters)


output_kf           <- filter_kf(world, universe, parameters, parameters_kf())
rms(world$X - output_kf$X_a)

param_var           <- parameters_var()
param_var$mean_B    <- mean(output_kf$B_a)
output_var          <- filter_var(world, universe, parameters, param_var)
rms(world$X - output_var$X_a)

param_henkf         <- parameters_henkf()
param_henkf$mean_B  <- mean(output_kf$B_a)
output_henkf        <- filter_henkf(world, universe, parameters, param_henkf)
rms(world$X - output_henkf$X_a)

param_enkf          <- parameters_enkf()
param_enkf$mean_A   <- mean(output_kf$A)
output_enkf         <- filter_enkf (world, universe, parameters, param_enkf)
rms(world$X - output_enkf$X_a)


# HBEF simplest

message(" ")
message("HBEF simplest")

param_hbef          <- parameters_hbef()
param_hbef$phi <- 30 #30 
param_hbef$chi <- 5 #5 
param_hbef$mean_A   <- mean(output_kf$A)
param_hbef$mean_Q   <- mean(universe$Q * parameters$distort_Q)
output_hbef_simplest  <- filter_hbef(world, universe, parameters, param_hbef) 
rms(world$X - output_hbef_simplest$X_a)


# RESULTS

RMSE <- matrix(0, nrow = 5, ncol = 1)

colnames(RMSE) <- 'RMSE'
rownames(RMSE) <- c('KF', 'Var', 'EnKF', 'HEnKF', 'HBEF_simplest')
RMSE[1,1] <- rms(world$X - output_kf$X_a)
RMSE[2,1] <- rms(world$X - output_var$X_a)
RMSE[3,1] <- rms(world$X - output_enkf$X_a)
RMSE[4,1] <- rms(world$X - output_henkf$X_a)
RMSE[5,1] <- rms(world$X - output_hbef_simplest$X_a)

# Bar Plot

RMSE_mx <- matrix(0, nrow=2, ncol=3)


RMSE_mx[1,1] <- RMSE[3,1] - RMSE[1,1] # enkf
RMSE_mx[2,1] <- 0

RMSE_mx[1,2] <- RMSE[5,1]  - RMSE[1,1] # HBEF simplest
RMSE_mx[2,2] <- 0

bars <- c('EnKF','HBEF-simplest') # ,'HBEF-MC')
barplot(RMSE_mx[,1:2], beside=FALSE, names.arg=bars[1:2], col=c("deepskyblue4","gold"), 
        main="Analysis RMSEs (relative to KF)")
        
print(RMSE) 

# Stats of differences from the KF

rms(output_var  $X_a         - output_kf$X_a)
rms(output_enkf $X_a         - output_kf$X_a)
rms(output_henkf$X_a         - output_kf$X_a)
rms(output_hbef_simplest$X_a - output_kf$X_a)

rms(world$X - output_kf$X_a)


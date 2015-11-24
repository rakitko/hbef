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
## 23 Nov 2015
################################################################################
library(mixAK)
library(MCMCpack)
source('functions.R')


parameters          <- create_parameters_universe_world()
universe            <- generate_universe(parameters)
world               <- generate_world(universe, parameters)

output_kf           <- filter_kf(world, universe, parameters, parameters_kf())

param_var           <- parameters_var()
param_var$mean_B    <- mean(output_kf$B_a)
output_var          <- filter_var(world, universe, parameters, param_var)

output_enkf         <- filter_enkf(world, universe, parameters, parameters_enkf())

param_henkf         <- parameters_henkf()
param_henkf$mean_B  <- mean(output_kf$B_a)
output_henkf        <- filter_henkf(world, universe, parameters, param_henkf)


# HBEF Full posterior M=500

message(" ")
message("HBEF Full posterior. M=500")

param_hbef          <- parameters_hbef()
param_hbef$use_L_o         <- TRUE
param_hbef$approximation   <- FALSE

message("use_L_o=", param_hbef$use_L_o)
message("approximation=", param_hbef$approximation)

param_hbef$mean_A   <- mean(output_kf$A)
param_hbef$mean_Q   <- mean(universe$Q)
output_hbef_full_500  <- filter_hbef(world, universe, parameters, param_hbef)



# HBEF Full posterior M=100

message(" ")
message("HBEF Full posterior. M=100")

param_hbef          <- parameters_hbef()
param_hbef$use_L_o         <- TRUE
param_hbef$approximation   <- FALSE
param_hbef$size_for_MC <- 100
  
message("use_L_o=", param_hbef$use_L_o)
message("approximation=", param_hbef$approximation)
message("size_for_MC=", param_hbef$size_for_MC)

param_hbef$mean_A   <- mean(output_kf$A)
param_hbef$mean_Q   <- mean(universe$Q)
output_hbef_full_100  <- filter_hbef(world, universe, parameters, param_hbef)



# HBEF Full posterior M=30

message(" ")
message("HBEF Full posterior. M=30")

param_hbef          <- parameters_hbef()
param_hbef$use_L_o         <- TRUE
param_hbef$approximation   <- FALSE
param_hbef$size_for_MC <- 30

message("use_L_o=", param_hbef$use_L_o)
message("approximation=", param_hbef$approximation)
message("size_for_MC=", param_hbef$size_for_MC)

param_hbef$mean_A   <- mean(output_kf$A)
param_hbef$mean_Q   <- mean(universe$Q)
output_hbef_full_30  <- filter_hbef(world, universe, parameters, param_hbef)



# HBEF Approximated posterior

message(" ")
message("HBEF Approximated posterior")

param_hbef          <- parameters_hbef()

param_hbef$use_L_o         <- TRUE
param_hbef$approximation   <- TRUE

message("use_L_o=", param_hbef$use_L_o)
message("approximation=", param_hbef$approximation)

param_hbef$mean_A   <- mean(output_kf$A)
param_hbef$mean_Q   <- mean(universe$Q)
output_hbef_approx  <- filter_hbef(world, universe, parameters, param_hbef)

# HBEF-simplest: No L_o

message(" ")
message("HBEF No L_o")

param_hbef          <- parameters_hbef()

param_hbef$use_L_o         <- FALSE
param_hbef$approximation   <- FALSE

message("use_L_o=", param_hbef$use_L_o)

param_hbef$mean_A   <- mean(output_kf$A)
param_hbef$mean_Q   <- mean(universe$Q)

output_hbef_simplest  <- filter_hbef(world, universe, parameters, param_hbef) 


RMSE <- matrix(NA, nrow = 9, ncol = 1)
colnames(RMSE) <- 'RMSE'
rownames(RMSE) <- c('KF', 'Var', 'EnKF', 'HEnKF', 'HBEF_full_500', 'HBEF_full_100', 'HBEF_full_30', 'HBEF_approx', 'HBEF_simplest')
RMSE[1,1] <- rms(world$X - output_kf$X_a)
RMSE[2,1] <- rms(world$X - output_var$X_a)
RMSE[3,1] <- rms(world$X - output_enkf$X_a)
RMSE[4,1] <- rms(world$X - output_henkf$X_a)
RMSE[5,1] <- rms(world$X - output_hbef_full_500$X_a)
RMSE[6,1] <- rms(world$X - output_hbef_full_100$X_a)
RMSE[7,1] <- rms(world$X - output_hbef_full_30$X_a)
RMSE[8,1] <- rms(world$X - output_hbef_approx$X_a)
RMSE[9,1] <- rms(world$X - output_hbef_simplest$X_a)

print(RMSE)  

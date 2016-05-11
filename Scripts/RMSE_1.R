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
## 6 July 2015
################################################################################
library(mixAK)
library(MCMCpack)
source('functions.R')

print("kf and hbef")

print("Generate unverse and world")

parameters          <- create_parameters_universe_world()
universe            <- generate_universe(parameters)
world               <- generate_world(universe, parameters)

print("kf")
output_kf           <- filter_kf(world, universe, parameters, parameters_kf())

print("hbef")
param_hbef          <- parameters_hbef()
param_hbef$mean_A   <- mean(output_kf$A)
param_hbef$mean_Q   <- mean(universe$Q)
output_hbef         <- filter_hbef(world, universe, parameters, param_hbef)

RMSE <- matrix(NA, nrow = 2, ncol = 1)
colnames(RMSE) <- 'RMSE'
rownames(RMSE) <- c('KF','HBEF')
RMSE[1,1] <- rms(world$X - output_kf$X_a)
RMSE[2,1] <- rms(world$X - output_hbef$X_a)

print(RMSE)  

################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## Calculate_data_for_B_evaluation.R
#
# This file contains the R script that computes data needed to estimate 
# the ``true'' background-error variances B_k for all filters
# as well as the variances of the ``truth'' V_k,
# as in the paper
# ``Hierarchical Bayes Ensemble Kalman Filtering''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
#
# This script creates a number of data files to be used in some other HBEF-package scripts,
# so it's best to run in first.
# 
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 10 July 2015
################################################################################
library(mixAK)
library(MCMCpack)
source('functions.R')

file.create("X_true")
file.create("X_f_kf")
file.create("X_f_var")
file.create("X_f_enkf")
file.create("X_f_henkf")
file.create("X_f_hbef")
file.create("B_a_kf")
file.create("B_a_var")
file.create("B_a_enkf")
file.create("B_a_henkf")
file.create("B_a_hbef")
file.create('Pi')
file.create('S_me')
file.create('S_pe')
file.create('Q_f')
file.create('Pi_f')
file.create('Q_tilde')
file.create('Pi_tilde')



parameters          <- create_parameters_universe_world()
parameters$time     <- 20000
parameters$L        <- 200
universe            <- generate_universe(parameters)



set.seed(1)
seeds_for_world   <- sample(1:10000, parameters$L)
seeds_for_filters <- sample(1:10000, parameters$L)

for(repeat_index in (1:parameters$L)){
  print(c('Repeat number:', repeat_index))
  parameters$seed_for_world   <- seeds_for_world[repeat_index]
  parameters$seed_for_filters <- seeds_for_filters[repeat_index]
  world                       <- generate_world(universe, parameters)
  output_kf                   <- filter_kf(world, universe, parameters, parameters_kf())
  param_var                   <- parameters_var()
  param_var$mean_B            <- mean(output_kf$B_a)
  output_var                  <- filter_var(world, universe, parameters, param_var)
  output_enkf                 <- filter_enkf(world, universe, parameters, parameters_enkf())
  param_henkf                 <- parameters_henkf()
  param_henkf$mean_B          <- mean(output_kf$B_a)
  output_henkf                <- filter_henkf(world, universe, parameters, param_henkf)
  param_hbef                  <- parameters_hbef()
  param_hbef$mean_A           <- mean(output_kf$A)
  param_hbef$mean_Q           <- mean(universe$Q)
  output_hbef                 <- filter_hbef(world, universe, parameters, param_hbef)

  write.table(world$X[1:parameters$time], file="X_true", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_kf$X_f[1:parameters$time], file="X_f_kf", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_var$X_f[1:parameters$time], file="X_f_var", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_enkf$X_f[1:parameters$time], file="X_f_enkf", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_henkf$X_f[1:parameters$time], file="X_f_henkf", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_hbef$X_f[1:parameters$time], file="X_f_hbef", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_kf$B_a[1:parameters$time], file="B_a_kf", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_var$B_a[1:parameters$time], file="B_a_var", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_enkf$B_a[1:parameters$time], file="B_a_enkf", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_henkf$B_a[1:parameters$time], file="B_a_henkf", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_hbef$B_a[1:parameters$time], file="B_a_hbef", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_hbef$Pi[1:parameters$time], file="Pi", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_hbef$S_me[1:parameters$time], file="S_me", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_hbef$S_pe[1:parameters$time], file="S_pe", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_hbef$Q_f[1:parameters$time], file="Q_f", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_hbef$Pi_f[1:parameters$time], file="Pi_f", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_hbef$Q_tilde[1:parameters$time], file="Q_tilde", row.names=FALSE, col.names=FALSE, append = TRUE)
  write.table(output_hbef$Pi_tilde[1:parameters$time], file="Pi_tilde", row.names=FALSE, col.names=FALSE, append = TRUE)
}

################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## RMSE_N_NoHEnKF_3HBversions.R
#
# HEnKF excluded, all 3 versions of HBEF are examined 
# (1) Full posterior
# (2) Approximated posterior.
# (3) No L_o term at all (the simplest version)
#
# This file contains the R script that computes the RMSEs for the state x
# for the 5 filters (KF, Var, EnKF,  and HBEF) -- as functions of the ensemble size N, 
# see Fig.8(top,left) in the paper
#
# ``Hierarchical Bayes Ensemble Kalman Filtering''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
# 
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 27 Sep 2015 Tsy updated
################################################################################
library(mixAK)
library(MCMCpack)
source('functions.R')


parameters      <- create_parameters_universe_world()
range           <- c(2:10) # range of N

# initialization of arrays
rms_kf    <- c(1:length(range))
rms_var   <- c(1:length(range))
rms_enkf  <- c(1:length(range))
rms_hbef  <- c(1:length(range))
rms_hbef_approx    <- c(1:length(range))
rms_hbef_simplest  <- c(1:length(range))

t1  <-  1
t2  <-  parameters$time
t2

for (i in (1:length(range))){
  print(c('Number of range point: ',i))
  parameters$N <- range[i]
  parameters          <- update_parameters(parameters)
  universe            <- generate_universe(parameters)
  world               <- generate_world(universe, parameters)
  output_kf           <- filter_kf(world, universe, parameters, parameters_kf())
  param_var           <- parameters_var()
  param_var$mean_B    <- mean(output_kf$B_a)
  output_var          <- filter_var(world, universe, parameters, param_var)
  output_enkf         <- filter_enkf(world, universe, parameters, parameters_enkf())

  # HBEF Full posterior
  
  message(" ")
  message("HBEF Full posterior")
  
  param_hbef          <- parameters_hbef()
  
  message("use_L_o=", param_hbef$use_L_o)
  message("approximation=", param_hbef$approximation)
  
  param_hbef$mean_A   <- mean(output_kf$A)
  param_hbef$mean_Q   <- mean(universe$Q)
  output_hbef_full  <- filter_hbef(world, universe, parameters, param_hbef)
  
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
  
  # HBEF No L_o
  
  message(" ")
  message("HBEF No L_o")
  
  param_hbef          <- parameters_hbef()
  
  param_hbef$use_L_o         <- FALSE
  
  message("use_L_o=", param_hbef$use_L_o)
  
  param_hbef$mean_A   <- mean(output_kf$A)
  param_hbef$mean_Q   <- mean(universe$Q)
  
  output_hbef_simplest  <- filter_hbef(world, universe, parameters, param_hbef) 
  
  
  rms_kf[i]    <- rms(world$X[t1:t2] - output_kf$X_a[t1:t2])
  rms_var[i]   <- rms(world$X[t1:t2] - output_var$X_a[t1:t2])
  rms_enkf[i]  <- rms(world$X[t1:t2] - output_enkf$X_a[t1:t2])
  rms_hbef[i]  <- rms(world$X[t1:t2] - output_hbef_full$X_a[t1:t2])
  rms_hbef_approx[i]  <-   rms(world$X[t1:t2] - output_hbef_approx$X_a[t1:t2])
  rms_hbef_simplest[i]  <- rms(world$X[t1:t2] - output_hbef_simplest$X_a[t1:t2])
}

rms_kf

rms_hbef


pdf("RMSE_N_NoHEnKF_3HBversions.pdf", width=5.51, height=5.51)
plot(range,rms_kf, type="n", ylim=c(0,0.6), ylab="RMSE minus KF-reference-RMSE",xlab="N",xaxt='n')
axis(1, at=range, labels=range)
abline(v=range,h = seq(0,4.5,0.1), col = "gray80", lwd=0.25, lty = 3)

lines(range, rms_kf-rms_kf, lwd=3, col="darkblue", lty=1)
lines(range, rms_var-rms_kf, lwd=4, col='navajowhite3', lty=1)
lines(range, rms_enkf-rms_kf, lwd=4, col="cadetblue", lty=1)
lines(range, rms_hbef_approx-rms_kf, lwd=4, col='darkmagenta', lty=1)
lines(range, rms_hbef_simplest-rms_kf, lwd=4, col='darkolivegreen3', lty=1)
lines(range, rms_hbef-rms_kf, lwd=4, col='lightcoral', lty=1)  

leg.txt<-(c("KF reference",'Var','EnKF','HBEF_approx','HBEF_simplest','HBEF_full'))
leg.col<-c("darkblue", "navajowhite3","cadetblue",'darkmagenta','darkolivegreen3',"lightcoral")
legend("topright", inset=0,leg.txt,col=leg.col, lwd=c(3,4,4,4,4,4), pt.lwd=3, cex=1.1, pt.cex=1, bg="white")
dev.off()

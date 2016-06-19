################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## RMSE_Q_distort.R
#
# Compute and plot  the RMSEs for the state x as functions of the coefficient of distortion
# of Q---for  KF, Var, EnKF, HEnKF, and three flavors of HBEF, specifically, 
#
# 1) HBEF with the non-approximated posterior and the Monte-Carlo size $M=500$,
#
# 2) HBEF with the approximated posterior (Inverse Wishart pdfs for the posterior 
# distributions of P and Q), and
#
# 3) HBEF with no feedback from observations to the covariances at all 
# (L_o=const, the posterior is defined to be the ``sub-posterior'' here).
#
# The results are presented in Fig.6 of the paper
#
# ``A Hierarchical Bayes Ensemble Kalman Filter''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
#   
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 14 July 2015
## 29 Nov 2015 M Tsy
## 29 May 2016 M Tsy
################################################################################
library(mixAK)
library(MCMCpack)
source('functions.R')

parameters      <- create_parameters_universe_world()
parameters$time <- 200000
parameters$N    <- 5
parameters$std_eta <- 1.0
parameters          <- update_parameters(parameters)

range <-c(1/16,1/8,0.25,0.5,1,2,4,8,16)

rms_kf      <- c(1:length(range))
rms_enkf    <- c(1:length(range))
rms_henkf   <- c(1:length(range))

rms_hbef_full     <- c(1:length(range))
rms_hbef_full_theta_infty     <- c(1:length(range))
rms_hbef_approx   <- c(1:length(range))
rms_hbef_simplest <- c(1:length(range))


t1  <-  1
t2  <-  parameters$time
for (i in (1:length(range))){
  print(c('Number of range point: ',i))
  parameters$distort_Q <- range[i]
  parameters          <- update_parameters(parameters)
  universe            <- generate_universe(parameters)
  world               <- generate_world(universe, parameters)
  
  output_kf           <- filter_kf(world, universe, parameters, parameters_kf())
  
  param_enkf          <- parameters_enkf()
  param_enkf$mean_A   <- mean(output_kf$A)
  output_enkf         <- filter_enkf(world, universe, parameters, param_enkf)
  
  param_henkf         <- parameters_henkf()
  param_henkf$mean_B  <- mean(output_kf$B_a)
  output_henkf        <- filter_henkf(world, universe, parameters, param_henkf)
  
  param_hbef          <- parameters_hbef()
  param_hbef$mean_A   <- mean(output_kf$A)
  param_hbef$mean_Q   <- mean(universe$Q)
  
  # FULL
  source('functions.R')
  param_hbef$approximation <- FALSE
  param_hbef$use_L_o       <- TRUE
  output_hbef_full         <- filter_hbef(world, universe, parameters, param_hbef)
  
  # FULL, theta=1000000
  
  param_hbef$approximation <- FALSE
  param_hbef$use_L_o       <- TRUE
  param_hbef$theta         <- 1000000
  output_hbef_full_theta_infty <- filter_hbef(world, universe, parameters, param_hbef)
  
  
  # APPROX
  
#  param_hbef$approximation <- TRUE 
#  param_hbef$use_L_o       <- TRUE
#  output_hbef_approx      <- filter_hbef(world, universe, parameters, param_hbef)
  
  # SIMPLEST
  
  param_hbef$approximation <- FALSE
  param_hbef$use_L_o       <- FALSE 
  output_hbef_simplest <- filter_hbef(world, universe, parameters, param_hbef)
  
  
  rms_kf[i]    <- rms(world$X[t1:t2] - output_kf$X_a[t1:t2])
  rms_enkf[i]  <- rms(world$X[t1:t2] - output_enkf$X_a[t1:t2])
  rms_henkf[i] <- rms(world$X[t1:t2] - output_henkf$X_a[t1:t2])
  
  rms_hbef_full[i]  <-     rms(world$X[t1:t2] - output_hbef_full$X_a[t1:t2])
  rms_hbef_full_theta_infty[i] <- rms(world$X[t1:t2] - output_hbef_full_theta_infty$X_a[t1:t2])
  rms_hbef_simplest[i]  <- rms(world$X[t1:t2] - output_hbef_simplest$X_a[t1:t2])
}
rms_hbef_full
rms_hbef_full_theta_infty
rms_hbef_simplest



pdf("RMSE_Q_distort.pdf", width=5.51, height=5.51)
plot(range,rms_kf, type="n", log="x",ylim=c(0.75,1.5), ylab="RMSE",xlab="Distortion coefficient of Q ",xaxt='n')
axis(1, at=range, labels=c("1/16","1/8","1/4","1/2","1","2","4","8","16"))
abline(v=range,h = seq(0,4.5,0.2), col = "gray80", lwd=0.25, lty = 3)

lines(range, rms_kf,    lwd=3, col="darkblue",       lty=1)
lines(range, rms_enkf,  lwd=4, col="cadetblue",      lty=1)
lines(range, rms_henkf, col='darkolivegreen3', lwd=4,lty=1)

lines(range, rms_hbef_simplest, col='lightcoral', lwd=3, lty=1)
points(range, rms_hbef_full_theta_infty, col='lightcoral', lwd=4, pch=3)
points(range, rms_hbef_full,    col='lightcoral', lwd=4, pch=8)

leg.txt<-(c("KF","EnKF", "HEnKF", "HBEF simplest","HBEF MCarlo, theta>>1",  "HBEF MCarlo"))
leg.col<-c("darkblue","cadetblue", "darkolivegreen3","lightcoral","lightcoral", "lightcoral")
legend("topright", inset=0,leg.txt,col=leg.col, lwd=c(3,3,3,3,NA,NA), pch=c(NA,NA,NA,NA,3,8), pt.lwd=3, cex=1.0, pt.cex=1, bg="white")
dev.off()



pdf("RMSE_Q_distort_bw.pdf", width=5.51, height=5.51)
plot(range,rms_kf, type="n", log="x",ylim=c(0.75,1.5), ylab="RMSE",xlab="Distortion coefficient of Q ",xaxt='n')
axis(1, at=range, labels=c("1/16","1/8","1/4","1/2","1","2","4","8","16"))
abline(v=range,h = seq(0,4.5,0.2), col = "gray80", lwd=0.25, lty = 3)

lines(range, rms_kf,   lwd=3, lty=1)
lines(range, rms_enkf, lwd=3, lty=5)
lines(range, rms_henkf,lwd=3, lty=6)
lines(range, rms_hbef_simplest, lwd=4, lty=2)
points(range, rms_hbef_full_theta_infty, lwd=4, pch=3)
points(range, rms_hbef_full,    lwd=4, pch=8)

leg.txt<-(c("KF","EnKF", "HEnKF", "HBEF simplest","HBEF MCarlo, theta>>1",  "HBEF MCarlo"))
legend("topright", inset=0,leg.txt, lwd=c(3,3,3,4,NA,NA), lty=c(1,5,6,2,NA,NA), pch=c(NA,NA,NA,NA,3,8), pt.lwd=3, cex=1.0, pt.cex=1, bg="white")
dev.off()


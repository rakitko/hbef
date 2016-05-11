################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## RMSE_Q_distort_NoHEnKFnoApproxHBEF.R
#
# Compute and plot  the RMSEs for the state x as functions of the coefficient of distortion
# of Q---for  KF, Var, EnKF, and two flavors of HBEF, specifically, 
#
# 1) the full version of HBEF: with the non-approximated posterior and the Monte-Carlo size $M=500$,
#
# 2) the simplest version of HBEF: with no feedback from observations to the covariances at all 
# (L_o=const, the posterior is defined to be the ``sub-posterior'' here).
#
# The results are presented in Fig.9 of the paper
#
# ``Hierarchical Bayes Ensemble Kalman Filtering''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
#   
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 21 Oct 2015
################################################################################
library(mixAK)
library(MCMCpack)
source('functions.R')

parameters      <- create_parameters_universe_world()
parameters$time <- 200000
parameters$N    <- 3
parameters$std_eta <- 1.5
parameters          <- update_parameters(parameters)

range <-c(1/16,1/8,0.25,0.5,1,2,4,8,16)

rms_kf      <- c(1:length(range))
rms_var     <- c(1:length(range))
rms_enkf    <- c(1:length(range))
rms_hbef    <- c(1:length(range))
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
  
  param_var           <- parameters_var()
  param_var$mean_B    <- mean(output_kf$B_a)
  output_var          <- filter_var(world, universe, parameters, param_var)
  
  output_enkf         <- filter_enkf(world, universe, parameters, parameters_enkf())

  param_hbef          <- parameters_hbef()
  param_hbef$mean_A   <- mean(output_kf$A)
  param_hbef$mean_Q   <- mean(universe$Q)
  output_hbef         <- filter_hbef(world, universe, parameters, param_hbef)
  
  param_hbef$use_L_o  <- FALSE 
  output_hbef_simplest <- filter_hbef(world, universe, parameters, param_hbef)
  
  
  rms_kf[i]    <- rms(world$X - output_kf$X_a[t1:t2])
  rms_var[i]   <- rms(world$X - output_var$X_a[t1:t2])
  rms_enkf[i]  <- rms(world$X - output_enkf$X_a[t1:t2])
  rms_hbef[i]  <- rms(world$X - output_hbef$X_a[t1:t2])
  rms_hbef_simplest[i]  <- rms(world$X - output_hbef_simplest$X_a[t1:t2])
}


pdf("RMSE_Q_distort_light.pdf", width=5.51, height=5.51)
plot(range,rms_kf, type="n", log="x",ylim=c(1.05,2.0), ylab="RMSE",xlab="Distortion coefficient of Q ",xaxt='n')
axis(1, at=range, labels=c("1/16","1/8","1/4","1/2","1","2","4","8","16"))
abline(v=range,h = seq(0,4.5,0.2), col = "gray80", lwd=0.25, lty = 3)

lines(range, rms_kf,    lwd=3, col="darkblue",       lty=1)
lines(range, rms_var,   lwd=4, col='navajowhite3',   lty=1)
lines(range, rms_enkf,  lwd=4, col="cadetblue",      lty=1)
lines(range, rms_hbef,  lwd=4, col='lightcoral',     lty=1)  
points(range, rms_hbef_simplest, col='lightcoral', pch=8, lwd=3)

leg.txt<-(c("KF",'Var','EnKF', 'HBEF full', 'HBEF simplest'))
leg.col<-c("darkblue", "navajowhite3","cadetblue", "lightcoral", "lightcoral")
legend("topright", inset=0,leg.txt,col=leg.col, lwd=c(3,3,3,3,NA), pch=c(NA,NA,NA,NA,8), pt.lwd=3, cex=1.0, pt.cex=1, bg="white")
dev.off()

pdf("RMSE_Q_distort_light_bw.pdf", width=5.51, height=5.51)
plot(range,rms_kf, type="n", log="x",ylim=c(1.05,2.0), ylab="RMSE",xlab="Distortion coefficient of Q ",xaxt='n')
axis(1, at=range, labels=c("1/16","1/8","1/4","1/2","1","2","4","8","16"))
abline(v=range,h = seq(0,4.5,0.2), col = "gray80", lwd=0.25, lty = 3)

lines(range, rms_kf,   lwd=3, lty=1)
lines(range, rms_var,  lwd=2, lty=4)
lines(range, rms_enkf, lwd=3, lty=5)
lines(range, rms_hbef, lwd=4, lty=2)  
points(range, rms_hbef_simplest, pch=8, lwd=2)

leg.txt<-(c("KF",'Var','EnKF', 'HBEF full', 'HBEF simplest'))
legend("topright", inset=0,leg.txt, lwd=c(3,2,3,4,NA), lty=c(1,4,5,2,NA), pch=c(NA,NA,NA,NA,8), pt.lwd=3, cex=1.0, pt.cex=1, bg="white")

dev.off()


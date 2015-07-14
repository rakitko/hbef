################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## RMSE_N.R
#
# This file contains the R script that computes the RMSEs for the state x
# for the 5 filters (KF, Var, EnKF, HEnKF, and HBEF) -- as functions of the ensemble size N, 
# see Fig.8(top,left) in the paper
#
# ``Hierarchical Bayes Ensemble Kalman Filtering''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
# 
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 13 July 2015
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
rms_henkf <- c(1:length(range))
rms_hbef  <- c(1:length(range))

t1  <-  1
t2  <-  parameters$time
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
  param_henkf         <- parameters_henkf()
  param_henkf$mean_B  <- mean(output_kf$B_a)
  output_henkf        <- filter_henkf(world, universe, parameters, param_henkf)
  param_hbef          <- parameters_hbef()
  param_hbef$mean_A   <- mean(output_kf$A)
  param_hbef$mean_Q   <- mean(universe$Q)
  output_hbef  <- filter_hbef(world, universe, parameters, param_hbef)
  
  rms_kf[i]    <- rms(world$X[t1:t2] - output_kf$X_a[t1:t2])
  rms_var[i]   <- rms(world$X[t1:t2] - output_var$X_a[t1:t2])
  rms_enkf[i]  <- rms(world$X[t1:t2] - output_enkf$X_a[t1:t2])
  rms_henkf[i] <- rms(world$X[t1:t2] - output_henkf$X_a[t1:t2])
  rms_hbef[i]  <- rms(world$X[t1:t2] - output_hbef$X_a[t1:t2])
}



pdf("RMSE_N.pdf", width=5.51, height=5.51)
plot(range,rms_kf, type="n", ylim=c(0,0.6), ylab="RMSE minus KF-reference-RMSE",xlab="N",xaxt='n')
axis(1, at=range, labels=range)
abline(v=range,h = seq(0,4.5,0.1), col = "gray80", lwd=0.25, lty = 3)

lines(range, rms_kf-rms_kf, lwd=3, col="darkblue", lty=1)
lines(range, rms_var-rms_kf, lwd=4, col='navajowhite3', lty=1)
lines(range, rms_enkf-rms_kf, lwd=4, col="cadetblue", lty=1)
lines(range, rms_henkf-rms_kf, lwd=4, col='darkolivegreen3', lty=1)
lines(range, rms_hbef-rms_kf, lwd=4, col='lightcoral', lty=1)  

leg.txt<-(c("KF reference",'Var','EnKF', 'HEnKF','HBEF'))
leg.col<-c("darkblue", "navajowhite3","cadetblue", "darkolivegreen3","lightcoral")
legend("topright", inset=0,leg.txt,col=leg.col, lwd=c(3,3,3,3,3), pt.lwd=3, cex=1.1, pt.cex=1, bg="white")
dev.off()


pdf("RMSE_N_bw.pdf", width=5.51, height=5.51)
plot(range,rms_kf, type="n", ylim=c(0,0.6), ylab="RMSE minus KF-reference-RMSE",xlab="N",xaxt='n')
axis(1, at=range)
abline(v=range,h = seq(0,4.5,0.1), col = "gray80", lwd=0.25, lty = 3)

lines(range, rms_kf-rms_kf,   lwd=3, lty=1)
lines(range, rms_var-rms_kf,  lwd=2, lty=4)
lines(range, rms_enkf-rms_kf, lwd=3, lty=5)
lines(range, rms_henkf-rms_kf,lwd=3, lty=6)
lines(range, rms_hbef-rms_kf, lwd=4, lty=2) 

leg.txt<-(c("KF reference",'Var','EnKF', 'HEnKF','HBEF'))
legend("topright", inset=0,leg.txt, lty=c(1,4,5,6,2), lwd=c(3,2,3,3,4), pt.lwd=c(NA,NA,NA,3,NA), cex=1.0, pt.cex=1, bg="white")
dev.off()

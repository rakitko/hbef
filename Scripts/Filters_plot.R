################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## Filters_plot.R
#
# This file contains the R script that computes the time series plot Fig.7 in the paper
#
# ``Hierarchical Bayes Ensemble Kalman Filtering''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
#
# 
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 13 July 2015
################################################################################
library(mixAK)
library(MCMCpack)
source('functions.R')


parameters      <- create_parameters_universe_world()
universe        <- generate_universe(parameters)
parameters$time <- 4000
world           <- generate_world(universe, parameters)
output_kf       <- filter_kf(world, universe, parameters, parameters_kf())
output_var      <- filter_var(world, universe, parameters, parameters_var())
output_enkf     <- filter_enkf(world, universe, parameters, parameters_enkf())
output_henkf    <- filter_henkf(world, universe, parameters, parameters_henkf())
output_hbef     <- filter_hbef(world, universe, parameters, parameters_hbef())


t1 <- 3000
t2 <- 3200
pdf("Filters.pdf", width=7.48, height=5.51)

par(mai=c(1.2,1.2,0.7,0.7))

plot(1, xlim=c(t1,t2), ylim=c(-5,5), type="n", xlab='Time', ylab='X')
#abline(v=seq(0, 5000, 50),h = seq(-15,15,5), col = "gray80", lwd=0.25, lty = 2)
points(c(t1:t2),world$X_obs[c(t1:t2)], col="blue", cex=0.8)
lines(c(t1:t2),world$X[c(t1:t2)]) 
lines(c(t1:t2),output_kf$X_a[c(t1:t2)],   col="darkblue",lwd=2)
lines(c(t1:t2),output_enkf$X_a[c(t1:t2)], col="cadetblue",lwd=2)
lines(c(t1:t2),output_hbef$X_a[c(t1:t2)], col="lightcoral",lwd=3)
leg.txt<-(c("Truth",'KF','EnKF','HBEF'))
leg.col<-c("black", "darkblue","cadetblue", "lightcoral")
legend("bottom", inset=0, horiz = TRUE,leg.txt,col=leg.col, lwd=c(1,2,2,2), cex=0.9, pt.cex=1, bg="white")
dev.off()

pdf("Filters_bw.pdf", width=7.48, height=5.51)

par(mai=c(1.2,1.2,0.7,0.7))

plot(1,xlim=c(t1,t2), ylim=c(-5,5), type="n", xlab='Time', ylab='X', lwd=1)
#abline(v=seq(0, 5000, 50),h = seq(-15,15,5), col = "gray80", lwd=0.25, lty = 2)
points(c(t1:t2),world$X_obs[c(t1:t2)], cex=0.8)
lines(c(t1:t2),world$X[c(t1:t2)],         lty=3,lwd=2) 
lines(c(t1:t2),output_kf$X_a[c(t1:t2)],   lty=1,lwd=1)
lines(c(t1:t2),output_enkf$X_a[c(t1:t2)], lty=5,lwd=1.5)
lines(c(t1:t2),output_hbef$X_a[c(t1:t2)], lty=2,lwd=2)
leg.txt<-(c("Truth",'KF','EnKF','HBEF'))
legend("bottom",  inset=0, horiz = TRUE,leg.txt, lty=c(3,1,5,2), lwd=c(2,1,1.5,2), cex=0.9, pt.cex=1, bg="white")
dev.off()

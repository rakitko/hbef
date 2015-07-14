################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## QfSme_timser.R
#
# This file contains the R script that plots a segment of the time series for:
# Q^f-Q and S^{me} - Q,
# see Fig.5(left) in the paper
#
# ``Hierarchical Bayes Ensemble Kalman Filtering''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
#
# 
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 12 July 2015
################################################################################
library(mixAK)
library(MCMCpack)
source('functions.R')

parameters          <- create_parameters_universe_world()
parameters$time     <- 5000
universe            <- generate_universe(parameters)
world               <- generate_world(universe, parameters)
output_hbef         <- filter_hbef(world, universe, parameters, parameters_hbef())

t1 <- 3000
t2 <- 3200

pdf("QfSme_timser.pdf", width=7.48, height=5.51)
plot(c(t1:t2), output_hbef$Q_a[(t1-1):(t2-1)] - universe$Q[t1:t2], type="l",ylim=c(-1.5,1.5), xlab="Time", ylab="Difference from Q", col="navajowhite3", lwd=1.0)
abline(h = seq(-1.5,1.5,0.5), v = seq(t1,t2,20), col = "gray80", lwd=0.25, lty = 2)

lines(c(t1:t2),output_hbef$S_me[t1:t2] - universe$Q[t1:t2], col="cadetblue", lwd=1.0)
lines(c(t1:t2),universe$Q[t1:t2] - universe$Q[t1:t2], col="black", lwd=1.5)
leg.txt<-(c(expression(paste('Q'^'f')),
            expression(paste('S'^'me')),
            expression(paste('Q'))))

leg.col<-c("navajowhite3",
           "cadetblue", 
           "black")
legend("topright", inset=0,leg.txt, col=leg.col, lty=1, lwd=c(1,1,1.5), cex=1.1, pt.cex=1, bg="white")
dev.off()


pdf("QfSme_timser_bw.pdf", width=7.48, height=5.51)
plot(c(t1:t2), output_hbef$Q_a[(t1-1):(t2-1)] - output_hbef$Q_a[t1:t2], ylim=c(-1.5,1.5), type="l", xlab="Time", ylab="Difference from Q", lty=2, lwd=1.0)
abline(h = seq(-1.5,1.5,0.5), v = seq(t1,t2,20), col = "gray80", lwd=0.25, lty = 2)
lines(c(t1:t2),output_hbef$S_me[t1:t2] - universe$Q[t1:t2], lty=1, lwd=1.0)
lines(c(t1:t2),universe$Q[t1:t2] - universe$Q[t1:t2], col="black", lwd=2.0)
leg.txt<-(c(expression(paste('Q'^'f')),
            expression(paste('S'^'me')),
            expression(paste('Q'))))

legend("topright", inset=0,leg.txt, lty=c(1,2,1), lwd=c(1,1,2.0), cex=1.1, pt.cex=1, bg="white")
dev.off()
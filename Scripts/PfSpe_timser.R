################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## PfSpe_timser.R
#
# This file contains the R script that plots a segment of the time series for:
# Pi^f-P and S^{pe} - P,
# see Fig.6(right) in the paper
#
# ``Hierarchical Bayes Ensemble Kalman Filtering''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
#
# Dependencies:
#
# Before running this script, you need to execute the script
# Calculate_data_for_B_evaluation.R,
# which writes the time series of X_true and X_f_hbef needed here.
# 
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 11 July 2015
################################################################################
library(mixAK)
library(MCMCpack)
source('functions.R')

parameters          <- create_parameters_universe_world()
parameters$time     <- 20000
universe            <- generate_universe(parameters)
world               <- generate_world(universe, parameters)
output_kf           <- filter_kf(world, universe, parameters, parameters_kf())
param_hbef          <- parameters_hbef()
param_hbef$mean_A   <- mean(output_kf$A)
param_hbef$mean_Q   <- mean(universe$Q)
output_hbef         <- filter_hbef(world, universe, parameters, parameters_hbef())

parameters$L        <- 200
L                   <- parameters$L


X_true_mat          <- matrix(ncol = parameters$time, nrow = L)
X_f_hbef_mat        <- matrix(ncol = parameters$time, nrow = L)

B_hbef              <- c(1:parameters$time)  

X_true              <- read.table("X_true")
X_f_hbef            <- read.table("X_f_hbef")


for(i in (1:L)){
  X_true_mat[i,]    <- X_true[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  X_f_hbef_mat[i,]  <- X_f_hbef[((1+parameters$time*(i-1)):(parameters$time*i)),1]
}

for(i in (1:parameters$time)){
  B_hbef[i]         <- var(X_f_hbef_mat[,i]-X_true_mat[,i])
}

P                   <- B_hbef - universe$Q

t1 <- 3000
t2 <- 3200

pdf("PfSpe_timser.pdf", width=7.48, height=5.51)
plot(c(t1:t2), output_hbef$Pi_f[t1:t2] - P[t1:t2], type="l", xlab="Time", ylab="Difference from P", col="cadetblue4", lwd=2.0)
abline(h = seq(-10,10,5), v = seq(t1,t2,20), col = "gray80", lwd=0.25, lty = 2)
lines(c(t1:t2),output_hbef$S_pe[t1:t2] - P[t1:t2], col="lightgoldenrod3", lwd=2.0)
lines(c(t1:t2),P[t1:t2] - P[t1:t2], col="black", lwd=1.5)
leg.txt<-(c(expression(paste('P'^'f')),
            expression(paste('S'^'pe')),
            expression(paste('P'))))

leg.col<-c("cadetblue4",
           "lightgoldenrod3", 
           "black")
legend("topright", inset=0,leg.txt, col=leg.col, lty=1, lwd=c(2,2,1.5), cex=1.1, pt.cex=1, bg="white")
dev.off()

pdf("PfSpe_timser_bw.pdf", width=7.48, height=5.51)
plot(c(t1:t2), output_hbef$Pi_f[t1:t2] - P[t1:t2], type="l", xlab="Time", ylab="Difference from P", lwd=1.5)
abline(h = seq(-10,10,5), v = seq(t1,t2,20), col = "gray80", lwd=0.25, lty = 2)
lines(c(t1:t2),output_hbef$S_pe[t1:t2] - P[t1:t2], lty=2, lwd=1.5)
lines(c(t1:t2),P[t1:t2] - P[t1:t2], lwd=1.0)
leg.txt<-(c(expression(paste('P'^'f')),
            expression(paste('S'^'pe')),
            expression(paste('P'))))

legend("topright", inset=0,leg.txt, lty=c(1,2,1), lwd=c(1.5,1.5,1.0), cex=1.1, pt.cex=1, bg="white")
dev.off()

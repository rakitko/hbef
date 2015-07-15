################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## PfSpe_RMSbias.R
#
# This file contains the R script that plots the dependencies of
# RMS(bias(Pi^f-P))  and RMS(bias(S^{pe} - P)), as well as
# RMS(bias(Pi^f-Pi)) and RMS(bias(S^{pe} - Pi))
# on the number of independent runs L,
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
# which writes the time series of Pi_f and S_pe needed here.
#
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 12 July 2015
################################################################################
library(mixAK)
library(MCMCpack)
source('functions.R')

parameters          <- create_parameters_universe_world()
parameters$time     <- 20000
universe            <- generate_universe(parameters)
parameters$L        <- 200
L                   <- parameters$L

Pi_f_mat         <- matrix(ncol = parameters$time, nrow = L)
S_pe_mat         <- matrix(ncol = parameters$time, nrow = L)
Pi_mat           <- matrix(ncol = parameters$time, nrow = L)
X_true_mat          <- matrix(ncol = parameters$time, nrow = L)
X_f_hbef_mat        <- matrix(ncol = parameters$time, nrow = L)
B_hbef              <- c(1:parameters$time)  

X_true              <- read.table("X_true")
X_f_hbef            <- read.table("X_f_hbef")
Pi_f                <- read.table('Pi_f')
S_pe                <- read.table('S_pe')
Pi                  <- read.table('Pi')


for(i in (1:L)){
  X_true_mat[i,]    <- X_true[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  X_f_hbef_mat[i,]  <- X_f_hbef[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  Pi_mat[i,]        <- Pi[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  Pi_f_mat[i,]      <- Pi_f[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  S_pe_mat[i,]      <- S_pe[((1+parameters$time*(i-1)):(parameters$time*i)),1]
}

for(i in (1:parameters$time)){
  B_hbef[i]         <- var(X_f_hbef_mat[,i]-X_true_mat[,i])
}

P                   <- B_hbef - universe$Q

L <- 200
mean_Pi_f         <- c(1:parameters$time)
mean_S_pe         <- c(1:parameters$time)
mean_exact_Pi_f   <- c(1:parameters$time)
mean_exact_S_pe   <- c(1:parameters$time)


L_range <- c(1,10,25,50,75,100,125,150,175,200)
arr_f <- c(1:length(L_range))
arr_s <- c(1:length(L_range))
arr_f_true <- c(1:length(L_range))
arr_s_true <- c(1:length(L_range))

t1 <- 1
t2 <- 2000

for(j in (1:length(L_range))){
  for( i in (2:parameters$time)){
    mean_Pi_f[i] <- mean(Pi_f_mat[(1:L_range[j]),i])
    mean_S_pe[i] <- mean(S_pe_mat[(1:L_range[j]),i])
    mean_exact_Pi_f[i] <- mean(Pi_f_mat[(1:L_range[j]),i] - Pi_mat[(1:L_range[j]),i])
    mean_exact_S_pe[i] <- mean(S_pe_mat[(1:L_range[j]),i] - Pi_mat[(1:L_range[j]),i])
  }
  arr_f_true[j] <- rms(mean_Pi_f[t1:t2]-P[t1:t2])
  arr_s_true[j] <- rms(mean_S_pe[t1:t2]-P[t1:t2])
  arr_f[j] <- rms(mean_exact_Pi_f[t1:t2])
  arr_s[j] <- rms(mean_exact_S_pe[t1:t2])
}

pdf("PfSpe_RMSbias.pdf", width=5.51, height=5.51)
plot(L_range,arr_f, type="n",ylim=c(0,5), ylab="RMS(bias)",xlab="Number of realizations (L)",xaxt='n')
axis(1, at=L_range)

abline(v=c(1,25,50,75,100,125,150,175,200),h = seq(0,5,0.5), col = "gray80", lwd=0.25, lty = 3)

lines(L_range, arr_f, lwd=3, col="navajowhite3", lty=1)
lines(L_range, arr_s, lwd=4, col='cadetblue', lty=1)
lines(L_range, arr_f_true, lwd=3, col="lightcoral", lty=1)
lines(L_range, arr_s_true, lwd=4, col='darkolivegreen3', lty=1)


leg.txt <- c(expression(paste(Pi^'f',' vs. ',P)), expression(paste('S'^'pe',' vs. ',P)), expression(paste(Pi^'f',' vs. ',Pi)), expression(paste('S'^'pe',' vs. ',Pi)))

leg.col<-c( "lightcoral", "darkolivegreen3", "navajowhite3","cadetblue")
legend("topright", inset=0,leg.txt,col=leg.col, lwd=c(3,3,3,3), pt.lwd=3, cex=1.1, pt.cex=1, bg="white")
dev.off()



pdf("PfSpe_RMSbias_bw.pdf", width=5.51, height=5.51)
plot(L_range,arr_f, type="n",ylim=c(0,5), ylab="RMS(bias)",xlab="Number of realizations (L)",xaxt='n')
axis(1, at=L_range)

abline(v=c(1,25,50,75,100,125,150,175,200),h = seq(0,5,0.5), col = "gray80", lwd=0.25, lty = 3)

lines(L_range, arr_f, lwd=3, lty=1)
lines(L_range, arr_s, lwd=3, lty=2)
lines(L_range, arr_f_true, lwd=1.5,  lty=1)
lines(L_range, arr_s_true, lwd=3,  lty=4)


leg.txt <- c(expression(paste(Pi^'f',' vs. ',P)), expression(paste('S'^'pe',' vs. ',P)), expression(paste(Pi^'f',' vs. ',Pi)), expression(paste('S'^'pe',' vs. ',Pi)))

legend("topright", inset=0,leg.txt, lty=c(1,4,1,2), lwd=c(1.5,3,3,3), pt.lwd=3, cex=1.1, pt.cex=1, bg="white")
dev.off()

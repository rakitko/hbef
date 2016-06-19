################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## ScndFltErrStats.R
#
# This file contains the R script that computes and plots time series 
# of the biases (along with their bootstrap confidence intervals) 
in the sample variances $S_k$ of the forecast-error ensemble 
# for the EnKF and the HBEF.
# The errors are computed w.r.t.the estimated true $B_k$.
#
# See Figs.3 and 5 in the paper
#
# ``A Hierarchical Bayes Ensemble Kalman Filter''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
#
# Dependencies:
#
# Before running this script, you need to execute the script
# Calculate_data_for_B_evaluation.R,
# which writes the time series of X_true, X_f_enkf, X_f_hbef as well as
# S_me, S_pe, B_a_enkf, B_a_hbef.
# needed here.
# 
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 11 July 2015
## 20 Dec 2015
## 31 May 2016 Tsy
################################################################################
library(mixAK)
library(MCMCpack)
source('functions.R')

parameters          <- create_parameters_universe_world()
parameters$time     <- 10000
parameters$L        <- 500
L                   <- parameters$L
ntime               <- parameters$time
N                   <- 5 # ensemble size


X_true_mat          <- matrix(ncol = parameters$time, nrow = L)
X_f_hbef_mat        <- matrix(ncol = parameters$time, nrow = L)
X_f_enkf_mat        <- matrix(ncol = parameters$time, nrow = L)

S_me_mat            <- matrix(ncol = parameters$time, nrow = L)
S_pe_mat            <- matrix(ncol = parameters$time, nrow = L)

S_enkf_mat          <- matrix(ncol = parameters$time, nrow = L)
S_hbef_mat          <- matrix(ncol = parameters$time, nrow = L)
B_a_hbef_mat        <- matrix(ncol = parameters$time, nrow = L)


# simulated:

S_enkf_mat_sim      <- matrix(ncol = parameters$time, nrow = L)
S_hbef_mat_sim      <- matrix(ncol = parameters$time, nrow = L)

#==


B_true_enkf         <- c(1:parameters$time)
B_true_hbef         <- c(1:parameters$time)

B_a_hbef_RMSE       <- c(1:parameters$time)
B_a_hbef_bias       <- c(1:parameters$time)

S_hbef_mean         <- c(1:parameters$time)
S_hbef_var          <- c(1:parameters$time)
S_hbef_var_sim      <- c(1:parameters$time)
S_hbef_bias         <- c(1:parameters$time)
S_hbef_RMSE         <- c(1:parameters$time)

S_enkf_mean         <- c(1:parameters$time)
S_enkf_var          <- c(1:parameters$time)
S_enkf_var_sim      <- c(1:parameters$time)
S_enkf_bias         <- c(1:parameters$time)
S_enkf_RMSE         <- c(1:parameters$time)


# READ from the files

X_true              <- read.table("X_true")
X_f_enkf            <- read.table("X_f_enkf")
X_f_hbef            <- read.table("X_f_hbef")

S_me                <- read.table("S_me")   
S_pe                <- read.table("S_pe")   

B_a_enkf            <- read.table("B_a_enkf")
B_a_hbef            <- read.table("B_a_hbef")


# Fill in the (i_world,time) arrays

for(i_world in (1:L)){
  X_true_mat  [i_world,]  <- X_true  [((1+parameters$time*(i_world-1)):(parameters$time*i_world)),1]
  X_f_enkf_mat[i_world,]  <- X_f_enkf[((1+parameters$time*(i_world-1)):(parameters$time*i_world)),1]
  X_f_hbef_mat[i_world,]  <- X_f_hbef[((1+parameters$time*(i_world-1)):(parameters$time*i_world)),1]
  
  S_me_mat    [i_world,]  <- S_me    [((1+parameters$time*(i_world-1)):(parameters$time*i_world)),1]
  S_pe_mat    [i_world,]  <- S_pe    [((1+parameters$time*(i_world-1)):(parameters$time*i_world)),1]
  
  S_enkf_mat[i_world,]    <- B_a_enkf[((1+parameters$time*(i_world-1)):(parameters$time*i_world)),1]
  B_a_hbef_mat[i_world,]  <- B_a_hbef[((1+parameters$time*(i_world-1)):(parameters$time*i_world)),1]
  S_hbef_mat[i_world,]    <- S_me_mat[i_world,] + S_pe_mat[i_world,]
}

# Estm B_true

for(i in (1:ntime)){
  B_true_enkf[i] <- var(X_f_enkf_mat[,i]-X_true_mat[,i])
  B_true_hbef[i] <- var(X_f_hbef_mat[,i]-X_true_mat[,i])
}

# Mean, Bias, & RMSE in S & B_a


for(i in (1:ntime)){
  S_enkf_mean[i]  <- mean(    S_enkf_mat [,i])
  S_enkf_var [i]  <- var (    S_enkf_mat [,i])
  S_enkf_bias[i] <- mean(    S_enkf_mat [,i] - B_true_enkf[i] )
  S_enkf_RMSE[i] <- sqrt( mean((S_enkf_mat[,i] - B_true_enkf[i])^2) )
  
  B_a_hbef_bias[i] <- mean(    B_a_hbef_mat [,i] - B_true_hbef[i] )
  B_a_hbef_RMSE[i] <- sqrt( mean((B_a_hbef_mat[,i] - B_true_hbef[i])^2) )
  
  S_hbef_mean[i] <- mean(    S_hbef_mat [,i])
  S_hbef_var [i] <- var (    S_hbef_mat [,i])
  S_hbef_bias[i] <- mean(    S_hbef_mat [,i] - B_true_hbef[i] )
  S_hbef_RMSE[i] <- sqrt( mean((S_hbef_mat[,i] - B_true_hbef[i])^2) )
}


# Plots

#-------------------------------------------------------------------
# BIAS in S w.r.t. B: Bootstrap confid intls

t1 <- 1 # 2000
t2 <- 10000 # 2400

nboot <- 500
conf <- 0.95

low_int_hbef <- t1:t2
upp_int_hbef <- t1:t2

low_int_enkf <- t1:t2
upp_int_enkf <- t1:t2

for(i in t1:t2){
  confid_intvl <- bootstrap(S_hbef_mat[,i] - B_true_hbef[i], nboot, conf)
  low_int_hbef[i] <- confid_intvl$ci_mean[1]
  upp_int_hbef[i] <- confid_intvl$ci_mean[2]
  
  confid_intvl <- bootstrap(S_enkf_mat[,i] - B_true_enkf[i], nboot, conf)
  low_int_enkf[i] <- confid_intvl$ci_mean[1]
  upp_int_enkf[i] <- confid_intvl$ci_mean[2]
}

# Portion of time the confidence interval for the BIAS in S_k did contain zero

time <- t2 - t1 +1
time0_enkf <- sum(low_int_enkf[t1:t2] < 0 & upp_int_enkf[t1:t2] >0)
portion0_enkf <- time0_enkf / time
portion0_enkf
time0_hbef <- sum(low_int_hbef[t1:t2] < 0 & upp_int_hbef[t1:t2] >0)
portion0_hbef <- time0_hbef / time
portion0_hbef

pdf("S_bias_EnKF_HBEF_timser.pdf",  width=7.48, height=5.51)
col <- rgb(159/255,182/255,205/255,0.4)
plot(c(t1:t2), S_enkf_bias[t1:t2], type="n", xlab="Time", 
     ylab="Error variances", ylim=c(-9,22))
polygon(c(t1:t2,t2:t1),c(low_int_enkf[t1:t2], upp_int_enkf[t2:t1]), 
        col=col, border=rgb(1,1,1))
polygon(c(t1:t2,t2:t1),c(low_int_hbef[t1:t2], upp_int_hbef[t2:t1]), 
        col=col, border=rgb(1,1,1))
lines(c(t1:t2),B_true_enkf[t1:t2], type="l", col="cadetblue",  lwd=1.5, lty=3)
lines(c(t1:t2),B_true_hbef[t1:t2], type="l", col="lightcoral", lwd=1.5)
lines(c(t1:t2),S_enkf_bias[t1:t2], type="l", col="cadetblue",  lwd=2.5, lty=3)
lines(c(t1:t2),S_hbef_bias[t1:t2], type="l", col="lightcoral", lwd=1.0)
abline(h=0)
leg.txt<-c("B_true: EnKF", "B_true: HBEF", "bias(S-B_true): EnKF", "bias(S-B_true): HBEF", "Confidence intervals")
leg.col<-c("cadetblue", "lightcoral", "cadetblue","lightcoral",col)
legend("topleft", inset=0,leg.txt, col=leg.col, lty=c(3,1,3,1,1), lwd=c(1.5,1.5,2.5,1,6), 
       #text.width=strwidth("10,000,000"), 
       cex=0.8, pt.cex=1, bg="white")
dev.off()


# b/w


pdf("S_bias_EnKF_HBEF_timser_bw.pdf",  width=7.48, height=5.51)
plot(c(t1:t2), S_enkf_bias[t1:t2], type="n", xlab="Time", 
     ylab="Error variances", ylim=c(-9,31))
col <- rgb(100/255,100/255,100/255,0.4)
polygon(c(t1:t2,t2:t1),c(low_int_enkf[t1:t2], upp_int_enkf[t2:t1]), 
        col=col, border=rgb(1,1,1))
polygon(c(t1:t2,t2:t1),c(low_int_hbef[t1:t2], upp_int_hbef[t2:t1]), 
        col=col, border=rgb(1,1,1))
lines(c(t1:t2),B_true_enkf[t1:t2], type="l",  lwd=1.5, lty=3)
lines(c(t1:t2),B_true_hbef[t1:t2], type="l",  lwd=1.5)
lines(c(t1:t2), S_enkf_bias[t1:t2], type="l", lwd=2.5, lty=3)
lines(c(t1:t2),S_hbef_bias[t1:t2], type="l",  lwd=1.0)
abline(h=0)
leg.txt<-c("B_true: EnKF", "B_true: HBEF", "bias(S-B_true): EnKF", "bias(S-B_true): HBEF", "Confidence intervals")
legend("topleft", inset=0,leg.txt, lty=c(3,1,3,1,1), lwd=c(1.5,1.5,2.5,1,6), col=c(1,1,1,1,col),
       #text.width=strwidth("10,000,000"), 
       cex=0.6, pt.cex=1, bg="white")
dev.off()




#-------------------------------------------------------------------
# RMSEs in B for HBEF & EnKF


t1 <- 2000 # 10000  #3000
t2 <- 2400 # 11000 #3200

plot(c(t1:t2), B_a_hbef_RMSE[t1:t2], type="l", xlab="Time", 
     ylab="RMSE(B_a - B_true)", col="lightcoral", lwd=2.0, ylim=c(0,20))
lines(c(t1:t2),S_enkf_RMSE[t1:t2], type="l",col="cadetblue", lwd=2.0)

leg.txt<-c("EnKF", "HBEF")
leg.col<-c("cadetblue", "lightcoral")
legend("topright", inset=0,leg.txt, col=leg.col, lty=1, lwd=c(2,2), cex=0.8, pt.cex=1, bg="white")

dev.copy(pdf, "B_RMSE_EnKF_HBEF_timser.pdf",  width=7.48, height=5.51)
dev.off()


# b/w


plot(c(t1:t2), B_a_hbef_RMSE[t1:t2], type="l", xlab="Time", 
     ylab="RMSE(B_a - B_true)",  lwd=2.0, lty=1, ylim=c(0,20))
lines(c(t1:t2),S_enkf_RMSE[t1:t2], type="l", lwd=2.0, lty=3)

leg.txt<-c("EnKF", "HBEF")
legend("topright", inset=0,leg.txt, lty=c(1,3), lwd=c(2,2), cex=0.8, pt.cex=1, bg="white")

dev.copy(pdf, "B_RMSE_EnKF_HBEF_timser_bw.pdf",  width=7.48, height=5.51)
dev.off()


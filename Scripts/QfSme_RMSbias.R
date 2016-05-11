################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## QfSme_RMSbias.R
#
# This file contains the R script that plots the dependencies of
# RMS(bias(Q^f-Q)) and RMS(bias(S^{me} - Q)) on the number of independent runs L,
# see Fig.6(left) in the paper
#
# ``Hierarchical Bayes Ensemble Kalman Filtering''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
#
# Dependencies:
#
# Before running this script, you need to execute the script
# Calculate_data_for_B_evaluation.R,
# which writes the time series of Q_f and S_me needed here.
#  
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 13 July 2015
################################################################################
source('functions.R')
parameters          <- create_parameters_universe_world()
parameters$time     <- 20000
universe            <- generate_universe(parameters)
parameters$L        <- 200
L                   <- parameters$L

Q_f_mat         <- matrix(ncol = parameters$time, nrow = L)
S_me_mat        <- matrix(ncol = parameters$time, nrow = L)

L                   <- parameters$L


mean_Q_f <- c(1:parameters$time)
mean_S_me <- c(1:parameters$time)

Q_f             <- read.table("Q_f")
S_me            <- read.table("S_me")

for(i in (1:L)){
  Q_f_mat[i,]     <- Q_f[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  S_me_mat[i,]    <- S_me[((1+parameters$time*(i-1)):(parameters$time*i)),1]
}

L_range <- c(1,10,25,50,75,100,125,150,175,200)

arr_f_true <- c(1:length(L_range))
arr_s_true <- c(1:length(L_range))

t1 <- 1
t2 <- 20000

for(j in (1:length(L_range))){
  for( i in (2:parameters$time)){
    mean_Q_f[i] <- mean(Q_f_mat[(1:L_range[j]),i])
    mean_S_me[i] <- mean(S_me_mat[(1:L_range[j]),i])
  }
  arr_f_true[j] <- rms(mean_Q_f[t1:t2]-universe$Q[t1:t2])
  arr_s_true[j] <- rms(mean_S_me[t1:t2]-universe$Q[t1:t2])
}



pdf("QfSme_RMSbias.pdf", width=5.51, height=5.51)
plot(L_range,L_range, type="n",ylim=c(0,2), ylab="RMS(bias)",xlab="Number of realizations (L)",xaxt='n')
axis(1, at=L_range)

abline(v=c(1,25,50,75,100,125,150,175,200),h = seq(0,5,0.5), col = "gray80", lwd=0.25, lty = 3)

lines(L_range, arr_f_true, lwd=3, col="navajowhite3", lty=1)
lines(L_range, arr_s_true, lwd=3, col='cadetblue', lty=1)


leg.txt <- c(expression(paste(Q^'f',' vs. ',Q)), expression(paste('S'^'me',' vs. ',Q)))

leg.col<-c("navajowhite3","cadetblue")
legend("topright", inset=0,leg.txt,col=leg.col, lwd=c(3,3,3,3), pt.lwd=3, cex=1.1, pt.cex=1, bg="white")
dev.off()



pdf("QfSme_RMSbias_bw.pdf", width=5.51, height=5.51)
plot(L_range,L_range, type="n",ylim=c(0,2), ylab="RMS(bias)",xlab="Number of realizations (L)",xaxt='n')
axis(1, at=L_range)

abline(v=c(1,25,50,75,100,125,150,175,200),h = seq(0,5,0.5), col = "gray80", lwd=0.25, lty = 3)

lines(L_range, arr_f_true, lwd=3, lty=1)
lines(L_range, arr_s_true, lwd=3, lty=2)


leg.txt <- c(expression(paste(Q^'f',' vs. ',Q)), expression(paste('S'^'me',' vs. ',Q)))

leg.col<-c("navajowhite3","cadetblue")
legend("topright", inset=0,leg.txt, lwd=c(3,3,3,3), lty=c(1,2), pt.lwd=3, cex=1.1, pt.cex=1, bg="white")
dev.off()


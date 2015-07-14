################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## qqplot_1.R
#
# This file contains the R script that computes the q-q (quantile-quantile) plot and 
# the Gaussian (normal) approximation for the unconditional distribution of x-m^f, 
# see Fig.4(right) in the paper
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
## 10 July 2015
################################################################################
source('functions.R')
parameters          <- create_parameters_universe_world()
parameters$time     <- 20000
parameters$L        <- 200
L                   <- parameters$L

X_true_mat          <- matrix(ncol = parameters$time, nrow = L)
X_f_hbef_mat        <- matrix(ncol = parameters$time, nrow = L)

X_true              <- read.table("X_true")
X_f_hbef            <- read.table("X_f_hbef")

for(i in (1:L)){
  X_true_mat[i,]    <- X_true[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  X_f_hbef_mat[i,]  <- X_f_hbef[((1+parameters$time*(i-1)):(parameters$time*i)),1]
}

M                   <- X_true_mat - X_f_hbef_mat
vec                 <- as.vector(M)
ind                 <- sample((1:length(vec)),5000)
y                   <- rnorm(5000)
vec_1               <- vec[ind]

t                   <- extRemes::qqplot(y,vec[ind], regress=FALSE, ylim=c(-20,20),xlim=c(-4.5,4.5), make.plot = FALSE,
                         xlab="Theoretical Quantiles", ylab="Sample Quantiles",
                         main = expression(paste('X'['f'],' - ','X',' | ','40',' < ','B'['true'],' < ','45')), pch=1,lwd=1, col="tomato")


pdf("qqplot1.pdf", width=5.51, height=5.51)
par(mai=c(1.2,1.2,0.7,0.7))
plot(1,type="n",ylim=c(-20,20),xlim=c(-4.5,4.5),xlab="Theoretical Quantiles", ylab="Sample Quantiles",
     main = expression(paste('X'^'f',' - ','X')))
abline(h = 0, v = 0, col = "gray80", lwd=0.25)
abline(h = c(-20,-10,10,20), v = c(-4,-2,2,4), col = "gray80", lwd=0.25, lty = 2)
points(t$qdata[,1], t$qdata[,2], pch=1,lwd=1, col="cadetblue")
col.line <- "gold1"
qqline(vec_1, col = col.line, lwd=2, lty=2)
dev.off()  

pdf("qqplot1_bw.pdf", width=5.51, height=5.51)
par(mai=c(1.2,1.2,0.7,0.7))
plot(1,type="n",ylim=c(-20,20),xlim=c(-4.5,4.5),xlab="Theoretical Quantiles", ylab="Sample Quantiles",
     main = expression(paste('X'^'f',' - ','X')))
abline(h = 0, v = 0, col = "gray80", lwd=0.25)
abline(h = c(-20,-10,10,20), v = c(-4,-2,2,4), col = "gray80", lwd=0.25, lty = 2)
points(t$qdata[,1], t$qdata[,2], pch=1,lwd=1)
qqline(vec_1, lwd=2, lty=2)
dev.off()  

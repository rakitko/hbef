################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
##  qqplot_2.R
#
# This file contains the R script that computes the q-q (quantile-quantile) plot and 
# the Gaussian (normal) approximation for the conditional distribution x-m^f|B, 
# see Fig.4(left) in the paper
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

ind_1               <- which((B_hbef > 4) & (B_hbef < 5))
ind_2               <- which((B_hbef > 19) & (B_hbef < 21))
ind_3               <- which((B_hbef > 40) & (B_hbef < 45))


vec_1               <- as.vector(X_f_hbef_mat[,ind_1]-X_true_mat[,ind_1])
vec_2               <- as.vector(X_f_hbef_mat[,ind_2]-X_true_mat[,ind_2])
vec_3               <- as.vector(X_f_hbef_mat[,ind_3]-X_true_mat[,ind_3])

ind_1               <- sample((1:length(vec_1)),5000)
ind_2               <- sample((1:length(vec_2)),5000)
ind_3               <- (1:length(vec_3))

vec_1               <- vec_1[ind_1]
vec_2               <- vec_2[ind_2]
vec_3               <- vec_3[ind_3]


t_1                 <- extRemes::qqplot(y,vec_1, regress=FALSE, make.plot = FALSE)
t_2                 <- extRemes::qqplot(y,vec_2, regress=FALSE, make.plot = FALSE)
t_3                 <- extRemes::qqplot(y,vec_3, regress=FALSE, make.plot = FALSE)


pdf("qqplot2.pdf", width=5.51, height=5.51)
par(mai=c(1.2,1.2,0.7,0.7))
plot(1,type="n",ylim=c(-20,20),xlim=c(-4.5,4.5),xlab="Theoretical Quantiles", ylab="Sample Quantiles",
     main = expression(paste('X'^'f',' - ','X',' | ','B')))
abline(h = 0, v = 0, col = "gray80", lwd=0.25)
abline(h = c(-20,-10,10,20), v = c(-4,-2,2,4), col = "gray80", lwd=0.25, lty = 2)
points(t_1$qdata[,1], t_1$qdata[,2], pch=1,lwd=1, col="cadetblue")
points(t_2$qdata[,1], t_2$qdata[,2], pch=8,lwd=1, col="lightcoral")
points(t_3$qdata[,1], t_3$qdata[,2], pch=2,lwd=1, col="navajowhite3")
col.line <- "gold1"
qqline(vec_1, col = col.line, lwd=2, lty=2)
qqline(vec_2, col = col.line, lwd=2, lty=2)
qqline(vec_3, col = col.line, lwd=2, lty=2)

leg.txt<-(c(expression(paste('4',' < ','B',' < ','5')),
            expression(paste('19',' < ','B',' < ','21')),
            expression(paste('34',' < ','B',' < ','39'))))
leg.col<-c("cadetblue", "lightcoral", "navajowhite3")
legend("bottomright",leg.txt, col=leg.col, pch=c(1,8,2), cex=1.2, bg="white")
dev.off()  


pdf("qqplot2_bw.pdf", width=5.51, height=5.51)
par(mai=c(1.2,1.2,0.7,0.7))
plot(1,type="n",ylim=c(-20,20),xlim=c(-4.5,4.5),xlab="Theoretical Quantiles", ylab="Sample Quantiles",
     main = expression(paste('X'^'f',' - ','X',' | ',"B")))
abline(h = 0, v = 0, col = "gray80", lwd=0.25)
abline(h = c(-20,-10,10,20), v = c(-4,-2,2,4), col = "gray80", lwd=0.25, lty = 2)
points(t_1$qdata[,1], t_1$qdata[,2], pch=1,lwd=1)
points(t_2$qdata[,1], t_2$qdata[,2], pch=8,lwd=1)
points(t_3$qdata[,1], t_3$qdata[,2], pch=2,lwd=1)
col.line <- "gold1"
qqline(vec_1, lwd=2, lty=2)
qqline(vec_2,  lwd=2, lty=2)
qqline(vec_3, lwd=2, lty=2)
leg.txt<-(c(expression(paste('4',' < ','B',' < ','5')),
            expression(paste('19',' < ','B',' < ','21')),
            expression(paste('34',' < ','B',' < ','39'))))

leg.col<-c("cadetblue", "lightcoral", "navajowhite3")
legend("bottomright",leg.txt, pch=c(1,8,2), cex=1.2, bg="white")
dev.off()
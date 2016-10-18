################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## Timeseries.R
# 
# This script computes a segment of the time series of the following variables:
# F_k
# Sigma_k
# V_k
# B_k for the HBEF
#
# see Fig.1 in the paper
#
# Before running this script, you need to execute the script Evaluate_B.R,
# which writes the time series of V_k and B_k_true for all the filters needed here.
#
# ``A Hierarchical Bayes Ensemble Kalman Filter''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
#
#
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 9 July 2015
## 4 Dec 2015
################################################################################
library(mixAK)
library(MCMCpack)
source('functions.R')

B_kf <- as.vector(t(as.matrix(read.table('B_kf.txt'))))
B_var <- as.vector(t(as.matrix(read.table('B_var.txt'))))
B_enkf <- as.vector(t(as.matrix(read.table('B_enkf.txt'))))
B_henkf <- as.vector(t(as.matrix(read.table('B_henkf.txt'))))
B_hbef <- as.vector(t(as.matrix(read.table('B_hbef.txt'))))
var_of_x_true <- as.vector(t(as.matrix(read.table('var_of_x_true.txt'))))

parameters      <- create_parameters_universe_world()
parameters$time <- 10000
universe        <- generate_universe(parameters)

t1 <- 2000
t2 <- 3000

pdf("Timeseries.pdf", width=7.48, height=11)
nf<-layout(matrix(c(1,2,3,4), ncol=1), c(19),c(lcm(5.588),lcm(5.08),lcm(5.08),lcm(5.08)),TRUE)
par(mai=c(0,0.8,0.2,0.4))
plot(1,xlim=c(t1,t2), ylim=c(0.6,1.2), type="n", ylab="F", axes = FALSE)
abline(h = 1, col = "gray50", lwd=2, lty = 2)
abline(v = seq(t1,t2,200), col = "gray80", lwd=0.25, lty = 2)
abline(h = seq(0.6,0.9,0.1), col = "gray80", lwd=0.25, lty = 2)
abline(h = seq(1.1,1.2,0.1), col = "gray80", lwd=0.25, lty = 2)
polygon(c(2155,2155,2173,2173),c(-100,100,100,-100), border=NA,col=rgb(240/255,128/255,128/255,0.4))
polygon(c(2760,2760,2815,2815),c(-100,100,100,-100), border=NA,col=rgb(240/255,128/255,128/255,0.4))
polygon(c(2020,2020,2040,2040),c(-100,100,100,-100), border=NA,col=rgb(95/255,158/255,160/255,0.4))
polygon(c(2380,2380,2430,2430),c(-100,100,100,-100), border=NA,col=rgb(95/255,158/255,160/255,0.4))
polygon(c(2560,2560,2590,2590),c(-100,100,100,-100), border=NA,col=rgb(95/255,158/255,160/255,0.4))  
polygon(c(2850,2850,2880,2880),c(-100,100,100,-100), border=NA,col=rgb(240/255,128/255,128/255,0.4))
lines(c(t1:t2), universe$F[t1:t2])
axis(side=1,seq(2000,3000,200))
axis(side=2,seq(0.6,1.2,0.1))
box()
mtext("(a)",side = 4, las=1, adj = -0.5)


par(mai=c(0,0.8,0,0.4))
plot(c(t1:t2),exp(universe$Sigma[t1:t2]), axes = FALSE, type="n", ylab=expression(sigma))
abline(v = seq(t1,t2,200), col = "gray80", lwd=0.25, lty = 2)
abline(h = seq(0.5,3.5,0.5), col = "gray80", lwd=0.25, lty = 2)
polygon(c(2155,2155,2173,2173),c(-100,100,100,-100), border=NA,col=rgb(240/255,128/255,128/255,0.4))
polygon(c(2760,2760,2815,2815),c(-100,100,100,-100), border=NA,col=rgb(240/255,128/255,128/255,0.4))
polygon(c(2020,2020,2040,2040),c(-100,100,100,-100), border=NA,col=rgb(95/255,158/255,160/255,0.4))
polygon(c(2380,2380,2430,2430),c(-100,100,100,-100), border=NA,col=rgb(95/255,158/255,160/255,0.4))
polygon(c(2560,2560,2590,2590),c(-100,100,100,-100), border=NA,col=rgb(95/255,158/255,160/255,0.4))  
polygon(c(2850,2850,2880,2880),c(-100,100,100,-100), border=NA,col=rgb(240/255,128/255,128/255,0.4))
lines(c(t1:t2),exp(universe$Sigma[t1:t2]))
axis(side=1,seq(2000,3000,200))
axis(side=2,seq(0.5,3.0,0.5))
box()
mtext("(b)",side = 4, las=1, adj=-0.5)

plot(c(t1:t2),var_of_x_true[t1:t2], type="n", ylab='V_true')
abline(v = seq(t1,t2,200), col = "gray80", lwd=0.25, lty = 2)
abline(h = seq(0,150,20), col = "gray80", lwd=0.25, lty = 2)
polygon(c(2155,2155,2173,2173),c(-100,150,150,-100), border=NA,col=rgb(240/255,128/255,128/255,0.4))
polygon(c(2760,2760,2815,2815),c(-100,150,150,-100), border=NA,col=rgb(240/255,128/255,128/255,0.4))
polygon(c(2020,2020,2040,2040),c(-100,150,150,-100), border=NA,col=rgb(95/255,158/255,160/255,0.4))
polygon(c(2380,2380,2430,2430),c(-100,150,150,-100), border=NA,col=rgb(95/255,158/255,160/255,0.4))
polygon(c(2560,2560,2590,2590),c(-100,150,150,-100), border=NA,col=rgb(95/255,158/255,160/255,0.4))  
polygon(c(2850,2850,2880,2880),c(-100,150,150,-100), border=NA,col=rgb(240/255,128/255,128/255,0.4))
lines(c(t1:t2),var_of_x_true[t1:t2])
mtext("(c)",side = 4, las=1, adj = -0.5)

plot(c(t1:t2),B_hbef[t1:t2], type="n", xlab="Time", ylab="B_true")
abline(v = seq(t1,t2,200), col = "gray80", lwd=0.25, lty = 2)
abline(h = seq(0,30,5), col = "gray80", lwd=0.25, lty = 2)
polygon(c(2155,2155,2173,2173),c(-100,100,100,-100), border=NA,col=rgb(240/255,128/255,128/255,0.4))
polygon(c(2760,2760,2815,2815),c(-100,100,100,-100), border=NA,col=rgb(240/255,128/255,128/255,0.4))
polygon(c(2020,2020,2040,2040),c(-100,100,100,-100), border=NA,col=rgb(95/255,158/255,160/255,0.4))
polygon(c(2380,2380,2430,2430),c(-100,100,100,-100), border=NA,col=rgb(95/255,158/255,160/255,0.4))
polygon(c(2560,2560,2590,2590),c(-100,100,100,-100), border=NA,col=rgb(95/255,158/255,160/255,0.4))  
polygon(c(2850,2850,2880,2880),c(-100,100,100,-100), border=NA,col=rgb(240/255,128/255,128/255,0.4))
lines(c(t1:t2),B_hbef[t1:t2])
mtext("(d)",side = 4, las=1, adj = -0.5)

mtext("Time", side=1, padj=3, outer=FALSE)

dev.off()



######################################################

pdf("Timeseries_bw.pdf", width=7.48, height=11)
nf<-layout(matrix(c(1,2,3,4), ncol=1), c(19),c(lcm(5.588),lcm(5.08),lcm(5.08),lcm(5.08)),TRUE)
par(mai=c(0,0.8,0.2,0.4))
plot(1,xlim=c(t1,t2), ylim=c(0.6,1.2), type="n", ylab="F", axes = FALSE)
abline(h = 1, col = "gray50", lwd=2, lty = 2)
abline(v = seq(t1,t2,200), col = "gray80", lwd=0.25, lty = 2)
abline(h = seq(0.6,0.9,0.1), col = "gray80", lwd=0.25, lty = 2)
abline(h = seq(1.1,1.2,0.1), col = "gray80", lwd=0.25, lty = 2)
polygon(c(2155,2155,2173,2173),c(-100,150,150,-100), border=NA,col=rgb(0.9,0.9,0.9,0.5))
polygon(c(2760,2760,2815,2815),c(-100,150,150,-100), border=NA,col=rgb(0.9,0.9,0.9,0.5))
polygon(c(2020,2020,2040,2040),c(-100,150,150,-100), border=NA,col=rgb(0.7,0.7,0.7,0.5))
polygon(c(2380,2380,2430,2430),c(-100,150,150,-100), border=NA,col=rgb(0.7,0.7,0.7,0.5))
polygon(c(2560,2560,2590,2590),c(-100,150,150,-100), border=NA,col=rgb(0.7,0.7,0.7,0.5))  
polygon(c(2850,2850,2880,2880),c(-100,150,150,-100), border=NA,col=rgb(0.9,0.9,0.9,0.5))
lines(c(t1:t2), universe$F[t1:t2])
axis(side=1,seq(2000,3000,200))
axis(side=2,seq(0.6,1.2,0.1))
box()
mtext("(a)",side = 4, las=1, adj = -0.5)


par(mai=c(0,0.8,0,0.4))
plot(c(t1:t2),exp(universe$Sigma[t1:t2]), axes = FALSE, type="n", ylab=expression(sigma))
abline(v = seq(t1,t2,200), col = "gray80", lwd=0.25, lty = 2)
abline(h = seq(0.5,3.5,0.5), col = "gray80", lwd=0.25, lty = 2)
polygon(c(2155,2155,2173,2173),c(-100,150,150,-100), border=NA,col=rgb(0.9,0.9,0.9,0.5))
polygon(c(2760,2760,2815,2815),c(-100,150,150,-100), border=NA,col=rgb(0.9,0.9,0.9,0.5))
polygon(c(2020,2020,2040,2040),c(-100,150,150,-100), border=NA,col=rgb(0.7,0.7,0.7,0.5))
polygon(c(2380,2380,2430,2430),c(-100,150,150,-100), border=NA,col=rgb(0.7,0.7,0.7,0.5))
polygon(c(2560,2560,2590,2590),c(-100,150,150,-100), border=NA,col=rgb(0.7,0.7,0.7,0.5))  
polygon(c(2850,2850,2880,2880),c(-100,150,150,-100), border=NA,col=rgb(0.9,0.9,0.9,0.5))
lines(c(t1:t2),exp(universe$Sigma[t1:t2]))
axis(side=1,seq(2000,3000,200))
axis(side=2,seq(0.5,3.0,0.5))
box()
mtext("(b)",side = 4, las=1, adj=-0.5)

plot(c(t1:t2),var_of_x_true[t1:t2], type="n", ylab='V_true')
abline(v = seq(t1,t2,200), col = "gray80", lwd=0.25, lty = 2)
abline(h = seq(0,100,20), col = "gray80", lwd=0.25, lty = 2)
polygon(c(2155,2155,2173,2173),c(-100,150,150,-100), border=NA,col=rgb(0.9,0.9,0.9,0.5))
polygon(c(2760,2760,2815,2815),c(-100,150,150,-100), border=NA,col=rgb(0.9,0.9,0.9,0.5))
polygon(c(2020,2020,2040,2040),c(-100,150,150,-100), border=NA,col=rgb(0.7,0.7,0.7,0.5))
polygon(c(2380,2380,2430,2430),c(-100,150,150,-100), border=NA,col=rgb(0.7,0.7,0.7,0.5))
polygon(c(2560,2560,2590,2590),c(-100,150,150,-100), border=NA,col=rgb(0.7,0.7,0.7,0.5))  
polygon(c(2850,2850,2880,2880),c(-100,150,150,-100), border=NA,col=rgb(0.9,0.9,0.9,0.5))
lines(c(t1:t2),var_of_x_true[t1:t2])
mtext("(c)",side = 4, las=1, adj = -0.5)

plot(c(t1:t2),B_hbef[t1:t2], type="n", ylab="B_true")
abline(v = seq(t1,t2,200), col = "gray80", lwd=0.25, lty = 2)
abline(h = seq(0,30,5), col = "gray80", lwd=0.25, lty = 2)
polygon(c(2155,2155,2173,2173),c(-100,100,100,-100), border=NA,col=rgb(0.9,0.9,0.9,0.5))
polygon(c(2760,2760,2815,2815),c(-100,100,100,-100), border=NA,col=rgb(0.9,0.9,0.9,0.5))
polygon(c(2020,2020,2040,2040),c(-100,100,100,-100), border=NA,col=rgb(0.7,0.7,0.7,0.5))
polygon(c(2380,2380,2430,2430),c(-100,100,100,-100), border=NA,col=rgb(0.7,0.7,0.7,0.5))
polygon(c(2560,2560,2590,2590),c(-100,100,100,-100), border=NA,col=rgb(0.7,0.7,0.7,0.5))  
polygon(c(2850,2850,2880,2880),c(-100,100,100,-100), border=NA,col=rgb(0.9,0.9,0.9,0.5))
lines(c(t1:t2),B_hbef[t1:t2])
mtext("(d)",side = 4, las=1, adj = -0.5)

mtext("Time", side=1, padj=3, outer=FALSE)

dev.off()

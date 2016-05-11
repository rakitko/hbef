################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## 
#
# 
#
#
# ``Hierarchical Bayes Ensemble Kalman Filtering''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
#
# 
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 11 July 2015
################################################################################
library(mixAK)
library(MCMCpack)
library(pscl)
source('functions.R')

parameters          <- create_parameters_universe_world()
parameters$time     <- 20000
universe            <- generate_universe(parameters)
parameters$L        <- 200
L                   <- parameters$L


X_true_mat          <- matrix(ncol = parameters$time, nrow = L)
X_f_hbef_mat        <- matrix(ncol = parameters$time, nrow = L)
Pi_mat              <- matrix(ncol = parameters$time, nrow = L)

B_hbef              <- c(1:parameters$time)  

X_true              <- read.table("X_true")
X_f_hbef            <- read.table("X_f_hbef")
Pi                  <- read.table("Pi")


for(i in (1:L)){
  X_true_mat[i,]    <- X_true[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  X_f_hbef_mat[i,]  <- X_f_hbef[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  Pi_mat[i,]        <- Pi[((1+parameters$time*(i-1)):(parameters$time*i)),1]
}

for(i in (1:parameters$time)){
  B_hbef[i]         <- var(X_f_hbef_mat[,i]-X_true_mat[,i])
}


num_of_plots <- 3
bot_bound <- c(1.95,5.95,9.75)
up_bound <- c(2.05,6.05,10.25)

P                   <- B_hbef - universe$Q

for(j in (1:num_of_plots)){
  
  P_mod  <- c()
  Pi_mod <- c()
  
  for(l in (1:parameters$time)){
    vec <- c()
    for(k in (1:parameters$L)){
      if((Pi_mat[k,l]>bot_bound[j])&&(Pi_mat[k,l]<up_bound[j])){
        vec <- c(vec, Pi_mat[k,l])
      }
    }
    if(length(vec)>0) P_mod <- c(P_mod, rep(P[l], length(vec)))
    Pi_mod <- c(Pi_mod, vec)
  }
  
  arg     <- seq(0.1,100,by=0.1)
  d       <- density.default(P_mod, from=0)
  val     <-c(1:length(arg))
  theta_P <- 4
    
  for(i in (1:length(arg))){
    val[i]<-pscl::densigamma(arg[i],theta_P/2+1, mean(Pi_mod)*theta_P/2)
  }
  
  d$y[1]<-0
  
  if (j==1) {
    d1   <- d
    P_1  <- P_mod
    val1 <- val
  }
  if (j==2) {
    d2   <- d
    val2 <- val
    P_2  <- P_mod
  }
  if (j==3) {
    d3   <- d
    val3 <- val
    P_3  <- P_mod
  }
}



pdf("dens_P.pdf", width=7.48, height=5.47)

nf<-layout(matrix(c(1,2,3), ncol=1), c(19),c(4.3,3.8,5.8),TRUE)
par(mai=c(0,0.8,0.2,0.2))
# 1
plot(1,type="n",ylim=c(0,0.6),xlim=c(0,30),axes=FALSE,xlab=expression(paste(Pi)), ylab=expression(paste('p','(',P,' | ','1.95',' < ', Pi,' < ','2.05',')')))
axis(side=1,seq(0,50,10))
axis(side=2,seq(0,0.6,0.1),labels=c("",seq(0.1,0.6,0.1)))
box()
abline(h = seq(0,0.6,0.1), v = seq(0,50,10), col = "gray80", lwd=0.25, lty = 2)
polygon(arg, val1, col=rgb(95/255,158/255,160/255,0.5), border="cadetblue4", lwd=1, lty=2)
polygon(d1, col=rgb(240/255,128/255,128/255,0.4), border="coral3", lwd=1, lty=1) 


leg.txt<-(c("Histogram","Inverse Gamma"))
leg.col<-c( "lightcoral","cadetblue")
legend("topright", inset=0,leg.txt,col=leg.col, lty=c(1,2), lwd=2, cex=1.3, pt.cex=1,bg="white")


par(mai=c(0,0.8,0,0.2))
# 2
plot(1,type="n",ylim=c(0,0.4),xlim=c(0,30),axes=FALSE,xlab=expression(paste(P)), ylab=expression(paste('p','(',P,' | ','5.95',' < ',Pi,' < ','6.05',')')))
axis(side=1,seq(0,50,10))
axis(side=2,seq(0,0.4,0.1),labels=c("",seq(0.1,0.4,0.1)))
box()
abline(h = seq(0,0.4,0.1), v = seq(0,50,10), col = "gray80", lwd=0.25, lty = 2)
polygon(arg, val2, col=rgb(95/255,158/255,160/255,0.5), border="cadetblue4", lwd=1, lty=2)


polygon(d2, col=rgb(240/255,128/255,128/255,0.4), border="coral3", lwd=1, lty=1) 

leg.txt<-(c("Histogram","Inverse Gamma"))
leg.col<-c( "lightcoral","cadetblue")
legend("topright", inset=0,leg.txt,col=leg.col, lty=c(1,2), lwd=2, cex=1.3, pt.cex=1, bg="white")

par(mai=c(0.8,0.8,0,0.2))

# 3
plot(1,type="n",ylim=c(0,0.4),xlim=c(0,30), axes=FALSE,xlab=expression(paste(P)), ylab=expression(paste('p','(',P,' | ','9.75',' < ',Pi,' < ','10.25',')')))
axis(side=1,seq(0,50,10))


axis(side=2,seq(0,0.4,0.1),labels=seq(0,0.4,0.1))
box()

abline(h = seq(0,0.4,0.1), v = seq(0,50,10), col = "gray80", lwd=0.25, lty = 2)


polygon(arg, val3, col=rgb(95/255,158/255,160/255,0.5), border="cadetblue4", lwd=1, lty=2)
polygon(d3, col=rgb(240/255,128/255,128/255,0.4), border="coral3", lwd=1, lty=1) 

leg.txt<-(c("Histogram","Inverse Gamma"))
leg.col<-c( "lightcoral","cadetblue")
legend("topright", inset=0,leg.txt,col=leg.col, lty=c(1,2), lwd=2, cex=1.3, pt.cex=1, bg="white")
dev.off()


## black and white #########################
pdf("dens_P_bw.pdf", width=7.48, height=5.47)
nf<-layout(matrix(c(1,2,3), ncol=1), c(19),c(4.3,3.8,5.8),TRUE)

par(mai=c(0,0.8,0.2,0.2))
# 1
plot(1,type="n",ylim=c(0,0.6),xlim=c(0,30),axes=FALSE,xlab=expression(paste(Pi)), ylab=expression(paste('p','(',P,' | ','1.95',' < ', Pi,' < ','2.05',')')))


axis(side=1,seq(0,50,10))
axis(side=2,seq(0,0.6,0.1),labels=c("",seq(0.1,0.6,0.1)))
box()
abline(h = seq(0,0.6,0.1), v = seq(0,50,10), col = "gray80", lwd=0.25, lty = 2)
lines(arg, val1, lty=2, lwd=2)
lines(d1,lty=1,lwd=2)

leg.txt<-(c("Histogram","Inverse Gamma"))
legend("topright", inset=0,leg.txt, lty=c(1,2), lwd=2, cex=1.3, pt.cex=1,bg="white")

par(mai=c(0,0.8,0,0.2))
# 2
plot(1,type="n",ylim=c(0,0.4),xlim=c(0,30),axes=FALSE,xlab=expression(paste(P)), ylab=expression(paste('p','(',P,' | ','5.95',' < ',Pi,' < ','6.05',')')))
axis(side=1,seq(0,50,10))
axis(side=2,seq(0,0.4,0.1),labels=c("",seq(0.1,0.4,0.1)))
box()
abline(h = seq(0,0.4,0.1), v = seq(0,50,10), col = "gray80", lwd=0.25, lty = 2)
lines(arg, val2, lty=2, lwd=2)
lines(d2,lty=1,lwd=2)

leg.txt<-(c("Histogram","Inverse Gamma"))
legend("topright", inset=0,leg.txt, lty=c(1,2), lwd=2, cex=1.3, pt.cex=1, bg="white")

par(mai=c(0.8,0.8,0,0.2))

# 3
plot(1,type="n",ylim=c(0,0.4),xlim=c(0,30), axes=FALSE,xlab=expression(paste(P)), ylab=expression(paste('p','(',P,' | ','9.75',' < ',Pi,' < ','10.25',')')))
axis(side=1,seq(0,50,10))


axis(side=2,seq(0,0.4,0.1),labels=seq(0,0.4,0.1))
box()

abline(h = seq(0,0.4,0.1), v = seq(0,50,10), col = "gray80", lwd=0.25, lty = 2)
lines(arg, val3, lty=2, lwd=2)
lines(d3,lty=1,lwd=2)
leg.txt<-(c("Histogram","Inverse Gamma"))
legend("topright", inset=0,leg.txt, lty=c(1,2), lwd=2, cex=1.3, pt.cex=1, bg="white")
dev.off()

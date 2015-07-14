################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## Evaluate_B.R
#
# Estimate the ``true'' background-error variances B_k for all filters
# as well as the variances of the ``truth'' V_k,
# as in the paper
#
# ``Hierarchical Bayes Ensemble Kalman Filtering''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
#
# Dependencies:
#
# Before running this script, you need to execute the script
# Calculate_data_for_B_evaluation.R,
# which writes the time series needed here (see the rows with read.table(.) below).
# 
# 
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 11 July 2015
################################################################################
source('functions.R')
parameters          <- create_parameters_universe_world()
parameters$time     <- 20000
parameters$L        <- 200

L                   <- parameters$L
X_true_mat          <- matrix(ncol = parameters$time, nrow = L)
X_f_kf_mat          <- matrix(ncol = parameters$time, nrow = L)
X_f_var_mat         <- matrix(ncol = parameters$time, nrow = L)
X_f_enkf_mat          <- matrix(ncol = parameters$time, nrow = L)
X_f_henkf_mat       <- matrix(ncol = parameters$time, nrow = L)
X_f_hbef_mat        <- matrix(ncol = parameters$time, nrow = L)

B_a_kf_mat          <- matrix(ncol = parameters$time, nrow = L)
B_a_var_mat         <- matrix(ncol = parameters$time, nrow = L)
B_a_enkf_mat        <- matrix(ncol = parameters$time, nrow = L)
B_a_henkf_mat       <- matrix(ncol = parameters$time, nrow = L)
B_a_hbef_mat        <- matrix(ncol = parameters$time, nrow = L)

X_true              <- read.table("X_true")
X_f_kf              <- read.table("X_f_kf")
X_f_var             <- read.table("X_f_var")
X_f_enkf            <- read.table("X_f_enkf")
X_f_henkf           <- read.table("X_f_henkf")
X_f_hbef            <- read.table("X_f_hbef")

B_a_kf              <- read.table("B_a_kf")
B_a_var             <- read.table("B_a_var")
B_a_enkf            <- read.table("B_a_enkf")
B_a_henkf           <- read.table("B_a_henkf")
B_a_hbef            <- read.table("B_a_hbef")

for(i in (1:L)){
  X_true_mat[i,]    <- X_true[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  X_f_kf_mat[i,]    <- X_f_kf[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  X_f_var_mat[i,]   <- X_f_var[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  X_f_enkf_mat[i,]  <- X_f_enkf[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  X_f_henkf_mat[i,] <- X_f_henkf[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  X_f_hbef_mat[i,]  <- X_f_hbef[((1+parameters$time*(i-1)):(parameters$time*i)),1]

  B_a_kf_mat[i,]    <- B_a_kf[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  B_a_var_mat[i,]   <- B_a_var[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  B_a_enkf_mat[i,]  <- B_a_enkf[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  B_a_henkf_mat[i,] <- B_a_henkf[((1+parameters$time*(i-1)):(parameters$time*i)),1]
  B_a_hbef_mat[i,]  <- B_a_hbef[((1+parameters$time*(i-1)):(parameters$time*i)),1]
}


B_kf                <- c(1:parameters$time)
B_var               <- c(1:parameters$time)
B_enkf              <- c(1:parameters$time)
B_henkf             <- c(1:parameters$time)
B_hbef              <- c(1:parameters$time)  
var_of_x_true       <- c(1:parameters$time)

for(i in (1:parameters$time)){
  B_kf[i]           <- var(X_f_kf_mat[,i]-X_true_mat[,i])
  B_var[i]          <- var(X_f_var_mat[,i]-X_true_mat[,i])
  B_enkf[i]         <- var(X_f_enkf_mat[,i]-X_true_mat[,i])
  B_henkf[i]        <- var(X_f_henkf_mat[,i]-X_true_mat[,i])
  B_hbef[i]         <- var(X_f_hbef_mat[,i]-X_true_mat[,i])
  var_of_x_true[i]  <- var(X_true_mat[,i])
}

write(var_of_x_true, 'var_of_x_true.txt')
write(B_kf, 'B_kf.txt')
write(B_var, 'B_var.txt')
write(B_enkf, 'B_enkf.txt')
write(B_henkf, 'B_henkf.txt')
write(B_hbef, 'B_hbef.txt')


mean_B_a_kf         <- 0
mean_B_a_var        <- 0
mean_B_a_enkf       <- 0
mean_B_a_henkf      <- 0
mean_B_a_hbef       <- 0

rms_B_a_kf          <- 0   
rms_B_a_var         <- 0
rms_B_a_enkf        <- 0
rms_B_a_henkf       <- 0
rms_B_a_hbef        <- 0

t1                  <- 1
t2                  <- parameters$time

for(t in 1:parameters$L){
  mean_B_a_kf       <- mean_B_a_kf + mean(B_a_kf_mat[t,(t1:t2)] - B_kf[(t1:t2)])
  mean_B_a_var      <- mean_B_a_var + mean(B_a_var_mat[t,(t1:t2)] - B_var[(t1:t2)])
  mean_B_a_enkf     <- mean_B_a_enkf + mean(B_a_enkf_mat[t,(t1:t2)] - B_enkf[(t1:t2)])
  mean_B_a_henkf    <- mean_B_a_henkf + mean(B_a_henkf_mat[t,(t1:t2)] - B_henkf[(t1:t2)])
  mean_B_a_hbef     <- mean_B_a_hbef + mean(B_a_hbef_mat[t,(t1:t2)] - B_hbef[(t1:t2)])
  
  rms_B_a_kf        <- rms_B_a_kf + rms(B_a_kf_mat[t,(t1:t2)] - B_kf[(t1:t2)])
  rms_B_a_var       <- rms_B_a_var + rms(B_a_var_mat[t,(t1:t2)] - B_var[(t1:t2)])
  rms_B_a_enkf      <- rms_B_a_enkf + rms(B_a_enkf_mat[t,(t1:t2)] - B_enkf[(t1:t2)])
  rms_B_a_henkf     <- rms_B_a_henkf + rms(B_a_henkf_mat[t,(t1:t2)] - B_henkf[(t1:t2)])
  rms_B_a_hbef      <- rms_B_a_hbef + rms(B_a_hbef_mat[t,(t1:t2)] - B_hbef[(t1:t2)])  
}


Statistics <- matrix(NA, nrow = 5, ncol = 3)
colnames(Statistics) <- c('Mean(B_hat-B)', 'RMSE(B_hat-B)', 'Mean(B)') 
rownames(Statistics) <- c('KF','Var', 'EnKF', 'HEnKF', 'HBEF')
Statistics[1,1] <- mean_B_a_kf/parameters$L
Statistics[2,1] <- mean_B_a_var/parameters$L
Statistics[3,1] <- mean_B_a_enkf/parameters$L 
Statistics[4,1] <- mean_B_a_henkf/parameters$L
Statistics[5,1] <- mean_B_a_hbef/parameters$L

Statistics[1,2] <- rms_B_a_kf/parameters$L
Statistics[2,2] <- rms_B_a_var/parameters$L
Statistics[3,2] <- rms_B_a_enkf/parameters$L 
Statistics[4,2] <- rms_B_a_henkf/parameters$L
Statistics[5,2] <- rms_B_a_hbef/parameters$L

Statistics[1,3] <- mean(B_kf)
Statistics[2,3] <- mean(B_var)
Statistics[3,3] <- mean(B_enkf) 
Statistics[4,3] <- mean(B_henkf)
Statistics[5,3] <- mean(B_hbef)

print(Statistics)  

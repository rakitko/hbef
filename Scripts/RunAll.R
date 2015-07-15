################################################################################
## This software is in the public domain, furnished "as is", without technical
## support, and with no warranty, express or implied, as to its usefulness for
## any purpose.
## 
## RunAll.R
#
# Launch all the scripts in the HBEF package
# to produce the results reported in the paper 
#
# ``Hierarchical Bayes Ensemble Kalman Filtering''
# by Michael Tsyrulnikov and Alexander Rakitko
# submitted to Physica D.
# 
# 
## Authors: Alexander Rakitko  (rakitko@gmail.com) and Michael Tsyrulnikov 
## 14 July 2015
################################################################################
source('Calculate_data_for_B_evaluation.R')
source('Evaluate_B.R')
source('Filters_plot.R')
source('qqplot_1.R')
source('qqplot_2.R')
source('PfSpe_RMSbias.R')
source('PfSpe_timser.R')
source('QfSme_RMSbias.R')
source('QfSme_timser.R')
source('RMSE_N.R')
source('RMSE_Q_distort.R')
source('RMSE_R.R')
source('RMSE_pi.R')
source('RMSE_sdSigma.R')
source('Timeseries.R')
source('dens_P.R')
source('dens_Pi.R')
source('dens_Q.R')
source('RMSE.R')

# hbef
The Hierarchical Bayes Ensemble Kalman Filter.

The R package HBEF allows one to carry out numerical experiments with the Hierarchical Bayes Ensemble Kalman Filter (HBEF).
The HBEF is proposed and described in the paper

Michael Tsyrulnikov and Alexander Rakitko,
A Hierarchical Bayes ensemble Kalman Filter. 
Physica D (Nonlinear Phenomena), 2017, v.338, 1-16, doi:10.1016/j.physd.2016.07.009.

With this package,  one can reproduce all the results reported in the paper.

The paper proposes a new kind of ensemble Kalman Filter that allows for the uncertainty in the 
model-error and predictability-error covariance matrices. These two matrices are treated
as random matrices and updated along with the state using observations and forecast ensemble data.
In this update, ensemble members are assimilated as generalized observations and 
ordinary observations are allowed to influence the covariances.

The main author of the package is Alexander Rakitko rakitko@gmail.com.

The corresponding author of the paper is Michael Tsyrulnikov mik.tsyrulnikov@gmail.com.

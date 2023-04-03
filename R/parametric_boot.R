# Parametric sampling in R
# Data from binomial(15, θ) for an unknown θ
x = c(3, 5, 7, 9, 11, 13)
binomSize = 15 # known size of binomial
n = length(x) # sample size
thetahat = mean(x)/binomSize # MLE for θ
nboot = 5000 # number of bootstrap samples to use
# nboot parametric samples of size n; organize in a matrix
tmpdata = rbinom(n*nboot, binomSize, thetahat)
bootstrapsample = matrix(tmpdata, nrow=n, ncol=nboot)
# Compute bootstrap means thetahat* and differences delta*
thetahatstar = colMeans(bootstrapsample)/binomSize
deltastar = thetahatstar - thetahat
# Find quantiles and make the bootstrap confidence interval
d = quantile(deltastar, c(.1,.9))
ci = thetahat - c(d[2], d[1])

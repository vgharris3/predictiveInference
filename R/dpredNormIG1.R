#' The Normal-Inverse Gamma Predictive Distribution
#'
#' dpredNormIG1 returns the predictive probability of future observations based on
#'
#' @param ypred vector of values for which predictive probability is desired
#' @param y vector of (observed) Normally-distributed values
#' @param mu0 mean of prior observations (set default)
#' @param k0 number of prior observations (default is 1)
#' @param sig20 variance of prior observations (set default)
#' @param nu0 number of prior observations (default is 1)
#' @param S size of random sample used for density approximation (default is 10^5)
#' @param Jeffreys flag for use of Jeffrey's prior
#'
#' @return random sample of size n Normal - Inverse Gamma predictive probability
#' @export
#'
#' @examples 1
dpredNormIG1 = function(ypred,y,mu0=0,k0=1,sig20=1,nu0=1,S = 100000,Jeffreys=FALSE){

  #ERROR HANDLING

  ##########################
  ##########################
  #NOT ENOUGH / TOO FEW PARAMETERS
  ##########################
  ##########################

  if(k0 <= 0){
    stop("k0 <= 0: k0 must be greater than 0")
    return(1)
  }

  if(sig20 <= 0){
    stop("sig20<=0:  sig20 must be greater than 0")
    return(1)
  }

  if(nu0 <= 0){
    stop("k0 <= 0: nu0 must be greater than 0")
    return(1)
  }

  #Why do I have this limit on the length of ypred??
  if(length(ypred) > 1000){
    stop("length of input vector ypred must not exceed 1000")
  }

  nobs = length(y);
  meanobs = mean(y);
  s2 = stats::var(y);

  if(Jeffreys){

    location = meanobs
    scale = sqrt(s2)*sqrt(1+1/nobs)

    xt = seq(1,2.5,len=100)
    yt = (1/scale) * stats::dt((xt - location)/scale,df = nobs-1)

    #Using t(n-1) distribution resulting from Jeffrey's prior:
    #(theta - ybar)/(s/sqrt(n))|y1,...,yn ~ t(n-1)
    #yields t distribution with location ybar and scale
    #std(y)*sqrt(1+1/n)

    xd = stats::dt((ypred-location)/scale,df=nobs-1)/scale

  } else {

    #Obtain a random sample from which to approximate the density
    rs = rpredNormIG1(S,y,mu0,k0,sig20,nu0,Jeffreys=FALSE)

    #Estimating density using R's density() function on random sample
    #(obtained from https://stackoverflow.com/questions/28077500/find-the-probability-density-of-a-new-data-point-using-density-function-in-r)
    d <- stats::density(rs)
    h = d$bw

    myKDE_vec <- function(tvec){
      kernelValues <- matrix(0, nrow = length(tvec), ncol = length(rs))
      rsmat_old = do.call("rbind",replicate(length(tvec),rs,simplify=FALSE))
      rsmat = t(matrix(replicate(length(tvec),rs),ncol = length(tvec)))
      transformed = (tvec - rsmat)/h
      kernelValues = stats::dnorm(transformed, mean = 0, sd = 1)/h
      return(apply(kernelValues,1,sum)/length(rs))
    }

    xd = myKDE_vec(ypred)
    #xd = 1

  }

  return(xd)
}

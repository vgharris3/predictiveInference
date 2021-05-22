#' The Normal-Inverse Gamma Predictive Distribution
#'
#' ppredNormIG returns the cumulative predictive probability of future observations based on
#'
#' @param x vector of values for which predictive probability is desired
#' @param obs vector of (observed) Normally-distributed values
#' @param mu0 mean of prior observations (set default)
#' @param k0 number of prior observations (default is 1)
#' @param sig20 variance of prior observations (set default)
#' @param nu0 number of prior observations (default is 1)
#' @param S size of random sample used for density approximation (default is 10^6)
#' @param Jeffreys flag for use of Jeffrey's prior
#'
#' @return random sample of size n Normal - Inverse Gamma predictive probability
#' @export
#'
#' @examples 1
ppredNormIG = function(x,obs,mu0=0,k0=1,sig20=1,nu0=1,S = 100000,Jeffreys=FALSE){

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

  if(length(x) > 1000){
    stop("length of input vector x must not exceed 1000")
  }

  nobs = length(obs);
  meanobs = mean(obs);
  s2 = stats::var(obs);

  if(Jeffreys){

    #Using t(n-1) distribution resulting from Jeffrey's prior:
    #(theta - ybar)/(s/sqrt(n))|y1,...,yn ~ t(n-1)
    #Finding t-distribution density for "shifted_scaled"
    #ss = (theta - obs)

    shifted_scaled = stats::dt(obs,df=length(obs)-1)

  } else {

    #Obtain a random sample from which to approximate the density
    rs = rpredNormIG(S,obs,mu0,k0,sig20,nu0,Jeffreys=FALSE)

    #Estimating density using R's density() function on random sample
    #(taken from https://stackoverflow.com/questions/28077500/find-the-probability-density-of-a-new-data-point-using-density-function-in-r)
    #d <- stats::density(rs)
    #h = d$bw

    #myKDE_vec <- function(tvec){
    #  kernelValues <- matrix(0, nrow = length(tvec), ncol = length(rs))
    #  rsmat_old = do.call("rbind",replicate(length(tvec),rs,simplify=FALSE))
    #  rsmat = t(matrix(replicate(length(tvec),rs),ncol = length(tvec)))
    #  transformed = (tvec - rsmat)/h
    #  kernelValues = stats::dnorm(transformed, mean = 0, sd = 1)/h
    #  return(apply(kernelValues,1,sum)/length(rs))
    #}

    #xd = myKDE_vec(x)
    #xd = 1

    sample_log = rep(0,length(rs))
    for(i in 1:length(rs)){
      if(rs[i]<=x){sample_log[i] = 1}
    }


  }

  return(sum(sample_log)/length(sample_log))
}

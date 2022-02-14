#' The Normal-Inverse Gamma Predictive Distribution
#'
#' ppredNormIG1 returns the cumulative predictive probability of future observations based on
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
ppredNormIG1 = function(ypred,y,mu0=0,k0=1,sig20=1,nu0=1,S = 100000,Jeffreys=FALSE){

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

  if(length(ypred) > 1000){
    stop("length of input vector x must not exceed 1000")
  }

  nobs = length(y);
  meanobs = mean(y);
  s2 = stats::var(y);

  if(Jeffreys){

    #Using t(n-1) distribution resulting from Jeffrey's prior:
    #(theta - ybar)/(s/sqrt(n))|y1,...,yn ~ t(n-1)
    #yields t distribution with location ybar and scale
    #std(y)*sqrt(1+1/n)

    location = meanobs
    scale = sqrt(s2)*sqrt(1+1/nobs)

    xp = stats::pt((ypred-location)/scale,df=nobs-1)

  } else {

  #Obtain a random sample from which to approximate the density
  rs = rpredNormIG1(S,y,mu0,k0,sig20,nu0,Jeffreys=FALSE)
  #rsmat = matrix(rpredNormIG(S*length(x),obs,mu0,k0,sig20,nu0,Jeffreys=FALSE),nrow = length(x))

  #xrsmat = cbind(x,rsmat)

  Frs = stats::ecdf(rs)
  xp = Frs(ypred)

  }

  #return(apply(xrsmat,1,function(x) length(x[x<=x[1]])-1)/ncol(rsmat))

  #return(Frs(x))
  return(xp)

}

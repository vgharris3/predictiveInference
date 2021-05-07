#' Title
#'
#' dpredPG returns a random sample of size n from the Poisson-Gamma predictive probability distribution given observations obs~Poi(theta)
#' and theta~Gamma(alpha,beta)
#'
#' @param n vector of integers for which predictive probability is desired
#' @param obs vector of (observed) Poisson-distributed counts
#' @param mu0 mean of prior observations (set default)
#' @param k0 number of prior observations (default is 1)
#' @param sig20 variance of prior observations (set default)
#' @param nu0 number of prior observations (default is 1)
#' @param Jeffreys flag for use of Jeffrey's prior
#'
#' @return random sample of size n Normal - Inverse Gamma predictive probability
#' @export
#'
#' @examples 1
dpredNormIG = function(n,obs,mu0=0,k0=1,sig20=1,nu0=1,Jeffreys=FALSE){

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

  nobs = length(obs);
  meanobs = mean(obs);
  s2 = var(obs);

  if(Jeffreys){

    #Using t(n-1) distribution resulting from Jeffrey's prior:
    #(theta - ybar)/(s/sqrt(n))|y1,...,yn ~ t(n-1)
    #Finding t-distribution density for "shifted_scaled"
    #ss = (theta - obs)

    shifted_scaled = dt(obs,df=length(obs)-1)

  } else {

    rs = rpredNormIG(n,obs,mu0,k0,sig20,nu0,Jeffreys=FALSE)

  }

  return(rs)
}

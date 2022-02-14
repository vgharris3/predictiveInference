#' The Normal-Inverse Gamma Predictive Distribution
#'
#' rpredNormIG1 returns a random sample of size n from the Normal-Inverse Gamma predictive probability distribution
#'
#' @param S desired random sample size
#' @param y vector of (observed) Normally-distributed values
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
rpredNormIG1 = function(S,y,mu0=0,k0=1,sig20=1,nu0=1,Jeffreys=FALSE){

  #ERROR HANDLING

  ##########################
  ##########################
  #NOT ENOUGH / TOO FEW PARASETERS
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


  nobs = length(y);
  meanobs = mean(y);
  s2 = stats::var(y);

  if(Jeffreys){
    location = meanobs
    scale = sqrt(s2)*sqrt(1+1/nobs)
    rs = location + scale * stats::rt(S,nobs-1)
  }

  else{

    kn = k0 + nobs; nun = nu0 + nobs
    mun = (k0*mu0+nobs*meanobs)/kn
    sig2n = (nu0*sig20 + (nobs-1)*s2 + k0*nobs*(meanobs-mu0)^2/kn)/nun

    sig2.postsample = 1/stats::rgamma(S,nun/2,sig2n*nun/2)
    theta.postsample = stats::rnorm(S,mun,sqrt(sig2.postsample/kn))

    rs = numeric(S)

    rs = stats::rnorm(S,theta.postsample,sqrt(sig2.postsample))

  }

  return(rs)
}

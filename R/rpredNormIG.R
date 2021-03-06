#' The Normal-Inverse Gamma Predictive Distribution
#'
#' rpredNormIG returns a random sample of size n from the Normal-Inverse Gamma predictive probability distribution
#'
#' @param n desired random sample size
#' @param obs vector of (observed) Normally-distributed values
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
rpredNormIG = function(n,obs,mu0=0,k0=1,sig20=1,nu0=1,Jeffreys=FALSE){

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
  s2 = stats::var(obs);

  if(Jeffreys){
    #s2 = var(y)
    #meanobs = mean(obs)
    #n = length(y)
    location = meanobs
    scale = sqrt(s2)*sqrt(1+1/nobs)
    #xt = seq(1,2.5,len=100)
    #yt = (1/scale) * dt((xt - location)/scale,df = nobs-1)

    #rs = metRology::rt.scaled(n,nobs-1,location,scale)
    rs = location + scale * stats::rt(n,nobs-1)
  }

  else{

    #nobs = length(obs);
    #meanobs = mean(obs);
    #s2 = stats::var(obs);

    kn = k0 + nobs; nun = nu0 + nobs
    mun = (k0*mu0+nobs*meanobs)/kn
    sig2n = (nu0*sig20 + (nobs-1)*s2 + k0*nobs*(meanobs-mu0)^2/kn)/nun

    sig2.postsample = 1/stats::rgamma(n,nun/2,sig2n*nun/2)
    theta.postsample = stats::rnorm(n,mun,sqrt(sig2.postsample/kn))

    #plot(theta.postsample,sig2.postsample,pch = 20)

    rs = numeric(n)

   # for(i in 1:n){
    #  rs[i] = rnorm(1,theta.postsample,sqrt(sig2.postsample))
    #}

    rs = stats::rnorm(n,theta.postsample,sqrt(sig2.postsample))

  }

  return(rs)
}

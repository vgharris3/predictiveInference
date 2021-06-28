#' The Normal-Inverse Gamma Predictive Distribution for Two Samples
#'
#' rpredNormIG2 returns a random sample of size n from the Normal-Inverse Gamma predictive probability distribution for 2 samples
#'
#' @param n desired random sample size (default is 1)
#' @param y1 vector of (observed) Normally-distributed values (sample 1)
#' @param y2 vector of (observed) Normally-distributed values (sample 2)
#' @param mu0 mean of pooled average from prior knowledge (default is 0)
#' @param g20 variance of pooled average from prior knowledge (default is 1)
#' @param d0 mean of delta from prior knowledge, where delta is half the population difference in means (default is 0)
#' @param t20 variance of delta from prior knowledge (default is 1)
#' @param s20 variance of sample from prior knowledge (default is 1)
#' @param nu0 number of prior observations (default is 1)
#'
#' @return random sample of size n Normal - Inverse Gamma predictive probability for each population sample
#' @export
#'
#' @examples 1
rpredNormIG2 = function(n=1,y1,y2,mu0=0,g20=1,d0=0,t20=1,s20=1,nu0=1){

  #ERROR HANDLING

  ##########################
  ##########################
  #NOT ENOUGH / TOO FEW PARAMETERS
  ##########################
  ##########################

  if(n <= 0){
    warning("n <= 0.  Setting n = 1.")
    n = 1
  }

  if(n!=round(n)){
    warning("n is not a whole number.  Setting n = round(n).")
    n = round(n)
  }

  if(g20 <= 0){
    stop("g20<=0:  g20 must be greater than 0")
    return(1)
  }

  if(t20 <= 0){
    stop("t20<=0:  t20 must be greater than 0")
    return(1)
  }

  if(s20 <= 0){
    stop("s20<=0:  s20 must be greater than 0")
    return(1)
  }

  if(nu0 <= 0){
    stop("k0 <= 0: nu0 must be greater than 0")
    return(1)
  }

  n1<-length(y1) ; n2<-length(y2)

  ## prior parameters
#  mu0<-50 ; g20<-625
#  d0<-0 ; t20<-625
#  s20<-100; nu0<-1

  ## starting values
  mu<- ( mean(y1) + mean(y2) )/2
  del<- ( mean(y1) - mean(y2) )/2

  ## Gibbs sampler
  MU<-DEL<-S2<-NULL
  Y12<-NULL
#  set.seed(1)
  for(s in 1:n)
  {

    ##update s2
    s2<-1/stats::rgamma(1,(nu0+n1+n2)/2,
                 (nu0*s20+sum((y1-mu-del)^2)+sum((y2-mu+del)^2) )/2)
    ##

    ##update mu
    var.mu<-  1/(1/g20+ (n1+n2)/s2 ) #this is gamma^2_n -- call it g2n
    mean.mu<- var.mu*( mu0/g20 + sum(y1-del)/s2 + sum(y2+del)/s2 )  #this is mu_n -- call it mun
    mu<-stats::rnorm(1,mean.mu,sqrt(var.mu))
    ##

    ##update del
    var.del<-  1/(1/t20+ (n1+n2)/s2 )
    mean.del<- var.del*( d0/t20 + sum(y1-mu)/s2 - sum(y2-mu)/s2 )
    del<-stats::rnorm(1,mean.del,sqrt(var.del))
    ##

    ##save parameter values
    MU<-c(MU,mu) ; DEL<-c(DEL,del) ; S2<-c(S2,s2)
    Y12<-rbind(Y12,c(stats::rnorm(2,mu+c(1,-1)*del,sqrt(s2) ) ) )
  }

  rs = Y12
  #HOW TO RETURN SEVERAL THINGS ACCESSIBLE VIA $??

  return(rs)
}

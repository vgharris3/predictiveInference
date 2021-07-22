#' The Normal-Inverse Gamma Predictive Distribution for k Samples
#'
#' rpredNormIGk returns a random sample of size n from the Normal-Inverse Gamma predictive probability distribution for 2 samples
#'
#' @param N desired random sample size
#' @param Y 2-column matrix of multiple samples of observed Normally-distributed values.  Column 1:  integer indices identifying the samples.  Column 2:  data values for the sample indicated by the corresponding index in column 1.
#' @param nu0 number of prior observations
#' @param s20 within-sample variance from prior knowledge
#' @param eta0 number of prior observations
#' @param t20 between-sample variance from prior knowledge
#' @param mu0 mean of pooled average from prior knowledge
#' @param g20 variance of pooled average from prior knowledge
#'
#' @return random sample of size n Normal - Inverse Gamma predictive probability for each population sample
#' @export
#'
#' @examples 1
rpredNormIGk = function(N=1,Y,nu0=1,s20=1,eta0=1,t20=1,mu0=0,g20=1){

  #ERROR HANDLING

  ##########################
  ##########################
  #NOT ENOUGH / TOO FEW PARAMETERS
  ##########################
  ##########################

  if(N <= 0){
    warning("N <= 0.  Setting N = 1.")
    N = 1
  }

  if(N!=round(N)){
    warning("N is not a whole number.  Setting N = round(n).")
    N = round(N)
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

  if(eta0 <= 0){
    stop("eta0 <= 0: eta0 must be greater than 0")
    return(1)
  }

#   n1<-length(y1) ; n2<-length(y2)
#
#   ## prior parameters
# #  mu0<-50 ; g20<-625
# #  d0<-0 ; t20<-625
# #  s20<-100; nu0<-1
#
#   ## starting values
#   mu<- ( mean(y1) + mean(y2) )/2
#   del<- ( mean(y1) - mean(y2) )/2
#
#   ## Gibbs sampler
#   MU<-DEL<-S2<-NULL
#   Y12<-NULL
# #  set.seed(1)
#   for(s in 1:n)
#   {
#
#     ##update s2
#     s2<-1/stats::rgamma(1,(nu0+n1+n2)/2,
#                  (nu0*s20+sum((y1-mu-del)^2)+sum((y2-mu+del)^2) )/2)
#     ##
#
#     ##update mu
#     var.mu<-  1/(1/g20+ (n1+n2)/s2 ) #this is gamma^2_n -- call it g2n
#     mean.mu<- var.mu*( mu0/g20 + sum(y1-del)/s2 + sum(y2+del)/s2 )  #this is mu_n -- call it mun
#     mu<-stats::rnorm(1,mean.mu,sqrt(var.mu))
#     ##
#
#     ##update del
#     var.del<-  1/(1/t20+ (n1+n2)/s2 )
#     mean.del<- var.del*( d0/t20 + sum(y1-mu)/s2 - sum(y2-mu)/s2 )
#     del<-stats::rnorm(1,mean.del,sqrt(var.del))
#     ##
#
#     ##save parameter values
#     MU<-c(MU,mu) ; DEL<-c(DEL,del) ; S2<-c(S2,s2)
#     Y12<-rbind(Y12,c(stats::rnorm(2,mu+c(1,-1)*del,sqrt(s2) ) ) )
#   }
#
#   rs = Y12
  #HOW TO RETURN SEVERAL THINGS ACCESSIBLE VIA $??
  #### MCMC approximation to posterior for the hierarchical normal model

  ## weakly informative priors
  # nu0<-1  ; s20<-100
  # eta0<-1 ; t20<-100
  # mu0<-50 ; g20<-25

  ## starting values
  m<-length(Y)
  ytilde<-n<-sv<-ybar<-rep(NA,m)
  for(j in 1:m)
  {
    ybar[j]<-mean(Y[[j]])
    sv[j]<-var(Y[[j]])
    n[j]<-length(Y[[j]])
  }
  theta<-ybar
  sigma2<-mean(sv)
  mu<-mean(theta)
  tau2<-var(theta)

  ## setup MCMC
  # set.seed(1)
  #S<-5000
  S<-N
  YTILDE<-matrix(nrow=S,ncol=m)
  THETA<-matrix( nrow=S,ncol=m)
  MST<-matrix( nrow=S,ncol=3)

  ## MCMC algorithm
  for(s in 1:S)
  {

    # sample new values of the thetas
    for(j in 1:m)
    {
      vtheta<-1/(n[j]/sigma2+1/tau2)
      etheta<-vtheta*(ybar[j]*n[j]/sigma2+mu/tau2)
      theta[j]<-rnorm(1,etheta,sqrt(vtheta))
    }

    #sample new value of sigma2
    nun<-nu0+sum(n)
    ss<-nu0*s20;for(j in 1:m){ss<-ss+sum((Y[[j]]-theta[j])^2)}
    sigma2<-1/rgamma(1,nun/2,ss/2)

    #predict ytilde_j for each theta_j
    for(j in 1:m){
      ytilde[j] = rnorm(1,mean=theta[j],sd=sqrt(sigma2))
    }

    #sample a new value of mu
    vmu<- 1/(m/tau2+1/g20)
    emu<- vmu*(m*mean(theta)/tau2 + mu0/g20)
    mu<-rnorm(1,emu,sqrt(vmu))

    # sample a new value of tau2
    etam<-eta0+m
    ss<- eta0*t20 + sum( (theta-mu)^2 )
    tau2<-1/rgamma(1,etam/2,ss/2)

    #store results
    YTILDE[s,]<-ytilde
    THETA[s,]<-theta
    MST[s,]<-c(mu,sigma2,tau2)

  }


  mcmc1<-list(YTILDE=YTILDE,THETA=THETA,MST=MST)

#Note: in mcmc1 each row contains a single prediction for each separate data set.
#The number of rows is the number of desired predictions (input N).
#The ith column therefore contains all the predictions for the ith data set.

  return(mcmc1)
}

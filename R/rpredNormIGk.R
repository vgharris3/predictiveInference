#' The Normal-Inverse Gamma Predictive Distribution for k Samples
#'
#' rpredNormIGk returns a random sample of size N from the Normal-Inverse Gamma predictive probability distribution for each of the k samples, along with N posterior draws of the within-group and between-groups unknown parameters.
#'
#' @param S desired random sample size
#' @param Y 2-column matrix of multiple samples of observed Normally-distributed values.  Column 1:  integer indices identifying the samples.  Column 2:  data values for the sample indicated by the corresponding index in column 1.
#' @param nu0 number of prior observations
#' @param s20 within-sample variance from prior knowledge
#' @param eta0 number of prior observations
#' @param t20 between-sample variance from prior knowledge
#' @param mu0 mean of pooled average from prior knowledge
#' @param g20 variance of pooled average from prior knowledge
#'
#' @return random sample of size S Normal - Inverse Gamma predictive probability for each population sample
#' @export
#'
#' @examples 1
rpredNormIGk = function(S=1,Y,nu0=1,s20=1,eta0=1,t20=1,mu0=0,g20=1){

  #ERROR HANDLING

  ##########################
  ##########################
  #NOT ENOUGH / TOO FEW PARAMETERS
  ##########################
  ##########################

  if(S <= 0){
    warning("S <= 0.  Setting S = 1.")
    S = 1
  }

  if(S!=round(S)){
    warning("S is not a whole number.  Setting S = round(S).")
    S = round(S)
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

  #### MCMC approximation to posterior for the hierarchical normal model

  ## Put data in list form
  Ylist<-list()
  J<-max(Y[,1]) #col. 1 of Y is an index identifying the data in col. 2
  n<-ybar<-ymed<-s2<-rep(0,J)
  for(j in 1:J) {
    Ylist[[j]]<-Y[ Y[,1]==j,2] #separate out data into separate vectors by index.  Y is a list of those vectors, which need not be all the same length
    ybar[j]<-mean(Ylist[[j]]) #record mean of each data set
    ymed[j]<-median(Ylist[[j]]) #record median of each data set
    n[j]<-length(Ylist[[j]]) #record number data points in each data set
    s2[j]<-var(Ylist[[j]]) #record variance of each data set
  }

  ## starting values
  m<-length(Ylist)
  ytilde<-n<-sv<-ybar<-rep(NA,m)
  for(j in 1:m)
  {
    ybar[j]<-mean(Ylist[[j]])
    sv[j]<-var(Ylist[[j]])
    n[j]<-length(Ylist[[j]])
  }
  theta<-ybar
  sigma2<-mean(sv)
  mu<-mean(theta)
  tau2<-var(theta)

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

  #YTILDE:  Nxm matrix of predictions (N predictions for each of the m data sets in Y)
  #THETA:  Nxm matrix of drawn within-group posterior means
  #MST:  Nx3 matrix of drawn posterior within-group variance and between-groups mean and variance

  return(mcmc1)
}

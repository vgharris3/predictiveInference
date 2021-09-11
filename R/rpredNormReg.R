#' Normal Regression Predictive Inference
#'
#' rpredNormReg returns a random sample of size ? from ? predictive probability distribution
#'
#' @param N desired random sample size
#' @param Xpred predictor vector(s) for which predictions are desired
#' @param X Explanatory variables matrix
#' @param y Vector of prior observations
#' @param nu0 number of prior observations
#' @param s20 ols estimate of variance (perform regression and compute this instead of user input?)
#' @param gprior flag to use Zellner's g-prior.  Default to TRUE
#'
#' @return random sample of size ? from ? predictive probability distribution
#' @export
#'
#' @examples 1
rpredNormReg = function(N=1,Xpred,X,y,nu0=1,s20=1,gprior = TRUE){

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

   if(nu0 <= 0){
    stop("k0 <= 0: nu0 must be greater than 0")
    return(1)
  }

  #Other checks:
  #X is singluar
  #Dimensions of X and y

  ###############################
  ###############################
  ###############################

  n = dim(X)[1] # = length(y)
  p = dim(X)[2] # = length(beta)
  g = length(y) #Zellner's g-prior

  if(gprior){

    Hg = (g/(g+1))*X%*%solve(t(X)%*%X)%*%t(X)
    SSRg = t(y)%*%( diag(1,nrow=n) - Hg)%*%y

    s2 = 1/rgamma(N, (nu0 + n)/2, (nu0*s20 + SSRg)/2 )

    Vb = g*solve(t(X)%*%X)/(g+1)
    Eb = Vb%*%t(X)%*%y

    E = matrix(rnorm(N*p,0,sqrt(s2)),N,p)
    beta = t( t(E%*%chol(Vb)) + c(Eb))

    result = Xpred%*%t(beta) # compute y prediction for input X vector for which prediction is desired.

  } else { #DO THIS PART TOO

    # fit.ls<-lm(y~-1+ X)
    # beta.0<-rep(0,p) ; Sigma.0<-diag(c(150,30,6,5)^2,p)
    # nu.0<-1 ; sigma2.0<- 15^2
    #
    # beta.0<-fit.ls$coef
    # nu.0<-1  ; sigma2.0<-sum(fit.ls$res^2)/(n-p)
    # Sigma.0<- solve(t(X)%*%X)*sigma2.0*n
    #
    #
    #
    # S<-5000
    #
    # rmvnorm<-function(n,mu,Sigma)
    # { # samples from the multivariate normal distribution
    #   E<-matrix(rnorm(n*length(mu)),n,length(mu))
    #   t(  t(E%*%chol(Sigma)) +c(mu))
    # }

    #MORE TO THIS

    #### Bayesian estimation via MCMC
    # n<-length(y)
    # X<-cbind(rep(1,n),x1,x2,x1*x2)
    # p<-dim(X)[2]
    #
    # fit.ls<-lm(y~-1+ X)
    # beta.0<-rep(0,p) ; Sigma.0<-diag(c(150,30,6,5)^2,p)
    # nu.0<-1 ; sigma2.0<- 15^2
    #
    # beta.0<-fit.ls$coef
    # nu.0<-1  ; sigma2.0<-sum(fit.ls$res^2)/(n-p)
    # Sigma.0<- solve(t(X)%*%X)*sigma2.0*n
    #
    #
    # S<-5000
    #
    # rmvnorm<-function(n,mu,Sigma)
    # { # samples from the multivariate normal distribution
    #   E<-matrix(rnorm(n*length(mu)),n,length(mu))
    #   t(  t(E%*%chol(Sigma)) +c(mu))
    # }
    #
    # ## some convenient quantites
    # n<-length(y)
    # p<-length(beta.0)
    # iSigma.0<-solve(Sigma.0)
    # XtX<-t(X)%*%X
    #
    # ## store mcmc samples in these objects
    # beta.post<-matrix(nrow=S,ncol=p)
    # sigma2.post<-rep(NA,S)
    #
    # ## starting value
    # set.seed(1)
    # sigma2<- var( residuals(lm(y~0+X)) )
    #
    # ## MCMC algorithm
    # for( scan in 1:S) {
    #
    #   #update beta
    #   V.beta<- solve(  iSigma.0 + XtX/sigma2 )
    #   E.beta<- V.beta%*%( iSigma.0%*%beta.0 + t(X)%*%y/sigma2 )
    #   beta<-t(rmvnorm(1, E.beta,V.beta) )
    #
    #   #update sigma2
    #   nu.n<- nu.0+n
    #   ss.n<-nu.0*sigma2.0 + sum(  (y-X%*%beta)^2 )
    #   sigma2<-1/rgamma(1,nu.n/2, ss.n/2)
    #
    #   #save results of this scan
    #   beta.post[scan,]<-beta
    #   sigma2.post[scan]<-sigma2
    # }
    #
    # round( apply(beta.post,2,mean), 3)

  }

  return(result)
}

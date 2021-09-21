#' Normal Regression Predictive Inference
#'
#' rpredNormReg returns a random sample of size ? from ? predictive probability distribution
#'
#' @param S desired random sample size
#' @param Xpred predictor vector(s) for which predictions are desired
#' @param X Explanatory variables matrix
#' @param y Vector of prior observations
#' @param nu0 number of prior observations
#' @param beta0 vector of prior regression coefficients
#' @param Sigma0 prior variance-covariance matrix
#' @param s20 prior variance
#' @param gprior flag to use Zellner's g-prior.  Default to TRUE
#'
#' @return random sample of size ? from ? predictive probability distribution
#' @export
#'
#' @examples 1
rpredNormReg = function(S=1,Xpred,X,y,beta0,Sigma0,nu0=1,s20=1,gprior = TRUE){

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
    warning("S is not a whole number.  Setting S = round(n).")
    S = round(S)
  }

   if(nu0 <= 0){
    stop("k0 <= 0: nu0 must be greater than 0")
    return(1)
  }

  #Other checks:
  #X is singluar
  #Dimensions of X and y
  #Dimensions of Sigma0

  ###############################
  ###############################
  ###############################

  n = dim(X)[1] # = length(y)
  p = dim(X)[2] # = length(beta)

  if(gprior){

    g = length(y) # Zellner's g-prior
    Hg = (g/(g+1))*X%*%solve(t(X)%*%X)%*%t(X)
    SSRg = t(y)%*%( diag(1,nrow=n) - Hg)%*%y

    s2 = 1/rgamma(S, (nu0 + n)/2, (nu0*s20 + SSRg)/2 )

    Vb = g*solve(t(X)%*%X)/(g+1)
    Eb = Vb%*%t(X)%*%y

    E = matrix(rnorm(S*p,0,sqrt(s2)),S,p)
    beta = t( t(E%*%chol(Vb)) + c(Eb))

    result$betas = beta
    result$predictions = Xpred%*%t(beta) # compute y prediction for input X vector for which prediction is desired.

  } else { #DO THIS PART TOO

    rmvnorm<-function(n,mu,Sigma)
    { # samples from the multivariate normal distribution
      E<-matrix(rnorm(n*length(mu)),n,length(mu))
      t(  t(E%*%chol(Sigma)) +c(mu))
    }

    ## some convenient quantities
    #n<-length(y)
    #p<-length(beta.0)
    iSigma.0<-solve(Sigma.0)                            #iSigma.0 = inverse of Sigma.0
    XtX<-t(X)%*%X

    ## store mcmc samples in these objects
    beta.post<-matrix(nrow=S,ncol=p)                    #storage for S instances for each coefficient of the p predictors
    sigma2.post<-rep(NA,S)                              #storage for S instances of variance corresponding to the betas

    ## starting value
    #set.seed(1)
    sigma2<- s20                #starting value for sigma2
    beta.0 = beta0              #starting value for beta, which gets reused in computing E(beta|y,X,sig2)

    ## MCMC algorithm
    for( scan in 1:S) {

      #update beta                                        #Formulas and steps Hoff p. 155
      V.beta<- solve(  iSigma.0 + XtX/sigma2 )            #Conditional variance of the regression coefficients
      E.beta<- V.beta%*%( iSigma.0%*%beta.0 + t(X)%*%y/sigma2 ) #Conditional mean of the variance coefficients
      beta<-t(rmvnorm(1, E.beta,V.beta) )                 #Gibbs sampler step:  update beta

      #update sigma2
      nu.n<- nu.0+n                                       #numerator of 1st term of inverse gamma sample
      ss.n<-nu.0*sigma2.0 + sum(  (y-X%*%beta)^2 )        #numerator of 2nd term of inverse gamma sample
      sigma2<-1/rgamma(1,nu.n/2, ss.n/2)                  #Gibbs ampler step:  update sigma^2

      #save results of this scan
      beta.post[scan,]<-beta                              #Store updated beta in current row of beta matrix
      sigma2.post[scan]<-sigma2                           #Store updated sigma^2 in current position in sigma^2 vector
    }

    round( apply(beta.post,2,mean), 3)                  #compute mean of Gibbs sampled betas (for check)

#    result = list(beta.post, Xpred%*%t(beta.post))

    result$betas = beta.post
    result$predictions = Xpred%*%t(beta.post)

  }

  return(result)
}

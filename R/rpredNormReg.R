#' Normal Regression Predictive Inference
#'
#' rpredNormReg returns a random sample of size ? from ? predictive probability distribution
#'
#' @param S desired random sample size
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
rpredNormReg = function(S=1,Xpred,X,y,nu0=1,s20=1,gprior = TRUE){

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

  ###############################
  ###############################
  ###############################

  n = dim(X)[1] # = length(y)
  p = dim(X)[2] # = length(beta)
  g = length(y) # Zellner's g-prior

  if(gprior){

    Hg = (g/(g+1))*X%*%solve(t(X)%*%X)%*%t(X)
    SSRg = t(y)%*%( diag(1,nrow=n) - Hg)%*%y

    s2 = 1/rgamma(S, (nu0 + n)/2, (nu0*s20 + SSRg)/2 )

    Vb = g*solve(t(X)%*%X)/(g+1)
    Eb = Vb%*%t(X)%*%y

    E = matrix(rnorm(S*p,0,sqrt(s2)),S,p)
    beta = t( t(E%*%chol(Vb)) + c(Eb))

    result = Xpred%*%t(beta) # compute y prediction for input X vector for which prediction is desired.
    #result = Xpred%*%beta # compute y prediction for input X vector for which prediction is desired.

  } else { #DO THIS PART TOO



  }

  return(result)
}

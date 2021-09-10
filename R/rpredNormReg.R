#' Normal Regression Predictive Inference
#'
#' rpredNormReg returns a random sample of size ? from ? predictive probability distribution
#'
#' @param N desired random sample size
#' @param X Explanatory variables matrix
#' @param y Vector of prior observations
#' @param nu0 number of prior observations
#' @param s20 ols estimate of variance [perform regression and compute this instead of user input?]
#'
#' @return random sample of size ? from ? predictive probability distribution
#' @export
#'
#' @examples 1
rpredNormReg = function(N=1,X,y,nu0=1,s20=1){

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

  n<-length(y)
  p<-dim(X)[2]

  fit.ls<-lm(y~-1+ X)
  beta.0<-rep(0,p) ; Sigma.0<-diag(c(150,30,6,5)^2,p)
  nu.0<-1 ; sigma2.0<- 15^2

  beta.0<-fit.ls$coef
  nu.0<-1  ; sigma2.0<-sum(fit.ls$res^2)/(n-p)
  Sigma.0<- solve(t(X)%*%X)*sigma2.0*n

  g = length(y) #Zellner's g-prior



  return(1)
}

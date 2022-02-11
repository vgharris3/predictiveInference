#' The Exponential-Gamma Predictive Distribution
#'
#' @param ypred survival time for which predictive density is desired
#' @param y ordered list of observed survival times and censored survival times (y1,...,yN), yi ~ exp(theta), theta ~ Gamma(dt, gm)
#' @param c censoring indicator (vector of 0s and 1s for censored and fully observed, respectively)
#' @param dt Gamma shape parameter for distribution of theta
#' @param gm Gamma rate parameter for distribution of theta
#'
#' @return Exponential-Gamma predictive probability
#' @export
#'
#' @examples 1
dpredEG = function(ypred, y, c, dt, gm){

  #dpredEG returns the predictive probability of surviving past [OR RATHER DYING AT] time x, given d observed and N-d censored copies out of XN total events
  #XN = (x1,...,xN), xi ~ exp(theta)
  #theta ~ Gamma(dt, gm)
  #d < n s.t. (x1,...,xd) are fully observed and (xd+1, ..., xN) are censored

  #p+pred+GE
  #p = cdf (like R distribution functions)
  #pred = prediction
  #GE = Gamma Exponential

  #x is a time past which dpredEG will return the predictive probability of survival Pr(X = x)
  #XN = (x1,...,xN), xi ~ exp(theta) where theta ~ Gamma(dt, gm)
  #d < n s.t. (x1,...,xd) are fully observed and (xd+1, ..., xN) are censored

  #ERROR HANDLING

  if(d > length(y)){
    stop("d > length(y):  The number of observed copies cannot exceed the total number of copies")
    return (1)
  }

  if(dt <= 0){
    stop("dt<=0:  dt must be greater than 0")
    return(1)
  }

  if(gm <= 0){
    stop("gm <= 0: gm must be greater than 0")
    return(1)
  }

  N = length(y)
  ybar = mean(y)

  d = sum(c) # number of fully observed copies

  numerator = log(d+dt) + (d+dt)*log(gm+N*ybar)
  denominator = (d+dt+1)*log(gm+N*ybar+ypred)

  f_y = exp(numerator - denominator)

  return(f_y)
}

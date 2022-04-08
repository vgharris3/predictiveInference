#' The Poisson-Gamma Predictive Distribution
#'
#' ppredPG returns the cumulative predictive probability of future observation having value {1,...,xmax} given observations obs~Poi(theta)
#' and theta~Gamma(alpha,beta)
#'
#' @param ypred non-negative integer count for which predictive probability is desired
#' @param y vector of (observed) Poisson-distributed counts
#' @param alpha sum of counts from beta prior observations for gamma prior distribution on theta
#' @param beta number of prior observations for gamma prior distribution on theta
#' @param lower logical; if TRUE (default), probabilities are P(R<=q) otherwise P(R>q)
#'
#' @return Poisson-Gamma cumulative predictive probability
#' @export
#'
#' @examples 1
ppredPG = function(ypred, y, alpha = 1, beta=1, lower=TRUE){

  #d+pred+PG
  #d = density (like R distribution functions)
  #pred = prediction
  #PG = Poisson-Gamma

  #ERROR HANDLING

  if(min(ypred) < 0){
    stop("ypred < 0:  values of ypred must be non-negative integers")
    return (1)
  }

  if(alpha <= 0){
    stop("a<=0:  a must be greater than 0")
    return(1)
  }

  if(beta <= 0){
    stop("b <= 0: b must be greater than 0")
    return(1)
  }

  if(min(y) < 0){
    stop("All observations must be non-negative")
  }

  evalVec = 0:max(ypred)

  f_y = dpredPG(0:max(ypred),y,alpha,beta)

  F_y = cumsum(f_y)

  returnVec = F_y[which(evalVec %in% ypred)]

  if(lower){return(returnVec)}
  else {return(1 - returnVec)}

}

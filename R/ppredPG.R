#' Title
#'
#' ppredPG returns the cumulative predictive probability of future observation having value {1,...,xmax} given observations obs~Poi(theta)
#' and theta~Gamma(alpha,beta)
#'
#' @param obs vector of (observed) Poisson-distributed counts
#' @param x non-negative integer count for which predictive probability is desired
#' @param alpha sum of counts from beta prior observations for gamma prior distribution on theta
#' @param beta number of prior observations for gamma prior distribution on theta
#' @param lower logical; if TRUE (default), probabilities are P(R<=q) otherwise P(R>q)
#'
#' @return Poisson-Gamma cumulative predictive probability
#' @export
#'
#' @examples 1
ppredPG = function(obs, x, alpha = 1, beta=1, lower=TRUE){

  #d+pred+PG
  #d = density (like R distribution functions)
  #pred = prediction
  #PG = Poisson-Gamma

  #ERROR HANDLING

  if(min(x) < 0){
    stop("x < 0:  xmax must be a non-negative integer")
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

  if(min(obs) <= 0){
    stop("All observations must be non-negative")
  }

  f_x = dpredPG(obs,x,alpha,beta)

  F_x = cumsum(f_x)

  if(lower){return(F_x)}
  else {return(1 - F_x)}

}

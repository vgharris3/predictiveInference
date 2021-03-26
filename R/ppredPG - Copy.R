#' Title
#'
#' ppredBB returns the predictive probability of <= x future successes out of M trials, given s observed successes out of a sample of
#size N, and user input shape parameters for Beta prior on Pr(success)
#'
#' @param obs vector of (observed) Poisson-distributed counts
#' @param xmax maximum integer for which predictive probability is desired
#' @param alpha sum of counts from beta prior observations for gamma prior distribution on theta
#' @param beta number of prior observations for gamma prior distribution on theta
#' @param lower logical; if TRUE (default), probabilities are P(R<=q) otherwise P(R>q)
#'
#' @return Poisson-Gamma cumulative predictive probability
#' @export
#'
#' @examples 1
ppredPG = function(obs, xmax, alpha = 1, beta=1, lower=TRUE){

  #d+pred+PG
  #d = density (like R distribution functions)
  #pred = prediction
  #PG = Poisson-Gamma

  #ERROR HANDLING

  if(xmax < 0){
    stop("xmax < 0:  xmax must be a non-negative integer")
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

  f_x = dpredPG(obs,xmax,alpha,beta)

  F_x = cumsum(f_x)

  if(lower){return(F_x[x])}
  else {return(1 - F_x[x])}

}

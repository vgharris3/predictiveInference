#' Title
#'
#' dpredPG returns the predictive probability of a future observation having value {1,...,xmax} give observations obs~Poi(theta)
#' and theta~Gamma(alpha,beta)
#'
#'
#' @param obs vector of (observed) Poisson-distributed counts
#' @param xmax maximum integer for which predictive probability
#' @param alpha sum of counts from beta prior observations for gamma prior distribution on theta
#' @param beta number of prior observations for gamma prior distribution on theta
#'
#' @return Poisson-Gamma predictive probability
#' @export
#'
#' @examples 1
dpredPG <- function(obs, xmax, alpha = 1, beta=1){

  #dpredPG returns the predictive probability of x future successes out of M trials, given s observed successes out of a sample of
  #size N, and user input shape parameters for Beta prior on Pr(success)

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

  #start here.

  #r[] is the number of successes in the M future observations

  #r = x;

  #numerator = lgamma(M+1) + lgamma(N+alpha+beta) + lgamma(r+s+alpha) + lgamma(M+N-r-s+beta);
  #denominator = lgamma(r+1) + lgamma(M-r+1) + lgamma(alpha+s) + lgamma(N-s+beta) + lgamma(M+N+alpha+beta);
  #f_x = exp(numerator - denominator)

  #return(f_x)
  return(sum(obs));
}

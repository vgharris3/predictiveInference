#' The Poisson-Gamma Predictive Distribution
#'
#' dpredPG returns the predictive probability of a future observation having value {1,...,x} given observations obs~Poi(theta)
#' and theta~Gamma(alpha,beta)
#'
#'
#' @param ypred vector of integers for which predictive probability is desired
#' @param y vector of (observed) Poisson-distributed counts
#' @param alpha sum of counts from beta prior observations for gamma prior distribution on theta
#' @param beta number of prior observations for gamma prior distribution on theta
#'
#' @return Poisson-Gamma predictive probability
#' @export
#'
#' @examples 1
dpredPG <- function(ypred, y, alpha = 1, beta=1){

  #d+pred+PG
  #d = density (like R distribution functions)
  #pred = prediction
  #PG = Poisson-Gamma

  #ERROR HANDLING

  if(min(ypred) < 0){
    stop("y < 0:  values of y must be non-negative integers")
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
    stop("All yervations must be non-negative")
  }

  a = alpha;
  b = beta;
  sobs = sum(y);
  n = length(y);

  f_y = stats::dnbinom(ypred,size = a + sobs, mu = (a + sobs)/(b + n));

  return(f_y);

}

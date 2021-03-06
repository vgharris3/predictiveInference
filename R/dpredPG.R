#' The Poisson-Gamma Predictive Distribution
#'
#' dpredPG returns the predictive probability of a future observation having value {1,...,x} given observations obs~Poi(theta)
#' and theta~Gamma(alpha,beta)
#'
#'
#' @param x vector of integers for which predictive probability is desired
#' @param obs vector of (observed) Poisson-distributed counts
#' @param alpha sum of counts from beta prior observations for gamma prior distribution on theta
#' @param beta number of prior observations for gamma prior distribution on theta
#'
#' @return Poisson-Gamma predictive probability
#' @export
#'
#' @examples 1
dpredPG <- function(x, obs, alpha = 1, beta=1){

  #d+pred+PG
  #d = density (like R distribution functions)
  #pred = prediction
  #PG = Poisson-Gamma

  #ERROR HANDLING

  if(min(x) < 0){
    stop("x < 0:  values of x must be non-negative integers")
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

  if(min(obs) < 0){
    stop("All observations must be non-negative")
  }

  a = alpha;
  b = beta;
  sobs = sum(obs);
  n = length(obs);

  f_x = stats::dnbinom(x,size = a + sobs, mu = (a + sobs)/(b + n));

  # checking with formula
  #num1 = lgamma(alpha + sobs + ytilde);
  #den1 = lgamma(alpha + sobs) + lgamma(ytilde+1);
  #fact1 = exp(num1 - den1);

  #fact2 = exp((alpha + sobs)*log((beta + n)/(beta+n+1)));

  #fact3 = exp(ytilde*log(1/(beta+n+1)));

  #f_x2 = fact1*fact2*fact3;

  return(f_x);

}

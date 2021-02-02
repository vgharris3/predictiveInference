#' Title
#'
#' @param n desired random sample size
#' @param XN (x1,...,xN), xi ~ exp(theta), theta ~ Gamma(dt, gm)
#' @param d < n s.t. (x1,...,xd) are fully observed and (xd+1, ..., xN) are censored
#' @param dt Gamma shape parameter for distribution of theta
#' @param gm Gamma rate parameter for distribution of theta
#'
#' @return random sample of size n from the Exponential-Gamma predictive probability distribution
#' @export
#'
#' @examples 1
rpredEG = function(n, XN, d, dt, gm){
  a = d + dt
  b = gm + sum(XN)
  theta = stats::rgamma(1,shape=a,rate=b)

  return(stats::rexp(n,theta))
}

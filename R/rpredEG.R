#' The Exponential-Gamma Predictive Distribution
#'
#' @param S desired random sample size
#' @param y list of observed and censored survival times (y1,...,yN), yi ~ Exp(theta), theta ~ Gamma(dt, gm)
#' @param c censoring indicator (vector of 0s and 1s for censored and fully observed, respectively)
#' @param dt Gamma shape parameter for distribution of theta
#' @param gm Gamma rate parameter for distribution of theta
#'
#' @return random sample of size S from the Exponential-Gamma predictive probability distribution
#' @export
#'
#' @examples 1
rpredEG = function(S, y, c, dt, gm){


  d = sum(c) # number of fully observed copies

  a = d + dt
  b = gm + sum(y)
  # theta = stats::rgamma(1,shape=a,rate=b)
  theta = stats::rgamma(S,shape=a,rate=b)

  return(stats::rexp(S,theta))
}

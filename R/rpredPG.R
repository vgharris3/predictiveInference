#' Title
#'
#' rpredPG returns a random sample of size n from the Poisson-Gamma predictive probability distribution given observations obs~Poi(theta)
#' and theta~Gamma(alpha,beta)
#'
#' @param n desired random sample size
#' @param obs vector of (observed) Poisson-distributed counts
#' @param x maximum integer for which predictive probability is desired
#' @param alpha sum of counts from beta prior observations for gamma prior distribution on theta
#' @param beta number of prior observations for gamma prior distribution on theta
#'
#' @return Poisson-Gamma cumulative predictive probability
#' @export
#'
#' @examples 1
rpredPG = function(n,obs, x, alpha = 1, beta=1){

  #ERROR HANDLING

  ##########################
  ##########################
  #NOT ENOUGH / TOO FEW PARAMETERS
  ##########################
  ##########################

  if(min(x) < 0){
    stop("x < 0:  x must be a non-negative integer")
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

  F_x = ppredPG(obs, x, alpha, beta)

  u = stats::runif(n)

  rank_list = numeric(n)

  for(i in 1:n) {
    rankF = which(abs(F_x - u[i]) == min(abs(F_x - u[i])))

    if(F_x[rankF] > u[i]){ rankF = rankF - 1 }

    rank_list[i] = rankF
  }

  return(rank_list)
}

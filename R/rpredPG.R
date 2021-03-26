#' Title
#'
#' rpredPG returns a random sample of size n from the Poisson-Gamma predictive probability distribution given observations obs~Poi(theta)
#' and theta~Gamma(alpha,beta)
#'
#' @param n desired random sample size
#' @param obs vector of (observed) Poisson-distributed counts
#' @param xmax maximum integer for which predictive probability is desired
#' @param alpha sum of counts from beta prior observations for gamma prior distribution on theta
#' @param beta number of prior observations for gamma prior distribution on theta
#'
#' @return Poisson-Gamma cumulative predictive probability
#' @export
#'
#' @examples 1
rpredPG = function(n,obs, xmax, alpha = 1, beta=1){

  #ERROR HANDLING

  ##########################
  ##########################
  #NOT ENOUGH / TOO FEW PARAMETERS
  ##########################
  ##########################

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

  F_x = ppredPG(obs, xmax, alpha, beta)

  u = stats::runif(n)

  rank_list = numeric(n)

  for(i in 1:n) { rank_list[i] = which(abs(F_x - u[i]) == min(abs(F_x - u[i]))) }

  for(i in 1:n) {
    rank = which(abs(F_x - u[i]) == min(abs(F_x - u[i])))

    if(F_x[rank] > u[i]){ rank = rank - 1 }

    rank_list[i] = rank
  }

  return(rank_list)
}

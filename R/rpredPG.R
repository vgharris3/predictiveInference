#' Title
#'
#' rpredPG returns a random sample of size n from the Poisson-Gamma predictive probability distribution given observations obs~Poi(theta)
#' and theta~Gamma(alpha,beta)
#'
#' @param n desired random sample size
#' @param obs vector of (observed) Poisson-distributed counts
#' @param alpha sum of counts from beta prior observations for gamma prior distribution on theta
#' @param beta number of prior observations for gamma prior distribution on theta
#'
#' @return Poisson-Gamma cumulative predictive probability
#' @export
#'
#' @examples 1
rpredPG = function(n,obs, alpha = 1, beta=1){

  #ERROR HANDLING

  ##########################
  ##########################
  #NOT ENOUGH / TOO FEW PARAMETERS
  ##########################
  ##########################

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

  #Need to establish upper bound of support for sampling (lower bound is 0)
  #Modified Bisection Method (rounding to integers)
    #Computed midpoint values will be fed into the function f_x(midpoint) - eps, which will feed midpoint into dnbinom, which requires integers
  #First reach to the right from the expected value E_x, until f_x(U) - eps < 0
  #Then employ Bisection Method, starting with E_x and U as the lower and upper ends of the first interval
  #Along the way, always round the middle value to the nearest integer

  #Set desired tolerance epsilon for root function
  eps = 0.0000001
  eps = .Machine$double.eps^0.5
  #Compute predictive expected value of x
  E_x = round((alpha + sum(obs))/(beta + length(obs)))
  Lower = E_x
  #Reach right to find upper bound of starting interval for Bisection Method
  fLower = dpredPG(E_x,obs,alpha,beta) - eps
  if(fLower <= 0){ stop("Density at Expected Value computed to be less than epsilon.")}

  fUpper = fLower
  reachExp = 0
  while(fUpper > 0){
    reachExp = reachExp + 1
    Upper = E_x + 10^reachExp
    fUpper = dpredPG(Upper,obs,alpha,beta) - eps
  }

  limit_found = 0
  right_end = E_x

  while(!limit_found){
    #Establish new interval
    Mid = round(mean(c(Lower,Upper)))
    fMid = dpredPG(Mid,obs,alpha,beta) - eps
    if(fMid > 0){
      Lower = Mid #Upper is still Upper
      if(dpredPG(Mid+1,obs,alpha,beta) - eps <= 0){ #if f(Mid) and f(Mid+1) are straddling 0
        limit_found = 1
        right_end = Mid + 1
      }
    } else { #fMid <= 0
      Upper = Mid #Lower is still Lower
      if(dpredPG(Mid-1,obs,alpha,beta) - eps > 0){ #if f(Mid) and f(Mid-1) are straddling 0
        limit_found = 1
        right_end = Mid
      }
    }
  }

  x = 0:right_end

  F_x = ppredPG(x, obs, alpha, beta)

  u = stats::runif(n)

  rank_list = numeric(n)

  for(i in 1:n) {
    rankF = which(abs(F_x - u[i]) == min(abs(F_x - u[i])))

    if(F_x[rankF] > u[i]){ rankF = rankF - 1 }

    rank_list[i] = rankF
  }

  return(rank_list)
}

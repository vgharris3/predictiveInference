#' Title
#'
#' rpredBB returns a random sample of size n from the Beta-Binomial predictive probability distribution computed from s observed successes out of a sample of size N, M future observations, and user input shape parameters for Beta prior on Pr(success)
#'
#' @param n desired random sample size
#' @param N number of observed independent binary variables
#' @param s number of observed successes out of N
#' @param M number of future observations for which prediction is desired (note x <= M)
#' @param alpha shape factors for prior Beta distribution on theta
#' @param beta rate factor for prior Beta distribution on theta
#'
#' @return random sample of size n from the Beta-Binomial predictive probability distribution
#' @export
#'
#' @examples 1
rpredBB = function(n, N, s, M, alpha = 1, beta = 1){
  #rpredBB returns a random sample of size n from the Beta-Binomial predictive probability distribution computed from s observed
  #successes out of a sample of size N, M future observations, and user input shape parameters for Beta prior on Pr(success)

  #r+pred+BB
  #r = random sample (like R distribution functions)
  #pred = prediction
  #BB = BetaBinomial
  #FF = Future Fraction

  #ERROR HANDLING

  ##########################
  ##########################
  #NOT ENOUGH / TOO FEW PARAMETERS
  ##########################
  ##########################

  if(s > N){
    stop("s > N:  The number of successes (s) cannot exceed the number of observations (N)")
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

  #############################################################
  #############################################################
  ##CHANGE THIS TO A CALL TO ppredBB() TO OBTAIN F_x()

  #r[] is the number of successes in the M future observations

  r = 1:M

  numerator = lgamma(M+1) + lgamma(N+alpha+beta) + lgamma(r+s+alpha) + lgamma(M+N-r-s+beta)
  denominator = lgamma(r+1) + lgamma(M-r+1) + lgamma(alpha+s) + lgamma(N-s+beta) + lgamma(M+N+alpha+beta)
  f_x = exp(numerator - denominator)

  F_x = cumsum(f_x)
  #############################################################
  #############################################################

  F_x = ppredBB(1:M,N,s,M,alpha,beta)

  u = stats::runif(n)

  rank_list = numeric(n)

  for(i in 1:n) {
    rankF = which(abs(F_x - u[i]) == min(abs(F_x - u[i])))

    if(F_x[rankF] > u[i]){ rankF = rankF - 1 }

    rank_list[i] = rankF
  }

  return(rank_list)
}

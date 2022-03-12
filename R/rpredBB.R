#' The Beta-Binomial Predictive Distribution
#'
#' rpredBB returns a random sample of size n from the Beta-Binomial predictive probability distribution computed
#' from s observed successes out of a sample of size N, M future observations, and user input shape parameters for Beta prior on Pr(success)
#'
#' @param S desired random sample size
#' @param N number of observed independent binary variables
#' @param t number of observed successes out of N observations
#' @param M number of future observations for which prediction is desired
#' @param alpha shape factor for prior Beta distribution on theta
#' @param beta rate factor for prior Beta distribution on theta
#'
#' @return random sample of size n from the Beta-Binomial predictive probability distribution
#' @export
#'
#' @examples 1
rpredBB = function(S, N, t, M, alpha = 1, beta = 1){
  #rpredBB returns a random sample of size N from the Beta-Binomial predictive probability distribution computed from t observed
  #successes out of a sample of size N, M future observations, and user input shape parameters for Beta prior on Pr(success)

  #The inverse transform method is used to generate the random sample

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

  if(t > N){
    stop("t > N:  The number of observed successes (t) cannot exceed the number of observations (N)")
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


  F_x = ppredBB(0:M,N,t,M,alpha,beta)

  u = stats::runif(S)

  rank_list = numeric(S)

  for(i in 1:S) {
    rankF = which(abs(F_x - u[i]) == min(abs(F_x - u[i])))

    if(F_x[rankF] > u[i]){ rankF = rankF - 1 }

    rank_list[i] = rankF
  }

  return(rank_list)
}

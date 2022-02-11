#' The Beta-Binomial Predictive Distribution
#'
#' ppredBB returns the predictive probability of <= x future successes out of M trials, given s observed successes out of a sample of
#size N, and user input shape parameters for Beta prior on Pr(success)
#'
#' @param tpred vector of integers: number(s) of successes for which prediction is desired
#' @param N number of observed independent binary variables
#' @param t number of observed successes out of N
#' @param M number of future observations for which prediction is desired (note x <= M)
#' @param alpha shape factors for prior Beta distribution on theta
#' @param beta rate factor for prior Beta distribution on theta
#' @param lower logical; if TRUE (default), probabilities are P(R<=q) otherwise P(R>q)
#'
#' @return Beta-Binomial cumulative predictive probability
#' @export
#'
#' @examples 1
ppredBB = function(tpred, N, t, M, alpha = 1, beta=1, lower=TRUE){

  #ppredBB returns the predictive probability of <= x future successes out of M trials, given s observed successes out of a sample of
  #size N, and user input shape parameters for Beta prior on Pr(success)

  #p+pred+BB+FF
  #p = cdf (like R distribution functions)
  #pred = prediction
  #BB = BetaBinomial
  #FF = Future Fraction

  #tpred is a vector of integers, for which the function will compute P(X<=x)
  #N is the number of observed independent binary variables
  #s is the number of observed successes out of N
  #M is the number of future observations for which prediction is desired (note x <= M)

  #alpha, beta:  shape factors for prior Beta distribution on theta, default is (1,1) (equiv. to Uniform dist.)
  #Note:  theta = Pr(success) for the N binary variables.  theta ~ Beta(alpha,beta).  User inputs the shape parameters for the Beta prior.
  #Interpretation:  alpha + beta = sample size, alpha/(alpha+beta) is fraction of sample that are successes
  #lower:  logical; if TRUE (default), probabilities are P[R<=q] otherwise P[R>q]

  #ERROR HANDLING

  if(t > N){
    stop("t > N:  The number of successes (t) cannot exceed the number of observations (N)")
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

  if(min(tpred) < 0){
    stop("min(x) < 0")
  }

  if(max(tpred) > M){
    stop("max(tpred) > M:  The number of successes (x) for which P(Tpred<=tpred) will be computed cannot exceed the number of future observations (M)")
    return (1)
  }

  #r[] is the number of successes in the M future observations

  f_t = dpredBB(0:max(tpred), N, t, M, alpha, beta)

  F_t = cumsum(f_t)

  return_index = which(tpred == tpred)

  if(lower){return(F_t[return_index])}
  else {return(1 - F_t[return_index])}


}

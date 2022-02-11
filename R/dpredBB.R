#' The Beta-Binomial Predictive Distribution
#'
#' dpredBB returns the predictive probability of x future successes out of M trials, given s observed successes out of a sample of size
#' N, and user input shape parameters for Beta prior on Pr(success)
#'
#' @param tpred vector of integers: number(s) of successes for which prediction is desired
#' @param N number of observed independent binary variables
#' @param t number of observed successes out of N observations
#' @param M number of future observations (note tpred <= M)
#' @param alpha first shape factor for prior Beta distribution on theta
#' @param beta second shape factor for prior Beta distribution on theta
#'
#' @return Beta-Binomial predictive probability
#' @export
#'
#' @examples 1
dpredBB <- function(tpred, N, t, M, alpha = 1, beta=1){

  #dpredBB returns the predictive probability of x future successes out of M trials, given s observed successes out of a sample of
  #size N, and user input shape parameters for Beta prior on Pr(success)

  #p+pred+BB
  #p = cdf (like R distribution functions)
  #pred = prediction
  #BB = BetaBinomial

  #tpred is a vector of integers, for which the function will compute P(X=x)
  #N is the number of observed independent binary variables
  #t is the number of observed successes out of N
  #M is the number of future observations for which prediction is desired (note x <= M)

  #alpha, beta:  shape factors for prior Beta distribution on theta, default is (1,1) (equiv. to Uniform dist.)
  #Note:  theta = Pr(success) for the N binary variables.  theta ~ Beta(alpha,beta).  User inputs the shape parameters for the Beta prior.
  #Interpretation:  alpha + beta = sample size, alpha/(alpha+beta) is fraction of sample that are successes
  #lower:  logical; if TRUE (default), probabilities are P[R<=q] otherwise P[R>q]

  #ERROR HANDLING

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

  if(min(tpred) < 0){
    stop("min(y) < 0")
  }

  if(max(tpred) > M){
    stop("max(y) > M:  The number of observations (tpred) for which P(Tpred<=tpred) will be computed cannot exceed the number of future observations (M)")
    return (1)
  }

  #r[] is the number of successes in the M future observations

  r = tpred;

  numerator = lgamma(M+1) + lgamma(N+alpha+beta) + lgamma(r+t+alpha) + lgamma(M+N-r-t+beta);
  denominator = lgamma(r+1) + lgamma(M-r+1) + lgamma(alpha+t) + lgamma(N-t+beta) + lgamma(M+N+alpha+beta);
  f_t = exp(numerator - denominator)

  return(f_t)

}

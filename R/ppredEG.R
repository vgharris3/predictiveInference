#' The Exponential-Gamma Predictive Distribution
#'
#' @param ypred survival time(s) for which predictive density is desired
#' @param y ordered list of observed survival times and censored survival times (y1,...,yN), yi ~ exp(theta), theta ~ Gamma(dt, gm)
#' @param c censoring indicator (vector of 0s and 1s for censored and fully observed, respectively)
#' @param dt Gamma shape parameter for distribution of theta
#' @param gm Gamma rate parameter for distribution of theta
#'
#' @return Exponential-Gamma cumulative predictive probability
#' @export
#'
#' @examples 1
ppredEG = function(ypred, y, c, dt, gm){

  #ppredEG returns the cumulative predictive probability of surviving to [OR RATHER DYING BY] any time up to max(x) (i.e. Pr(X <= x)),
    #given d observed and N-d censored copies out of XN total events

  #p+pred+GE
  #p = cdf (like R distribution functions)
  #pred = prediction
  #GE = Gamma Exponential

  #x is a time past which dpredEG will return the predictive probability of survival Pr(X = x)
  #XN = (x1,...,xN), xi ~ exp(theta) where theta ~ Gamma(dt, gm)
  #d < n s.t. (x1,...,xd) are fully observed and (xd+1, ..., xN) are censored

  #ERROR HANDLING

  if(sum(c) > length(y)){
    stop("sum(c) > length(y):  The number of observed copies cannot exceed the total number of copies")
    return (1)
  }

  if(dt <= 0){
    stop("dt<=0:  dt must be greater than 0")
    return(1)
  }

  if(gm <= 0){
    stop("gm <= 0: gm must be greater than 0")
    return(1)
  }


  d = sum(c) # number of fully observed copies

  Freturn = numeric(length(ypred))

  for (i in 1:length(ypred)){
    Freturn[i] = stats::integrate(dpredEG,lower = 0,upper = ypred[i],y = y, c = c, dt = dt, gm = gm)
  }

  return(Freturn)

}

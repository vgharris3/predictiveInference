#' The Exponential-Gamma Predictive Distribution
#'
#' @param x survival time
#' @param XN (x1,...,xN), xi ~ exp(theta), theta ~ Gamma(dt, gm)
#' @param d < n s.t. (x1,...,xd) are fully observed and (xd+1, ..., xN) are censored
#' @param dt Gamma shape parameter for distribution of theta
#' @param gm Gamma rate parameter for distribution of theta
#'
#' @return Exponential-Gamma cumulative predictive probability
#' @export
#'
#' @examples 1
ppredEG = function(x, XN, d, dt, gm){

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

  if(d > length(XN)){
    stop("d > length(XN):  The number of observed copies cannot exceed the total number of copies")
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

#  if(min(x) <= 0){
#    stop("x_i < - for some i:  x_i must be greater than or equal to zero for all i")
#  }


  Freturn = numeric(length(x))

  for (i in 1:length(x)){
    Freturn[i] = stats::integrate(dpredEG,lower = 0,upper = x[i],XN = XN, d = d, dt = dt, gm = gm)
  }

  return(Freturn)

}

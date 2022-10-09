#' Random Generation for Truncated Univariate Normal 
#'
#' Draws from truncated univariate normal distribution within an interval.
#'  
#' @param mean vector of means. The length is the number of observations. 
#' @param sd standard deviation. Defaults to 1.
#' @param lower a scalar of lower bound for truncation, or a vector of 
#'   lower bounds with the same length as \code{mean}. 
#' @param upper a scalar of upper bound for truncation, or a vector of 
#'   upper bounds with the same length as \code{mean}. 
#' @param n number of random samples when \code{mean} is a scalar.
#'
#' @return Returns a vector of random numbers following the specified 
#'   truncated univariate normal distribution. 
#'
#' @examples
#' set.seed(1203)
#' x = rtnorm(mean=rep(1,1000), sd=2, lower=-2, upper=3)
#' summary(x)
#' 
#' # use the alternative form of input
#' set.seed(1203)
#' x = rtnorm(mean=1, sd=2, lower=-2, upper=3, n=1000)
#' summary(x)
#' 
#' @export
rtnorm <- function(mean, sd=1, lower, upper, n=NULL) {
  if (is.null(n)) {
    n = length(mean);
  } else if (n >= 1) {
    if (length(mean) > 1 || length(lower) > 1 || length(upper) > 1) {
      stop("mean, lower, and upper must be a scalar when n >= 1");
    }
    n = round(n);
    mean = rep(mean, n);
  } else {
    stop("Invalid input of n");
  }
  
  n1 = length(lower); 
  n2 = length(upper); 
  if (!(n==n1 || n1==1)) {
    stop("lower must have the same length as mean");
  }
  if (!(n1==n2)) {
    stop("lower and upper must have the same length");
  }
  
  as.vector(rtnormcpp(mean, sd, lower, upper))
}


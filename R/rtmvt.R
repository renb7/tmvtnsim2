#' Random Generation for Truncated Multivariate t 
#'
#' Draws from truncated multivariate t distribution subject to 
#' linear inequality constraints represented by a matrix. 
#' 
#' @param mean \code{n x p} matrix of means. The number of rows is the number 
#'   of observations. The number of columns is the dimension of the problem.
#' @param sigma \code{p x p} covariance matrix. 
#' @param nu degrees of freedom for Student-t distribution.
#' @param blc \code{m x p} matrix of coefficients for linear inequality 
#'   constraints. If \code{NULL}, the \code{p x p} identity matrix will be used.
#' @param lower \code{n x m} or \code{1 x m} matrix of lower bounds for 
#'   truncation. 
#' @param upper \code{n x m} or \code{1 x m} matrix of upper bounds for 
#'   truncation. 
#' @param init \code{n x p} or \code{1 x p} matrix of initial values. 
#'   If \code{NULL}, default initial values will be generated. 
#' @param burn number of burn-in iterations. Defaults to 10.
#' @param n number of random samples when \code{mean} is a vector.
#' 
#' @return Returns an \code{n x p} matrix of random numbers following the 
#'   specified truncated multivariate t distribution. 
#'
#' @examples
#' # Example 1: full rank blc
#' d = 3;
#' rho = 0.5;
#' sigma = matrix(0, d, d);
#' sigma = rho^abs(row(sigma) - col(sigma));
#' nu = 10;
#' blc = diag(1,d);
#' n = 1000;
#' mean = matrix(rep(1:d,n), nrow=n, ncol=d, byrow=TRUE);
#' set.seed(1203)
#' result = rtmvt(mean, sigma, nu, blc, -1, 1, burn=50)
#' apply(result, 2, summary)
#'
#' # Example 2: use the alternative form of input
#' set.seed(1203)
#' result = rtmvt(mean=1:d, sigma, nu, blc, -1, 1, burn=50, n)
#' apply(result, 2, summary)
#'
#' # Example 3: non-full rank blc, different means
#' d = 3;
#' rho = 0.5;
#' sigma = matrix(0, d, d);
#' sigma = rho^abs(row(sigma) - col(sigma));
#' nu = 10;
#' blc = matrix(c(1,0,1,1,1,0), nrow=d-1, ncol=d, byrow=TRUE)
#' n = 100;
#' set.seed(3084)
#' mean = matrix(runif(n*d), nrow=n, ncol=d);
#' result = rtmvt(mean, sigma, nu, blc, -1, 1, burn=50)
#' apply(result, 2, summary)
#'
#' @export
rtmvt <- function(mean, sigma, nu, blc=NULL, lower, 
                  upper, init=NULL, burn=10, n=NULL) {
  if (is.matrix(sigma)) {
    p = ncol(sigma);
  } else {
    p = 1;
    sigma = as.matrix(sigma);
  }
  
  mean = matrix(mean, ncol=p);
  
  if (is.null(blc)) {
    blc = diag(p);
  } else {
    blc = matrix(blc, ncol=p);
  }
  
  m = nrow(blc);

  lower = matrix(lower, ncol=m);
  upper = matrix(upper, ncol=m);
  
  
  if (is.null(n)) {
    n = nrow(mean);
  } else if (n >= 1) {
    if (nrow(mean) > 1 || nrow(lower) > 1 || nrow(upper) > 1) {
      stop("mean, lower, and upper must be a vector when n >= 1");
    }
    n = round(n);
    mean = matrix(mean, nrow=n, ncol=p, byrow=TRUE);
  } else {
    stop("Invalid input of n");
  }
  
  
  n1 = nrow(lower); 
  n2 = nrow(upper); 
  if (!(n==n1 || n1==1)) {
    stop("lower must be a row vector or have the same number of rows as mean");
  }
  if (!(n1==n2)) {
    stop("lower and upper must have the same number of rows");
  }
  
  if (is.null(init)) {
    init = matrix(0, n1, p);
  } else {
    init = matrix(init, ncol=p);
  }
  
  n3 = nrow(init);
  if (!(n1==n3)) {
    stop("init must have the same number of rows as lower and upper");
  }
  
  rtmvtcpp(mean, sigma, nu, blc, lower, upper, init, burn);
}


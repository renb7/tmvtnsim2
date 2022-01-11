#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' Random Generation for Truncated Multivariate Normal 
//'
//' Draws from truncated multivariate normal distribution subject to 
//' linear inequality constraints represented by a matrix. 
//' 
//' @param mean \code{n x p} matrix of means. The number of rows is the number 
//'   of observations. The number of columns is the dimension of the problem.
//' @param sigma \code{p x p} covariance matrix. 
//' @param blc \code{m x p} matrix of coefficients for linear inequality 
//'   constraints.
//' @param lower \code{n x m} or \code{1 x m} matrix of lower bounds for 
//'   truncation. 
//' @param upper \code{n x m} or \code{1 x m} matrix of upper bounds for 
//'   truncation. 
//' @param init \code{n x p} or \code{1 x p} matrix of initial values.
//' @param burn Number of burn-in iterations. Defaults to 10.
//' 
//' @return Returns a \code{n x p} matrix of random numbers following the 
//'   specified truncated multivariate normal distribution. 
//'
//' @examples
//' # Example 1: full rank
//' d = 3;
//' rho = 0.9;
//' sigma = matrix(0, d, d);
//' sigma = rho^abs(row(sigma) - col(sigma));
//' blc = diag(1,d);
//' n = 1000;
//' mean = matrix(rep(1:d,n), nrow=n, ncol=d, byrow=TRUE);
//' lower = matrix(rep(-1,d), nrow=1);
//' upper = matrix(rep(1,d), nrow=1);
//' init = matrix(rep(0,d), nrow=1);
//' set.seed(1203)
//' result = rtmvnorm(mean, sigma, blc, lower, upper, init, burn=50)
//' apply(result, 2, summary)
//' 
//' # Example 2: non-full rank, invalid initial values
//' d = 3;
//' rho = 0.5;
//' sigma = matrix(0, d, d);
//' sigma = rho^abs(row(sigma) - col(sigma));
//' blc = matrix(c(1,1,1,0,1,0,1,0,1),ncol=d);
//' n = 100;
//' mean = matrix(rep(1:d,n), nrow=n, ncol=d, byrow=TRUE);
//' lower = matrix(rep(-1,d), nrow=1);
//' upper = matrix(rep(1,d), nrow=1);
//' init = matrix(rep(0.8,d), nrow=1);
//' set.seed(1228)
//' result = rtmvnorm(mean, sigma, blc, lower, upper, init, burn=50)
//' apply(result, 2, summary)
//' 
//' # Example 3: means, lower, or upper bounds differ across samples
//' d = 3;
//' rho = 0.5;
//' sigma = matrix(0, d, d);
//' sigma = rho^abs(row(sigma) - col(sigma));
//' blc = matrix(c(1,0,1,1,1,0),nrow=d-1,ncol=d,byrow=TRUE)
//' n = 100;
//' set.seed(3084)
//' mean = matrix(runif(n*d), nrow=n, ncol=d);
//' lower = matrix(rep(-1,d-1), nrow=1);
//' upper = matrix(rep(1,d-1), nrow=1);
//' init = matrix(rep(0.8,d), nrow=1);
//' result = rtmvnorm(mean, sigma, blc, lower, upper, init, burn=50)
//' apply(result, 2, summary)
//' 
//' @export
// [[Rcpp::export]]
arma::mat rtmvnorm(const arma::mat& mean, 
                   const arma::mat& sigma, 
                   const arma::mat& blc,
                   const arma::mat& lower, 
                   const arma::mat& upper,
                   arma::mat& init,
                   const arma::uword burn = 10) {
  const unsigned int n=mean.n_rows, p=mean.n_cols;
  arma::mat x(n,p); // output samples
  
  // draw from truncated univariate normal
  if (p==1) {
    if (blc(0,0) > 0.0) { 
      x.col(0) = rtnorm(mean.col(0), sqrt(sigma(0,0)), 
            lower.col(0)/blc(0,0), 
            upper.col(0)/blc(0,0));
    } else if (blc(0,0) < 0.0) {
      x.col(0) = rtnorm(mean.col(0), sqrt(sigma(0,0)), 
            upper.col(0)/blc(0,0), 
            lower.col(0)/blc(0,0));
    } else {
      arma::vec lower1 {R_NegInf};
      arma::vec upper1 {R_PosInf};
      x.col(0) = rtnorm(mean.col(0), sqrt(sigma(0,0)), 
            lower1, upper1);
    }
    return x;
  }
  
  // check dimensions
  const unsigned int n1=lower.n_rows, n2=upper.n_rows, n3=init.n_rows;
  const unsigned int m=blc.n_rows, m1=lower.n_cols, m2=upper.n_cols;
  
  if (!(n==n1 || n1==1)) {
    Rcpp::stop("lower must be a row vector or have the same number of rows as mean");
  }
  
  if (!(n1==n2 && n1==n3)) {
    Rcpp::stop("lower, upper, and init must have the same number of rows");
  }
  
  if (!(m==m1)) {
    Rcpp::stop("lower must have as many columns as the number of constraints");
  }
  
  if (!(m1==m2)) {
    Rcpp::stop("lower and upper must have the same number of columns");
  }
  
  if (!(blc.n_cols==p)) {
    Rcpp::stop("blc must have the same number of columns as mean");
  }
  
  if (!(init.n_cols==p)) {
    Rcpp::stop("init must have the same number of columns as mean");
  }
  
  
  // check boundary conditions
  if (sum(any(lower >= upper)) > 0) {
    Rcpp::stop("lower bound must be smaller than upper bound");
  }
  
  // check initial values, and generate automatically if needed
  arma::mat blct = blc.t();
  arma::mat initc = init*blct;
  if (sum(all(initc >= lower + 1e-8 && initc <= upper - 1e-8)) < m) {
    arma::mat estimate(n1,m);
    for (arma::uword i=0; i<n1; i++) {
      for (arma::uword j=0; j<m; j++) {
        if (std::isinf(lower(i,j)) && std::isinf(upper(i,j))) {
          estimate(i,j) = 0;  
        } else if (std::isinf(lower(i,j))) {
          estimate(i,j) = upper(i,j) - 1e-8;  
        } else if (std::isinf(upper(i,j))) {
          estimate(i,j) = lower(i,j) + 1e-8;  
        } else {
          estimate(i,j) = 0.5*(lower(i,j) + upper(i,j));  
        }
      }
    }
    init = estimate * arma::pinv(blct);
  }
  
  // check whether a matrix has identical rows
  auto f = [](const arma::mat& y) {
    const unsigned int n=y.n_rows, p=y.n_cols; 
    bool identical=1;
    if (n>1) {
      for (arma::uword i=1; i<n; i++) {
        for (arma::uword j=0; j<p; j++) {
          if (y(i,j) != y(i-1,j)) {
            identical = 0;
            break;
          }
        }
        if (identical==0) break;
      }
    }
    return identical;
  };
  
  
  // transform to the problem with identity covariance matrix
  arma::mat cholt = arma::trans(arma::chol(sigma));
  arma::mat R = blc*cholt;
  arma::uvec js(p);
  std::iota(js.begin(), js.end(), 0);
  Rcpp::NumericVector mu {0}; 
  
  // Gibbs step for sampling truncated multivariate normal
  auto g = [p, m, R, js, mu](arma::vec a, arma::vec b, arma::vec& z) {
    for (arma::uword j=0; j<p; j++) {
      // set up linear inequality constraints for z(j)
      arma::vec rj = R.col(j);
      arma::uvec j2 = arma::find(js != j);
      arma::mat Rj = R.cols(j2);
      arma::vec zj = z(j2);
      arma::vec atemp = a - Rj*zj;
      arma::vec btemp = b - Rj*zj;
      
      // determine lower and upper bounds for z(j)
      double lowerj=R_NegInf, upperj=R_PosInf;
      for (arma::uword k=0; k<m; k++) {
        if (rj(k) != 0) {
          double ak=atemp(k)/rj(k), bk=btemp(k)/rj(k);
          if (rj(k) > 0) {
            if (ak > lowerj) lowerj = ak;
            if (bk < upperj) upperj = bk;
          } else {
            if (bk > lowerj) lowerj = bk;
            if (ak < upperj) upperj = ak;
          }
        }
      }
      
      // generate z(j) for truncated univariate normal
      arma::vec lowerj1(1); lowerj1(0) = lowerj;
      arma::vec upperj1(1); upperj1(0) = upperj;
      z(j) = rtnorm(mu, 1, lowerj1, upperj1)(0);
    }
  };
  
  
  arma::vec mean1(p), lower1(m), upper1(m), init1(p);
  arma::vec a(m), b(m), z(p);
  
  // obtain burn + n samples in case of identical means, lower and upper bounds 
  if (f(mean) && f(lower) && f(upper)) {
    mean1 = arma::trans(mean.row(0));
    lower1 = arma::trans(lower.row(0));
    upper1 = arma::trans(upper.row(0));
    init1 = arma::trans(init.row(0));
    
    a = lower1 - blc*mean1;
    b = upper1 - blc*mean1;
    z = solve(cholt, init1-mean1);
    
    for (arma::uword i=0; i<burn+n; i++) {
      g(a, b, z);
      if (i>=burn) x.row(i-burn) = arma::trans(cholt*z + mean1);
    }
    return x;
  }
  
  // obtain (burn+1)*n samples for non-identical means, lower, or upper bounds
  if (n1==1) {
    lower1 = arma::trans(lower.row(0));
    upper1 = arma::trans(upper.row(0));
    init1 = arma::trans(init.row(0));
  }
  
  for (arma::uword i=0; i<n; i++) {
    mean1 = arma::trans(mean.row(i));
    if (n1==n) {
      lower1 = arma::trans(lower.row(i));
      upper1 = arma::trans(upper.row(i));
      init1 = arma::trans(init.row(i));
    }
    
    a = lower1 - blc*mean1;
    b = upper1 - blc*mean1;
    z = solve(cholt, init1-mean1);
    
    for (arma::uword i2=0; i2<burn+1; i2++) {
      g(a, b, z);
    }
    x.row(i) = arma::trans(cholt*z + mean1);
  }
  return x;
}


//' Random Generation for Truncated Multivariate t 
//'
//' Draws from truncated multivariate t distribution subject to 
//' linear inequality constraints represented by a matrix. 
//' 
//' @param mean \code{n x p} matrix of means. The number of rows is the number 
//'   of observations. The number of columns is the dimension of the problem.
//' @param sigma \code{p x p} covariance matrix. 
//' @param nu degrees of freedom for Student-t distribution.
//' @param blc \code{m x p} matrix of coefficients for linear inequality 
//'   constraints.
//' @param lower \code{n x m} or \code{1 x m} matrix of lower bounds for 
//'   truncation. 
//' @param upper \code{n x m} or \code{1 x m} matrix of upper bounds for 
//'   truncation. 
//' @param init \code{n x p} or \code{1 x p} matrix of initial values.
//' @param burn Number of burn-in iterations. Defaults to 10.
//' 
//' @return Returns a \code{n x p} matrix of random numbers following the 
//'   specified truncated multivariate t distribution. 
//'
//' @examples
//' # Example 1: full rank
//' d = 3;
//' rho = 0.5;
//' nu = 10;
//' sigma = matrix(0, d, d);
//' sigma = rho^abs(row(sigma) - col(sigma));
//' blc = diag(1,d);
//' n = 1000;
//' mean = matrix(rep(1:d,n), nrow=n, ncol=d, byrow=TRUE);
//' lower = matrix(rep(-1,d), nrow=1);
//' upper = matrix(rep(1,d), nrow=1);
//' init = matrix(rep(0.8,d), nrow=1);
//' set.seed(1203)
//' result = rtmvt(mean, sigma, nu, blc, lower, upper, init, burn=50)
//' apply(result, 2, summary)
//' 
//' # Example 2: non-full rank, different means
//' d = 3;
//' rho = 0.5;
//' sigma = matrix(0, d, d);
//' sigma = rho^abs(row(sigma) - col(sigma));
//' nu = 10;
//' blc = matrix(c(1,0,1,1,1,0),nrow=d-1,ncol=d,byrow=TRUE)
//' n = 100;
//' set.seed(3084)
//' mean = matrix(runif(n*d), nrow=n, ncol=d);
//' lower = matrix(rep(-1,d-1), nrow=1);
//' upper = matrix(rep(1,d-1), nrow=1);
//' init = matrix(rep(0.8,d), nrow=1);
//' result = rtmvt(mean, sigma, nu, blc, lower, upper, init, burn=50)
//' apply(result, 2, summary)
//'
//' @export
// [[Rcpp::export]]
arma::mat rtmvt(const arma::mat& mean, 
                const arma::mat& sigma, 
                const double nu,
                const arma::mat& blc,
                const arma::mat& lower, 
                const arma::mat& upper,
                arma::mat& init,
                const arma::uword burn = 10) {
  const unsigned int n=mean.n_rows, p=mean.n_cols;
  arma::mat x(n,p); // output samples
  
  // draw from truncated univariate normal
  if (p==1) {
    if (blc(0,0) > 0.0) { 
      x.col(0) = rtnorm(mean.col(0), sqrt(sigma(0,0)), 
            lower.col(0)/blc(0,0), 
            upper.col(0)/blc(0,0));
    } else if (blc(0,0) < 0.0) {
      x.col(0) = rtnorm(mean.col(0), sqrt(sigma(0,0)), 
            upper.col(0)/blc(0,0), 
            lower.col(0)/blc(0,0));
    } else {
      arma::vec lower1 {R_NegInf};
      arma::vec upper1 {R_PosInf};
      x.col(0) = rtnorm(mean.col(0), sqrt(sigma(0,0)), 
            lower1, upper1);
    }
    return x;
  }
  
  // check dimensions
  const unsigned int n1=lower.n_rows, n2=upper.n_rows, n3=init.n_rows;
  const unsigned int m=blc.n_rows, m1=lower.n_cols, m2=upper.n_cols;
  
  if (!(n==n1 || n1==1)) {
    Rcpp::stop("lower must be a row vector or have the same number of rows as mean");
  }
  
  if (!(n1==n2 && n1==n3)) {
    Rcpp::stop("lower, upper, and init must have the same number of rows");
  }
  
  if (!(m==m1)) {
    Rcpp::stop("lower must have as many columns as the number of constraints");
  }
  
  if (!(m1==m2)) {
    Rcpp::stop("lower and upper must have the same number of columns");
  }
  
  if (!(blc.n_cols==p)) {
    Rcpp::stop("blc must have the same number of columns as mean");
  }
  
  if (!(init.n_cols==p)) {
    Rcpp::stop("init must have the same number of columns as mean");
  }
  
  
  // check boundary conditions
  if (sum(any(lower >= upper)) > 0) {
    Rcpp::stop("lower bound must be smaller than upper bound");
  }
  
  // check initial values, and generate automatically if needed
  arma::mat blct = blc.t();
  arma::mat initc = init*blct;
  if (sum(all(initc >= lower + 1e-8 && initc <= upper - 1e-8)) < m) {
    arma::mat estimate(n1,m);
    for (arma::uword i=0; i<n1; i++) {
      for (arma::uword j=0; j<m; j++) {
        if (std::isinf(lower(i,j)) && std::isinf(upper(i,j))) {
          estimate(i,j) = 0;  
        } else if (std::isinf(lower(i,j))) {
          estimate(i,j) = upper(i,j) - 1e-8;  
        } else if (std::isinf(upper(i,j))) {
          estimate(i,j) = lower(i,j) + 1e-8;  
        } else {
          estimate(i,j) = 0.5*(lower(i,j) + upper(i,j));  
        }
      }
    }
    init = estimate * arma::pinv(blct);
  }
  
  // check whether a matrix has identical rows
  auto f = [](const arma::mat& y) {
    const unsigned int n=y.n_rows, p=y.n_cols; 
    bool identical=1;
    if (n>1) {
      for (arma::uword i=1; i<n; i++) {
        for (arma::uword j=0; j<p; j++) {
          if (y(i,j) != y(i-1,j)) {
            identical = 0;
            break;
          }
        }
        if (identical==0) break;
      }
    }
    return identical;
  };
  
  
  // transform to the problem with identity covariance matrix
  arma::mat cholt = arma::trans(arma::chol(sigma));
  arma::mat R = blc*cholt;
  arma::uvec js(p);
  std::iota(js.begin(), js.end(), 0);
  Rcpp::NumericVector mu {0};
  
  // Gibbs step for sampling truncated multivariate t
  auto g = [p, m, R, js, mu, nu](arma::vec a, arma::vec b, double u, 
                                 arma::vec& z) {
    double denom = sqrt(u/nu);
    arma::vec y = z*denom;  // initial value for tmvnorm for this step
    arma::vec ay = a*denom, by = b*denom;
    
    for (arma::uword j=0; j<p; j++) {
      // set up linear inequality constraints for z(j)
      arma::vec rj = R.col(j);
      arma::uvec j2 = arma::find(js != j);
      arma::mat Rj = R.cols(j2);
      arma::vec yj = y(j2);
      arma::vec atempy = ay - Rj*yj;
      arma::vec btempy = by - Rj*yj;
      
      // determine lower and upper bounds for z(j)
      double lowerj=R_NegInf, upperj=R_PosInf;
      for (arma::uword k=0; k<m; k++) {
        if (rj(k) != 0) {
          double ak=atempy(k)/rj(k), bk=btempy(k)/rj(k);
          if (rj(k) > 0) {
            if (ak > lowerj) lowerj = ak;
            if (bk < upperj) upperj = bk;
          } else {
            if (bk > lowerj) lowerj = bk;
            if (ak < upperj) upperj = ak;
          }
        }
      }
      
      // generate y(j) for truncated univariate normal
      arma::vec lowerj1(1); lowerj1(0) = lowerj;
      arma::vec upperj1(1); upperj1(0) = upperj;
      y(j) = rtnorm(mu, 1, lowerj1, upperj1)(0);
    }
    
    z = y/denom; // update tmvt for this step;
  };
  
  
  arma::vec mean1(p), lower1(m), upper1(m), init1(p);
  arma::vec a(m), b(m), z(p);
  
  // obtain burn + n samples in case of identical means, lower and upper bounds 
  if (f(mean) && f(lower) && f(upper)) {
    mean1 = arma::trans(mean.row(0));
    lower1 = arma::trans(lower.row(0));
    upper1 = arma::trans(upper.row(0));
    init1 = arma::trans(init.row(0));
    
    a = lower1 - blc*mean1;
    b = upper1 - blc*mean1;
    z = solve(cholt, init1-mean1);
    
    for (arma::uword i=0; i<burn+n; i++) {
      double u = R::rchisq(nu);
      g(a, b, u, z);
      if (i>=burn) x.row(i-burn) = arma::trans(cholt*z + mean1);
    }
    return x;
  }
  
  // obtain (burn+1)*n samples for non-identical means, lower, or upper bounds
  if (n1==1) {
    lower1 = arma::trans(lower.row(0));
    upper1 = arma::trans(upper.row(0));
    init1 = arma::trans(init.row(0));
  }
  
  for (arma::uword i=0; i<n; i++) {
    mean1 = arma::trans(mean.row(i));
    if (n1==n) {
      lower1 = arma::trans(lower.row(i));
      upper1 = arma::trans(upper.row(i));
      init1 = arma::trans(init.row(i));
    }
    
    a = lower1 - blc*mean1;
    b = upper1 - blc*mean1;
    z = solve(cholt, init1-mean1);
    
    for (arma::uword i2=0; i2<burn+1; i2++) {
      double u = R::rchisq(nu);
      g(a, b, u, z);
    }
    x.row(i) = arma::trans(cholt*z + mean1);
  }
  return x;
}



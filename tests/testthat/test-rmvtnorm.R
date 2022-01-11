test_that("full-rank, identical means, lower and upper bounds", {
  d <- 3
  rho <- 0.9
  sigma <- matrix(0, nrow=d, ncol=d)
  sigma <- rho^abs(row(sigma) - col(sigma))
  D1 <- diag(1,d) # Full rank
  
  set.seed(1203)
  ans.1 <- tmvmixnorm::rtmvn(n=1000, Mean=1:d, sigma, D=D1, lower=rep(-1,d), upper=rep(1,d),
                             int=rep(0,d), burn=50, thin=0)
  chk.1 <- apply(ans.1, 2, summary)
  
  set.seed(1203)
  n = 1000;
  mean = matrix(rep(1:d,n), nrow=n, ncol=d, byrow=TRUE)
  lower = matrix(rep(-1,d), nrow=1);
  upper = matrix(rep(1,d), nrow=1);
  init = matrix(rep(0,d), nrow=1);
  ans.1b = rtmvnorm(mean, sigma, D1, lower, upper, init, burn=50)
  chk.1b <- apply(ans.1b, 2, summary)
  
  expect_equal(chk.1, chk.1b)
})


test_that("non-full rank, identical means, lower and upper bounds", {
  d <- 3
  rho <- 0.5
  sigma <- matrix(0, nrow=d, ncol=d)
  sigma <- rho^abs(row(sigma) - col(sigma))
  D2 <- matrix(c(1,1,1,0,1,0,1,0,1),ncol=d)
  
  set.seed(1228)
  ans.2 <- tmvmixnorm::rtmvn(n=100, Mean=1:d, sigma, D=D2, lower=rep(-1,d),
                 upper=rep(1,d), burn=10, thin=0)
  chk.2 <- apply(ans.2, 2, summary)
  
  set.seed(1228)
  n = 100;
  mean = matrix(rep(1:d,n), nrow=n, ncol=d, byrow=TRUE)
  lower = matrix(rep(-1,d), nrow=1);
  upper = matrix(rep(1,d), nrow=1);
  init = matrix(rep(0.8,d), nrow=1);
  ans.2b = rtmvnorm(mean, sigma, D2, lower, upper, init, burn=10)
  chk.2b <- apply(ans.2b, 2, summary)
  
  expect_equal(chk.2, chk.2b)
})


test_that("non-full rank, different means", {
  d <- 3
  rho <- 0.5
  sigma <- matrix(0, nrow=d, ncol=d)
  sigma <- rho^abs(row(sigma) - col(sigma))
  D3 <- matrix(c(1,0,1,1,1,0),nrow=d-1,ncol=d,byrow=TRUE)

  set.seed(3084)
  n = 100;
  mean = matrix(runif(n*d), nrow=n, ncol=d);
  lower = matrix(rep(-1,d-1), nrow=1);
  upper = matrix(rep(1,d-1), nrow=1);
  init = matrix(rep(0.8,d), nrow=1);
  ans.3 = rtmvnorm(mean, sigma, D3, lower, upper, init, burn=10)
  chk.3 <- apply(ans.3, 2, summary)

  chk.3b = matrix(c(-0.95373733, -1.28766084, -1.38358884,
                    -0.28714981, -0.35548892, -0.27641684,
                    0.05976233,  0.08396816,  0.09040197,
                    0.05838389,  0.04167013,  0.10280794,
                    0.36374043,  0.41420365,  0.51658866,
                    1.77116523,  0.97481806,  1.26606795), 
                  nrow=6, ncol=3, byrow=TRUE);
  dimnames(chk.3b) = dimnames(chk.3);
  
  expect_equal(chk.3, chk.3b)
})

test_that("truncated mvt, non-full rank, different means", {
  d = 3;
  rho = 0.5;
  sigma = matrix(0, d, d);
  sigma = rho^abs(row(sigma) - col(sigma));
  nu = 10;
  blc = matrix(c(1,0,1,1,1,0),nrow=d-1,ncol=d,byrow=TRUE)
  n = 100;
  set.seed(3084)
  mean = matrix(runif(n*d), nrow=n, ncol=d);
  lower = matrix(rep(-1,d-1), nrow=1);
  upper = matrix(rep(1,d-1), nrow=1);
  init = matrix(rep(0.8,d), nrow=1);
  result = rtmvt(mean, sigma, nu, blc, lower, upper, init, burn=50)
  chk.4 <- apply(result, 2, summary)
  
  chk.4b = matrix(c(-1.40342755, -1.74138329, -1.68153340,
                    -0.39575796, -0.48564615, -0.50035150,
                    -0.04051560, -0.09808149,  0.02762573,
                    -0.03028324, -0.04621082,  0.04633071,
                    0.42064496,  0.36668839,  0.59122084,
                    1.11797524,  1.98752407,  1.76118335), 
                  nrow=6, ncol=3, byrow=TRUE);
  dimnames(chk.4b) = dimnames(chk.4);
  
  expect_equal(chk.4, chk.4b)
})


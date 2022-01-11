#' @name tmvtnsim-package
#' @aliases tmvtnsim-package
#' @docType package
#'
#' @title Truncated Multivariate Normal and t Distribution Simulation
#'
#' @description Simulation of random vectors from truncated multivariate 
#' normal and t distributions based on the algorithms proposed by 
#' Yifang Li and Sujit K. Ghosh (2015) <doi:10.1080/15598608.2014.996690>. 
#' We allow the mean, lower and upper bounds to differ across samples 
#' to accommodate regression problems. The algorithms are implemented 
#' in C++ and hence are highly efficient. 
#' 
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Yifang Li and Sujit K. Ghosh. Efficient Sampling Methods for Truncated
#' Multivariate Normal and Student-t Distributions Subject to Linear 
#' Inequality Constraints. Journal of Statistical Theory and Practice. 
#' 2015;9:712-732. \doi{10.1080/15598608.2014.996690}
#'
#' @useDynLib tmvtnsim, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
NULL


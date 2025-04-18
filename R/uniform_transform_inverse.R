#' @title inverse transformation function
#' @description
#' This function uses a reverts the transformation for uniform distributions
#' @param theta_trans vector of parameters to transformed back
#' @param args a list of arguments as returned by \link[SMCfeatures]{define_args_logistic_growth}
#' @return vector of transformed parameters
#' @export
#'
uniform_transform_inverse <- function(theta_trans,args){
  # This function inverses the transform for uniform distributions.
  lower_mat <- matrix(args$lower, nrow=nrow(theta_trans), ncol=ncol(theta_trans),
                      byrow=TRUE)
  upper_mat <- matrix(args$upper, nrow=nrow(theta_trans), ncol=ncol(theta_trans),
                      byrow=TRUE)

  theta <- (lower_mat+upper_mat*exp(theta_trans))/(1+exp(theta_trans))
  return(theta)
}

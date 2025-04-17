#' @title transformation function
#' @description
#' This function uses a transform for uniform distributions, so all values are within the uniform bounds.
#'
#' @param theta vector of parameters to transform
#' @param args a list of arguments as returned by \link[SMCfeatures]{define_args_logistic_growth}
#' @return vector of transformed parameters
#' @export
#'
uniform_transform <- function(theta,args){
  # This function uses a transform for uniform distributions, so that all
  # values are within the uniform bounds.
  lower_mat <- matrix(args$lower, nrow=nrow(theta), ncol=ncol(theta),
                      byrow=TRUE)
  upper_mat <- matrix(args$upper, nrow=nrow(theta), ncol=ncol(theta),
                      byrow=TRUE)
  return(log((theta-lower_mat)/(upper_mat-theta)))

}

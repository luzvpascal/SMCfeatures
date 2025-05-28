#' @title Uniform sampling of parameters
#' @description
#' This function uniformly samples parameters between lower and upper bounds defined in args
#'
#' @param args a list of arguments as returned by \link[SMCfeatures]{define_args_logistic_growth}
#' @return a vector of parameter samples
#' @export
#' @import stats

uniform_sampler <- function(args){
  stats::runif(args$n_params, args$lower,args$upper)
} #uniform distribution

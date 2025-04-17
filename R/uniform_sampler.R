#' @title Uniform sampling of parameters
#' @description
#' Uniform sampling of parameters between lower and upper bounds defined in args
#'
#' @param args a list of arguments as returned by \link[SMCfeatures]{args_function}
#' @return a vector of sampled parameters
#' @export
#' @import stats

uniform_sampler <- function(args){
  stats::runif(args$n_params, args$lower,args$upper)
} #uniform distribution

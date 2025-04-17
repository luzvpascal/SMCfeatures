#' @title logLikelihood_logistic_growth
#' @description
#' Gaussian likelihood function
#' @param parameters parameters
#' @param simulation list of simulation arguments
#' @param args a list of arguments as returned by \link[SMCfeatures]{args_function}
#' @return log likelihood
#' @export
log_likelihood_logistic_growth <- function(parameters,
                                         simulation,
                                         args){

  n_data <- length(args$output_data)
  sigma <- parameters[length(parameters)]

  if (args$include_data_constraints){#if the user includes data
    return( -n_data*log(sigma)-sum(sum((simulation$output_data-args$output_data)^2/(2*sigma^2))))
  } else{
    return(0)
  }
}

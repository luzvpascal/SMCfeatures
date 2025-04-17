#' @title Log-Likelihood function for the logistic growth case study
#' @description
#' Gaussian likelihood function
#' @param parameters vector of parameters
#' @param simulation list of simulation arguments
#' @param args a list of arguments as returned by \link[SMCfeatures]{define_args_logistic_growth}
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

#' @title Log-Likelihood function for the logistic growth case study
#' @description
#' Gaussian likelihood function
#' @param parameters vector of parameters, with last parameter the standard deviation
#' @param simulation list like object that contains output_data (output of simulation)
#' @param args a list like object of arguments that includes a boolean include_data_constraints, that indicates if data constraints should be included, and output_data (output of observed data)
#' @return log likelihood of simulation outputs against empirical data if args$include_data_constraints set to TRUE, 0 otherwise.
#' @export
log_likelihood_Gaussian <- function(parameters,
                                   simulation,
                                   args){

  n_data <- length(args$output_data)
  sigma <- parameters[length(parameters)]

  if (args$include_data_constraints){#if the user includes data
    return( -n_data*log(sigma)-sum(sum((simulation$y_data-args$output_data)^2/(2*sigma^2))))
  } else{
    return(0)
  }
}

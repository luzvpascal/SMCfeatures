#' @title Model on data for logistic growth case study
#' @description
#' Simulate outputs of calibrated model for the logistic growth case study
#' @param parameters vector of parameters (growth rate, capacity, initial abundance)
#' @param simulation list of simulation arguments, as returned by \link[SMCfeatures]{model_logistic_growth_constraint}
#' @param args a list of arguments as returned by \link[SMCfeatures]{define_args_logistic_growth}
#' @return list like object of simulation with added time, and y (output of simulation)
#' @export
model_logistic_growth_plots <- function(parameters,
                                       simulation,
                                       args){
  #Inputs parameters and outputs a simulation

  # Extract information
  r <- parameters[1]
  K <- parameters[2]
  y0 <- parameters[3]

  # Simulation details
  simulation$t <- seq(0, max(args$input_data), length.out=100)
  simulation$y <- (K*y0)/(y0+(K-y0)*exp(-simulation$t*r))
  return(simulation)
}


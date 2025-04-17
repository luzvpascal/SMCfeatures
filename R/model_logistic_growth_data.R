#' @title model_logistic_growth_data
#' @description
#' Simulate
#' @param parameters vector of parameters (growth rate, capacity, initial abundance)
#' @param simulation list of simulation arguments
#' @param args a list of arguments as returned by \link[SMCfeatures]{args_function}
#' @return simulation with added t_data, and y_data (output of simulation on input data)
#' @export
model_logistic_growth_data <- function(parameters,
                                       simulation,
                                       args){
  # likelihood needs data at all time points
  #Inputs parameters and outputs a simulation

  # Extract information
  r <- parameters[1]
  K <- parameters[2]
  y0 <- parameters[3]

  # Simulation details
  simulation$t_data <- args$input_data
  simulation$y_data <- (K*y0)/(y0+(K-y0)*exp(-simulation$t_data*r))
  return(simulation)
}


#' @title Model simulation on input data
#' @description
#' Simulate outputs of calibrated model, on data points only
#'
#' @param parameters vector of parameters (growth rate, capacity, initial abundance)
#' @param simulation list of simulation arguments, as returned by \link[SMCfeatures]{model_logistic_growth_constraint}
#' @param args a list of arguments as defined in \link[SMCfeatures]{SMCfeatures}
#' @return list like object representing simulation outputs
#' t_data: time steps for data (args$input_data)
#' y_data: output of simulation on input data
#' @export
model_on_data <- function(parameters,
                         simulation,
                         args){
  # likelihood needs data at all time points
  #Inputs parameters and outputs a simulation

  # Simulation details
  simulation$t_data <- args$input_data
  simulation$y_data <- args$simulate_model(parameters,
                                             simulation$t_data)
  return(simulation)
}

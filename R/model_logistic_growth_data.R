#' @title Model on data for logistic growth case study, on data points only
#' @description
#' Simulate outputs of calibrated model, on data points only
#' @param parameters vector of parameters (growth rate, capacity, initial abundance)
#' @param simulation list of simulation arguments, as returned by \link[SMCfeatures]{model_logistic_growth_constraint}
#' @param args a list of arguments as returned by \link[SMCfeatures]{define_args_logistic_growth}
#' @return ist like object representing simulation outputs
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
  simulation$y_data <- args$model_simulation(parameters,
                                             simulation$t_data)
  return(simulation)
}

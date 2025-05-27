#' @title Model simulation on input data
#' @description
#' Simulate outputs of calibrated model, on data points only
#' @param parameters vector of parameters (growth rate, capacity, initial abundance)
#' @param args a list of arguments as defined in \link[SMCfeatures]{SMCfeatures}
#' @return ist like object representing simulation outputs
#' t_data: time steps for data (args$input_data)
#' y_data: output of simulation on input data
#' @export
model_for_plots_logistic_growth <- function(parameters,
                                            args){
  # likelihood needs data at all time points
  #Inputs parameters and outputs a simulation

  # Simulation details
  simulation <- list()

  simulation$t_plot <- seq(0,50)
  simulation$y_plot <- args$simulate_model(parameters,
                                           simulation$t_plot)
  return(simulation)
}

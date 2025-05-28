#' @title Model simulation on input data
#' @description
#' Simulate outputs of model for plotting purposes
#'
#' @param parameters vector of parameters (growth rate, capacity, initial abundance)
#' @param simulation list of simulation arguments, as returned by \link[SMCfeatures]{model_logistic_growth_constraint}
#' @param args a list of arguments as defined in \link[SMCfeatures]{SMCfeatures}
#' @return list like object representing simulation outputs
#' t_data: time steps for data (args$input_data)
#' y_data: output of simulation on input data
#' @export
model_for_plots_logistic_growth <- function(parameters,
                                            simulation,
                                            args){
  # likelihood needs data at all time points
  #Inputs parameters and outputs a simulation

  # Simulation details
  simulation$t_plot <- seq(0,50)
  simulation$y_plot <- args$simulate_model(parameters,
                                           simulation$t_plot)
  return(simulation)
}

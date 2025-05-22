#' @title model_logistic_growth_constraint
#' @description
#' Simulate
#' @param parameters vector of parameters (growth rate, capacity, initial abundance, sigma)
#' @return list like object representing simulation outputs
#' t_const: vector of time constraints (time step 5 and 50)
#' y_const: output of simulation on input t_const
#' K: carrying capacity
#' @export
model_logistic_growth_constraint <- function(parameters){
  #Inputs parameters and outputs a simulation

  simulation <- list()

  # Simulation details
  simulation$t_const <- c(5, 50)
  simulation$y_const <- SMCfeatures::model_logistic_growth(parameters,
                                                           simulation$t_const)
  simulation$K <- K

  return(simulation)
}

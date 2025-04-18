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
  # Extract information
  r <- parameters[1]
  K <- parameters[2]
  y0 <- parameters[3]

  # Simulation details
  simulation$t_const <- c(5, 50)
  simulation$y_const <- (K*y0)/(y0+(K-y0)*exp(-simulation$t_const*r))
  simulation$K <- K

  return(simulation)
}

#' @title Simulate constraints and calculate discrepancy function
#' @description
#' Simulate constraints and calculate discrepancy function. In this example, there are two non-empirical constraints: at year 5, the coral cover must be lower that 10%, and at year 50, the coral cover must be higher that K-1%.
#' @param parameters a vector of parameters. Here parameters must be input in this order r, K, y0 and sigma
#' @param args a list of arguments as defined in the function \link[SMCfeatures]{SMCfeatures}
#' @return a list with
#' simulation list like object with t_const, y_const, K.
#' a summary statistic (discrepancy measure).
#' @export
#'
discrepancy_logistic_growth <- function(parameters,
                                        args){

  simulation <- list()
  # Simulation details
  simulation$t_const <- c(5, 50)
  simulation$y_const <- SMCfeatures::model_logistic_growth(parameters,
                                                           simulation$t_const)
  simulation$K <- parameters[2]
  if (args$include_expert_constraints){

    ## 5yr constraint
    # max coral population in first 5 years
    max_population_in_5 <- simulation$y_const[1]

    # the max allowable coral population in first 5 years
    tolerance_in_5 <- 10

    #measure constraint discrepancy
    discrepancy_5 <- max(max_population_in_5 - tolerance_in_5,0)

    ## 50yr constraint
    # difference between carrying capacity and final coral population at 50 years
    difference_from_K <- simulation$K - simulation$y_const[2]

    # max allowable difference
    tolerance_in_50 <- 1 #within 1# of carrying cap

    #measure constraint discrepancy
    discrepancy_50 <- max(difference_from_K - tolerance_in_50,0)

    ## Discrepancy measure
    discrepancy <- discrepancy_5 + discrepancy_50

    return(list(simulation,discrepancy))
  } else {
    return(list(simulation,args$discrepancy_final))
  }
}

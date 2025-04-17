#' @title Discrepancy function
#' @description
#' Discrepancy measure
#' @param simulation list of simulation arguments
#' @param args a list of arguments as returned by \link[SMCfeatures]{args_function}
#' @return summary statistic (discrepancy measure).
#' @export
#'
discrepancy_logistic_growth <- function(simulation,
                                        args){

  if (args$include_expert_constraints){
    ## 5yr constraint
    # max coral population in first 5 years
    max_population_in_5 <- max(simulation$y_const[simulation$t_const==5])

    # the max allowable coral population in first 5 years
    tolerance_in_5 <- 10

    #measure constraint discrepancy
    discrepancy_5 <- max(max_population_in_5 - tolerance_in_5,0)

    ## 50yr constraint
    # difference between carrying capacity and final coral population at 50 years
    difference_from_K <- simulation$K - simulation$y_const[simulation$t_const==50]

    # max allowable difference
    tolerance_in_50 <- 1 #within 1# of carrying cap

    #measure constraint discrepancy
    discrepancy_50 <- max(difference_from_K - tolerance_in_50,0);

    ## Discrepancy measure
    discrepancy <- discrepancy_5 + discrepancy_50

    return(discrepancy)
  } else {
    return(args$discrepancy_final)
  }
}

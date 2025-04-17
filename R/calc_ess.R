#' @title Calculate effective sample size
#' @description
#' Calculates the effective sample size
#' @param newgamma new temperature
#' @param oldgamma old temperature
#' @param part_w particule weights
#' @param part_loglike particule log likelihood
#' @return Effective sample size
#' @export
calc_ess <- function(newgamma,oldgamma,part_w,part_loglike){
  #define global arguments
  log_part_w = log(part_w) + (newgamma-oldgamma)*part_loglike

  if (newgamma == oldgamma){
   log_part_w <- log(part_w)
  }

  #numerically stabilize before exponentiating
  log_part_w <- log_part_w - max(log_part_w)
  part_w <- exp(log_part_w)

  #normalize weights
  part_w <- part_w/sum(part_w)

  # compute effective sample size
  ess <- 1/sum(part_w^2)

  return(ess)
}

#' @title Monte Carlo Markov Chain
#' @description
#' Runs MCMC
#' @param args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @param mcmc_steps number of MCMC steps
#' @param params_transformed vector of values of parameters for current particle
#' @param cov_matrix covariance matrix
#' @param disc current discrepancy of the particle
#' @param dist_next next objective distance
#' @param acc_count number of times the particle has been moved
#' @param loglike log-likelihood of parameters
#' @param newgamma new temperature
#' @return list
#' params_transformed: vector of moved parameters
#' acc_count: number of times particle has been moved
#' sims: simulation outputs on expert constraint
#' loglike log likelihood of moved parameters
#' @export
#' @import MASS

MCMC_combined <- function(args,
                          mcmc_steps,
                          params_transformed,
                          cov_matrix,
                          disc,
                          dist_next,
                          acc_count,
                          loglike,
                          newgamma){
  for (r in seq(mcmc_steps)){
    # Gaussian random walk
    prop_transformed <- MASS::mvrnorm(n=1,mu=params_transformed,Sigma=cov_matrix)
    # Transform back to calculate prior probs and discrepancy
    prop <- args$trans_finv(matrix(prop_transformed, nrow=1),args)

    # Calculate prior probabilities
    prior_curr <- args$prior_pdf(matrix(params_transformed, nrow=1))
    prior_prop <- args$prior_pdf(matrix(prop_transformed, nrow=1))

    # early rejection (assumes symmetric proposals)
    if(((is.nan(prior_prop/prior_curr))| (runif(1) > prior_prop/prior_curr))){
      next
    }

    #find proposal discrepancy
    sims_prop <- args$simulate_constraint(prop)
    disc_prop <- args$calculate_discrepancy(sims_prop, args)

    # ABC part of the acceptance probability
    # reject if discrepancy is above threshold
    if (disc_prop > dist_next) {
      next
    }

    # then the metropolis-hastings ratio is satisfied
    sims_prop <- args$simulate_data(prop,sims_prop,args)
    log_like_prop  <- args$calculate_log_likelihood(prop,sims_prop,args)
    mh  <- exp(newgamma*(log_like_prop - loglike))

    if (runif(1) < mh){
      # then the metropolis-hastings ratio is satisfied
      params_transformed <- prop_transformed
      disc <- disc_prop
      acc_count <- acc_count + 1
      sims <- sims_prop
      loglike <- log_like_prop
    }
  }
  return(list(params_transformed=params_transformed,
              disc=disc,
              acc_count=acc_count,
              sims=sims,
              loglike=loglike))
}

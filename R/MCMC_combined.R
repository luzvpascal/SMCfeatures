#' @title Monte Carlo Markov Chain
#' @description
#' Runs MCMC
#' @param i line number of particle to be moved
#' @param sim_args a list of arguments as returned by \link[EEMtoolbox]{args_function}
#' @param mcmc_trials number of MCMC steps to try before selecting appropriate number.
#' @param dist_next next objective distance
#' @param part_vals matrix of current values of parameters of each particle (each particle represented by a row)
#' @param part_s vector of discrepancy measures of each particle (each particle represented by a row)
#' @param part_sim matrix of current simulation values: equilibriums and eigenvalues of jacobian (each particle represented by a row)
#' @param cov_matrix covariance matrix
#' @param summ_func function calculating equilibrium points and real parts of the Jacobians eigenvalues to summarise ecosystem features. Default =summarise_ecosystem_features_GLV. Options include summarise_ecosystem_features_Baker (automatically chosen if model="Bimler-Baker") and summarise_ecosystem_features_Gompertz, (automatically chosen if model="Gompertz"). Needs to be defined if model="customized" chosen.
#' @param disc_func summary statistic (discrepancy measure).
#' @param trans_finv inverse of trans_f function.
#' @param pdf joint probability density function.
#' @param loglike log-likelihood
#' @param newgamma new temperature
#' @return list
#' part_vals: updated value of parameter values for particle i
#' part_s: discrepancy measure for particle i
#' part_sim: summary of ecosystem features for particle i
#' i_acc: number of times particle movement accepted
#' sims_mcmc: number of successful walks
#' @export
#' @import MASS

MCMC_combined <- function(
                 mcmc_steps,
                 params_transformed,
                 cov_matrix,
                 args,
                 disc,
                 acc_count,
                 dist_next,
                 sims,
                 loglike,
                 newgamma){
  # documentation
  sims_mcmc <-  0
  i_acc <- 0
  for (r in seq(mcmc_steps)){
    # Gaussian random walk
    prop_transformed <- MASS::mvrnorm(n=1,mu=params_transformed,Sigma=cov_matrix)
    # Transform back to calculate prior probs and discrepancy
    prop <- args$trans_finv(prop_transformed,args)

    # Calculate prior probabilities
    prior_curr <- args$pdf(params_transformed)
    prior_prop <- args$pdf(prop_transformed)

    # early rejection (assumes symmetric proposals)
    if(((is.nan(prior_prop/prior_curr))| (runif(1) > prior_prop/prior_curr))){
      next
    }

    #find proposal discrepancy
    sims_prop <- args$simulate_constraint(prop)
    disc_prop <- args$calculate_discrepancy(sims_prop, args)

    sims_mcmc=sims_mcmc+1
    # ABC part of the acceptance probability

    # reject if discrepancy is above threshold
    if (dist_prop > dist_next) {
      next
    }

    # then the metropolis-hastings ratio is satisfied
    sims_prop <- args$simulate_data(prop,args,sims_prop)
    log_like_prop  <- args$calculate_log_likelihood(prop,sims_prop,args)
    mh  <- exp(newgamma*(log_like_prop - loglike))

    if (runif(1) < mh){
      # then the metropolis-hastings ratio is satisfied
      params_transformed <- prop_transformed
      disc <- disc_prop
      sims <- sims_prop
      loglike <- log_like_prop
      acc_count <- acc_count + 1
    }
  }
  return(list(params_transformed=params_transformed,
              disc=disc,
              acc_count=acc_count,
              sims=sims,
              loglike=loglike))
}

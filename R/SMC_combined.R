#' @title Generation of model ensembles Sequential Monte Carlo - Approximate Bayesian Computation.
#' @description
#' Generation of model ensembles using Adaptive sequential Monte Carlo for approximate Bayesian computation
#' @param args a list of arguments as returned by \link[SMCfeatures]{define_args_logistic_growth}
#' @param n_particles Number of desired ensemble members. Default to 10000
#' @param mcmc_trials number of MCMC steps to try before selecting appropriate number. Default 10
#' @param discrepancy_final target discrepancy threshold. Default 0. If zero, p_acc_min is used to determine stopping criteria.
#' @param a_disc tuning parameter for adaptive selection of discrepancy threshold sequence.
#' @param a_like tuning parameter for adaptive selection of likelihood ESS sequence
#' @param c tuning parameter for choosing the number of MCMC iterations in move step. Default 0.01
#' @param p_acc_min minimum acceptable acceptance rate in the MCMC interations before exit. Default 0.0001
#' @param n_cores Number of cores desired to be used for sampling. Default set to 1 core (sequential sampling).
#' @return list: sims=number of simulations
#' part_vals=parameter values
#' part_s=discrepancy value
#' prior_sample=prior distribution
#' @export
#' @import parallel
#' @import doParallel
#' @import foreach
#'
SMC_combined <- function(args,
                         n_particles,
                         mcmc_trials,
                         a_disc,
                         a_like,
                         c,
                         p_acc_min,
                         n_cores){

  ESS <- a_like*n_particles # target effective sample size based on retaining a_like% of particles

  ###########################################################
  # Prior sample -- Initialise ##############################
  ###########################################################

  print('Begin sampling from prior distribution')

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  #sample priors
  param_vals <- foreach::foreach(i = 1:n_particles, .combine="rbind") %dopar% {
    args$sampler(args)
  }
  param_vals <- unname(param_vals)

  #simulate model for the constraint
  param_sims_const <- foreach::foreach(i = 1:n_particles) %dopar% {
    args$simulate_constraint(param_vals[i,])
  }

  #calculate discrepancy
  param_disc <- foreach::foreach(i = 1:n_particles, .combine="rbind") %dopar% {
    args$calculate_discrepancy(param_sims_const[[i]],args)
  }
  param_disc <- c(unname(param_disc))

  #simulate model for the likelihood
  param_sims <- foreach::foreach(i = 1:n_particles) %dopar% {
    args$simulate_data(param_vals[i,],param_sims_const[[i]],args)
  }

  #calculate likelihood
  param_loglike <- foreach::foreach(i = 1:n_particles) %dopar% {
    args$calculate_log_likelihood(param_vals[i,],param_sims[[i]],args)
  }
  param_loglike <- unlist(param_loglike)
  #stop cluster
  parallel::stopCluster(cl)
  rm(cl)

  #save prior sample
  prior <- param_vals

  param_w <- rep(1,n_particles)/n_particles # particle weightings

  print(paste0('Initial constraint acceptance: ',
               sum(param_disc<=args$discrepancy_final),
               " of ",
               n_particles,
               " particles"))

  ###########################################################
  # SMC #####################################################
  ###########################################################
  #Set up to run the sequence of distributions
  gamma_t <- 0 #initial temperature for likelihood annealing
  dist_max <- max(param_disc) #initial maximum discrepancy
  SMC_count <- 0 #track sequential steps

  ## Sequential distributions
  #Run until likelihood is fully integrated and constraints are met

  if (gamma_t < 1 | dist_max > args$discrepancy_final){
    print('Begin SMC-ABC search method')
  }

  while(gamma_t < 1 | dist_max > args$discrepancy_final){
    SMC_count  <- SMC_count+1 #track number of SMC steps
    print(paste0("SMC iteration: ", SMC_count))

    #########################################################################################################
    # -----REWEIGHT----- Select next discrepancy and likelihood thresholds, and weight particles accordingly#
    #########################################################################################################

    #UPDATE DISCREPANCY AND REWEIGHT ###################################
    #select discrepancy threshold
    dist_next  <- max(quantile(param_disc,a_disc),args$discrepancy_final) #retaining only the top a_disc percent

    # weight according to new discrepancy threshold -- set weights for
    # particles with discrepancy above threshold to 0
    param_loglike_weighted_by_discrepancy_threshold <- param_loglike
    param_loglike_weighted_by_discrepancy_threshold[param_disc>dist_next] = -Inf #check here

    param_w_weighted_by_discrepancy_threshold <- param_w
    param_w_weighted_by_discrepancy_threshold[param_disc>dist_next] <- 0


    #UPDATE LIKELIHOOD AND REWEIGHT ###################################
    #check effective sample size (ESS) for final posterior (when gamma_t=1)
    ess_posterior <- SMCfeatures::calc_ess(1,
                                           gamma_t,
                                           param_w_weighted_by_discrepancy_threshold,
                                           param_loglike_weighted_by_discrepancy_threshold)
    #if ESS at posterior is acceptable (above the minimum amount ESS)
    if (ess_posterior > ESS){
      newgamma <- 1 #sample the posterior
    } else{
      #find the temperature (gamma_t) at which the effective sample size is equal to the threshold ESS
      solution <- pracma::fzero(function(x){SMCfeatures::calc_ess(x,
                                                     gamma_t,
                                                     param_w_weighted_by_discrepancy_threshold,
                                                     param_loglike_weighted_by_discrepancy_threshold)-ESS},
                                gamma_t)
      newgamma <- solution$x
    }

    #calculate particle likelihood weights for new temperature (gamma)
    log_particle_w <- log(param_w) +(newgamma - gamma_t)*param_loglike
    log_particle_w <- log_particle_w - max(log_particle_w) #numerically stabilise
    param_w <- exp(log_particle_w) #return from log space
    param_w <- param_w/sum(param_w) #normalise weights

    #set weights to 0 if discrepancy is above tolerance
    param_w[param_disc>dist_next] <- 0
    param_w <- param_w/sum(param_w) #normalise weights

    #PREPARE TO RESAMPLE  ###################################

    #print the new discrepancy threshold and likelihood annealing temperature
    print(paste0('Discrepancy threshold is ',  round(dist_next, digits = 5), ' (aiming for ', args$discrepancy_final, ")"))
    print(paste0("Annealing temperature is ", round(newgamma, digits = 5), " (aiming for 1.00)"))

    #apply transform
    param_vals_transformed <- args$trans_f(param_vals, args)
    param_vals <- matrix(0,ncol=ncol(param_vals), nrow=nrow(param_vals))


    ###########################################################################
    ## -----RESAMPLE-----    ##################################################
    ###########################################################################
    #duplicate good particles
    r <- sample(seq(n_particles),n_particles, replace=TRUE, prob=param_w) #duplicate good ones

    param_vals_transformed <- param_vals_transformed[r,]
    param_loglike <- param_loglike[r]
    param_sims <- param_sims[r]
    param_disc <- param_disc[r]

    #reset weights
    param_w = rep(1,n_particles)/n_particles # particle weightings

    #calc cov for random walk MH-MCMC step #check where this goes?
    cov_matrix = (2.38^2)*cov(param_vals_transformed)/dim(param_vals_transformed)[2]

    ###########################################################################
    ## -----MOVE-----    ##################################################
    ###########################################################################

    # MCMC MOVE STEPS
    acc_counter = rep(0, n_particles) #track particle move acceptance

    #initial trial of MCMC iterations####
    #parallel for loop
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    mcmc_outcome <- foreach::foreach(i=1:n_particles,
                                     .packages =c('SMCfeatures')) %dopar% {
                                       SMCfeatures::MCMC_combined(args,
                                                                  mcmc_trials,
                                                                  param_vals_transformed[i,],
                                                                  cov_matrix,
                                                                  param_disc[i],
                                                                  dist_next,
                                                                  param_sims[i],
                                                                  param_loglike[i],
                                                                  acc_counter[i],
                                                                  newgamma)
                                     }
    parallel::stopCluster(cl)#stop cluster
    rm(cl)

    ## check dimensions here####
    param_vals_transformed <- matrix(unlist(lapply(mcmc_outcome, `[[`, 1)),
                                     nrow=length(mcmc_outcome), byrow = TRUE)

    param_disc <- (unlist(lapply(mcmc_outcome, `[[`, 2)))

    acc_counter <- (unlist(lapply(mcmc_outcome, `[[`, 3)))

    param_sims <- lapply(mcmc_outcome, `[[`, 4)

    param_loglike <- (unlist(lapply(mcmc_outcome, `[[`, 5)))

    # determine number of MCMC iterations to perform
    acc_rate <- sum(acc_counter)/(mcmc_trials*n_particles) #calculate acceptance rate
    mcmc_iters <- floor(log(c)/log(1-acc_rate)+1) #predict how many mcmc steps are needed to get approx unchanged paricles c#

    #do the remaining MCMC iterations####
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    mcmc_outcome <- foreach::foreach(i=1:n_particles,
                                     .packages =c('SMCfeatures')) %dopar% {
                                       SMCfeatures::MCMC_combined(args,
                                                                  mcmc_iters-mcmc_trials,
                                                                  param_vals_transformed[i,],
                                                                  cov_matrix,
                                                                  param_disc[i],
                                                                  dist_next,
                                                                  param_sims[i],
                                                                  param_loglike[i],
                                                                  acc_counter[i],
                                                                  newgamma)
                                     }
    parallel::stopCluster(cl)#stop cluster
    rm(cl)

    ## check dimensions here####
    param_vals_transformed <- matrix(unlist(lapply(mcmc_outcome, `[[`, 1)),
                                     nrow=length(mcmc_outcome), byrow = TRUE)

    param_disc <- (unlist(lapply(mcmc_outcome, `[[`, 2)))

    acc_counter <- (unlist(lapply(mcmc_outcome, `[[`, 3)))

    param_sims <- lapply(mcmc_outcome, `[[`, 4)

    param_loglike <- (unlist(lapply(mcmc_outcome, `[[`, 5)))


    #update number of mcmc trials required in next SMC step
    num_mcmc_iters = max(0, mcmc_iters - mcmc_trials) + mcmc_trials
    mcmc_trials = ceiling(mcmc_iters/2)

    # PREPARE FOR REWEIGHTING

    #convert back to theta for reweighting
    param_vals <- args$trans_finv(param_vals_transformed,args)

    #update the temperature and discrepancy max
    gamma_t = newgamma
    dist_max = max(param_disc)

    # calc percentage acceptance
    p_acc <- sum(acc_counter)/(num_mcmc_iters*n_particles)

    if (p_acc < p_acc_min){ # if acceptance is too low, then bail
      print('Getting out as MCMC acceptance rate is below acceptable threshold')
      break
    }
  }

  print('SMC algorithm complete')

  ## simulate
  posterior <- param_vals

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  param_sims_posterior <- foreach::foreach(i = 1:n_particles) %dopar% {
    args$simulate_plots(posterior[i,],param_sims[[i]],args)
  }
  param_sims_prior <- foreach::foreach(i = 1:n_particles) %dopar% {
    args$simulate_plots(prior[i,],param_sims[[i]],args)
  }
  parallel::stopCluster(cl)#stop cluster
  rm(cl)


  return(list(posterior=posterior,
              prior=prior,
              param_sims_posterior=param_sims_posterior,
              param_sims_prior=param_sims_prior))
}

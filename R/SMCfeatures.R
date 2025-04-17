#' @title SMCfeatures
#' @description
#' Main function for any parameterisation of expert elicited constraints
#' @param n_particles Number of desired ensemble members. Default to 10000
#' @param define_args summary statistic (discrepancy measure). Default discrepancy_logistic_growth
#' @param mcmc_trials number of MCMC steps to try before selecting appropriate number. Default 10
#' @param dist_final target discrepancy threshold. Default 0. If zero, p_acc_min is used to determine stopping criteria.
#' @param a_disc tuning parameter for adaptive selection of discrepancy threshold sequence.
#' @param a_like tuning parameter for adaptive selection of likelihood ESS sequence
#' @param c tuning parameter for choosing the number of MCMC iterations in move step. Default 0.01
#' @param p_acc_min minimum acceptable acceptance rate in the MCMC interations before exit. Default 0.0001
#' @param n_cores Number of cores desired to be used for sampling. Default set to 1 core (sequential sampling).
#' @param output_prior logical. If set to TRUE, algorithm returns prior distributions of parameters ensemble of parameters. Default FALSE
#' @param output_args logical. If set to TRUE, algorithm returns output from SMCfeatures::args_function for this problem
#' @param output_discrepancy logical. If set to TRUE, algorithm returns discrepancy values
#' @param output_matrix logical. If set to TRUE, algorithm returns interaction matrix and growthrates
#' @return vector of transformed parameters
#' @export
SMCfeatures <- function(define_args=SMCfeatures::define_args_logistic_growth,
                        include_expert_constraints=TRUE,
                        include_data_constraints=TRUE,
                        n_particles=10000,
                        mcmc_trials=10,
                        discrepancy_final=0,
                        a_disc=0.6,
                        a_like=0.3,
                        c=0.01,
                        p_acc_min=0.0001,
                        n_cores = 1L){

  if (n_cores == 1L){
    print("The code will run on 1 cluster only (sequential).")
    print("Change the parameter 'n_cores' to parallelise code")
  }
  # Tuning parameters
  if(a_like >= a_disc){
    print('WARNING: a_like must be smaller than a_disc.')
    print('Code will break because you are requesting to retain %f percent of particles with low discrepancy')
    print('and then to retain percent of particles with low discrepany and high likelihood')
  }

  # Defining special arguments ####

  #PRIOR - define the prior distribution for the modelling scenario
  # This function returns "args" a list that contains:
  # lower - an array of lower parameter bounds
  # upper - an array of upper parameter bounds
  # prior_sampler - a function that samples a parameter set from the prior
  # trans_f - a function that transforms the prior into a more manageable space
  # trans_f_inv - a function that back transforms the prior into its original space
  # prior_pdf - the probability density of the transform space
  # param_labels - the names of the parameters in latex style
  # n_params - the number of model parameters
  # simulate_constraint - a function that simulates the model outputs needed to assess if constraints are met
  # simulate_data - a function that simulates the model outputs needed to assess likelihood of producing data
  # simulate_plots - a function that simulates the model outputs needed to plot results
  args <- define_args()

  # #all other relevant parameters
  args$include_expert_constraints <- include_expert_constraints
  args$include_data_constraints <- include_data_constraints
  args$discrepancy_final <- discrepancy_final

  # args$n_particles <- n_particles
  # args$mcmc_trials <- mcmc_trials
  # args$a_disc <- a_disc
  # args$a_like <- a_like
  # args$c <- c
  # args$p_acc_min <- p_acc_min
  # args$n_cores <- n_cores

  ### RUN SMC ####
  outputs <- SMCfeatures::SMC_combined(args,
                                       n_particles,
                                       mcmc_trials,
                                       discrepancy_final,
                                       a_disc,
                                       a_like,
                                       c,
                                       p_acc_min,
                                       n_cores)

  return(outputs)
}

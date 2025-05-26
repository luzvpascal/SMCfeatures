#' @title SMCfeatures
#' @description
#' Main function for any parameterisation of expert elicited constraints
#' @param lower a vector of lower parameter bounds. If using the Gaussian log-likelihood function (\link[SMCfeatures]{log_likelihood_Gaussian}), the last number of this array describes the standard deviation.
#' @param upper a vector of upper parameter bounds. If using the Gaussian log-likelihood function (\link[SMCfeatures]{log_likelihood_Gaussian}), the last number of this array describes the standard deviation.
#' @param simulate_model a function that simulates the model outputs needed to assess if constraints are met
#' @param calculate_discrepancy a function that calculates discrepancy score. See  \link[SMCfeatures]{discrepancy_logistic_growth} for an example.
#' @param simulate_data a function that simulates the model outputs needed to assess likelihood of producing data
#' @param input_data a vector of input data. Default NA.
#' @param output_data a vector of output data. Default NA.
#' @param calculate_log_likelihood log likelihood function. Default \link[SMCfeatures]{log_likelihood_Gaussian}.
#' @param include_expert_constraints boolean indicating whether expert ellicited constraints should be included. Default TRUE
#' @param include_data_constraints boolean indicating whether data should be included. Default TRUE. If input_data or output_data is NA, then include_data_constraints is set to FALSE
#' @param prior_sampler sampling function that generates random vectors from the joint prior distribution. Default \link[SMCfeatures]{uniform_sampler}.
#' @param trans_f transform of prior parameter space to ensure unbounded support for MCMC sampling. Default \link[SMCfeatures]{uniform_transform}
#' @param trans_finv inverse of trans_f function. Default \link[SMCfeatures]{uniform_transform_inverse}
#' @param prior_pdf joint probability density function. Default \link[SMCfeatures]{uniform_pdf_transformed}
#' @param n_particles Number of desired ensemble members. Default to 10000
#' @param mcmc_trials number of MCMC steps to try before selecting appropriate number. Default 10
#' @param discrepancy_final target discrepancy threshold. Default 0. If zero, p_acc_min is used to determine stopping criteria.
#' @param a_disc tuning parameter for adaptive selection of discrepancy threshold sequence. Default 0.6
#' @param a_like tuning parameter for adaptive selection of likelihood ESS sequence. Default 0.3
#' @param c tuning parameter for choosing the number of MCMC iterations in move step. Default 0.01
#' @param p_acc_min minimum acceptable acceptance rate in the MCMC interations before exit. Default 0.0001
#' @param n_cores Number of cores desired to be used for sampling. Default set to 1 core (sequential sampling).
#' @return vector of transformed parameters
#' @export
SMCfeatures <- function(upper,
                        lower,
                        simulate_model,
                        calculate_discrepancy,
                        simulate_data = SMCfeatures::model_on_data,
                        input_data=NA,
                        output_data=NA,
                        calculate_log_likelihood=SMCfeatures::log_likelihood_Gaussian,
                        include_expert_constraints=TRUE,
                        include_data_constraints=TRUE,
                        sampler=SMCfeatures::uniform_sampler,
                        trans_f=SMCfeatures::uniform_transform,
                        trans_finv=SMCfeatures::uniform_transform_inverse,
                        prior_pdf=SMCfeatures::uniform_pdf_transformed,
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
  args <- list()

  ## global arguments ####
  args$upper <- upper
  args$lower <- lower
  args$n_params <- length(args$upper)

  # prior functions
  args$sampler <- sampler
  args$trans_f <- trans_f
  args$trans_finv <- trans_finv
  args$prior_pdf <- prior_pdf

  #model simulation
  args$simulate_model <- simulate_model #inputs (parameters,inputs), returns [simulation]
  args$simulate_data <- SMCfeatures::model_on_data #inputs (parameters,args,simulation), returns [simulation]

  #discrepancy
  args$calculate_discrepancy <- calculate_discrepancy

  #data and log likelihood
  args$input_data <- input_data#times
  args$output_data <- output_data#coral cover observation
  args$calculate_log_likelihood  <- calculate_log_likelihood


  # #all other relevant parameters
  args$include_expert_constraints <- include_expert_constraints
  if (is.na(input_data[1])){
    print("WARNING: No data input")
    print("The argument include_data_constraints will be set to FALSE")
    args$include_data_constraints <- FALSE
  }
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
                                       a_disc,
                                       a_like,
                                       c,
                                       p_acc_min,
                                       n_cores)

  return(outputs)
}

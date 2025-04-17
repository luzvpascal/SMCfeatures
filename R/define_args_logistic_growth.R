#' @title Arguments for SMCfeatures
#' @description
#' Extract arguments necessary to run SMCfeatures
#' @return A list of arguments defining the problem.
#' @export
define_args_logistic_growth <- function(){
  #define global arguments
  args <- list()
  args$parameter_labels <- c("r","K","y0","sigma")
  args$upper <- c(0.5, 80, 5, 20)
  args$lower <- c(0, 60, 0, 0)
  args$n_params <- length(args$upper)

  # prior functions
  args$sampler <- SMCfeatures::uniform_sampler
  args$trans_f <- SMCfeatures::uniform_transform
  args$trans_finv <- SMCfeatures::uniform_transform_inverse
  args$prior_pdf <- SMCfeatures::uniform_pdf_transformed

  #model simulation
  args$simulate_constraint <- SMCfeatures::model_logistic_growth_constraint #inputs (parameters,args), returns [simulation]
  args$simulate_data <- SMCfeatures::model_logistic_growth_data #inputs (parameters,args,simulation), returns [simulation]
  args$simulate_plots <- SMCfeatures::model_logistic_growth_plots #inputs (parameters,args,simulation), returns [simulation]

  #data and log likelihood
  args$input_data <- c(1, 3, 10)#times
  args$output_data <- c(4, 4, 10)#coral cover observation
  args$calculate_log_likelihood  <- SMCfeatures::log_likelihood_logistic_growth

  #discrepancy
  args$calculate_discrepancy <- SMCfeatures::discrepancy_logistic_growth

  return(args)
}

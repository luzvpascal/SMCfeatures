# SMCfeatures

# Installation

```r
packages_to_install <- c("parallel", "doParallel", "foreach","dplyr", "foreach", "ggplot2", "magrittr", "MASS", "parallel", "parallelly","stats","tidyr")

# Install packages if not already installed
install_if_not_installed <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
}

# Apply the function to install packages
invisible(lapply(packages_to_install, install_if_not_installed))

#Install from github
library(devtools)
install_github("luzvpascal/SMCfeatures",
               host = "https://api.github.com")
```

# Getting started with a simple example

## Mathematical model
We consider a simple logistic growth example, representing the evolution of coral growth after a disturbance. The logistic growth model has the following shape:

```math
\frac{\rm{d}y}{\rm{d}t} = r.y(t)\left(1-\frac{y(t)}{K}\right),\quad y(0) = y_0,
```
where $y (t)$ represents the coral cover at year $t$ (% area), $t$ is the time (years), $r$ is the instrinsic growth rate ($year^{-1}$), $K$ is the carrying capacity (% area) and $y_0$ is the initial coral cover (% area).

The analytical solution of this ODE is :
```math
y(t) = \frac{Ky_0}{y_0+(K-y_0)e^{-rt}}.
```

Three parameters are uncertain: $r$, $K$ and $y_0$.

## Step 1: Defining the model function into R

The first step of our approach is to define a function that simulates the model for any parameters $r$, $K$ and $y_0$. For the logistic growth example, this function is already implemented in the package, as `model_logistic_growth`.
```r
#' @title Model of logistic growth case study
#' @description
#' Simulate outputs of logistic growth model
#' @param parameters vector of parameters (growth rate, capacity, initial abundance)
#' @param input input data. here it is a vector of time steps of interest
#' @return
#' y_data: output of simulation on input data as a vector
#' @export
model_logistic_growth <- function(parameters,
                                  input){
  # likelihood needs data at all time points
  #Inputs parameters and outputs a simulation

  # Extract information
  r <- parameters[1]
  K <- parameters[2]
  y0 <- parameters[3]

  # Simulation details
  y_data <- (K*y0)/(y0+(K-y0)*exp(-input*r))
  return(y_data)
}
```

This 
```

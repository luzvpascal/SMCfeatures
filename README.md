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

#Getting started with a simple example

We consider a simple logistic growth example, representing the evolution of coral growth after a disturbance. The logistic growth model has the following shape:

```math
\frac{ds}{dt} = r.s\left(1-\frac{s}{K}\right)
```

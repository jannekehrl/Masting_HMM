## Simple R script to test Rstan installation
## E.g. for new computer

##First steps:
# Install R and Rstudio
# Install matching Rtools 
# Install packages from stan universe directly (not within Rstudio)
    # Note if you did install within Rstudio, remove those first

##Here is command to install packages
install.packages(
  c("rstan", "StanHeaders"),
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)

##If you need to know what versions of stan packages?
packageVersion("rstan")
packageVersion("StanHeaders")
packageVersion("Rcpp")

##Does R find Rtools? And what version?
pkgbuild::has_rtools()

#check which g++ and make version
system("g++ --version") #might return G__~1.EXE (GCC) version #, this is OK
Sys.which("make") #did it find make? also check path for rtools version

#check whether a program can be compiled in C++
pkgbuild::has_build_tools(debug = TRUE) #you are in trouble if it returns false

#Show the makevars.win. You may need this (although I didn't for 2nd setup)
writeLines(readLines(file.path(Sys.getenv("HOME"), ".R", "Makevars.win")))

#Load Rstan
library(rstan)

# Optional: Configure RStan to use multiple cores and save compiled models
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Run the built-in example model (from stan) to test installation
# If installation correct, you will get many lines of text including compiling
# and ending in reports on the four chains (ieration, sampling, warmup and time)
example(stan_model, package = "rstan", run.dontrun = TRUE)
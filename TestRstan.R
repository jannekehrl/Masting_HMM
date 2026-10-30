## Simple R script to test Rstan installation

##What versions of stan packages?
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

#Show the makevars
writeLines(readLines(file.path(Sys.getenv("HOME"), ".R", "Makevars.win")))

#Load Rstan
library(rstan)

#Run a simple model
stancode <- 'data {real y_mean;} parameters {real y;} model {y ~ normal(y_mean,1);}'
mod <- stan_model(model_code = stancode, verbose = TRUE)


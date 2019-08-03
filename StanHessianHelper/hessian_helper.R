library(rstan)
library(Rcpp)

load_model <- function(model_file, stan_data) {
  hamiltonian_env = new.env()
  hamiltonian_env$data = stan_data

  data_file = tempfile()
  if(length(stan_data) == 0) {
    file.create(data_file)
  } else {
    stan_rdump(names(stan_data), file = data_file, envir = list2env(stan_data))
  }

  # create linear_regression_model.cpp from Stan file
  if(!file.exists(model_file)) {
    stop("model '", model_file,"' does not exist")
  }
  code = stanc(model_file, model_name = 'log_density', obfuscate_model_name = FALSE)
  write(code$cppcode, "./StanHessianHelper/model_log_density.hpp")

  # compile with RCPP
  # system(paste0('cp ', paste0(getSrcDirectory(function(dummy) {dummy})), 'helper.cpp helper.cpp'))
  #sourceCpp(paste0(system.file(package = "RNUTS"), "/StanHessianHelper/helper.cpp"))
  sourceCpp("./StanHessianHelper/helper.cpp")
  set_data(data_file)
  file.remove(data_file)

  hamiltonian_env$U = function(q) -jacobian(q)$u
  hamiltonian_env$GradU = function(q) -jacobian(q)$jac
  hamiltonian_env$HessU = function(q) -hessian(q)$hess
  hamiltonian_env$HessVecProd = function(q, vec) -hessian_vector(q, vec)$hessv

  return(hamiltonian_env)
}

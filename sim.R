###############################################################################
#
# Here we perform a small simulation study to compare the computational time 
# needed by the function `rpkbd()` in the `Directional` package and the function
# `rpkb()` from the `QuadratiK` package for generating random samples from the 
# Poisson kernel-based distribution (PKBD)
# 
# #############################################################################

# Load required packages
library("snowfall") # for parallel computing
library("QuadratiK")
library("Directional")

## Divide the combinations of parameters in different files to be run 
## simultaneously on different cores
## 
replace_special_chars <- function(x) {
  pattern <- "[^A-Za-z0-9.,]"
  cleaned_str <- gsub(pattern, "", x)
  return(cleaned_str)
}

args <- commandArgs(trailingOnly = TRUE)
args <- strsplit(args, ",")
args <- unlist(lapply(args, replace_special_chars))
print(args)

## Function for: 
##  - generate random samples with all the considered methods;
##  - perform parameters estimation on the data sets generated with `rpkbd()` 
##    and `rpkb()` with method `rejacg`;
##  - return execution times for random generation and estimation process, 
##    together with the obtained estimates.
OneSimulation<-function(size, d, rho){

  mu_vec <- c(1, rep(0, d-1))
  time_dir <- system.time({
    x_dir <- rpkbd(n = size, mu = mu_vec, rho = rho)
  }) 
  
  time_acg <- system.time({
    x_q_acg <- rpkb(n = size, mu = mu_vec, rho = rho, method = "rejacg")
  })
  
  if(rho < 0.9){
    time_vmf <- system.time({
      x_q_vmf <- rpkb(n = size, mu = mu_vec, rho = rho, method = "rejvmf")
    })
  }
  
  time_psaw <- system.time({
    x_q_psaw <- rpkb(n = size, mu = mu_vec, rho = rho, method = "rejpsaw")
  })
  
  time_est_dir1 <- system.time({
    est_dir1 <- pkbd.mle(x = x_dir)
  })
  
  time_est_dir2 <- system.time({
    est_dir2 <- pkbd.mle(x = x_q_acg$x)
  })
  
  time_est_q1 <- system.time({
    est_q1 <- pkbc(dat = x_dir, nClust = 1)
  })
  
  time_est_q2 <- system.time({
    est_q2 <- pkbc(dat = x_q_acg$x, nClust = 1)
  })
  
  if (rho < 0.9){
    res <- c(time_dir[3], time_acg[3], time_vmf[3], time_psaw[3],
             est_dir1$mu, est_dir1$rho, est_dir2$mu, est_dir2$rho,
             est_q1@res_k[[1]]$params$mu, est_q1@res_k[[1]]$params$rho,
             est_q2@res_k[[1]]$params$mu, est_q2@res_k[[1]]$params$rho,
             time_est_dir1[3], time_est_dir2[3], time_est_q1[3], time_est_q2[3])
  } else {
    res <- c(time_dir[3], time_acg[3], 0, time_psaw[3],
             est_dir1$mu, est_dir1$rho, est_dir2$mu, est_dir2$rho,
             est_q1@res_k[[1]]$params$mu, est_q1@res_k[[1]]$params$rho,
             est_q2@res_k[[1]]$params$mu, est_q2@res_k[[1]]$params$rho,
             time_est_dir1[3], time_est_dir2[3], time_est_q1[3], time_est_q2[3])
  }
  
  return(res)
}

## Function for: 
##  - setting the parameters with correct format
##  - perform random sample generation
##  - write obtained results
simula.internal <- function(line_num, nrep, size, d, rho){

  line_num <- as.numeric(line_num)
  nrep <- as.numeric(nrep)
  size <- as.numeric(size)
  d <- as.numeric(d)
  rho <- as.numeric(rho)
  
  set.seed(1234 + nrep)
  res <- OneSimulation(size=size, d=d, rho=rho)

  res <- c(line_num, nrep, size, d, rho, res)
  res <- paste(res, collapse = " ")
  res <- paste("\n", res, sep = " ")
  write(res, file = paste("results/RESULTS_pkbd_d",d,"_", line_num %/% 33, ".txt", sep = ""), append = TRUE)
}

simula.internal(args[1], args[2], args[3], args[4], args[5])

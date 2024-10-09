### Set parameters for Simulation ###
N.simulations <- 100

dim_values <- c(2,3,5,10)
P_n_obs <- c(100,500, 1000, 2000, 5000, 10000) 
rho_values <- c(0.6, 0.7, 0.8, 0.9, 0.95)

params <- as.matrix(expand.grid(1:N.simulations,P_n_obs, dim_values, rho_values))
params <- as.data.frame(params)

# Run the following from the secoond time to check only the missing combinations
par_done <- read.table("results/res_done.txt")
par_done <- par_done[,2:5]
colnames(par_done) <- colnames(params)
par_done <- split(par_done, seq(nrow(par_done)))
params <- split(params, seq(nrow(params)))

params <- setdiff(params, par_done)
params <- matrix(unlist(params), ncol=4,byrow = TRUE)

write.table(params, 'parameters.txt', row.names = TRUE, col.names = FALSE, sep = ',', quote = F)

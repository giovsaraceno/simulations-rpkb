##########################################################################
##
## Code for merge the obtained results for random generation from PKBD
##
##########################################################################

rm(list=ls())
# Set working directory
setwd("/Users/sarac/Dropbox/OnGoing/rpkb_jmva/")

# set the array used in the simulations
# merge all the results in a list with an element for each dimension considered
# (since the length of the estimate of the location vector changes accordingly)
arrays <- seq(0,100) 
d_values <- c(2, 3, 5, 10)
results <- list()
for(j in seq_along(d_values)){
   res <- matrix(ncol = 17 + 4*d_values[j])
   for(i in arrays){
      file.name <- paste0("results/RESULTS_pkbd_d",d_values[j],"_",i,".txt")
      if(file.exists(file.name)){
         dat_i <- read.table(paste0("results/RESULTS_pkbd_d",d_values[j],"_",i,".txt"))
         res <- rbind(res, dat_i)
      } else{
         cat(paste0("Array ", i, " is missing","\n"))
      }
   }
   res <- res[-1,]
   results[[j]] <- res
}

# Save data
save(results, file="res_rpkbd.RData")

# Save only the combination of parameters executed in order to verify if there 
# are missing combinations and then run them.
res_check <- matrix(ncol = 5)
for(j in seq_along(d_values)){
   res_check <- rbind(res_check, results[[j]][,1:5])
}
res_check <- res_check[-1,]
res_check <- res_check[order(res_check$V2,res_check$V3,res_check$V4, res_check$V5),]
which(duplicated(res_check[,2:5]))
#res_tot <- res_tot[-which(duplicated(res_tot[,2:5])),]
res_check[1:30,]
write.table(res_check, "results/res_done.txt")



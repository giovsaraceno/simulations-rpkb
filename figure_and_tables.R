############################################################################
##
## Investigation of the performance of the functions for generating from the 
## Poisson kernel-based distribution using QuadratiK and Directional packages
##
############################################################################

rm(list=ls())

# Load required packages
library(ggplot2)
library(tidyverse)

# Load results
load("res_rpkbd.RData")

############################################################################
## Investigate time for generating samples from pkbd
## Select only the execution time for the random generation part 
res_time <- matrix(ncol = 8)
colnames(res_time) <- colnames(results[[1]][,2:9])
for(j in 1:4){
  res_time <- rbind(res_time, results[[j]][,2:9])
}
res_time <- as.data.frame(res_time[-1,])
colnames(res_time) <- c("rep", "size", "d", "rho","rpkbd", "rejacg", "rejvmf", "rejpsaw")
str(res_time)

# Compute the average computational time over the 100 replications
# and re-organize the data set such that there is the variable "Time" with the
# obtained values and the variable "Method" reporting the corresponding method
# used for the generation
res_mean <- res_time %>%
  group_by(size, d, rho) %>%
  summarise(
    rpkbd = mean(rpkbd),
    rejacg = mean(rejacg),
    rejvmf = mean(rejvmf),
    rejpsaw = mean(rejpsaw),
  ) %>%
  ungroup() %>%
  gather(key = "Method", value = "Time", rpkbd:rejpsaw) 

res_mean$size <- factor(res_mean$size)
res_mean$d <- factor(res_mean$d)
res_mean$Method <- factor(res_mean$Method)


## Plots for fixed dimension
for(d in unique(res_mean$d)){
  pl <- ggplot(res_mean[which(res_mean$d==d),], aes(x=rho, y=Time, color=Method))+
    geom_line(aes(color = Method, linetype = size), linewidth = 1.1, alpha=.6) +
    labs(y="Time", x = expression(rho))+
    theme_minimal()+
    theme_light()+
    theme(legend.title = element_text(size=16),
          legend.text = element_text(size = 18),
          plot.title = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          strip.text = element_text(size = 14)) +
    scale_color_brewer(palette='Set1')
  ggsave(pl,filename = paste0("figures/time_rho_d",d,".pdf"),width = 10,height = 6)
}

# The 'rejvmf' method requires an elevated computational time.
# We exclude this method from the following plots for sake of a better 
# representation.
res_mean <- res_time %>%
  group_by(size, d, rho) %>%
  summarise(
    rpkbd = mean(rpkbd),
    rejacg = mean(rejacg),
    rejpsaw = mean(rejpsaw),
  ) %>%
  ungroup() %>%
  gather(key = "Method", value = "Time", c(rpkbd,rejacg,rejpsaw)) 

res_mean$size <- factor(res_mean$size)
res_mean$d <- factor(res_mean$d)
res_mean$Method <- factor(res_mean$Method)

## Plot for lower sample sizes 
pl <- ggplot(res_mean[which(res_mean$size %in% c(500,1000,2000) & res_mean$d %in% c(2,3,10)),], aes(x=rho, y=Time, color=Method))+
  geom_line(aes(color = Method), linewidth = 1.1, alpha=.6) +
  facet_grid( d ~ size, scales = "free")+
  labs(y="Time", x = expression(rho))+
  theme_minimal()+
  theme_light()+
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  scale_color_brewer(palette='Set1')
ggsave(pl,filename = paste0("figures/time_rho_low.pdf"),width = 6,height = 6)

## Re-do figure without legend to insert it into the manuscript
pl <- ggplot(res_mean[which(res_mean$size %in% c(500,1000,2000) & res_mean$d %in% c(2,3,10)),], aes(x=rho, y=Time, color=Method))+
  geom_line(aes(color = Method), linewidth = 1.1, alpha=.6) +
  facet_grid( d ~ size, scales = "free")+
  labs(y="Time", x = expression(rho))+
  theme_minimal()+
  theme_light()+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  scale_color_brewer(palette='Set1')
ggsave(pl,filename = paste0("figures/time_rho_low_noLeg.pdf"),width = 6,height = 6)

# Plot for larger sample sizes
pl <- ggplot(res_mean[which(res_mean$size %in% c(5000,10000) & res_mean$d %in% c(2,3,10)),], aes(x=rho, y=Time, color=Method))+
  geom_line(aes(color = Method), linewidth = 1.1, alpha=.6) +
  facet_grid( d ~ size, scales = "free")+
  labs(y="Time", x = expression(rho))+
  theme_minimal()+
  theme_light()+
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  scale_color_brewer(palette='Set1')
ggsave(pl,filename = paste0("figures/time_rho_high.pdf"),width = 6,height = 6)

############################################################################
## We now investigate the execution time for estimating the parameters of PKBD
## using the function `pkbd.mle()` from the `Directional` package and the 
## function `pkbc` from `QuadratiK`. Notice that  `pkbc` is originally 
## implemented for estimating the model parameters of a mixture of PKBDs for 
## clustering. 
est_time <- matrix(ncol = 8)
for(j in 1:4){
  n_r <- ncol(results[[j]])
  est_time <- rbind(est_time, as.matrix(results[[j]][,c(2:5,(n_r-3):n_r)]))
}
est_time <- as.data.frame(est_time[-1,])
colnames(est_time) <- c("rep", "size", "d", "rho","rpkbd_1","rpkbd_2", "rejacg_1", "rejacg_2")
str(est_time)

## Legend: 
##  - 1: data are generated using `rpkbd`
##  - 2: data are generated using `rpkb` with method `rejacg`

est_time_mean <- est_time %>%
  group_by(size, d, rho) %>%
  summarise(
    Directional_1 = mean(rpkbd_1),
    QuadratiK_1 = mean(rejacg_1),
    Directional_2 = mean(rpkbd_2),
    QuadratiK_2 = mean(rejacg_2),
  ) %>%
  ungroup() %>%
  gather(key = "Method", value = "Time", Directional_1:QuadratiK_2) 

est_time_mean$size <- factor(est_time_mean$size)
est_time_mean$d <- factor(est_time_mean$d)
est_time_mean$Method <- factor(est_time_mean$Method)

## Time (for estimation) vs rho for lower sample sizes 
pl <- ggplot(est_time_mean[which(est_time_mean$d %in% c(2,3,10) & est_time_mean$size %in% c(500,1000,2000)),], aes(x=rho, y=Time, color=Method))+
geom_line(aes(color = Method), linewidth = 1.1, alpha=.6) +
facet_grid(d ~ size, scales = "free")+
labs(y="Time", x = expression(rho))+
theme_minimal()+
theme_light()+
theme(legend.title = element_text(size=16),
      legend.text = element_text(size = 18),
      plot.title = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      strip.text = element_text(size = 14)) +
scale_color_brewer(palette='Set1')
ggsave(pl,filename = paste0("figures/time_rho_est_low_both.pdf"),width = 6,height = 6)

pl <- ggplot(est_time_mean[which(est_time_mean$d %in% c(2,3,10) & est_time_mean$size %in% c(5000,10000)),], aes(x=rho, y=Time, color=Method))+
  geom_line(aes(color = Method), linewidth = 1.1, alpha=.6) +
  facet_grid(d ~ size, scales = "free")+
  labs(y="Time", x = expression(rho))+
  theme_minimal()+
  theme_light()+
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  scale_color_brewer(palette='Set1')
ggsave(pl,filename = paste0("figures/time_rho_est_high_both.pdf"),width = 6,height = 6)

## The computational time does not change if the data are generated with the 
## `rpkbd` function or by `rpkb`. We then select only one case.
est_time_mean <- est_time %>%
  group_by(size, d, rho) %>%
  summarise(
    Directional = mean(rpkbd_2),
    QuadratiK = mean(rejacg_2),
  ) %>%
  ungroup() %>%
  gather(key = "Method", value = "Time", Directional:QuadratiK) 

est_time_mean$d <- factor(est_time_mean$d)
est_time_mean$Method <- factor(est_time_mean$Method)

## Plots 
pl <-  ggplot(est_time_mean[which(est_time_mean$d %in% c(2,3,10) & est_time_mean$size %in% c(500,1000,2000)),], aes(x=rho, y=Time, color=Method))+
  geom_line(aes(color = Method), linewidth = 1.1, alpha=.6) +
  facet_grid(d ~ size, scales = "free")+
  labs(y="Time", x = expression(rho))+
  theme_minimal()+
  theme_light()+
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  scale_color_brewer(palette='Set1')
ggsave(pl,filename = paste0("figures/time_rho_est_low.pdf"),width = 6,height = 6)

pl <- ggplot(est_time_mean[which(est_time_mean$d %in% c(2,3,10) & est_time_mean$size %in% c(5000,10000)),], aes(x=rho, y=Time, color=Method))+
  geom_line(aes(color = Method), linewidth = 1.1, alpha=.6) +
  facet_grid(d ~ size, scales = "free")+
  labs(y="Time", x = expression(rho))+
  theme_minimal()+
  theme_light()+
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  scale_color_brewer(palette='Set1')
ggsave(pl,filename = paste0("figures/time_rho_est_high.pdf"),width = 6,height = 6)

## In both cases, the function `pkbc` takes more time in estimating the 
## parameters. This is more evident for increasing sample size. We create an 
## additional plot of Time vs n. 
## 
est_time_mean$rho <- factor(est_time_mean$rho)
pl <- ggplot(est_time_mean[which(est_time_mean$d %in% c(2,3,10) & est_time_mean$rho %in% c(0.7,0.8,0.9,0.95)),], aes(x=size, y=Time, color=Method))+
  geom_line(aes(color = Method), linewidth = 1.1, alpha=.6) +
  facet_grid(d ~ rho, scales = "free")+
  labs(y="Time", x = "n")+
  theme_minimal()+
  theme_light()+
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  scale_color_brewer(palette='Set1')
ggsave(pl,filename = paste0("figures/time_n_est.pdf"),width = 6,height = 6)
### The clustering algorithm used for estimating the parameters of 
### a single PKBD is slower than the MLE implemented in Directional
### 

###############################################################################
## We now investigate the performance of the clustering algorithm and the 
## MLE in estimating the model parameters 

## - Select the estimates of rho and mu 
## - Compute the difference between the estimated rho and the "true" value 
##   used for generating the data set.
## - Compute the MSE (Mean squared error) between the estimated mu and the 
##   value used for generating the data set.
d_values <- c(2,3,5,10)
est_res <- matrix(ncol = 16)
for(j in 1:4){
  res <- results[[j]]
  nr <- ncol(res)
  d <- d_values[j]
  rho <- res[,5]
  mu_vec <- c(1, rep(0,d-1))
  
  mu_d1 <- res[,10:(10+d-1)]
  rho_d1 <- res[,10+d]
  mu_d2 <- res[,(10+d+1):(10+2*d)]
  rho_d2 <- res[,10+2*d+1]
  
  mu_q1 <- res[,(10+2*d+2):(10+3*d+1)]
  rho_q1 <- res[,10+3*d+2]
  mu_q2 <- res[,(10+3*d+3):(10+4*d+2)]
  rho_q2 <- res[,10+4*d+3]
  
  est_res <- rbind(est_res, 
                    cbind(as.matrix(res[,c(2:5)]),
                          rho_d1 - rho,
                          rho_d2 - rho,
                          rho_q1 - rho,
                          rho_q2 - rho,
                          sqrt(rowSums(t(apply(mu_d1, 1, function(x) x-mu_vec)^2))),
                          sqrt(rowSums(t(apply(mu_d2, 1, function(x) x-mu_vec)^2))),
                          sqrt(rowSums(t(apply(mu_q1, 1, function(x) x-mu_vec)^2))),
                          sqrt(rowSums(t(apply(mu_q2, 1, function(x) x-mu_vec)^2))),
                          apply(mu_d1, 1, function(x) acos(x%*%mu_vec)),
                          apply(mu_d2, 1, function(x) acos(x%*%mu_vec)),
                          apply(mu_q1, 1, function(x) acos(x%*%mu_vec)),
                          apply(mu_q2, 1, function(x) acos(x%*%mu_vec))))
}

est_res <- as.data.frame(est_res[-1,])
colnames(est_res) <- c("rep", "size", "d", "rho","rhod_1","rhod_2", "rhoq_1", "rhoq_2","MSEd_1","MSEd_2", "MSEq_1", "MSEq_2","thetad_1","thetad_2", "thetaq_1", "thetaq_2")
str(est_res)

# Note: the labels 1 and 2 are as before!

### Estimation of rho
##
est_res_mean <- est_res %>%
  group_by(size, d, rho) %>%
  summarise(
    Directional_1 = mean(rhod_1),
    QuadratiK_1 = mean(rhoq_1),
    Directional_2 = mean(rhod_2),
    QuadratiK_2 = mean(rhoq_2)
  ) %>%
  ungroup() %>%
  gather(key = "Method", value = "delta_rho", Directional_1:QuadratiK_2) 

est_res_mean$d <- factor(est_res_mean$d)
est_res_mean$Method <- factor(est_res_mean$Method)

pl <- ggplot(est_res_mean[which(est_res_mean$d %in% c(2,3,5,10) & est_res_mean$size %in% c(500,1000,2000)),], aes(x=rho, y=delta_rho, color=Method))+
  geom_line(aes(color = Method), linewidth = 1.1, alpha=.6) +
  facet_grid(d ~ size, scales = "free")+
  labs(y=expression(hat(rho) - rho), x = expression(rho))+
  theme_minimal()+
  theme_light()+
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  scale_color_brewer(palette='Set1')
ggsave(pl,filename = paste0("figures/rho_hat_low_both.pdf"),width = 6,height = 6)

pl <- ggplot(est_res_mean[which(est_res_mean$d %in% c(2,3,5,10) & est_res_mean$size %in% c(5000,10000)),], aes(x=rho, y=delta_rho, color=Method))+
  geom_line(aes(color = Method), linewidth = 1.1, alpha=.6) +
  facet_grid(d ~ size, scales = "free")+
  #labs(y=expression(hat(rho) - rho), x = expression(rho))+
  theme_minimal()+
  theme_light()+
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  scale_color_brewer(palette='Set1')
ggsave(pl,filename = paste0("figures/rho_hat_high_both.pdf"),width = 6,height = 6)

## Results depend if data points are generated using the Directional package
## or the QuadratiK package. However the performance in estimating rho is 
## the same.
## In particular, if data are generated using "rpkbd" the difference between 
## the estimated rho and the "true" value has greater variability, while if 
## the function "rpkb" is used the difference has less variation and it is 
## closer to zero.


### Investigation of estimation of mu using MSE
### 
est_res_mean <- est_res %>%
  group_by(size, d, rho) %>%
  summarise(
    Directional_1 = mean(MSEd_1),
    QuadratiK_1 = mean(MSEq_1),
    Directional_2 = mean(MSEd_2),
    QuadratiK_2 = mean(MSEq_2)
  ) %>%
  ungroup() %>%
  gather(key = "Method", value = "MSE", Directional_1:QuadratiK_2) 

est_res_mean$d <- factor(est_res_mean$d)
est_res_mean$Method <- factor(est_res_mean$Method)

## Plots
pl <- ggplot(est_res_mean[which(est_res_mean$d %in% c(2,3,5,10) & est_res_mean$size %in% c(500,1000,2000)),], aes(x=rho, y=MSE, color=Method))+
  geom_line(aes(color = Method), linewidth = 1.1, alpha=.6) +
  facet_grid(d ~ size, scales = "free")+
  #labs(y=expression(hat(rho) - rho), x = expression(rho))+
  theme_minimal()+
  theme_light()+
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  scale_color_brewer(palette='Set1')
ggsave(pl,filename = paste0("figures/mse_rho_low_both.pdf"),width = 6,height = 6)

pl <- ggplot(est_res_mean[which(est_res_mean$d %in% c(2,3,5,10) & est_res_mean$size %in% c(5000,10000)),], aes(x=rho, y=MSE, color=Method))+
  geom_line(aes(color = Method), linewidth = 1.1, alpha=.6) +
  facet_grid(d ~ size, scales = "free")+
  #labs(y=expression(hat(rho) - rho), x = expression(rho))+
  theme_minimal()+
  theme_light()+
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  scale_color_brewer(palette='Set1')
ggsave(pl,filename = paste0("figures/mse_rho_high_both.pdf"),width = 6,height = 6)
## For estimating the location vector mu, the considered functions perform 
## almost identically, independently of the method used for generate the data

### For curiosity we compute also the angle between vectors (mu_hat and mu) as
### performance measure for the estimation of mu.
est_res_mean <- est_res %>%
  group_by(size, d, rho) %>%
  summarise(
    Directional_1 = mean(thetad_1),
    QuadratiK_1 = mean(thetaq_1),
    Directional_2 = mean(thetad_2),
    QuadratiK_2 = mean(thetaq_2)
  ) %>%
  ungroup() %>%
  gather(key = "Method", value = "theta", Directional_1:QuadratiK_2) 

est_res_mean$d <- factor(est_res_mean$d)
est_res_mean$Method <- factor(est_res_mean$Method)


##
ggplot(est_res_mean[which(est_res_mean$d %in% c(2,3,5,10) & est_res_mean$size %in% c(500,1000,2000)),], aes(x=rho, y=theta, color=Method))+
  geom_line(aes(color = Method), linewidth = 1.1, alpha=.6) +
  facet_grid(d ~ size, scales = "free")+
  #labs(y=expression(hat(rho) - rho), x = expression(rho))+
  theme_minimal()+
  theme_light()+
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  scale_color_brewer(palette='Set1')

ggplot(est_res_mean[which(est_res_mean$d %in% c(2,3,5,10) & est_res_mean$size %in% c(5000,10000)),], aes(x=rho, y=theta, color=Method))+
  geom_line(aes(color = Method), linewidth = 1.1, alpha=.6) +
  facet_grid(d ~ size, scales = "free")+
  #labs(y=expression(hat(rho) - rho), x = expression(rho))+
  theme_minimal()+
  theme_light()+
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  scale_color_brewer(palette='Set1')
## We obtain identical results.
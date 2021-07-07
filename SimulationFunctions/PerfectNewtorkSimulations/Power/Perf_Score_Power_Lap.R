###############################
##
## Project: MetaboGuru
##
## Purpose: Power of Score Tests
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2020-11-11
##
## ---------------------------
## Notes:
##   
##
## ---------------------------

# Helpful Functions -------------------------------------------------------

`%nin%` <- Negate(`%in%`)
source('PerfectNetworkSimFunctions.R')

nsim <- 10000  ## number of simulations
n <- 160 ## sample size
# mX <- matrix(1, n)
# b0 <- 0.2644 ## intercept term
sd.y <- 1.3688 ## standard deviation of Y
tau <- 1 ## Tuning parameter for regularization kernel of normalized laplacian
set.seed(4)
X <- data.frame(X1 = rep(0:1, each = n/2),
                X2 = runif(n, 0, 5))
b0 <- c(0.2644, 0.5, 0.25)
H0.form <- formula(Y~X1+X2)

# Same Size ---------------------------------------------------------------

# * 15 --------------------------------------------------------------------
p <- 15 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

perf_davie_L_ss_pw_15 <- plyr::ldply(graph.list, Perf_Davie_SameSize, include.network = "L",
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(perf_davie_L_ss_pw_15$pos_def); unique(perf_davie_L_ss_pw_15$Power)

# * 30 --------------------------------------------------------------------
p <- 30 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

perf_davie_L_ss_pw_30 <- plyr::ldply(graph.list, Perf_Davie_SameSize, include.network = "L",
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(perf_davie_L_ss_pw_30$pos_def); unique(perf_davie_L_ss_pw_30$Power)

# * 45 --------------------------------------------------------------------
p <- 45 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

perf_davie_L_ss_pw_45 <- plyr::ldply(graph.list, Perf_Davie_SameSize, include.network = "L",
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(perf_davie_L_ss_pw_45$pos_def); unique(perf_davie_L_ss_pw_45$Power)

# Small Graph -------------------------------------------------------------

# * 15 --------------------------------------------------------------------
p <- 15 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

perf_davie_L_sm_pw_15 <- plyr::ldply(graph.list, Perf_Davie_SmallGraph, include.network = "L",
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(perf_davie_L_sm_pw_15$pos_def); unique(perf_davie_L_sm_pw_15$Power)

# * 30 --------------------------------------------------------------------
p <- 30 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

perf_davie_L_sm_pw_30 <- plyr::ldply(graph.list, Perf_Davie_SmallGraph, include.network = "L",
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(perf_davie_L_sm_pw_30$pos_def); unique(perf_davie_L_sm_pw_30$Power)

# * 45 --------------------------------------------------------------------
p <- 45 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

perf_davie_L_sm_pw_45 <- plyr::ldply(graph.list, Perf_Davie_SmallGraph, include.network = "L",
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(perf_davie_L_sm_pw_45$pos_def); unique(perf_davie_L_sm_pw_45$Power)

# Diff Density ------------------------------------------------------------

# * 15 --------------------------------------------------------------------
p <- 15 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

perf_davie_L_dm_pw_15 <- plyr::ldply(graph.list, Perf_Davie_DiffDens, include.network = "L",
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau,
                                     new.edge.prob=0.05) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(perf_davie_L_dm_pw_15$pos_def); unique(perf_davie_L_dm_pw_15$Power)

## ## ## ##
perf_davie_L_dh_pw_15 <- plyr::ldply(graph.list, Perf_Davie_DiffDens, include.network = "L",
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau,
                                     new.edge.prob=0.15) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(perf_davie_L_dh_pw_15$pos_def); unique(perf_davie_L_dh_pw_15$Power)

# * 30 --------------------------------------------------------------------
p <- 30 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

perf_davie_L_dm_pw_30 <- plyr::ldply(graph.list, Perf_Davie_DiffDens, include.network = "L",
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau,
                                     new.edge.prob=0.05) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(perf_davie_L_dm_pw_30$pos_def); unique(perf_davie_L_dm_pw_30$Power)

## ## ## ##
perf_davie_L_dh_pw_30 <- plyr::ldply(graph.list, Perf_Davie_DiffDens, include.network = "L",
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau,
                                     new.edge.prob=0.15) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(perf_davie_L_dh_pw_30$pos_def); unique(perf_davie_L_dh_pw_30$Power)

# * 45 --------------------------------------------------------------------
p <- 45 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

perf_davie_L_dm_pw_45 <- plyr::ldply(graph.list, Perf_Davie_DiffDens, include.network = "L",
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau,
                                     new.edge.prob=0.05) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(perf_davie_L_dm_pw_45$pos_def); unique(perf_davie_L_dm_pw_45$Power)

## ## ## ##
perf_davie_L_dh_pw_45 <- plyr::ldply(graph.list, Perf_Davie_DiffDens, include.network = "L",
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau,
                                     new.edge.prob=0.15) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(perf_davie_L_dh_pw_45$pos_def); unique(perf_davie_L_dh_pw_45$Power)



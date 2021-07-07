###############################
##
## Project: MetaboGuru
##
## Purpose: Partial mismatch for Score test
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2020-11-13
##
## ---------------------------
## Notes:
##   
##
## ---------------------------

# Helpful Functions -------------------------------------------------------

`%nin%` <- Negate(`%in%`)
source('PartialNetworkSimFunctions.R')

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
pp <- 0.40

# Same Size ---------------------------------------------------------------

# * 15 --------------------------------------------------------------------
p <- 15 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

part40_davie_ss_pw_15 <- plyr::ldply(graph.list, Part_Davie_SameSize, perc.perm = pp,
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(part40_davie_ss_pw_15$pos_def); unique(part40_davie_ss_pw_15$Power)

# * 30 --------------------------------------------------------------------
p <- 30 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

part40_davie_ss_pw_30 <- plyr::ldply(graph.list, Part_Davie_SameSize, perc.perm = pp,
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(part40_davie_ss_pw_30$pos_def); unique(part40_davie_ss_pw_30$Power)

# * 45 --------------------------------------------------------------------
p <- 45 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

part40_davie_ss_pw_45 <- plyr::ldply(graph.list, Part_Davie_SameSize, perc.perm = pp,
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(part40_davie_ss_pw_45$pos_def); unique(part40_davie_ss_pw_45$Power)

# Small Graph -------------------------------------------------------------

# * 15 --------------------------------------------------------------------
p <- 15 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

part40_davie_sm_pw_15 <- plyr::ldply(graph.list, Part_Davie_SmallGraph, perc.perm = pp,
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(part40_davie_sm_pw_15$pos_def); unique(part40_davie_sm_pw_15$Power)

# * 30 --------------------------------------------------------------------
p <- 30 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

part40_davie_sm_pw_30 <- plyr::ldply(graph.list, Part_Davie_SmallGraph, perc.perm = pp,
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(part40_davie_sm_pw_30$pos_def); unique(part40_davie_sm_pw_30$Power)

# * 45 --------------------------------------------------------------------
p <- 45 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

part40_davie_sm_pw_45 <- plyr::ldply(graph.list, Part_Davie_SmallGraph, perc.perm = pp,
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(part40_davie_sm_pw_45$pos_def); unique(part40_davie_sm_pw_45$Power)

# Diff Density ------------------------------------------------------------

# * 15 --------------------------------------------------------------------
p <- 15 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

part40_davie_dm_pw_15 <- plyr::ldply(graph.list, Part_Davie_DiffDens, perc.perm = pp,
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau,
                                     new.edge.prob=0.05) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(part40_davie_dm_pw_15$pos_def); unique(part40_davie_dm_pw_15$Power)

## ## ## ##
part40_davie_dh_pw_15 <- plyr::ldply(graph.list, Part_Davie_DiffDens, perc.perm = pp,
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau,
                                     new.edge.prob=0.15) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(part40_davie_dh_pw_15$pos_def); unique(part40_davie_dh_pw_15$Power)

# * 30 --------------------------------------------------------------------
p <- 30 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

part40_davie_dm_pw_30 <- plyr::ldply(graph.list, Part_Davie_DiffDens, perc.perm = pp,
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau,
                                     new.edge.prob=0.05) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(part40_davie_dm_pw_30$pos_def); unique(part40_davie_dm_pw_30$Power)

## ## ## ##
part40_davie_dh_pw_30 <- plyr::ldply(graph.list, Part_Davie_DiffDens, perc.perm = pp,
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau,
                                     new.edge.prob=0.15) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(part40_davie_dh_pw_30$pos_def); unique(part40_davie_dh_pw_30$Power)

# * 45 --------------------------------------------------------------------
p <- 45 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

part40_davie_dm_pw_45 <- plyr::ldply(graph.list, Part_Davie_DiffDens, perc.perm = pp,
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau,
                                     new.edge.prob=0.05) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(part40_davie_dm_pw_45$pos_def); unique(part40_davie_dm_pw_45$Power)

## ## ## ##
part40_davie_dh_pw_45 <- plyr::ldply(graph.list, Part_Davie_DiffDens, perc.perm = pp,
                                     H0.form=H0.form, data = X, b0=b0,
                                     sd.y=sd.y, zz=zz, tau=tau,
                                     new.edge.prob=0.15) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(part40_davie_dh_pw_45$pos_def); unique(part40_davie_dh_pw_45$Power)


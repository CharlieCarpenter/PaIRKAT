###############################
##
## Project: MetaboGuru
##
## Purpose: Power from Simes method
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2021-06-14
##
## ---------------------------
## Notes:
##   
##
## ---------------------------

# Helpful functions
`%nin%` <- Negate(`%in%`)
library(tidyverse); library(magrittr)

source('~/Documents/Research/Current/MetaboGuru/Carpenter/RCode/DaviesSims/CompetitorSims/competitorFuns.R')

## Data Read In ----

nsim <- 10000  ## number of simulations
nperm <- 1000 ## number of permutation for score test
n <- 160 ## sample size
# mX <- matrix(1, n)
# b0 <- 0.2644 ## intercept term
sd.y <- 1.3688 ## standard deviation of Y
tau <- 1 ## Tuning parameter for regularization kernel of normalized laplacian
set.seed(4)
X <- data.frame(X1 = factor(rep(0:1, each = n/2)),
                X2 = runif(n, 0, 5))
b0 <- c(0.2644, 0.5, 0.25)
mX <- model.matrix(~X1+X2, data = X)

# Same Size ---------------------------------------------------------------

# * 15 --------------------------------------------------------------------
p <- 15 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

sms_ss_pw_15 <- plyr::ldply(graph.list, Simes_SameSize, include.graph = F,
                            mX=mX, b0=b0, sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(sms_ss_pw_15$pos_def); unique(sms_ss_pw_15$Power)

# * 30 --------------------------------------------------------------------
p <- 30 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

sms_ss_pw_30 <- plyr::ldply(graph.list, Simes_SameSize, include.graph = F,
                            mX=mX, b0=b0, sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(sms_ss_pw_30$pos_def); unique(sms_ss_pw_30$Power)

# * 45 --------------------------------------------------------------------
p <- 45 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

sms_ss_pw_45 <- plyr::ldply(graph.list, Simes_SameSize, include.graph = F,
                            mX=mX, b0=b0, sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(sms_ss_pw_45$pos_def); unique(sms_ss_pw_45$Power)

# Small Graph -------------------------------------------------------------

# * 15 --------------------------------------------------------------------
p <- 15 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

sms_sm_pw_15 <- plyr::ldply(graph.list, Simes_SmallGraph, include.graph = F,
                            mX=mX, b0=b0, sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(sms_sm_pw_15$pos_def); unique(sms_sm_pw_15$Power)

# * 30 --------------------------------------------------------------------
p <- 30 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

sms_sm_pw_30 <- plyr::ldply(graph.list, Simes_SmallGraph, include.graph = F,
                            mX=mX, b0=b0, sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(sms_sm_pw_30$pos_def); unique(sms_sm_pw_30$Power)

# * 45 --------------------------------------------------------------------
p <- 45 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

sms_sm_pw_45 <- plyr::ldply(graph.list, Simes_SmallGraph, include.graph = F,
                            mX=mX, b0=b0, sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(sms_sm_pw_45$pos_def); unique(sms_sm_pw_45$Power)

# Diff Density ------------------------------------------------------------

# * 15 --------------------------------------------------------------------
p <- 15 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

sms_dm_pw_15 <- plyr::ldply(graph.list, Simes_DiffDens, include.graph = F,
                            mX=mX, b0=b0, sd.y=sd.y, zz=zz, tau=tau,
                            new.edge.prob=0.05) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(sms_dm_pw_15$pos_def); unique(sms_dm_pw_15$Power)

## ## ## ##
sms_dh_pw_15 <- plyr::ldply(graph.list, Simes_DiffDens, include.graph = F,
                            mX=mX, b0=b0, sd.y=sd.y, zz=zz, tau=tau,
                            new.edge.prob=0.15) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(sms_dh_pw_15$pos_def); unique(sms_dh_pw_15$Power)

# * 30 --------------------------------------------------------------------
p <- 30 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

sms_dm_pw_30 <- plyr::ldply(graph.list, Simes_DiffDens, include.graph = F,
                            mX=mX, b0=b0, sd.y=sd.y, zz=zz, tau=tau,
                            new.edge.prob=0.05) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(sms_dm_pw_30$pos_def); unique(sms_dm_pw_30$Power)

## ## ## ##
sms_dh_pw_30 <- plyr::ldply(graph.list, Simes_DiffDens, include.graph = F,
                            mX=mX, b0=b0, sd.y=sd.y, zz=zz, tau=tau,
                            new.edge.prob=0.15) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(sms_dh_pw_30$pos_def); unique(sms_dh_pw_30$Power)

# * 45 --------------------------------------------------------------------
p <- 45 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0.1, p)

sms_dm_pw_45 <- plyr::ldply(graph.list, Simes_DiffDens, include.graph = F,
                            mX=mX, b0=b0, sd.y=sd.y, zz=zz, tau=tau,
                            new.edge.prob=0.05) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(sms_dm_pw_45$pos_def); unique(sms_dm_pw_45$Power)

## ## ## ##
sms_dh_pw_45 <- plyr::ldply(graph.list, Simes_DiffDens, include.graph = F,
                            mX=mX, b0=b0, sd.y=sd.y, zz=zz, tau=tau,
                            new.edge.prob=0.15) %>% 
  mutate(Power = sum(pVal < 0.05, na.rm = T)/n())
sum(sms_dh_pw_45$pos_def); unique(sms_dh_pw_45$Power)




###############################
##
## Project: MetaboGuru
##
## Purpose: TypeI Error rates of Cumu Tests
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
source('CompleteMismatchSimFunctions.R')

nsim <- 1000  ## number of simulations
nperm <- 1000 ## number of permutation for score test
n <- 160 ## sample size
mX <- matrix(1, n)
b0 <- 0.2644 ## intercept term
sd.y <- 1.3688 ## standard deviation of Y
delta <- 1 ## Tuning parameter for regularization kernel of normalized laplacian
set.seed(4)
# X <- data.frame(X1 = factor(rep(0:1, each = n/2)),
#                 X2 = runif(n, 0, 5))
# b0 <- c(0.2644, 0.5, 0.25)
# mX <- model.matrix(~X1+X2, data = X)

# Same Size ---------------------------------------------------------------

# * 15 --------------------------------------------------------------------
p <- 15 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_scor_ss_t1_15 <- plyr::ldply(graph.list, Comp_Score_SameSize,
                                  mX=mX, b0=b0, sd.y=sd.y, zz=zz, delta=delta) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_scor_ss_t1_15$pos_def); unique(comp_scor_ss_t1_15$TypeI)

# * 30 --------------------------------------------------------------------
p <- 30 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_scor_ss_t1_30 <- plyr::ldply(graph.list, Comp_Score_SameSize,
                                  mX=mX, b0=b0, sd.y=sd.y, zz=zz, delta=delta) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_scor_ss_t1_30$pos_def); unique(comp_scor_ss_t1_30$TypeI)

# * 45 --------------------------------------------------------------------
p <- 45 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_scor_ss_t1_45 <- plyr::ldply(graph.list, Comp_Score_SameSize,
                                  mX=mX, b0=b0, sd.y=sd.y, zz=zz, delta=delta) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_scor_ss_t1_45$pos_def); unique(comp_scor_ss_t1_45$TypeI)

# Small Graph -------------------------------------------------------------

# * 15 --------------------------------------------------------------------
p <- 15 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_scor_sm_t1_15 <- plyr::ldply(graph.list, Comp_Score_SmallGraph,
                                  mX=mX, b0=b0, sd.y=sd.y, zz=zz, delta=delta) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_scor_sm_t1_15$pos_def); unique(comp_scor_sm_t1_15$TypeI)

# * 30 --------------------------------------------------------------------
p <- 30 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_scor_sm_t1_30 <- plyr::ldply(graph.list, Comp_Score_SmallGraph,
                                  mX=mX, b0=b0, sd.y=sd.y, zz=zz, delta=delta) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_scor_sm_t1_30$pos_def); unique(comp_scor_sm_t1_30$TypeI)

# * 45 --------------------------------------------------------------------
p <- 45 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_scor_sm_t1_45 <- plyr::ldply(graph.list, Comp_Score_SmallGraph,
                                  mX=mX, b0=b0, sd.y=sd.y, zz=zz, delta=delta) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_scor_sm_t1_45$pos_def); unique(comp_scor_sm_t1_45$TypeI)

# Diff Density ------------------------------------------------------------

# * 15 --------------------------------------------------------------------
p <- 15 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_scor_dm_t1_15 <- plyr::ldply(graph.list, Comp_Score_DiffDens,
                                  mX=mX, b0=b0, sd.y=sd.y, zz=zz, delta=delta,
                                  new.edge.prob=0.05) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_scor_dm_t1_15$pos_def); unique(comp_scor_dm_t1_15$TypeI)

## ## ## ##
comp_scor_dh_t1_15 <- plyr::ldply(graph.list, Comp_Score_DiffDens,
                                  mX=mX, b0=b0, sd.y=sd.y, zz=zz, delta=delta,
                                  new.edge.prob=0.15) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_scor_dh_t1_15$pos_def); unique(comp_scor_dh_t1_15$TypeI)

# * 30 --------------------------------------------------------------------
p <- 30 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_scor_dm_t1_30 <- plyr::ldply(graph.list, Comp_Score_DiffDens,
                                  mX=mX, b0=b0, sd.y=sd.y, zz=zz, delta=delta,
                                  new.edge.prob=0.05) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_scor_dm_t1_30$pos_def); unique(comp_scor_dm_t1_30$TypeI)

## ## ## ##
comp_scor_dh_t1_30 <- plyr::ldply(graph.list, Comp_Score_DiffDens,
                                  mX=mX, b0=b0, sd.y=sd.y, zz=zz, delta=delta,
                                  new.edge.prob=0.15) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_scor_dh_t1_30$pos_def); unique(comp_scor_dh_t1_30$TypeI)

# * 45 --------------------------------------------------------------------
p <- 45 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_scor_dm_t1_45 <- plyr::ldply(graph.list, Comp_Score_DiffDens,
                                  mX=mX, b0=b0, sd.y=sd.y, zz=zz, delta=delta,
                                  new.edge.prob=0.05) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_scor_dm_t1_45$pos_def); unique(comp_scor_dm_t1_45$TypeI)

## ## ## ##
comp_scor_dh_t1_45 <- plyr::ldply(graph.list, Comp_Score_DiffDens,
                                  mX=mX, b0=b0, sd.y=sd.y, zz=zz, delta=delta,
                                  new.edge.prob=0.15) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_scor_dh_t1_45$pos_def); unique(comp_scor_dh_t1_45$TypeI)



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
zz <- rep(0, p)

comp_davie_ss_t1_15 <- plyr::ldply(graph.list, Comp_Davie_SameSize,
                                   H0.form=H0.form, data = X, b0=b0,
                                   sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_davie_ss_t1_15$pos_def); unique(comp_davie_ss_t1_15$TypeI)

# * 30 --------------------------------------------------------------------
p <- 30 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_davie_ss_t1_30 <- plyr::ldply(graph.list, Comp_Davie_SameSize,
                                   H0.form=H0.form, data = X, b0=b0,
                                   sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_davie_ss_t1_30$pos_def); unique(comp_davie_ss_t1_30$TypeI)

# * 45 --------------------------------------------------------------------
p <- 45 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_davie_ss_t1_45 <- plyr::ldply(graph.list, Comp_Davie_SameSize,
                                   H0.form=H0.form, data = X, b0=b0,
                                   sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_davie_ss_t1_45$pos_def); unique(comp_davie_ss_t1_45$TypeI)

# Small Graph -------------------------------------------------------------

# * 15 --------------------------------------------------------------------
p <- 15 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_davie_sm_t1_15 <- plyr::ldply(graph.list, Comp_Davie_SmallGraph,
                                   H0.form=H0.form, data = X, b0=b0,
                                   sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_davie_sm_t1_15$pos_def); unique(comp_davie_sm_t1_15$TypeI)

# * 30 --------------------------------------------------------------------
p <- 30 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_davie_sm_t1_30 <- plyr::ldply(graph.list, Comp_Davie_SmallGraph,
                                   H0.form=H0.form, data = X, b0=b0,
                                   sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_davie_sm_t1_30$pos_def); unique(comp_davie_sm_t1_30$TypeI)

# * 45 --------------------------------------------------------------------
p <- 45 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_davie_sm_t1_45 <- plyr::ldply(graph.list, Comp_Davie_SmallGraph,
                                   H0.form=H0.form, data = X, b0=b0,
                                   sd.y=sd.y, zz=zz, tau=tau) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_davie_sm_t1_45$pos_def); unique(comp_davie_sm_t1_45$TypeI)

# Diff Density ------------------------------------------------------------

# * 15 --------------------------------------------------------------------
p <- 15 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_davie_dm_t1_15 <- plyr::ldply(graph.list, Comp_Davie_DiffDens,
                                   H0.form=H0.form, data = X, b0=b0,
                                   sd.y=sd.y, zz=zz, tau=tau,
                                   new.edge.prob=0.05) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_davie_dm_t1_15$pos_def); unique(comp_davie_dm_t1_15$TypeI)

## ## ## ##
comp_davie_dh_t1_15 <- plyr::ldply(graph.list, Comp_Davie_DiffDens,
                                   H0.form=H0.form, data = X, b0=b0,
                                   sd.y=sd.y, zz=zz, tau=tau,
                                   new.edge.prob=0.15) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_davie_dh_t1_15$pos_def); unique(comp_davie_dh_t1_15$TypeI)

# * 30 --------------------------------------------------------------------
p <- 30 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_davie_dm_t1_30 <- plyr::ldply(graph.list, Comp_Davie_DiffDens,
                                   H0.form=H0.form, data = X, b0=b0,
                                   sd.y=sd.y, zz=zz, tau=tau,
                                   new.edge.prob=0.05) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_davie_dm_t1_30$pos_def); unique(comp_davie_dm_t1_30$TypeI)

## ## ## ##
comp_davie_dh_t1_30 <- plyr::ldply(graph.list, Comp_Davie_DiffDens,
                                   H0.form=H0.form, data = X, b0=b0,
                                   sd.y=sd.y, zz=zz, tau=tau,
                                   new.edge.prob=0.15) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_davie_dh_t1_30$pos_def); unique(comp_davie_dh_t1_30$TypeI)

# * 45 --------------------------------------------------------------------
p <- 45 ## size of network
set.seed(2)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
zz <- rep(0, p)

comp_davie_dm_t1_45 <- plyr::ldply(graph.list, Comp_Davie_DiffDens,
                                   H0.form=H0.form, data = X, b0=b0,
                                   sd.y=sd.y, zz=zz, tau=tau,
                                   new.edge.prob=0.05) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_davie_dm_t1_45$pos_def); unique(comp_davie_dm_t1_45$TypeI)

## ## ## ##
comp_davie_dh_t1_45 <- plyr::ldply(graph.list, Comp_Davie_DiffDens,
                                   H0.form=H0.form, data = X, b0=b0,
                                   sd.y=sd.y, zz=zz, tau=tau,
                                   new.edge.prob=0.15) %>% 
  mutate(TypeI = sum(pVal < 0.05, na.rm = T)/n())
sum(comp_davie_dh_t1_45$pos_def); unique(comp_davie_dh_t1_45$TypeI)

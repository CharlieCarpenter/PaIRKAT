###############################
##
## Project: MetaoGuru Aim 3
##
## Purpose: Changing Signal to Noise ratio ( Var(Beta*Z)/Var(Y-Beta*Z) )
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2020-04-16
##
## ---------------------------
## Notes:
##   
##
## ---------------------------

## Data Read In ----

# Helpful functions
`%nin%` <- Negate(`%in%`)
source("PaIRKAT/helpers.R")

## network simulation functions
source('CompleteMismatchSimFunctions.R')
source('PartialNetworkSimFunctions.R')
source('PerfectNetworkSimFunctions.R')
source('NoNetworkSimFunctions.R')

## Parameters ----
nsim <- 1000 ## number of simulations
nperm <- 1000 ## number of permutation for score test
n <- 160 ## sample size
X <- matrix(1, n) ## Null design matrix for permutation test
b0 <- 0.2644 ## intercept term
sd.y <- 1.3688 ## standard deviation of Y
delta <- 1 ## Tuning parameter for regularization kernel of normalized laplacian
p <- 30 ## 30 metabolites for everything

## Generating Data ----

set.seed(2345)
graph.list <- lapply(1:nsim, function(x) sample_pa(n = p, directed = F))
betas <- lapply(seq(0.1,2, length.out = 10), ## Tuning parameter to multiply betas
                function(a) a*rnorm(p, mean = 0.1, sd = 0.1))

sig.noise.data <- plyr::llply(betas, function(b){
  
  ## Generating data set from every graph
  Y.Z <- plyr::llply(graph.list, function(graph){
    ## Converts a graph into an adjacency matrix
    A <- as.matrix( get.adjacency(graph) )
    total_edge_n <- sum(A==1)
    
    ## Convert adjacency matrix into a "starter" precision matrix
    test <- A + diag(p)
    
    ## Make test positive using approach of Danaher et al (2014)
    ## Function in SimulationFuncitons.R
    Omega1 <- Danaher_pos_def(test)
    Z <- mvrnorm(n=n, rep(0,p), solve(Omega1))
    model <- model.matrix(~ 1 , data.frame(Z))
    Y <- rnorm(n, b0+Z%*%b, sd.y)
    V.Y <- var(b0+Z%*%b)
    V.e <- var(Y - (b0+Z%*%b))
    
    list(Y = Y, Z = Z, V.Y = V.Y, V.e = V.e, Ratio = V.Y/V.e)
  })
  ## Pulling Ratio from every data set 
  sig.noise <- Y.Z %>% map_dbl("Ratio") %>% mean
  
  list(Y.Z = Y.Z, Sig.Noise = sig.noise)
})
(sn <- unlist(map(sig.noise.data, "Sig.Noise")))

## No Network ----

set.seed(2)
no.net.sig.noise <- plyr::ldply(sig.noise.data, function(sig.dat){
  
  Y.Z <- sig.dat[["Y.Z"]]
  mapply(NoNet_Score_SigNoise, Y.Z, graph.list,
         MoreArgs = list(.X = X))
  }, .progress = "time")

## Complete Mismatch ----

set.seed(2)
comp.mis.sig.noise <- plyr::ldply(sig.noise.data, function(sig.dat){
  
  Y.Z <- sig.dat[["Y.Z"]]
  mapply(Comp_Score_SigNoise, Y.Z, graph.list, 
         MoreArgs = list(.delta = delta, .X = X))
  }, .progress = "time")

comp.mis.sig.noise.L <- plyr::ldply(sig.noise.data, function(sig.dat){
  
  Y.Z <- sig.dat[["Y.Z"]]
  mapply(Comp_Score_SigNoise, Y.Z, graph.list, 
         MoreArgs = list(include.network = "L",
                         .delta = delta, .X = X))
}, .progress = "time")

## Partial Mismatch 10----
perc.perm <- 0.1

set.seed(2)
part10.mis.sig.noise <- plyr::ldply(sig.noise.data, function(sig.dat){
  
  Y.Z <- sig.dat[["Y.Z"]]
  mapply(Part_Score_SigNoise, Y.Z, graph.list, 
         MoreArgs = list(.delta = delta, .X = X, 
                         .perc.perm = perc.perm))
  }, .progress = "time")

part10.mis.sig.noise.L <- plyr::ldply(sig.noise.data, function(sig.dat){
  
  Y.Z <- sig.dat[["Y.Z"]]
  mapply(Part_Score_SigNoise, Y.Z, graph.list, 
         MoreArgs = list(.delta = delta, .X = X, 
                         include.network = "L",
                         .perc.perm = perc.perm))
}, .progress = "time")

## Partial Mismatch 40----
perc.perm <- 0.4

set.seed(2)
part40.mis.sig.noise <- plyr::ldply(sig.noise.data, function(sig.dat){
  
  Y.Z <- sig.dat[["Y.Z"]]
  mapply(Part_Score_SigNoise, Y.Z, graph.list, 
         MoreArgs = list(.delta = delta, .X = X,
                         .perc.perm = perc.perm))
}, .progress = "time")

part40.mis.sig.noise.L <- plyr::ldply(sig.noise.data, function(sig.dat){
  
  Y.Z <- sig.dat[["Y.Z"]]
  mapply(Part_Score_SigNoise, Y.Z, graph.list, 
         MoreArgs = list(.delta = delta, .X = X,
                         include.network = "L",
                         .perc.perm = perc.perm))
}, .progress = "time")

## Partial Mismatch 70----
perc.perm <- 0.7

set.seed(2)
part70.mis.sig.noise <- plyr::ldply(sig.noise.data, function(sig.dat){
  
  Y.Z <- sig.dat[["Y.Z"]]
  mapply(Part_Score_SigNoise, Y.Z, graph.list, 
         MoreArgs = list(.delta = delta, .X = X,
                         .perc.perm = perc.perm))
}, .progress = "time")

part70.mis.sig.noise.L <- plyr::ldply(sig.noise.data, function(sig.dat){
  
  Y.Z <- sig.dat[["Y.Z"]]
  mapply(Part_Score_SigNoise, Y.Z, graph.list, 
         MoreArgs = list(.delta = delta, .X = X,
                         include.network = "L",
                         .perc.perm = perc.perm))
}, .progress = "time")

## Perfect Network ----

set.seed(2)
perf.net.sig.noise <- plyr::ldply(sig.noise.data, function(sig.dat){
  
  Y.Z <- sig.dat[["Y.Z"]]
  mapply(Perf_Score_SigNoise, Y.Z, graph.list, 
         MoreArgs = list(.delta = delta, .X = X))
  }, .progress = "time")

perf.net.sig.noise.L <- plyr::ldply(sig.noise.data, function(sig.dat){
  
  Y.Z <- sig.dat[["Y.Z"]]
  mapply(Perf_Score_SigNoise, Y.Z, graph.list, 
         MoreArgs = list(include.network = "L",
                         .delta = delta, .X = X))
}, .progress = "time")

## Result ----

st <- c("Perfect", "Partial Mismatch (10)", "Partial Mismatch (40)",
        "Partial Mismatch (70)", "Complete Mismatch", "No Network")

pow <- function(p) sum(p < 0.05)/length(p)

pp <- c(apply(perf.net.sig.noise, 1,  pow),
        apply(perf.net.sig.noise.L, 1,  pow),
        apply(part10.mis.sig.noise, 1, pow ),
        apply(part10.mis.sig.noise.L, 1, pow ),
        apply(part40.mis.sig.noise, 1, pow ),
        apply(part40.mis.sig.noise.L, 1, pow ),
        apply(part70.mis.sig.noise, 1, pow ),
        apply(part70.mis.sig.noise.L, 1, pow ),
        apply(comp.mis.sig.noise, 1, pow ),
        apply(comp.mis.sig.noise.L, 1, pow ),
        apply(no.net.sig.noise, 1, pow ))

sig.noise.result <- data.frame(Sig.Noise = rep(sn,11), 
                               Structure = factor(c(rep(st[-length(st)],
                                                        each = 2*length(sn)),
                                                    rep(st[length(st)], 
                                                        length(sn))),
                                                  levels = st),
                               Test = c(rep(c("RL", "L"), each = 10, 
                                            times = (length(st)-1)),
                                        rep("NN", 10)),
                               Power = pp)

## Go to "Generating Data" above for "sn"
bold <- element_text(face = "bold")

lbs <- lapply(c("\\textbf{No Network}", 
                "$\\widetilde{\\textbf{L}}$",
                "$\\widetilde{\\textbf{L}}_R$"), TeX)

ggplot(sig.noise.result, aes(x = Sig.Noise, y = Power,
                             linetype = Test, col = Structure)) +
  geom_line() + geom_point() + theme_bw() +
  scale_x_continuous(breaks = round(sn,3)[-c(2,4)]) +
  scale_linetype_manual(values = c('solid', 'dotted', 'twodash'),
                        labels = lbs) +
  theme(axis.text.x = element_text(face = "bold", angle = 65, hjust = 1),
        legend.title = bold, legend.text = bold,
        axis.text.y = bold, axis.title.x = bold, axis.title.y = bold) +
  labs(x = TeX("$\\frac{Var(\\beta Z)}{Var(\\epsilon)}$"))


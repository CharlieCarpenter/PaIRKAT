###############################
##
## Project: MetaboGuru: PaIRKAT
##
## Purpose: PCA testing for pathways
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2021-06-11
##
## ---------------------------
## Notes:
##   
##
## ---------------------------

# Helpful functions
`%nin%` <- Negate(`%in%`)
library(tidyverse); library(magrittr)

# PCA Functions -----------------------------------------------------------

## Same Size
PCA_SameSize <- function(graph, mX, b0, sd.y, zz, tau, m,
                         include.graph = T){
  
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  p <- nrow(A)
  
  ## Convert adjacency matrix into a "starter" precision matrix
  test <- A + diag(p)
  total_edge_n <- sum(A==1)
  
  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)
  if(pd){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=nrow(mX), rep(0,p), solve(Omega1))
    Y <- rnorm(n=nrow(mX), mX%*%b0+Z%*%zz, sd.y)
    V.Y <- var(mX%*%b0+Z%*%zz); V.e <- var(Y - mX%*%b0+Z%*%zz)
    
    L <- graph.laplacian(graph, normalized = TRUE)
    RL <- solve(diag(nrow(L)) +tau*L)
    
    Z <- scale(Z); if(include.graph) Z <- Z%*%RL
    pca_pval <- CCpcaTest(.Y=Y, .mX=mX, .Z=Z, .m=m)
    
  } else{
    pca_pval <- V.Y <- V.e <- NA
  }
  c(pVal = pca_pval, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

## Small 
PCA_SmallGraph <- function(graph, mX, b0, sd.y, zz, tau, m,
                           include.graph = T){
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  p <- nrow(A); total_edge_n <- sum(A==1)
  
  ## Convert adjacency matrix into a "starter" precision matrix
  test <- A + diag(p)
  
  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)
  if(pd){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=nrow(mX), rep(0,p), solve(Omega1))
    Y <- rnorm(n=nrow(mX), mX%*%b0+Z%*%zz, sd.y)
    V.Y <- var(mX%*%b0+Z%*%zz); V.e <- var(Y - mX%*%b0+Z%*%zz)
    
    ## Shrinking graphs
    nd <- new.degs(graph)
    As <- A[nd > 0, nd > 0]; Zs <- Z[ , nd > 0]
    gs <- graph_from_adjacency_matrix(As, mode = "undirected")
    
    L <- graph.laplacian(gs, normalized = TRUE)
    RL <- solve(diag(nrow(L)) +tau*L)
    
    Zs <- scale(Zs); if(include.graph) Zs <- Zs%*%RL
    pca_pval <- CCpcaTest(.Y=Y, .mX=mX, .Z=Zs, .m=ncol(Zs))
    
  } else{
    pca_pval <- V.Y <- V.e <- NA
  }
  c(pVal = pca_pval, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

## Diff Dens
PCA_DiffDens <- function(graph, mX, b0, sd.y, zz, tau, 
                         new.edge.prob, m, include.graph = T){
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  p <- nrow(A); total_edge_n <- sum(A==1)
  A_big <- increase.edge.density(A, edge.prob = new.edge.prob)
  
  ## Convert adjacency matrix into a "starter" precision matrix
  test <- A_big + diag(p)
  
  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)
  if(pd){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=nrow(mX), rep(0,p), solve(Omega1))
    Y <- rnorm(n=nrow(mX), mX%*%b0+Z%*%zz, sd.y)
    V.Y <- var(mX%*%b0+Z%*%zz); V.e <- var(Y - mX%*%b0+Z%*%zz)
    
    L <- graph.laplacian(graph, normalized = TRUE)
    RL <- solve(diag(nrow(L)) +tau*L)

    Z <- scale(Z); if(include.graph) Z <- Z%*%RL
    pca_pval <- CCpcaTest(.Y=Y, .mX=mX, .Z=Z, .m=m)
    
  } else{
    pca_pval <- V.Y <- V.e <- NA
  }
  c(pVal = pca_pval, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

PCA_SigNoise <- function(YZ, graph, H0.form, data, m,
                         include.graph = T, .tau){
  
  L <- graph.laplacian(graph, normalized = TRUE)
  RL <- solve(diag(nrow(L)) +tau*L)
  
  Z <- scale(YZ$Z)
  if(include.graph) Z <- Z%*%RL
  
  dd <- cbind(Y=YZ$Y, data)
  mX <- model.matrix(H0.form, data = dd)
  
  pca_pval <- CCpcaTest(.Y=dd$Y, .mX=mX, .Z=Z, .m=m)
}

# Simes Functions ---------------------------------------------------------

## Same Size
Simes_SameSize <- function(graph, mX, b0, sd.y, zz, tau, include.graph = T){
  
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  p <- nrow(A)
  
  ## Convert adjacency matrix into a "starter" precision matrix
  test <- A + diag(p)
  total_edge_n <- sum(A==1)
  
  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)
  if(pd){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=nrow(mX), rep(0,p), solve(Omega1))
    Y <- rnorm(n=nrow(mX), mX%*%b0+Z%*%zz, sd.y)
    V.Y <- var(mX%*%b0+Z%*%zz); V.e <- var(Y - mX%*%b0+Z%*%zz)
    
    L <- graph.laplacian(graph, normalized = TRUE)
    RL <- solve(diag(nrow(L)) +tau*L)
    
    Z <- scale(Z); if(include.graph) Z <- Z%*%RL
    sime_pval <- Simes(.Y=Y, .mX=mX, .Z=Z)
    
  } else{
    sime_pval <- V.Y <- V.e <- NA
  }
  c(pVal = sime_pval, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

## Small 
Simes_SmallGraph <- function(graph, mX, b0, sd.y, zz, tau, include.graph = T){
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  p <- nrow(A); total_edge_n <- sum(A==1)
  
  ## Convert adjacency matrix into a "starter" precision matrix
  test <- A + diag(p)
  
  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)
  if(pd){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=nrow(mX), rep(0,p), solve(Omega1))
    Y <- rnorm(n=nrow(mX), mX%*%b0+Z%*%zz, sd.y)
    V.Y <- var(mX%*%b0+Z%*%zz); V.e <- var(Y - mX%*%b0+Z%*%zz)
    
    ## Shrinking graphs
    nd <- new.degs(graph)
    As <- A[nd > 0, nd > 0]; Zs <- Z[ , nd > 0]
    
    gs <- graph_from_adjacency_matrix(As, mode = 'undirected')
    L <- graph.laplacian(gs, normalized = TRUE)
    RL <- solve(diag(nrow(L)) +tau*L)
    
    Zs <- scale(Zs); if(include.graph) Zs <- Zs%*%RL
    sime_pval <- Simes(.Y=Y, .mX=mX, .Z=Zs)
    
  } else{
    sime_pval <- V.Y <- V.e <- NA
  }
  c(pVal = sime_pval, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

## Diff Dens
Simes_DiffDens <- function(graph, mX, b0, sd.y, zz, tau, 
                           new.edge.prob, include.graph = T){
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  p <- nrow(A); total_edge_n <- sum(A==1)
  A_big <- increase.edge.density(A, edge.prob = new.edge.prob)
  
  ## Changing graph to new density from adj matrix
  # graph <- A_big %>% 
  #   graph_from_adjacency_matrix(mode = "undirected")
  # new.edge.density <- edge_density(graph)
  
  ## Convert adjacency matrix into a "starter" precision matrix
  test <- A_big + diag(p)
  
  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)
  if(pd){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=nrow(mX), rep(0,p), solve(Omega1))
    Y <- rnorm(n=nrow(mX), mX%*%b0+Z%*%zz, sd.y)
    V.Y <- var(mX%*%b0+Z%*%zz); V.e <- var(Y - mX%*%b0+Z%*%zz)
    
    L <- graph.laplacian(graph, normalized = TRUE)
    RL <- solve(diag(nrow(L)) +tau*L)
    
    Z <- scale(Z); if(include.graph) Z <- Z%*%RL
    sime_pval <- Simes(.Y=Y, .mX=mX, .Z=Z)
    
  } else{
    sime_pval <- V.Y <- V.e <- NA
  }
  c(pVal = sime_pval, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

Simes_SigNoise <- function(YZ, graph, H0.form, data,
                           include.graph = T, .tau){
  
  L <- graph.laplacian(graph, normalized = TRUE)
  RL <- solve(diag(nrow(L)) +tau*L)
  
  Z <- scale(YZ$Z)
  if(include.graph) Z <- Z%*%RL
  
  dd <- cbind(Y=YZ$Y, data)
  mX <- model.matrix(H0.form, data = dd)
  
  sime_pval <- Simes(.Y=dd$Y, .mX=mX, .Z=Z)
}

# Core Functions ----------------------------------------------------------

CCpcaTest <- function(.Y, .mX, .Z, .m){
  SVD <- svd(.Z)
  scores <- SVD$u %*% diag(SVD$d)
  scrs <- scores[,1:.m]
  
  ll1 <- lm(.Y ~ .mX)
  mXc <- cbind(.mX, scrs)
  ll2 <- lm(.Y ~ mXc)
  
  anova(ll1, ll2)$`Pr(>F)`[2]
}

Simes <- function(.Y, .mX, .Z){
  m <- ncol(.Z)
  ps <- numeric(m)
  
  ## Univariate tests
  for(i in 1:m){
    dd <- data.frame(Y = .Y, .mX, Z = .Z[,i])
    ll <- lm(Y ~ mX + Z, data = dd)
    
    ps[i] <- coef(summary(ll))[ncol(dd)-1, 4]
  }
  
  min(sort(ps)*(m)/(1:m)) # simes p-value adjustment
}

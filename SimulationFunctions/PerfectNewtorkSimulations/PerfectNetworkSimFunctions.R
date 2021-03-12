###############################
##
## Project: MetaoGuru Aim 3
##
## Purpose: Perfect Network Simulation Functions
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2020-02-11
##
## ---------------------------
## Notes:
## Doesn't use A_p
## Uses L_tilda <- normalize.laplacian(A) for K_regL
## i.e. perfect match
## ---------------------------

# source("PaIRKAT/helpers.R")

Perf_Scor_SameSize <- function(graph, mX, b0, sd.y, zz, delta,
                               rho = "median.pairwise", include.network = "R"){
  
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
    
    L_tilda <- graph.laplacian(graph, normalized = T) %>% as.matrix
    
    if(include.network == "R") K_regL <- solve(diag(p) + delta*L_tilda)
    if(include.network == "L") K_regL <- L_tilda
    
    Z <- scale(Z); ZL <- Z %*% K_regL
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    K <- Gaussian_kernel(rho, ZL)
    sc <- CC_Chisq_Score(K, Y, mX)
    } else{
      sc <- V.Y <- V.e <- NA
      }
  c(sc, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

Perf_Scor_SmallGraph <- function(graph, mX, b0, sd.y, zz, delta,
                                 rho = "median.pairwise", include.network = "R"){
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
    
    L_tilda <- As %>% 
      graph_from_adjacency_matrix(mode = "undirected") %>% 
      graph.laplacian(normalized = T) %>% as.matrix
    
    if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda)) + delta*L_tilda)
    if(include.network == "L") K_regL <- L_tilda
    
    Zs <- scale(Zs); ZL <- Zs %*% K_regL
    if(rho == "median.pairwise") rho <- median(dist(Zs))
    
    K <- Gaussian_kernel(rho, ZL)
    sc <- CC_Chisq_Score(K, Y, mX)
  } else{
    sc <- V.Y <- V.e <- NA
  }
  c(sc, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

Perf_Scor_DiffDens <- function(graph, mX, b0, sd.y, zz, delta, new.edge.prob,
                               rho = "median.pairwise", include.network = "R"){
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  p <- nrow(A); total_edge_n <- sum(A==1)
  A_big <- increase.edge.density(A, edge.prob = new.edge.prob)
  
  ## Changing graph to new density from adj matrix
  graph <- A_big %>% 
    graph_from_adjacency_matrix(mode = "undirected")
  new.edge.density <- edge_density(graph)
  
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
    
    L_tilda <- graph.laplacian(graph, normalized = T) %>% as.matrix
    
    if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda)) + delta*L_tilda)
    if(include.network == "L") K_regL <- L_tilda
    
    Z <- scale(Z); ZL <- Z %*% K_regL
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    K <- Gaussian_kernel(rho, ZL)
    sc <- CC_Chisq_Score(K, Y, mX)
  } else{
    sc <- V.Y <- V.e <- NA
  }
  c(sc, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

## A function built for power analysis on prebuilt data sets
Perf_Score_SigNoise <- function(YZ, graph, include.network = "R",
                                .delta, .X){
  
  L_tilda <- graph.laplacian(graph, normalized = T) %>% as.matrix
  
  if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda)) + .delta*L_tilda)
  if(include.network == "L") K_regL <- L_tilda
  
  Z <- scale(YZ$Z); rho <- median(dist(Z))
  ZL <- Z %*% K_regL
  
  K <- Gaussian_kernel(rho, ZL)
  sc <- CC_Chisq_Score(K, YZ$Y, .X)
  
  sc["pVal"]
}


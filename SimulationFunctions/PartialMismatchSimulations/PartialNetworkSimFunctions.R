###############################
##
## Project: MetaoGuru Aim 3
##
## Purpose: Functions for Partial mismatch network match
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2020-03-04
##
## ---------------------------
## Notes:
##
##
## ---------------------------

source("PaIRKAT/helpers.R")

## forces a change in specified number of edges (n_diff_edge)
Part_Score_SameSize <- function(graph, mX, b0, sd.y, zz, delta, perc.perm, # Decimals
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
    
    ## Graph with permuted edges
    A_p <- edge.perm(A, perc.perm)
    nse <- count_shared_edges(A, A_p)
    
    L_tilda_p <- A_p %>% 
      graph_from_adjacency_matrix(mode = "undirected") %>% 
      graph.laplacian(normalized = T) %>% as.matrix
    
    if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda_p)) + delta*L_tilda_p)
    if(include.network == "L") K_regL <- L_tilda_p
    
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


## performs same tests, but changes density of testing graph
Part_Score_SmallGraph <- function(graph, mX, b0, sd.y, zz, delta, perc.perm, # Decimals
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
    
    ## Shrinking graphs and data
    nd <- new.degs(graph)
    As <- A[nd > 0, nd > 0]; Zs <- Z[ , nd > 0]
    
    ## Graph with permuted edges
    A_p <- edge.perm(As, perc.perm)
    
    L_tilda_p <- A_p %>% 
      graph_from_adjacency_matrix(mode = "undirected") %>% 
      graph.laplacian(normalized = T) %>% as.matrix
    
    if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda_p)) + delta*L_tilda_p)
    if(include.network == "L") K_regL <- L_tilda_p
    
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

## Changes network density
Part_Score_DiffDens <- function(graph, mX, b0, sd.y, zz, delta, perc.perm, # Decimals
                                new.edge.prob,
                                rho = "median.pairwise", include.network = "R"){
  ## Converts a graph into an adjacency matrix
  A <- get.adjacency(graph) %>% as.matrix %>% 
    increase.edge.density(edge.prob = new.edge.prob)
  total_edge_n <- sum(A==1)
  graph <- graph_from_adjacency_matrix(A, mode = "undirected")
  new.edge.dens <- edge_density(graph)
  
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
    
    ## Graph with permuted edges
    A_p <- edge.perm(A, perc.perm)
    nse <- count_shared_edges(A, A_p)
    
    L_tilda_p <- A_p %>% 
      graph_from_adjacency_matrix(mode = "undirected") %>% 
      graph.laplacian(normalized = T) %>% as.matrix
    
    if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda_p)) + delta*L_tilda_p)
    if(include.network == "L") K_regL <- L_tilda_p
    
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
Part_Score_SigNoise <- function(YZ, graph, include.network = "R",
                                .perc.perm, .X, .delta){
  
  A <- as.matrix(get.adjacency(graph))
  A_p <- edge.perm(A, .perc.perm)
  
  L_tilda <- A_p %>% 
    graph_from_adjacency_matrix(mode = "undirected") %>% 
    graph.laplacian(normalized = T) %>% as.matrix
  
  if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda)) + .delta*L_tilda)
  if(include.network == "L") K_regL <- L_tilda
  
  Z <- scale(YZ$Z); rho <- median(dist(Z))
  Z <- Z %*% K_regL
  
  K <- Gaussian_kernel(rho, Z)
  sc <- CC_Chisq_Score(K, YZ$Y, .X)
  
  sc["pVal"]
}


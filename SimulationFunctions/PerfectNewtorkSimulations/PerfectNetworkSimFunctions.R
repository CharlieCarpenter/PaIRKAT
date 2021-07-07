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

Perf_Davie_SameSize <- function(graph, H0.form, data,
                                b0, sd.y, zz, tau, kernel = "G",
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
    mX <- as.matrix(cbind(1, data))
    Z <- mvrnorm(n=nrow(mX), rep(0,p), solve(Omega1))
    Y <- rnorm(n=nrow(mX), mX%*%b0+Z%*%zz, sd.y)
    V.Y <- var(mX%*%b0+Z%*%zz); V.e <- var(Y - mX%*%b0+Z%*%zz)
    
    dd <- cbind(Y, data)
    L_tilda <- graph.laplacian(graph, normalized = T) %>% as.matrix
    
    if(include.network == "R") K_regL <- solve(diag(p) + tau*L_tilda)
    if(include.network == "L") K_regL <- L_tilda
    
    Z <- scale(Z); ZL <- Z %*% K_regL
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    if(kernel == "G") K <- Gaussian_kernel(rho, ZL)
    if(kernel == "L") K <- plyKern(ZL, pow = 1, rho = 0)
    
    dav_pval <- SKAT.c(H0.form, data = dd, K=K)$Qq
    
  } else{
    dav_pval <- V.Y <- V.e <- NA
  }
  c(pVal = dav_pval, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

Perf_Davie_SmallGraph <- function(graph, H0.form, data,
                                  b0, sd.y, zz, tau, kernel = "G",
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
    mX <- as.matrix(cbind(1, data))
    Z <- mvrnorm(n=nrow(mX), rep(0,p), solve(Omega1))
    Y <- rnorm(n=nrow(mX), mX%*%b0+Z%*%zz, sd.y)
    V.Y <- var(mX%*%b0+Z%*%zz); V.e <- var(Y - mX%*%b0+Z%*%zz)
    
    dd <- cbind(Y, data)
    ## Shrinking graphs
    nd <- new.degs(graph)
    As <- A[nd > 0, nd > 0]; Zs <- Z[ , nd > 0]
    
    L_tilda <- As %>% 
      graph_from_adjacency_matrix(mode = "undirected") %>% 
      graph.laplacian(normalized = T) %>% as.matrix
    
    if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda)) + tau*L_tilda)
    if(include.network == "L") K_regL <- L_tilda
    
    Zs <- scale(Zs); ZL <- Zs %*% K_regL
    if(rho == "median.pairwise") rho <- median(dist(Zs))
    
    if(kernel == "G") K <- Gaussian_kernel(rho, ZL)
    if(kernel == "L") K <- plyKern(ZL, pow = 1, rho = 0)
    
    dav_pval <- SKAT.c(H0.form, data = dd, K=K)$Qq
    
  } else{
    dav_pval <- V.Y <- V.e <- NA
  }
  c(pVal = dav_pval, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

Perf_Davie_DiffDens <- function(graph, H0.form, data, b0, 
                                sd.y, zz, tau, new.edge.prob,
                                kernel = "G",
                                rho = "median.pairwise", include.network = "R"){
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
    mX <- as.matrix(cbind(1, data))
    Z <- mvrnorm(n=nrow(mX), rep(0,p), solve(Omega1))
    Y <- rnorm(n=nrow(mX), mX%*%b0+Z%*%zz, sd.y)
    V.Y <- var(mX%*%b0+Z%*%zz); V.e <- var(Y - mX%*%b0+Z%*%zz)
    
    dd <- cbind(Y, data)
    L_tilda <- graph.laplacian(graph, normalized = T) %>% as.matrix
    
    if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda)) + tau*L_tilda)
    if(include.network == "L") K_regL <- L_tilda
    
    Z <- scale(Z); ZL <- Z %*% K_regL
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    if(kernel == "G") K <- Gaussian_kernel(rho, ZL)
    if(kernel == "L") K <- plyKern(ZL, pow = 1, rho = 0)
    
    dav_pval <- SKAT.c(H0.form, data = dd, K=K)$Qq
    
  } else{
    dav_pval <- V.Y <- V.e <- NA
  }
  c(pVal = dav_pval, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

## A function built for power analysis on prebuilt data sets
Perf_Davie_SigNoise <- function(YZ, graph, H0.form, data, 
                                include.network = "R", .tau){
  
  L_tilda <- graph.laplacian(graph, normalized = T) %>% as.matrix
  
  if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda)) + .tau*L_tilda)
  if(include.network == "L") K_regL <- L_tilda
  
  Z <- scale(YZ$Z); rho <- median(dist(Z))
  ZL <- Z %*% K_regL
  K <- Gaussian_kernel(rho, ZL)
  
  dd <- cbind(Y = YZ$Y, data)
  SKAT.c(H0.form, data = dd, K=K)$Qq
}


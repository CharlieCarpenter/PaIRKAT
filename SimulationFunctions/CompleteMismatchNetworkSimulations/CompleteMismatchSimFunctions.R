###############################
##
## Project: MetaoGuru Aim 3
##
## Purpose: Simulation functions for complete network mismatches
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2020-03-04
##
## ---------------------------
## Notes:
## Forces no shared edges
##
## ---------------------------

source("PaIRKAT/helpers.R")

## Forces no shared edges
Comp_Davie_SameSize <- function(graph, H0.form, data, 
                                b0, sd.y, zz, tau, kernel = "G",
                                rho = "median.pairwise", include.network = "R"){
  
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  p <- nrow(A)
  mX <- as.matrix(cbind(1, data))
  
  ## Convert adjacency matrix into a "starter" precision matrix
  ## Adding transpose to make symmetric? (undirected?)
  test <- A + diag(p)
  
  ## New adjacency matrix from new sample graph
  gg2 <- sample_pa(p, directed = F)
  A_p <- as.matrix( get.adjacency(gg2) )
  
  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)
  if(pd){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=nrow(mX), rep(0,p), solve(Omega1))
    Y <- rnorm(n=nrow(mX), mX%*%b0+Z%*%zz, sd.y)
    V.Y <- var(mX%*%b0+Z%*%zz); V.e <- var(Y - mX%*%b0+Z%*%zz)
    dd <- cbind(Y, data)
    
    ## making symmetric adjacency matrix? could result in 2s?
    total_edge_n <- sum(A==1)
    
    nse <- count_shared_edges(A, A_p)
    L_tilda_p <- nomatch(A, A_p, nse) %>% ## Forcing no match
      graph_from_adjacency_matrix(mode = "undirected") %>%
      graph.laplacian(normalized = T) %>% as.matrix
    
    if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda_p)) + tau*L_tilda_p)
    if(include.network == "L") K_regL <- L_tilda_p
    
    Z <- scale(Z); ZL <- Z %*% K_regL
    
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    if(kernel == "G") K <- Gaussian_kernel(rho, ZL)
    if(kernel == "L") K <- plyKern(ZL, pow = 1, rho = 0)
    
    pVal <- SKAT.c(H0.form, data = dd, K=K)$Qq
    
  } else{
    pVal <- V.Y <- V.e <- NA
  }
  c(pVal=pVal, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

Comp_Davie_SmallGraph <- function(graph, H0.form, data, 
                                  b0, sd.y, zz, tau, kernel = "G",
                                  rho = "median.pairwise", include.network = "R"){
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  p <- nrow(A)
  mX <- as.matrix(cbind(1, data))
  
  ## Convert adjacency matrix into a "starter" precision matrix
  ## Adding transpose to make symmetric? (undirected?)
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
    dd <- cbind(Y, data)
    
    ## making symmetric adjacency matrix? could result in 2s?
    total_edge_n <- sum(A==1)
    
    ## Shrinking graphs
    nd <- new.degs(graph)
    As <- A[nd > 0, nd > 0]; Zs <- Z[ , nd > 0]
    ## New adjacency matrix from new sample graph
    As_ng <- sample_pa(sum(nd>0), directed = F) %>%
      get.adjacency %>% as.matrix
    
    ## shared edges
    nse <- count_shared_edges(As, As_ng)
    L_tilda_s <- nomatch(As, As_ng, nse) %>% ## Forcing no matching edges
      graph_from_adjacency_matrix(mode = "undirected") %>%
      graph.laplacian(normalized = T) %>% as.matrix
    
    ## Redoing until we get a non-singular new matrix
    ## (fully connected graph)
    i <- 1
    while(is.singular.matrix(diag(nrow(L_tilda_s)) + tau*L_tilda_s) &
          i < 101){ ## only gonna try 50 new plots
      
      L_tilda_s <- nomatch(As, As_ng, nse) %>% ## Forcing no matching edges
        graph_from_adjacency_matrix(mode = "undirected") %>%
        graph.laplacian(normalized = T) %>% as.matrix
      
      i <- i+1
    }
    
    non.singular <- !is.singular.matrix(diag(nrow(L_tilda_s)) + tau*L_tilda_s)
    if(non.singular){
      
      if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda_s)) + tau*L_tilda_s)
      if(include.network == "L") K_regL <- L_tilda_s
      
      Zs <- scale(Zs); ZL <- Zs %*% K_regL
      
      if(rho == "median.pairwise") rho <- median(dist(Z))
      
      if(kernel == "G") K <- Gaussian_kernel(rho, ZL)
      if(kernel == "L") K <- plyKern(ZL, pow = 1, rho = 0)
      
      pVal <- SKAT.c(H0.form, data = dd, K=K)$Qq
      
    } else pVal <- V.Y <- V.e <- NA
    
  }else{
    pVal <- V.Y <- V.e <- NA
  }
  c(pVal=pVal, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

Comp_Davie_DiffDens <- function(graph, H0.form, data,
                                b0, sd.y, zz, tau, kernel = "G",
                                new.edge.prob,
                                rho = "median.pairwise", include.network = "R"){
  
  edge_density1 <- edge_density(graph)
  mX <- as.matrix(cbind(1, data))
  
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  A_d <- increase.edge.density(A, edge.prob = new.edge.prob)
  ## Convert adjacency matrix into a "starter" precision matrix
  ## Adding transpose to make symmetric? (undirected?)
  test <- A_d + diag(p)
  
  ## New adjacency matrix from new sample graph
  gg2 <- sample_pa(p, directed = F)
  A_p <- as.matrix( get.adjacency(gg2) )
  
  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)
  if(pd){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=nrow(mX), rep(0,p), solve(Omega1))
    Y <- rnorm(n=nrow(mX), mX%*%b0+Z%*%zz, sd.y)
    V.Y <- var(mX%*%b0+Z%*%zz); V.e <- var(Y - mX%*%b0+Z%*%zz)
    dd <- cbind(Y, data)
    
    ## making symmetric adjacency matrix? could result in 2s?
    total_edge_n <- sum(A==1)
    
    nse <- count_shared_edges(A_d, A_p)
    
    g_p <- nomatch(A, A_p, nse) %>% ## Forcing no match
      graph_from_adjacency_matrix(mode = "undirected") 
    
    L_tilda_p <- graph.laplacian(g_p, normalized = T) %>% as.matrix
    
    if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda_p)) + tau*L_tilda_p)
    if(include.network == "L") K_regL <- L_tilda_p
    
    Z <- scale(Z); ZL <- Z %*% K_regL
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    if(kernel == "G") K <- Gaussian_kernel(rho, ZL)
    if(kernel == "L") K <- plyKern(ZL, pow = 1, rho = 0)
    
    pVal <- SKAT.c(H0.form, data = dd, K=K)$Qq
    
  } else{
    pVal <- V.Y <- V.e <- NA
  }
  c(pVal=pVal, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

## A function built for power analysis on prebuilt data sets
Comp_Davie_SigNoise <- function(YZ, graph, H0.form, data, 
                                include.network = "R", .tau){
  
  A <- as.matrix(get.adjacency(graph))
  gg2 <- sample_pa(ncol(A), directed = F)
  A_p <- as.matrix(get.adjacency(gg2))
  nse <- count_shared_edges(A, A_p)
  
  L_tilda_p <- nomatch(A, A_p, nse) %>% ## Forcing no match
    graph_from_adjacency_matrix(mode = "undirected") %>%
    graph.laplacian(normalized = T) %>% as.matrix
  
  if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda_p)) + .tau*L_tilda_p)
  if(include.network == "L") K_regL <- L_tilda_p
  
  Z <- scale(YZ$Z); rho <- median(dist(Z))
  Z <- Z %*% K_regL
  K <- Gaussian_kernel(rho, Z)
  
  dd <- cbind(Y = YZ$Y, data)
  SKAT.c(H0.form, data = dd, K=K)$Qq
}



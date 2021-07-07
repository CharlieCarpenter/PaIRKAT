###############################
##
## Project:
##
## Purpose:
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2020-03-04
##
## ---------------------------
## Notes:
##
## Doesn't include network with Z
## ---------------------------

source("PaIRKAT/helpers.R")

NoNet_Davie_SameSize <- function(graph, H0.form, data, b0,
                                 sd.y, zz, tau, kernel = "G",
                                 rho = "median.pairwise"){
  
  ## Converts a graph into an adjacency matrix
  adj <- as.matrix( get.adjacency(graph) )
  p <- nrow(adj); total_edge_n <- sum(adj == 1)
  ## Convert adjacency matrix into a "starter" precision matrix
  test <- adj + diag(p)
  
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
    
    # no network information
    # L <- graph.laplacian(graph, normalize = T)
    # K_regL <- solve(diag(p) + tau*L)
    
    Z <- scale(Z)
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    # no network information
    # ZL <- Z %*% K_regL
    
    if(kernel == "G") K <- Gaussian_kernel(rho, Z)
    if(kernel == "L") K <- plyKern(Z, pow = 1, rho = 0)
    
    dd <- cbind(Y, data)
    dav_pval <- SKAT.c(H0.form, data = dd, K=K)$Qq
    
  } else{
    dav_pval <- V.Y <- V.e <- NA
  }
  c(pVal = dav_pval, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

NoNet_Davie_SmallGraph <- function(graph, H0.form, data, b0,
                                   sd.y, zz, tau, kernel = "G",
                                   rho = "median.pairwise"){
  
  ## Converts a graph into an adjacency matrix
  adj <- as.matrix( get.adjacency(graph) )
  p <- nrow(adj); total_edge_n <- sum(adj == 1)
  ## Convert adjacency matrix into a "starter" precision matrix
  test <- adj + diag(p)
  
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
    
    ## Shrinking graphs / what we believe our pathway is
    nd <- new.degs(graph)
    #As <- A[nd > 0, nd > 0]
    Zs <- Z[ , nd > 0]
    
    # no network information
    # K_regL <- solve(diag(p) + s2*L_tilda)
    # Lc <- chol(K_regL)
    
    Z <- scale(Z)
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    # no network information
    # Z <- Z %*% t(Lc)
    
    if(kernel == "G") K <- Gaussian_kernel(rho, Z)
    if(kernel == "L") K <- plyKern(Z, pow = 1, rho = 0)
    
    dd <- cbind(Y, data)
    dav_pval <- SKAT.c(H0.form, data = dd, K=K)$Qq
    
  } else{
    dav_pval <- V.Y <- V.e <- NA
  }
  c(pVal = dav_pval, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

NoNet_Davie_DiffDens <- function(graph, H0.form, data, b0,
                                 sd.y, zz, tau, kernel = "G",
                                 new.edge.prob, rho = "median.pairwise"){
  ## Converts a graph into an adjacency matrix
  adj <- as.matrix( get.adjacency(graph) )
  p <- nrow(adj); total_edge_n <- sum(adj == 1)
  adj <- increase.edge.density(adj, edge.prob = new.edge.prob)
  ## Convert adjacency matrix into a "starter" precision matrix
  test <- adj + diag(p)
  
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
    
    # no network information
    # K_regL <- solve(diag(p) + s2*L_tilda)
    # Lc <- chol(K_regL)
    
    Z <- scale(Z)
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    if(kernel == "G") K <- Gaussian_kernel(rho, ZL)
    if(kernel == "L") K <- plyKern(ZL, pow = 1, rho = 0)
    
    dd <- cbind(Y, data)
    dav_pval <- SKAT.c(H0.form, data = dd, K=K)$Qq
    
  } else{
    dav_pval <- V.Y <- V.e <- NA
  }
  c(pVal = dav_pval, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

## A function built for power analysis on prebuilt data sets
NoNet_Davie_SigNoise <- function(YZ, graph, H0.form, data){
  Z <- scale(YZ$Z); rho <- median(dist(Z))
  K <- Gaussian_kernel(rho, Z)
  # print(head(YZ$Y)); print(class(YZ$Y))
  
  dd <- cbind(Y = YZ$Y, data)
  SKAT.c(H0.form, data = dd, K=K)$Qq
}


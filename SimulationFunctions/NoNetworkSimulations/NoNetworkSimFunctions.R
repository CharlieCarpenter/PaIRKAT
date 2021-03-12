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

# Score -------------------------------------------------------------------
NoNet_Score_SameSize <- function(graph, mX, b0, sd.y, zz, delta,
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
    Z <- mvrnorm(n=nrow(mX), rep(0,p), solve(Omega1))
    Y <- rnorm(n=nrow(mX), mX%*%b0+Z%*%zz, sd.y)
    V.Y <- var(mX%*%b0+Z%*%zz); V.e <- var(Y - mX%*%b0+Z%*%zz)
    
    # no network information
    # K_regL <- solve(diag(p) + s2*L_tilda)
    # Lc <- chol(K_regL)
    
    Z <- scale(Z)
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    # no network information
    # Z <- Z %*% t(Lc)
    
    K <- Gaussian_kernel(rho, Z)
    sc <- CC_Chisq_Score(K, Y, mX)
  } else{
    sc <- V.Y <- V.e <- NA
  }
  c(sc, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

NoNet_Score_SmallGraph <- function(graph, mX, b0, sd.y, zz, delta,
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
    
    K <- Gaussian_kernel(rho, Z)
    sc <- CC_Chisq_Score(K, Y, mX)
  } else{
    sc <- V.Y <- V.e <- NA
  }
  c(sc, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

NoNet_Score_DiffDens <- function(graph, mX, b0, sd.y, zz, delta,
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
    Z <- mvrnorm(n=nrow(mX), rep(0,p), solve(Omega1))
    Y <- rnorm(n=nrow(mX), mX%*%b0+Z%*%zz, sd.y)
    V.Y <- var(mX%*%b0+Z%*%zz); V.e <- var(Y - mX%*%b0+Z%*%zz)
    
    # no network information
    # K_regL <- solve(diag(p) + s2*L_tilda)
    # Lc <- chol(K_regL)
    
    Z <- scale(Z)
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    K <- Gaussian_kernel(rho, Z)
    sc <- CC_Chisq_Score(K, Y, mX)
  } else{
    sc <- V.Y <- V.e <- NA
  }
  c(sc, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

## A function built for power analysis on prebuilt data sets
NoNet_Score_SigNoise <- function(YZ, graph, .X){
  Z <- scale(YZ$Z); rho <- median(dist(Z))
  
  K <- Gaussian_kernel(rho, Z)
  sc <- CC_Chisq_Score(K, YZ$Y, .X)
  
  sc["pVal"]
}

# Permutations ------------------------------------------------------------
NoNetwork_SameSize <- function(graph, n, X, p, b0, sd.y, betas, delta, rho = "median.pairwise"){

  ## Converts a graph into an adjacency matrix
  adj <- as.matrix( get.adjacency(graph) )

  ## Convert adjacency matrix into a "starter" precision matrix
  test <- adj + diag(p)

  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)

  if(pd){
    ## Outcome and model.matrix
    Sigma <- solve(Omega1)
    Z <- mvrnorm(n=n, rep(0,p), Sigma)
    model <- model.matrix(~ 1 , data.frame(Z))
    Y <- rnorm(n, b0+Z%*%betas, sd.y)

    # no network information
    # K_regL <- solve(diag(p) + s2*L_tilda)
    # Lc <- chol(K_regL)

    Z <- scale(Z)
    if(rho == "median.pairwise") rho <- median(dist(Z))

    # no network information
    # Z <- Z %*% t(Lc)

    V.Y <- var(Y); V.e <- var(Y - b0+Z%*%betas)
    K <- Gaussian_kernel(rho, Z)

    perm.test <- Permutation_test(K, Y, X, nperm)

    Sc_stat <- perm.test["Stat.sc"]
    Pval <- perm.test["p_value"]
  }else{
    Sc_stat <- Pval <- total_edge_n <- shared_edge_n <- V.Y <- V.e <-  NA
  }

  c(Pval = Pval, Sc_stat = Sc_stat, edge_density1 = edge_density(graph),
    pos_def = pd, V.Y = V.Y, V.e = V.e)
}

##
NoNetwork_SmallGraph <- function(graph, n, X, p, b0, sd.y, betas, delta, rho = "median.pairwise"){

  ## Converts a graph into an adjacency matrix
  adj <- as.matrix( get.adjacency(graph) )

  ## Convert adjacency matrix into a "starter" precision matrix
  test <- adj + diag(p)

  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)

  if(pd){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=n, rep(0,p), solve(Omega1))
    model <- model.matrix(~ 1 , data.frame(Z))
    Y <- rnorm(n, b0+Z%*%betas, sd.y)

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

    K <- Gaussian_kernel(rho, Z)

    perm.test <- Permutation_test(K, Y, X, nperm)

    Sc_stat <- perm.test["Stat.sc"]
    Pval <- perm.test["p_value"]
  }else{
    Sc_stat <- Pval <- total_edge_n <- shared_edge_n <- NA
  }

  c(Pval = Pval, Sc_stat = Sc_stat, edge_density1 = edge_density(graph),
    pos_def = pd)
}

## Changing network desnities
NoNetwork_DiffDense <- function(graph, n, X, p, b0, sd.y, betas, new.edge.prob,
                               delta, rho = "median.pairwise"){
  
  ## Converts a graph into an adjacency matrix
  adj <- as.matrix( get.adjacency(graph) )
  adj <- increase.edge.density(adj, edge.prob = new.edge.prob)
  
  ## Convert adjacency matrix into a "starter" precision matrix
  test <- adj + diag(p)
  
  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)
  
  if(pd){
    ## Outcome and model.matrix
    Sigma <- solve(Omega1)
    Z <- mvrnorm(n=n, rep(0,p), Sigma)
    model <- model.matrix(~ 1 , data.frame(Z))
    Y <- rnorm(n, b0+Z%*%betas, sd.y)
    
    # no network information
    # K_regL <- solve(diag(p) + s2*L_tilda)
    # Lc <- chol(K_regL)
    
    Z <- scale(Z)
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    # no network information
    # Z <- Z %*% t(Lc)
    
    K <- Gaussian_kernel(rho, Z)
    
    perm.test <- Permutation_test(K, Y, X, nperm)
    
    Sc_stat <- perm.test["Stat.sc"]
    Pval <- perm.test["p_value"]
  }else{
    Sc_stat <- Pval <- total_edge_n <- shared_edge_n <- NA
  }
  
  c(Pval = Pval, Sc_stat = Sc_stat, edge_density1 = edge_density(graph),
    pos_def = pd)
}

## A function built for power analysis on prebuilt data sets
No.Network.Sig.Noise <- function(YZ, graph, .X, .nperm){
  Z <- scale(YZ$Z); rho <- median(dist(Z))
  
  K <- Gaussian_kernel(rho, Z)
  perm.test <- Permutation_test(K, YZ$Y, .X, .nperm)
  
  perm.test["p_value"]
}


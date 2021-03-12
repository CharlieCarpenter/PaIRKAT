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

source('~/Documents/Research/Current/MetaboGuru/Carpenter/RCode/ScoreStat/ScoreSimFunctions.R')

# Score -------------------------------------------------------------------
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

# Permutations ------------------------------------------------------------
PerfectNetwork_SameSize <- function(graph, n, X, p, b0, sd.y, betas, delta,
                                    rho = "median.pairwise", include.network = "full"){

  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )

  ## Convert adjacency matrix into a "starter" precision matrix
  test <- A + diag(p)

  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)
  if(pd){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=n, rep(0,p), solve(Omega1))
    model <- model.matrix(~ 1 , data.frame(Z))
    Y <- rnorm(n, X%*%betas+Z%*%zz, sd.y)
    V.Y <- var(X%*%betas+Z%*%zz); V.e <- var(Y - X%*%betas+Z%*%zz)
    
    total_edge_n <- sum(A==1)
    
    L_tilda <- graph.laplacian(graph, normalized = T) %>% as.matrix
    
    if(include.network == "R") K_regL <- solve(diag(p) + delta*L_tilda)
    if(include.network == "L") K_regL <- L_tilda
    
    Z <- scale(Z); ZL <- Z %*% K_regL
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    K <- Gaussian_kernel(rho, ZL)
    perm.test <- Permutation_test(K, Y, X, nperm)

    V.e.net <- var(Y - b0+Z%*%betas)
    Sc_stat <- perm.test["Stat.sc"]
    Pval <- perm.test["p_value"]
  }else{
    Sc_stat <- Pval <- total_edge_n <- shared_edge_n <- V.Y <- V.e <- V.e.net <- NA
  }

  c(Pval = Pval, Sc_stat = Sc_stat, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e, V.e.net = V.e.net,
    total_edge_n = total_edge_n, pos_def = pd)
}

## Smaller testing (missing nodes)
PerfectNetwork_SmallGraph <- function(graph, n, X, p, b0, sd.y, betas, delta,
                                      rho = "median.pairwise", include.network = "full"){

  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )

  ## Convert adjacency matrix into a "starter" precision matrix
  test <- A + diag(p)

  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)

  if(pd){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=n, rep(0,p), solve(Omega1))
    model <- model.matrix(~ 1 , data.frame(Z))
    Y <- rnorm(n, b0+Z%*%betas, sd.y)

    ## Shrinking graphs
    nd <- new.degs(graph)
    As <- A[nd > 0, nd > 0]
    Zs <- Z[ , nd > 0]

    total_edge_n <- sum(A==1)

    L_tilda <- As %>% 
      graph_from_adjacency_matrix(mode = "undirected") %>% 
      graph.laplacian(normalized = T) %>% as.matrix
    
    if(include.network == "R") K_regL <- solve(diag(p) + delta*L_tilda)
    if(include.network == "L") K_regL <- L_tilda
    
    Z <- scale(Z); ZL <- Z %*% K_regL
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    K <- Gaussian_kernel(rho, Zs)
    perm.test <- Permutation_test(K, Y, X, nperm)
  
    Sc_stat <- perm.test["Stat.sc"]
    Pval <- perm.test["p_value"]
  }else{
    Sc_stat <- Pval <- total_edge_n <- NA
  }

  c(Pval = Pval, Sc_stat = Sc_stat, edge_density1 = edge_density(graph),
    total_edge_n = total_edge_n, pos_def = pd)
}

## Changing Network density
PerfectNetwork_DiffDens <- function(graph, n, X, p, b0, sd.y, betas, delta, new.edge.prob,
                                    rho = "median.pairwise", include.network = "full"){
  
  edge_density1 <- edge_density(graph)
  
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  
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
    Z <- mvrnorm(n=n, rep(0,p), solve(Omega1))
    model <- model.matrix(~ 1 , data.frame(Z))
    Y <- rnorm(n, b0+Z%*%betas, sd.y)
    
    total_edge_n <- sum(A==1)
    
    L_tilda <- graph.laplacian(graph, normalized = T) %>% as.matrix
    
    if(include.network == "R") K_regL <- solve(diag(p) + delta*L_tilda)
    if(include.network == "L") K_regL <- L_tilda
    
    Z <- scale(Z); ZL <- Z %*% K_regL
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    K <- Gaussian_kernel(rho, Z)
    perm.test <- Permutation_test(K, Y, X, nperm)
    
    Sc_stat <- perm.test["Stat.sc"]
    Pval <- perm.test["p_value"]
  }else{
    Sc_stat <- Pval <- total_edge_n <- shared_edge_n <- new.edge.density <-  NA
  }
  
  c(Pval = Pval, Sc_stat = Sc_stat, edge_density1 = edge_density1,
    total_edge_n = total_edge_n, pos_def = pd, new.edge.density = new.edge.density)
}

## Changing Network density: Testing on low density graph
PerfectNetwork_DiffDens_Sm <- function(graph, n, X, p, b0, sd.y, betas, delta, new.edge.prob,
                                       rho = "median.pairwise", include.network = "full"){
  
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  
  A_big <- increase.edge.density(A, edge.prob = new.edge.prob)
  
  new.edge.density <- A_big %>% 
    graph_from_adjacency_matrix(mode = "undirected") %>% 
    edge_density
  
  ## Convert adjacency matrix into a "starter" precision matrix
  test <- A_big + diag(p)
  
  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)
  if(pd){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=n, rep(0,p), solve(Omega1))
    model <- model.matrix(~ 1 , data.frame(Z))
    Y <- rnorm(n, b0+Z%*%betas, sd.y)
    
    total_edge_n <- sum(A==1)
    
    L_tilda <- graph.laplacian(graph, normalized = T) %>% as.matrix
    
    if(include.network == "R") K_regL <- solve(diag(p) + delta*L_tilda)
    if(include.network == "L") K_regL <- L_tilda
    
    Z <- scale(Z); ZL <- Z %*% K_regL
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    K <- Gaussian_kernel(rho, Z)
    perm.test <- Permutation_test(K, Y, X, nperm)
    
    Sc_stat <- perm.test["Stat.sc"]
    Pval <- perm.test["p_value"]
  }else{
    Sc_stat <- Pval <- total_edge_n <- shared_edge_n <- new.edge.density <-  NA
  }
  
  c(Pval = Pval, Sc_stat = Sc_stat, edge_density = edge_density(graph),
    total_edge_n = total_edge_n, pos_def = pd, new.edge.density = new.edge.density)
}

## A function built for power analysis on prebuilt data sets
Perfect.Network.Sig.Noise <- function(YZ, graph, .delta, .X, .nperm){
  
  L_tilda <- graph.laplacian(graph, normalized = T) %>% as.matrix
  Z <- scale(YZ$Z); rho <- median(dist(Z))
  K_regL <- solve(diag(ncol(L_tilda)) + delta*L_tilda); Z <- Z %*% K_regL
  
  K <- Gaussian_kernel(rho, Z)
  perm.test <- Permutation_test(K, YZ$Y, .X, .nperm)
  
  perm.test["p_value"]
}

## Didn't really work. TypeI error rates were way off
# Cumulants ---------------------------------------------------------------
Perf_Cumu_SameSize <- function(graph, mX, b0, sd.y, zz, delta,
                               rho = "median.pairwise", include.network = "full"){
  
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
    K_regL <- solve(diag(p) + delta*L_tilda)
    Z <- scale(Z)
    
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    ## How to include smoothed network
    if(include.network == "full") ZL <- Z %*% K_regL
    if(include.network == "chol.decomp") {Lc <- chol(K_regL); ZL <- Z %*% t(Lc)}
    
    K <- Gaussian_kernel(rho, ZL)
    cu <- CC_cumuMatch(K, Y, mX)
  } else{
    cu <- V.Y <- V.e <- NA
  }
  c(cu, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

Perf_Cumu_SmallGraph <- function(graph, mX, b0, sd.y, zz, delta,
                                 rho = "median.pairwise", include.network = "full"){
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
    
    K_regL <- solve(diag(nrow(As)) + delta*L_tilda)
    Zs <- scale(Zs)
    
    if(rho == "median.pairwise") rho <- median(dist(Zs))
    
    ## How to include smoothed network
    if(include.network == "full") Zs <- Zs %*% K_regL
    if(include.network == "lower.triangle") {Lc <- chol(K_regL); Zs <- Zs %*% t(Lc)}
    
    K <- Gaussian_kernel(rho, Zs)
    cu <- CC_cumuMatch(K, Y, mX)
  } else{
    cu <- V.Y <- V.e <- NA
  }
  c(cu, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

Perf_Cumu_DiffDens <- function(graph, mX, b0, sd.y, zz, delta, new.edge.prob,
                               rho = "median.pairwise", include.network = "full"){
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
    K_regL <- solve(diag(p) + delta*L_tilda)
    Z <- scale(Z)
    
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    ## How to include smoothed network
    if(include.network == "full") Z <- Z %*% K_regL
    if(include.network == "lower.triangle") {Lc <- chol(K_regL); Z <- Z %*% t(Lc)}
    
    K <- Gaussian_kernel(rho, Z)
    cu <- CC_cumuMatch(K, Y, mX)
  } else{
    cu <- V.Y <- V.e <- NA
  }
  c(cu, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}





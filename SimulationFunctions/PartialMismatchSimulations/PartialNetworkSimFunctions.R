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

# Score -------------------------------------------------------------------

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

# Permutation -------------------------------------------------------------

## forces a change in specified number of edges (n_diff_edge)
partialNetwork_SameSize <- function(graph, n, X, p, perc.perm, # Decimals
                                    b0, sd.y, betas, delta, rho = "median.pairwise",
                                    include.network = "full"){

  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )

  ## Convert adjacency matrix into a "starter" precision matrix
  test <- A + diag(p)

  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)

  if(is.positive.definite(Omega1)){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=n, rep(0,p), solve(Omega1))
    model <- model.matrix(~ 1 , data.frame(Z))
    Y <- rnorm(n, b0+Z%*%betas, sd.y)
    V.Y <- var(Y); V.e <- var(Y - b0+Z%*%betas)
    
    ## making symmetric adjacency matrix? could result in 2s?
    total_edge_n <- sum(A==1)

    ## Graph with permuted edges
    A_p <- edge.perm(A, perc.perm)
    nse <- count_shared_edges(A, A_p)
    
    L_tilda_p <- A_p %>% 
      graph_from_adjacency_matrix(mode = "undirected") %>% 
      graph.laplacian(normalized = T) %>% as.matrix
    
    K_regL <- solve(diag(p) + delta*L_tilda_p)
    Z <- scale(Z)
    
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    ## How to include smoothed network
    if(include.network == "full") Z <- Z %*% K_regL
    if(include.network == "lower.triangle") {Lc <- chol(K_regL); Z <- Z %*% t(Lc)}
    
    K <- Gaussian_kernel(rho, Z)
    perm.test <- Permutation_test(K, Y, X, nperm)

    V.e.net <- var(Y - b0+Z%*%betas)
    Sc_stat <- perm.test["Stat.sc"]
    Pval <- perm.test["p_value"]
  }else{
    Sc_stat <- Pval <- total_edge_n <- nse <- V.Y <- V.e <- V.e.net <- NA
  }

  c(Pval = Pval, Sc_stat = Sc_stat, edge_density1 = edge_density(graph),
    total_edge_n = total_edge_n, nse = nse,
    V.Y = V.Y, V.e = V.e, V.e.net = V.e.net,
    pos_def = is.positive.definite(Omega1))
}


## performs same tests, but changes density of testing graph
partialNetwork_SmallGraph <- function(graph, n, X, p, perc.perm,
                                      b0, sd.y, betas, delta, rho = "median.pairwise",
                                      include.network = "full"){

  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )

  ## Convert adjacency matrix into a "starter" precision matrix
  test <- A + diag(p)

  ## New adjacency matrix from new sample graph
  # gg2 <- sample_pa(pn)
  # A_p <- as.matrix( get.adjacency(gg2) )

  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)

  if(pd){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=n, rep(0,p), solve(Omega1))
    model <- model.matrix(~ 1, data.frame(Z))
    Y <- rnorm(n, b0+Z%*%betas, sd.y)

    ## Shrinking graphs
    nd <- new.degs(graph)
    As <- A[nd > 0, nd > 0]
    Zs <- Z[ , nd > 0]

    ## Graph with permuted edges
    A_p <- edge.perm(As, perc.perm)

    L_tilda_p <- A_p %>% 
      graph_from_adjacency_matrix(mode = "undirected") %>% 
      graph.laplacian(normalized = T) %>% as.matrix
    
    K_regL <- solve(diag(nrow(L_tilda_p)) + delta*L_tilda_p)
    Zs <- scale(Zs)
    
    if(rho == "median.pairwise") rho <- median(dist(Zs))
    
    ## How to include smoothed network
    if(include.network == "full") Zs <- Zs %*% K_regL
    if(include.network == "lower.triangle") {Lc <- chol(K_regL); Zs <- Zs %*% t(Lc)}
    
    K <- Gaussian_kernel(rho, Zs)
    perm.test <- Permutation_test(K, Y, X, nperm)

    Sc_stat <- perm.test["Stat.sc"]
    Pval <- perm.test["p_value"]
  }else{
    Sc_stat <- Pval <- total_edge_n <- nse <- NA
  }

  c(Pval = Pval, Sc_stat = Sc_stat, edge_density1 = edge_density(graph), pos_def = pd)
}


## Changes network density
partialNetwork_DiffDens <- function(graph, n, X, p, perc.perm, # Decimals
                                    b0, sd.y, betas, delta, rho = "median.pairwise",
                                    new.edge.prob,
                                    include.network = "full"){
  
  edge_density1 <- edge_density(graph)

  ## Converts a graph into an adjacency matrix
  A <- get.adjacency(graph) %>% as.matrix %>% 
    increase.edge.density(edge.prob = new.edge.prob)
  
  graph <- graph_from_adjacency_matrix(A, mode = "undirected")
  
  new.edge.dens <- edge_density(graph)
  
  ## Convert adjacency matrix into a "starter" precision matrix
  test <- A + diag(p)
  
  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  
  if(is.positive.definite(Omega1)){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=n, rep(0,p), solve(Omega1))
    model <- model.matrix(~ 1 , data.frame(Z))
    Y <- rnorm(n, b0+Z%*%betas, sd.y)
    
    ## making symmetric adjacency matrix? could result in 2s?
    total_edge_n <- sum(A==1)
    
    ## Graph with permuted edges
    A_p <- edge.perm(A, perc.perm)
    nse <- count_shared_edges(A, A_p)
    
    L_tilda_p <- A_p %>% 
      graph_from_adjacency_matrix(mode = "undirected") %>% 
      graph.laplacian(normalized = T) %>% as.matrix
    
    K_regL <- solve(diag(p) + delta*L_tilda_p)
    Z <- scale(Z)
    
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    ## How to include smoothed network
    if(include.network == "full") Z <- Z %*% K_regL
    if(include.network == "lower.triangle") {Lc <- chol(K_regL); Z <- Z %*% t(Lc)}
    
    K <- Gaussian_kernel(rho, Z)
    perm.test <- Permutation_test(K, Y, X, nperm)
    
    Sc_stat <- perm.test["Stat.sc"]
    Pval <- perm.test["p_value"]
  }else{
    Sc_stat <- Pval <- total_edge_n <- nse <- NA
  }
  
  c(Pval = Pval, Sc_stat = Sc_stat, edge_density1 = edge_density1,
    new.edge.dens = new.edge.dens,
    total_edge_n = total_edge_n, nse = nse,
    pos_def = is.positive.definite(Omega1))
}


## Creates Z from higher density and tests on lower density graph
partialNetwork_DiffDens_Sm <- function(graph, n, X, p, perc.perm, # Decimals
                                    b0, sd.y, betas, delta, rho = "median.pairwise",
                                    new.edge.prob,
                                    include.network = "full"){
 
   ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  A_big <- increase.edge.density(A, edge.prob = new.edge.prob)
  
  new.edge.dens <- A_big %>% 
    graph_from_adjacency_matrix(mode = "undirected") %>% 
    edge_density()
  
  ## Convert adjacency matrix into a "starter" precision matrix
  test <- A_big + diag(p)
  
  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  
  if(is.positive.definite(Omega1)){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=n, rep(0,p), solve(Omega1))
    model <- model.matrix(~ 1 , data.frame(Z))
    Y <- rnorm(n, b0+Z%*%betas, sd.y)
    
    ## making symmetric adjacency matrix? could result in 2s?
    total_edge_n <- sum(A==1)
    
    ## Graph with permuted edges
    A_p <- edge.perm(A, perc.perm)
    nse <- count_shared_edges(A, A_p)
    
    L_tilda_p <- A_p %>% 
      graph_from_adjacency_matrix(mode = "undirected") %>% 
      graph.laplacian(normalized = T) %>% as.matrix
    
    K_regL <- solve(diag(p) + delta*L_tilda_p)
    Z <- scale(Z)
    
    if(rho == "median.pairwise") rho <- median(dist(Z))
    
    ## How to include smoothed network
    if(include.network == "full") Z <- Z %*% K_regL
    if(include.network == "lower.triangle") {Lc <- chol(K_regL); Z <- Z %*% t(Lc)}
    
    K <- Gaussian_kernel(rho, Z)
    perm.test <- Permutation_test(K, Y, X, nperm)
    
    Sc_stat <- perm.test["Stat.sc"]
    Pval <- perm.test["p_value"]
  }else{
    Sc_stat <- Pval <- total_edge_n <- nse <- NA
  }
  
  c(Pval = Pval, Sc_stat = Sc_stat, edge_density1 = edge_density(graph),
    total_edge_n = total_edge_n, nse = nse,
    new.edge.dens = new.edge.dens,
    pos_def = is.positive.definite(Omega1))
}

## A function built for power analysis on prebuilt data sets
Partial.Mismatch.Sig.Noise <- function(YZ, graph, .perc.perm, .X, .delta, .nperm){
  
  A <- as.matrix(get.adjacency(graph))
  A_p <- edge.perm(A, .perc.perm)
  
  L_tilda <- A_p %>% 
    graph_from_adjacency_matrix(mode = "undirected") %>% 
    graph.laplacian(normalized = T) %>% as.matrix
  
  Z <- scale(YZ$Z); rho <- median(dist(Z))
  K_regL <- solve(diag(ncol(L_tilda)) + .delta*L_tilda); Z <- Z %*% K_regL
  
  K <- Gaussian_kernel(rho, Z)
  perm.test <- Permutation_test(K, YZ$Y, .X, .nperm)
  
  perm.test["p_value"]
}


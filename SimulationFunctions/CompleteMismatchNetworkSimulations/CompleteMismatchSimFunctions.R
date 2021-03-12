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

# Score -------------------------------------------------------------------
## Forces no shared edges
Comp_Score_SameSize <- function(graph, mX, b0, sd.y, zz, delta,
                                  rho = "median.pairwise", include.network = "R"){
  
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  p <- nrow(A)
  
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
    
    ## making symmetric adjacency matrix? could result in 2s?
    total_edge_n <- sum(A==1)
    
    nse <- count_shared_edges(A, A_p)
    L_tilda_p <- nomatch(A, A_p, nse) %>% ## Forcing no match
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

Comp_Score_SmallGraph <- function(graph, mX, b0, sd.y, zz, delta,
                                    rho = "median.pairwise", include.network = "R"){
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  p <- nrow(A)
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
    while(is.singular.matrix(diag(nrow(L_tilda_s)) + delta*L_tilda_s) &
          i < 101){ ## only gonna try 50 new plots
      
      L_tilda_s <- nomatch(As, As_ng, nse) %>% ## Forcing no matching edges
        graph_from_adjacency_matrix(mode = "undirected") %>%
        graph.laplacian(normalized = T) %>% as.matrix
      
      i <- i+1
    }
    
    non.singular <- !is.singular.matrix(diag(nrow(L_tilda_s)) + delta*L_tilda_s)
    if(non.singular){
      
      if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda_s)) + delta*L_tilda_s)
      if(include.network == "L") K_regL <- L_tilda_s
      
      Zs <- scale(Zs); ZL <- Zs %*% K_regL
      
      if(rho == "median.pairwise") rho <- median(dist(Z))
      
      K <- Gaussian_kernel(rho, ZL)
      sc <- CC_Chisq_Score(K, Y, mX)
    } else sc <- V.Y <- V.e <- NA
    
  }else{
    sc <- V.Y <- V.e <- NA
  }
  c(sc, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e,
    total_edge_n = total_edge_n, pos_def = pd)
}

Comp_Score_DiffDens <- function(graph, mX, b0, sd.y, zz, delta, new.edge.prob,
                                rho = "median.pairwise", include.network = "R"){
  
  edge_density1 <- edge_density(graph)
  
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  A <- increase.edge.density(A, edge.prob = new.edge.prob)
  ## Convert adjacency matrix into a "starter" precision matrix
  ## Adding transpose to make symmetric? (undirected?)
  test <- A + diag(p)
  
  ## New adjacency matrix from new sample graph
  gg2 <- sample_pa(p, directed = F)
  A_p <- as.matrix( get.adjacency(gg2) ) %>% 
    increase.edge.density(edge.prob = new.edge.prob)
  
  ## Make test positive using approach of Danaher et al (2014)
  ## Function in SimulationFuncitons.R
  Omega1 <- Danaher_pos_def(test)
  pd <- is.positive.definite(Omega1)
  if(pd){
    ## Outcome and model.matrix
    Z <- mvrnorm(n=nrow(mX), rep(0,p), solve(Omega1))
    Y <- rnorm(n=nrow(mX), mX%*%b0+Z%*%zz, sd.y)
    V.Y <- var(mX%*%b0+Z%*%zz); V.e <- var(Y - mX%*%b0+Z%*%zz)
    
    ## making symmetric adjacency matrix? could result in 2s?
    total_edge_n <- sum(A==1)
    
    nse <- count_shared_edges(A, A_p)
    
    g_p <- nomatch(A, A_p, nse) %>% ## Forcing no match
      graph_from_adjacency_matrix(mode = "undirected") 
    
    L_tilda_p <- graph.laplacian(g_p, normalized = T) %>% as.matrix
    
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
Comp_Score_SigNoise <- function(YZ, graph, include.network = "R", .X, .delta){
  
  A <- as.matrix(get.adjacency(graph))
  gg2 <- sample_pa(ncol(A), directed = F)
  A_p <- as.matrix(get.adjacency(gg2))
  nse <- count_shared_edges(A, A_p)
  
  L_tilda_p <- nomatch(A, A_p, nse) %>% ## Forcing no match
    graph_from_adjacency_matrix(mode = "undirected") %>%
    graph.laplacian(normalized = T) %>% as.matrix

  if(include.network == "R") K_regL <- solve(diag(nrow(L_tilda_p)) + .delta*L_tilda_p)
  if(include.network == "L") K_regL <- L_tilda_p
  
  Z <- scale(YZ$Z); rho <- median(dist(Z))
  Z <- Z %*% K_regL
  
  K <- Gaussian_kernel(rho, Z)
  sc <- CC_Chisq_Score(K, YZ$Y, .X)
  
  sc["pVal"]
}

# Permutation -------------------------------------------------------------

## Changing network density
CompMismatch_DiffDense <- function(graph, n, X, p, b0, sd.y, betas, delta, new.edge.prob, 
                                  rho = "median.pairwise", include.network = "full"){
  
  edge_density1 <- edge_density(graph)
  
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  A <- increase.edge.density(A, edge.prob = new.edge.prob)
  ## Convert adjacency matrix into a "starter" precision matrix
  ## Adding transpose to make symmetric? (undirected?)
  test <- A + diag(p)
  
  ## New adjacency matrix from new sample graph
  gg2 <- sample_pa(p, directed = F)
  A_p <- as.matrix( get.adjacency(gg2) ) %>% 
    increase.edge.density(edge.prob = new.edge.prob)
  
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
    
    nse <- count_shared_edges(A, A_p)
    
    g_p <- nomatch(A, A_p, nse) %>% ## Forcing no match
      graph_from_adjacency_matrix(mode = "undirected") 
    
    L_tilda_p <- graph.laplacian(g_p, normalized = T) %>% as.matrix
    
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
    Sc_stat <- Pval <- total_edge_n <- nse <- V.Y <- V.e <- V.e.net <- NA
  }
  
  c(Pval = Pval, Sc_stat = Sc_stat, edge_density1 = edge_density1,
    new.edge.dens = edge_density(g_p),
    total_edge_n = total_edge_n, nse = nse, 
    pos_def = is.positive.definite(Omega1))
}

## Tests on original (lower density graph)
CompMismatch_DiffDense_Sm <- function(graph, n, X, p, b0, sd.y, betas, delta, new.edge.prob,
                                   rho = "median.pairwise", include.network = "full"){

  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )
  A_big <- increase.edge.density(A, edge.prob = new.edge.prob)
  ## Convert adjacency matrix into a "starter" precision matrix
  ## Adding transpose to make symmetric? (undirected?)
  test <- A_big + diag(p)

  ## New adjacency matrix from new sample graph
  gg2 <- sample_pa(p, directed = F)
  A_p <- as.matrix( get.adjacency(gg2) )

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

    nse <- count_shared_edges(A, A_p)

    ## Forcing no match between original low density graph and new low density graph
    L_tilda_p <- nomatch(A, A_p, nse) %>% 
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
    pos_def = is.positive.definite(Omega1))
}

## Forces no shared edges
CompMismatch_SameSize <- function(graph, n, X, p, b0, sd.y, betas, delta,
                                  rho = "median.pairwise", include.network = "full"){

  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )

  ## Convert adjacency matrix into a "starter" precision matrix
  ## Adding transpose to make symmetric? (undirected?)
  test <- A + diag(p)

  ## New adjacency matrix from new sample graph
  gg2 <- sample_pa(p, directed = F)
  A_p <- as.matrix( get.adjacency(gg2) )

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

    nse <- count_shared_edges(A, A_p)

    L_tilda_p <- nomatch(A, A_p, nse) %>% ## Forcing no match
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
    Sc_stat <- Pval <- total_edge_n <- nse <- NA
  }

  c(Pval = Pval, Sc_stat = Sc_stat, edge_density1 = edge_density(graph),
    V.Y = V.Y, V.e = V.e, V.e.net = V.e.net,
    total_edge_n = total_edge_n, nse = nse,
    pos_def = is.positive.definite(Omega1))
}

CompMismatch_SmallGraph <- function(graph, n, X, p, b0, sd.y, betas, delta,
                                    rho = "median.pairwise", include.network = "full"){
  ## Converts a graph into an adjacency matrix
  A <- as.matrix( get.adjacency(graph) )

  ## Convert adjacency matrix into a "starter" precision matrix
  ## Adding transpose to make symmetric? (undirected?)
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

    ## Shrinking graphs
    nd <- new.degs(graph)
    As <- A[nd > 0, nd > 0]
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
    while(is.singular.matrix(diag(nrow(L_tilda_s)) + delta*L_tilda_s) &
          i < 101){ ## only gonna try 50 new plots

      L_tilda_s <- nomatch(As, As_ng, nse) %>% ## Forcing no matching edges
        graph_from_adjacency_matrix(mode = "undirected") %>%
        graph.laplacian(normalized = T) %>% as.matrix

      i <- i+1
    }

    Zs <- Z[ , nd > 0]

    non.singular <- !is.singular.matrix(diag(nrow(L_tilda_s)) + delta*L_tilda_s)
    if(non.singular){

      K_regL <- solve(diag(nrow(L_tilda_s)) + delta*L_tilda_s)
      Zs <- scale(Zs)

      if(rho == "median.pairwise") rho <- median(dist(Zs))

      ## How to include smoothed network
      if(include.network == "full") Zs <- Zs %*% K_regL
      if(include.network == "lower.triangle") {Lc <- chol(K_regL); Zs <- Zs %*% t(Lc)}

      K <- Gaussian_kernel(rho, Zs)
      perm.test <- Permutation_test(K, Y, X, nperm)

      Sc_stat <- perm.test["Stat.sc"]
      Pval <- perm.test["p_value"]
    } else Sc_stat <- Pval <- total_edge_n <- nse <- NA

  }else{
    Sc_stat <- Pval <- total_edge_n <- nse <- NA
  }

  c(Pval = Pval, Sc_stat = Sc_stat, edge_density1 = edge_density(graph),
    non.singular = non.singular,
    total_edge_n = total_edge_n,  pos_def = is.positive.definite(Omega1))
}

## A function built for power analysis on prebuilt data sets
Complete.Mismatch.Sig.Noise <- function(YZ, graph, .X, .delta, .nperm){
  
  A <- as.matrix(get.adjacency(graph))
  gg2 <- sample_pa(ncol(A), directed = F)
  A_p <- as.matrix(get.adjacency(gg2))
  nse <- count_shared_edges(A, A_p)
  
  L_tilda_p <- nomatch(A, A_p, nse) %>% ## Forcing no match
    graph_from_adjacency_matrix(mode = "undirected") %>%
    graph.laplacian(normalized = T) %>% as.matrix
  
  Z <- scale(YZ$Z); rho <- median(dist(Z))
  K_regL <- solve(diag(ncol(L_tilda_p)) + .delta*L_tilda_p); Z <- Z %*% K_regL
  
  K <- Gaussian_kernel(rho, Z)
  perm.test <- Permutation_test(K, YZ$Y, .X, .nperm)
  
  perm.test["p_value"]
}


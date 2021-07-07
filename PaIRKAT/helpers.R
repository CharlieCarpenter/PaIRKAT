###############################
##
## Project: MetaboGuru
##
## Purpose: Shiny app helper functions
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2020-08-11
##
## ---------------------------
## Notes:
##   Get pdat from pathList
##   sapply getNetwork
##   sapply kernelTest
## ---------------------------

# Helpful functions
`%nin%` <- Negate(`%in%`)
pkgs <- c("tidyverse", "magrittr", "igraph", "matrixcalc",
          "MASS", "diffusr", "Matrix", "KEGGREST")

suppressMessages(lapply(pkgs, library, character.only = T))

########## Data Set Up #############

# pathVar = "KEGG",
getNetworks <- function(pathDat, metab, pdat, 
                        pathCol, pathID){
  
  # pathVar <- as.character(pathVar)  pathVar,
  # pdat <- pathList(pathDat, database, min.size)
  
  # pathVar, 

  networks <- lapply(pdat$testPaths$keggPath,
                     getNetwork, .comps = pdat$comps, 
                     .metab = metab, .pathDat = pathDat,
                     compoundReaction = pdat$compoundReaction,
                     .pathCol = pathCol, .pathID = pathID)
  
  ## Naming list
  names(networks) <- pdat$testPaths$pathwayNames[!is.na(pdat$testPaths$pathwayNames)]

  ## Removing networks without connections (getNetwork returns -1)
  keepNet <- !sapply(networks, is.double)
  networks <- networks[keepNet]
  pdat$testPaths <- pdat$testPaths[keepNet, ]
  pdat$comps <- pdat$comps[keepNet]
  
  return(list(networks = networks, pdat = pdat))
}

# pathVar,
pathList <- function(pathDat, min.size, pathID){
  
  .cr <- keggLink("compound", "reaction")
  hsapath <- unique(keggLink("pathway", "hsa"))
  
  cr <- substr(.cr, 5, nchar(.cr))
  reactions <- names(cr)
  
  compId <- as.character(pathDat[, pathID])
  compId <- unlist(strsplit(compId[!is.na(compId)], "[,]"))
  
  results <- data.frame(keggPath = hsapath, 
                        stringsAsFactors=FALSE)
  
  comps <- sapply(hsapath, function(p) keggGet(p)[[1]])
  comp <- sapply(comps, function(p) names(p$COMPOUND))
  compNames <- sapply(comps, function(p) p$NAME)
  results$inpathway <- sapply(comp, function(co) sum(compId %in% co))
  
  testPaths <- results[results$inpathway >= min.size, ]
  co <- sub(" - Homo sapiens (human)", "",
            compNames[names(compNames) %in% testPaths$keggPath],
            fixed = T)
  
  testPaths <- merge(testPaths,
                     data.frame(keggPath = names(co), 
                                pathwayNames = co),
                     by = "keggPath")
  
  return(list(testPaths = testPaths, 
              comps = comps[names(comps) %in% testPaths$keggPath],
              pathDat = pathDat[!is.na(pathDat[, pathID]) & !duplicated(pathDat[, pathID]), ],
              compoundReaction = cr)
  )
}

## Calculates Laplacian of metabolite pathway ## pathVar,
getNetwork <- function(pathId, .comps, .metab, .pathDat, 
                       compoundReaction, .pathCol, .pathID){
  
  target_compound <- names(.comps[[pathId]]$COMPOUND)
  
  cvnames <- .pathDat[.pathDat[[.pathID]] %in% target_compound, ]
  varnames <- cvnames[[.pathCol]]
 
  path.v <- varnames[varnames %in% names(.metab)]
  path.v <- path.v[path.v %in% names(.metab)[names(.metab) %in% varnames ]]
 
  cnames <- cvnames[varnames %in% names(.metab), .pathID, drop = TRUE]
  
  ncomp <- length(cnames)
  reactions <- names(compoundReaction)
  
  A <- matrix(-1, length(cnames), length(cnames))
  diag(A) <- 0
  for(i in 1:ncomp){
    r1 <- reactions[which(compoundReaction==cnames[i])]
    j <- 1
    while (j < i){
      r2 <- reactions[which(compoundReaction==cnames[j])]
      common <- intersect(r1, r2)
      
      if(length(common) > 0) {
        A[i, j] <- A[j, i] <- 1
      } else A[i, j] <- A[j, i] <- 0
      
      j <- j + 1
    }
  }
  
  if(sum(A)==0) return(-1)
  G <- igraph::graph_from_adjacency_matrix(A, mode="undirected")
  V(G)$label <- path.v
  return(G)
}

########## Model Functions #############

## Kernel test including network information through laplacian
PaIRKAT <- function(G, H0.form, data, out.type, tau = 1, metab){
  
  varnames <- V(G)$label
  ZZ <- scale(metab[, varnames[varnames %in% names(metab)]] )
  
  ## normalized Laplacian
  L <- graph.laplacian(G, normalized = T)
  rho <- median(dist(ZZ))
  Z <- ZZ %*% solve(diag(nrow(L)) + tau*L)
  K <- Gaussian_kernel(rho, Z)
  
  if(out.type == "C"){
    pp <- SKAT.c(H0.form, data = data, K=K)
  }
  
  if(out.type == "D"){
    pp <- SKAT.b(H0.form, data = data, K=K)
  }
  
  pp
}

## Function for making formula from "..." in functions
formula_fun <- function(covs){
  cc <- character(0)
  for(i in 2:length(covs)) cc <- paste(cc, covs[i], sep = "+")
  
  formula( paste0("~ ", covs[1], cc) )  ## pasting for final formula
}

Gaussian_kernel <- function(rho, Z){
  exp(-(1/rho)*as.matrix(dist(Z, method = "euclidean", upper = T)^2))
}

## Linear and poly Kern
plyKern <- function(Z, pow, rho=0){
  (Z%*%t(Z) + rho)^pow
}

## Calculates scale param "ka" (kappa) and df "nu"
scaleChi <- function(P, K){ 
  ## Pieces of Itilde
  Itt <- (tr(P%*%K%*%P%*%K))/2
  Its <- tr(P%*%K%*%P)/2
  Iss <- tr(P)/2
  
  Itilde <- Itt - Its %*% solve(Iss) %*% t(Its)
  e <- tr(P%*%K)/2
  
  c(ka = Itilde/(2*e), ## Scale
    nu = (2*e^2)/Itilde) ## Degrees of Freedom
  
  ## Original paper (Liu, Lin, Gosh 2008)
  ## has Iss = tr(P^2)/2 but P is idempotent
}

## Trace of a matrix
tr <- function(x) sum(diag(x))

## Functions for simulation

## Making pos def matrix for sim functions
Danaher_pos_def <- function(m, cc = 3.4, dd = 0.95){
  AA <- m*cc + t(m*cc)
  AA[AA>0] <- 1
  AA <- AA-diag(diag(AA))+diag(nrow(m))*dd
  AA <- AA/( as.matrix(rep(1,nrow(m))) ) %*% t(1.4*rowSums(abs(AA)))
  AA <- (AA+t(AA))/2
  AA <- AA-diag(diag(AA))+diag(nrow(m))*dd
  AA
}

# Pieces for davies -------------------------------------------------

#Compute the tail probability of 1-DF chi-square mixtures
KAT.pval <- function(Q.all, lambda, acc=1e-9,lim=1e6){
  pval = rep(0, length(Q.all))
  i1 = which(is.finite(Q.all))
  for(i in i1){
    tmp <- davies(Q.all[i],lambda,acc=acc,lim=lim)
    pval[i] = tmp$Qq
    
    if(tmp$ifault>0) warning(paste("ifault =", tmp$ifault))
    # pval[i] = Sadd.pval(Q.all[i],lambda)
  }
  return(pval)
}

SKAT.c <- function(formula.H0, data = NULL, K, adjusted = T,
                   acc = 0.00001, lim = 10000, tol = 1e-10) {
  
  m0 <- lm(formula.H0, data)
  mX <- model.matrix(formula.H0, data)
  
  res <- resid(m0); df <- nrow(mX)-ncol(mX)
  s2 <- sum(res^2)/df
  
  P0  <- diag(nrow(mX)) - mX %*% (solve(t(mX) %*% mX) %*% t(mX))
  PKP <- P0 %*% K %*% P0
  
  if(adjusted){
    q <- as.numeric(res %*% K %*% res /(s2*df))
    ee <- eigen(PKP - q * P0, symmetric = T)
    q <- 0 ## Redefining for adjusted stat
  } else{
    q <- as.numeric(res %*% K %*% res / s2)
    ee <- eigen(PKP, symmetric = T) 
  }
  
  lambda <- ee$values[abs(ee$values) >= tol]
  dav <- davies(q, lambda = sort(lambda, decreasing=T),
                acc = acc, lim = lim)
  
  c(dav, Q.adj=q)
}

SKAT.b <- function(formula.H0, data = NULL, K, adjusted = T,
                   acc = 0.00001, lim = 10000, tol = 1e-10) {
  
  X1 <- model.matrix(formula.H0, data)
  lhs <- formula.H0[[2]]
  y <- eval(lhs, data)
  
  y <- factor(y)
  
  
  if (nlevels(y) != 2) {
    stop('The phenotype is not binary!\n')
  } else {
    y <- as.numeric(y) - 1
  }
  
  glmfit <- glm(y ~ X1 - 1, family = binomial)
  
  betas <- glmfit$coef
  mu  <- glmfit$fitted.values
  eta <- glmfit$linear.predictors
  res.wk <- glmfit$residuals
  res <- y - mu
  
  w   <- mu * (1-mu)
  sqrtw <- sqrt(w)
  
  adj <- sum((sqrtw * res.wk)^2) 
  
  DX12 <- sqrtw * X1
  
  qrX <- qr(DX12)
  Q <- qr.Q(qrX)
  Q <- Q[, 1:qrX$rank, drop=FALSE]
  
  P0 <- diag(nrow(X1)) - Q %*% t(Q)
  
  DKD <- tcrossprod(sqrtw) * K
  tQK <- t(Q) %*% DKD
  QtQK <- Q %*% tQK 
  PKP <- DKD - QtQK - t(QtQK) + Q %*% (tQK %*% Q) %*% t(Q)
  q <- as.numeric(res %*% K %*% res) / adj
  ee <- eigen(PKP - q * P0, symmetric = T, only.values=T)  		
  lambda <- ee$values[abs(ee$values) >= tol]
  
  p.value <- KAT.pval(0, lambda=sort(lambda, decreasing = T), acc = acc, lim = lim) 
  
  return(list(p.value=p.value, Q.adj = q))
}



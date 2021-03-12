###############################
##
## Project: PaIRKAT GitHub Files
##
## Purpose: Primary PaIRKAT functions
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2021-03-11
##
## ---------------------------
## Notes:
##   
##
## ---------------------------

source('helpers.R')

# Example -----------------------------------------------------------------

# metabDat <- normalized, imputed, metabolite abundances with
#             metabolites as columns and subjects as rows.
#             Subject IDs should be the first column.
#           
# pathDat <- pathway metadata with a column of unique database
#            (e.g. KEGG) IDs and a column of metabolite names
#            matching colnames of metabDat.
#            
# clinDat <- Clinical variables with demographics, etc. that are
#            important model covariates. Should have subject IDs
#            matchin metabDat's first column.

## Gathering pathway information from KEGG of size >= 10
# pdat <- pathList(pathDat, min.size = 10, 
#                  pathID = "KEGG") 

## Compile pathway information into networks
# nets <- getNetworks(pathDat, metabDat, pdat,
#                     # Name of column with metabolite names
#                     pathCol = "BIOCHEMICAL", 
#                     # Name of column with compound IDs
#                     pathID = "KEGG") 

## Outcome
# Y = clinDat$Y 

## Model matrix
# X <- formula(~ X1 + X2 + X3)
# mm <- model.matrix(X, data = clinDat)

## Results set up
# npath <- nrow(nets$pdat$testPaths)
# pKat.rslt <- data.frame(Pathway = character(npath),
#                         `Pathway Size` = numeric(npath),
#                         `Score Statistic` = numeric(npath),
#                         pValue = numeric(npath))
# 
# ## Running PaIRKAT
# for (i in 1:npath) {
#   z <- PaIRKAT(G = nets$networks[[i]],
#                out.type = "C", # for continuous (or D for dichotomous)
#                Y = Y, model = mm,
#                tau = 1, metab = metabDat)
# 
#   pKat.rslt[i,] <- c(nets$pdat$testPaths$pathwayNames[i],
#                      nets$pdat$testPaths$inpathway[i],
#                      z['Q'], z['pVal'])
# }
# 
# pKat.rslt$pValueFDR <- p.adjust(pKat.rslt$pValue, method = "BH")



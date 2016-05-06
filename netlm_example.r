##################################################################################
### Dependencies
##################################################################################
library(statnet)
set.seed(1)

##################################################################################
### Functions
##################################################################################
createMatrix <- function(file,nlines,nrows,ncols,log1pflag=FALSE){
  matrix <- matrix(0,nrow=nrows, ncol=ncols)
  for(i in 1:nlines){
    c1_id <- file$V1[i]
    c2_id <- file$V2[i]
    matrix[c1_id,c2_id] <- file$V3[i]
    matrix[c2_id,c1_id] <- file$V3[i] 
  }
  if (log1pflag)
    log1p(matrix)
  else
    matrix
}

mrqap <- function(y, x, nperm, mode, diag, nullhyp, tstat, path){
  ### MRQAP
  ptm <- proc.time()
  M <- netlm(y , x , reps = nperm , mode = mode , diag = diag, nullhyp = nullhyp, test.statistic = tstat)
  t <- proc.time() - ptm  
  
  ### PRINT
  nlLabeled <- list()
  nlLabeled <- summary(M)
  nlLabeled$names <- c("Intercept", "Log Trade", "Log Distance", "Colonial")
  nlLabeled$coefficients = round(nlLabeled$coefficients, 4)
  
  ### TO FILE
  capture.output(cat(paste(print(nlLabeled , prmsd = TRUE),print(M$tstat , prmsd = TRUE),print(t, prmsd = TRUE), sep='\n')) , file = paste(path,paste(paste('summary',tstat,mode,nullhyp,nperm, sep='_'),'txt',sep='.'),sep='/') )
  t
}

##################################################################################
### Loadind data files
##################################################################################
file_trade <- read.table('data/country_trade_index.txt',stringsAsFactors=FALSE)
file_relig <- read.table('data/country_religion_index.txt',stringsAsFactors=FALSE)
file_dist <- read.table('data/country_distance_index.txt',stringsAsFactors=FALSE)
file_colon <- read.table('data/country_colonial_index.txt',stringsAsFactors=FALSE)
file_lang <- read.table('data/country_lang_index.txt',stringsAsFactors=FALSE)
country_num <- 249

##################################################################################
### Creating Matrices
##################################################################################
trade_matrix <- createMatrix(file_trade, 11639, country_num, country_num, TRUE)
distance_matrix <- createMatrix(file_dist, 30381, country_num, country_num, TRUE)
colonial_matrix <- createMatrix(file_colon, 134, country_num, country_num)
lang_matrix <- createMatrix(file_lang, 8105, country_num, country_num)
religion_matrix <- createMatrix(file_relig, 10667, country_num, country_num)

##################################################################################
### MODEL
##################################################################################
x <- array(dim=c(3, country_num, country_num))
x[1,,] <- trade_matrix
x[2,,] <- distance_matrix
x[3,,] <- colonial_matrix

##################################################################################
### MRQAP (BETA COEFFICIENTS)
### nullhyp: "qap", "qapspp", "qapy", "qapx", "qapallx"
### test.statistic: "t-value", "beta"
##################################################################################
nperm <- 10
path <- 'results'
M1 <- mrqap(lang_matrix, x, nperm, 'graph', FALSE, 'qap', 'beta', path)
M2 <- mrqap(lang_matrix, x, nperm, 'digraph', FALSE, 'qap', 'beta', path)
M3 <- mrqap(lang_matrix, x, nperm, 'graph', FALSE, 'qapy', 'beta', path)
M4 <- mrqap(lang_matrix, x, nperm, 'digraph', FALSE, 'qapy', 'beta', path)
M5 <- mrqap(lang_matrix, x, nperm, 'graph', FALSE, 'qapx', 'beta', path)
M6 <- mrqap(lang_matrix, x, nperm, 'digraph', FALSE, 'qapx', 'beta', path)
M7 <- mrqap(lang_matrix, x, nperm, 'graph', FALSE, 'qapallx', 'beta', path)
M8 <- mrqap(lang_matrix, x, nperm, 'digraph', FALSE, 'qapallx', 'beta', path)

##################################################################################
### MRQAP (BETA COEFFICIENTS)
### nullhyp: "qap", "qapspp", "qapy", "qapx", "qapallx"
### test.statistic: "t-value", "beta"
##################################################################################
M9 <- mrqap(lang_matrix, x, nperm, 'graph', FALSE, 'qap', 't-value', path)
M10 <- mrqap(lang_matrix, x, nperm, 'digraph', FALSE, 'qap', 't-value', path)
M11 <- mrqap(lang_matrix, x, nperm, 'graph', FALSE, 'qapy', 't-value', path)
M12 <- mrqap(lang_matrix, x, nperm, 'digraph', FALSE, 'qapy', 't-value', path)
M13 <- mrqap(lang_matrix, x, nperm, 'graph', FALSE, 'qapx', 't-value', path)
M14 <- mrqap(lang_matrix, x, nperm, 'digraph', FALSE, 'qapx', 't-value', path)
M15 <- mrqap(lang_matrix, x, nperm, 'graph', FALSE, 'qapallx', 't-value', path)
M16 <- mrqap(lang_matrix, x, nperm, 'digraph', FALSE, 'qapallx', 't-value', path)

##################################################################################
### Printing timings
##################################################################################
capture.output(cat(paste(print(M1,prmsd=TRUE),print(M2,prmsd=TRUE),print(M3,prmsd=TRUE),print(M4,prmsd=TRUE),print(M5,prmsd=TRUE),print(M6,prmsd=TRUE),print(M7,prmsd=TRUE),print(M8,prmsd=TRUE),print(M9,prmsd=TRUE),print(M10,prmsd=TRUE),print(M11,prmsd=TRUE),print(M12,prmsd=TRUE),print(M13,prmsd=TRUE),print(M14,prmsd=TRUE),print(M15,prmsd=TRUE),print(M16,prmsd=TRUE), sep='\n')) , file = 'timing.txt' )

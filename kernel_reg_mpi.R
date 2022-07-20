###
###
### Supervised classification for replicated point patterns with kernel regression
### Pawlasová, Dvořák, 2022
###


### ---------------------------------------------------------
### 1. libraries
suppressMessages(library("spatstat"))
suppressMessages(library("pracma"))
suppressMessages(library(pbdIO, quiet = TRUE))

### ---------------------------------------------------------
### 2. supporting functions

### computing the distance between two estimated pair correlation functions
g_dist <- function(g_X, g_Y, r){
  part = r[2]-r[1]
  j = 2:length(r)
  g_diff = g_X - g_Y
  g_sq = g_diff^2
  g_sq = g_sq[j]*part
  res = sum(g_sq)
  return(res)
}

### computing dissimilarity matrix
dissim_matrix <- function(pcf1, pcf2, r = "NULL"){
  n = dim(pcf1)[1]
  m = dim(pcf2)[1]
  D = matrix(nrow = n, ncol = m)
  
  for(i in 1:n){
    for(j in 1:m){
      D[i, j] = g_dist(pcf1[i,], pcf2[j,], r)
    }
  }
  
  return(D)
}


### function pcf.ppp from spatstat, some parts of the code deleted to make it faster, now it works for point patterns
### on rectangular wondows, no defaults available, translation correction only 

my_pcf = function(X, r){
  win <- Window(X)
  areaW <- area(win)
  npts <- npoints(X)
  lambda <- npts/areaW
  lambda2area <- areaW * lambda^2
  kernel <- "epanechnikov"
  correction <- "translate"
  divisor <- "d"
  stoyan = 0.15
  h <- stoyan/sqrt(lambda)
  hmax <- h
  bw <- h/sqrt(5)
  dr <- r[2L] - r[1L]
  b <- c(-dr, r)
  breaks = as.breakpts(b)
  rmax <- breaks$max
  what <- "all"
  close <- closepairs(X, rmax + hmax, what = what)
  dIJ <- close$d
  
  ### computing translation edge correction
  dx <- close$dx
  dy <- close$dy
  wide <- diff(win$xrange)
  high <- diff(win$yrange)
  edgewt <- wide * high/((wide - abs(dx)) * (high - abs(dy)))

  ### computing the kernel density estimation
  w <- edgewt/dIJ
  wtot <- sum(w)
  weights <- w/wtot
  totMass <- sum(weights)
  n.user <- length(r)
  n <- max(n.user, 512)
  if (n > 512) 
    n <- 2^ceiling(log2(n))
  from = 0
  to = rmax
  lo <- from - 4 * bw
  up <- to + 4 * bw
  y <- .Call(stats:::C_BinDist, dIJ, weights, lo, up, n) * totMass
  kords <- seq.int(0, 2 * (up - lo), length.out = 2L * n)
  kords[(n + 2):(2 * n)] <- -kords[n:2]
  a <- bw * sqrt(5)
  ax <- abs(kords)
  kords = ifelse(ax < a, 3/4 * (1 - (ax/a)^2)/a, 0)
  kords <- fft(fft(y) * Conj(fft(kords)), inverse = TRUE)
  kords <- pmax.int(0, Re(kords)[1L:n]/length(y))
  xords <- seq.int(lo, up, length.out = n)
  x <- seq.int(from, to, length.out = n.user)
  y <- approx(xords, kords, x)$y * wtot
  gT <- y/(2 * pi * lambda2area)
  
  return(gT)
}

### computing pcf for data
compute_pcf <- function(data, r){
  n <- length(data)
  nr <- length(r)
  pcf_data = matrix(nrow = n, ncol = nr)
  
    for (i in 1:n){
      pcf_data[i,] <- my_pcf(data[[i]], r = r)
    }
  return(pcf_data)   
}


### code for the kernel regresion, Ferraty and Vieu, https://www.math.univ-toulouse.fr/~ferraty/online-resources.html
### only slight changes were exectuted

class_knn_lcv <- function(Classes, Dtrain, Dtest, neigh_min=2, step="NULL", kind_of_kernel = "quadratic"){
  
  Classes <- as.vector(Classes)
  if(is.vector(Dtest)) 
    Dd <- as.matrix(t(Dtest))
  testfordim <- sum(dim(Dtrain)==dim(Dtest))==2           

  if(testfordim){
    twodatasets <- sum(Dtrain==Dtest)!=prod(dim(Dtrain))   
  } else
    twodatasets <- TRUE
  
  kernel <- get(kind_of_kernel)
  
  n1 <- ncol(Dtrain)
  if(step=="NULL")
    step <- ceiling(n1/100) 
  if(step == 0)
    step <- 1
  
  Knearest <- seq(from = neigh_min, to = n1 %/% 2, by = step)           
  kmax <- max(Knearest)
  
  Classes_estimated <- 0
  Bandwidth_opt <- 0
  nbclass <- max(Classes)
  BINARY <- matrix(0, n1, nbclass)
  for(g in 1:nbclass)
    BINARY[, g] <- as.numeric(Classes == g)                           
 
  HAT_PROB <- matrix(0, nrow = nbclass, ncol = length(Knearest))  
  
  Knn1 <- 0
  for(i in 1:n1) {
    Norm_diff <- Dtrain[, i] 
    Norm_order <- order(Norm_diff)                              
    zz <- sort(Norm_diff)[2:(kmax + 2)]                         
    Bandwidth <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])             
    z <- zz[ - (kmax + 1)]
    ZMAT <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
    UMAT <- ZMAT/Bandwidth
    KMAT <- kernel(UMAT)
    KMAT[col(KMAT) > row(KMAT)] <- 0                            
    Ind_curves <- Norm_order[2:(kmax + 1)]
    
    for(g in 1:nbclass) {                                               
      Ind_resp <- BINARY[Ind_curves, g]
      YMAT <- matrix(rep(Ind_resp, kmax), nrow = kmax, byrow = T)
      HAT_PROB[g, ] <- apply(YMAT[Knearest,  ] * KMAT[Knearest, ], 1, sum)
    }
    
    Kmatsumbyrow <- apply(KMAT[Knearest,  ], 1, sum)                 
    HAT_PROB <- HAT_PROB/matrix(Kmatsumbyrow, nrow(HAT_PROB), ncol(HAT_PROB), byrow=T)
    Criterium <- t(rep(1, nbclass)) %*% (HAT_PROB - BINARY[i,])^2     
    index <- order(as.vector(Criterium))[1]                           
    
    # results of LCV
    Knn1[i] <- Knearest[index]                                        
    Classes_estimated[i] <- order(HAT_PROB[, index])[nbclass]         
    Bandwidth_opt[i] <- Bandwidth[index]                             
  }
  

  Misclas_estimated <- sum(Classes_estimated != Classes)/n1                    
  

  if(twodatasets) {                                                           
    Bandwidth2 <- 0
    n2 <- nrow(Dtest)
    
    for(k in 1:n2) {
      Sm2k <- Dtest[k,]
      Sm2k_ord <- order(Dtest[k, ])                               
      knn <- Knn1[Sm2k_ord[1]]                                    
      Bandwidth2[k] <- sum(sort(Sm2k)[knn:(knn+1)])*0.5
    }
    
    KERNEL <- kernel(Dtest/Bandwidth2)
    KERNEL[KERNEL < 0] <- 0                                       
    KERNEL[KERNEL > 1] <- 0
    Denom <- apply(as.matrix(KERNEL), 1, sum)
    PROB_PREDICTED <- matrix(0, nrow = n2, ncol = nbclass)
    
    # compute posterior probabilities
    for(g in 1:nbclass) {                                         
      PROBKERNEL <- KERNEL %*% BINARY[, g]
      PROB_PREDICTED[, g] <- PROBKERNEL/Denom
    }
    
    Classes_predicted <- t(apply(PROB_PREDICTED, 1, order))
    Classes_predicted <- Classes_predicted[,nbclass]
    
    return(list(Estimated_classnumber = Classes_estimated, 
                Predicted_classnumber = Classes_predicted, 
                Bandwidths = Bandwidth_opt, 
                Misclas = Misclas_estimated))
  }
  else {
    return(list(Estimated_classnumber = Classes_estimated, 
                Bandwidths = Bandwidth_opt, 
                Misclas = Misclas_estimated))
  }
}


### -------------------------------------------------------------------
### 3. Hyperparameters

Ntrain <- as.numeric(commandArgs(TRUE)[1]) ### Train data - number of patterns in one group (read from keyboard)
Ntest <- as.numeric(commandArgs(TRUE)[2]) ### Test data - number of patterns in one group (read from keyboard) 
NRep <- as.numeric(commandArgs(TRUE)[3]) ### Number of replications (read from keyboard)

Mhalf <- Ntest + 1
Mwhole <- 2*Ntest

r <- seq(0, 0.25, 0.001) ### argument for the pcf, observation window is unit square 

### binary classification
classes <- c(rep(1, Ntrain), rep(2, Ntrain))
incomeclasses <- c(rep(1, Ntest), rep(2, Ntest))

comm.set.seed(seed = 42, diff = FALSE)

### -------------------------------------------------------------------
### 4. Actual computation part - MPI setting for possible replications computed in parallel 

### Epanechnikov kernel used in the classification
quadratic = function(u){
  res = (3/4)*(1 - u^2)
  res[res < 0] = 0
  return(res)
}


### binary classification Thomas vs Thomas (different values of parameter sigma)

doTheJob <- function(i){ 
    m <- rThomas(kappa = 20, scale = 0.1, mu = 6, win = owin(c(0,1),c(0,1)), saveparents=FALSE, nsim = Ntrain)
    mdef <- rThomas(kappa = 20, scale = 0.05, mu = 6, win = owin(c(0,1),c(0,1)), saveparents=FALSE, nsim = Ntrain)
    tdata <- c(m, mdef) ### training data
    
    income <- rThomas(kappa = 20, scale = 0.1, mu = 6, win = owin(c(0, 1),c(0,1)), saveparents=FALSE, nsim = Ntest)
    incomedef <- rThomas(kappa = 20, scale = 0.05, mu = 6, win = owin(c(0,1),c(0,1)), saveparents=FALSE, nsim = Ntest)
    incomedata <- c(income, incomedef) ### testing data 
    
    pcf_train <- compute_pcf(tdata, r=r)
    pcf_income <- compute_pcf(incomedata, r=r)
    
    Dtrain <- dissim_matrix(pcf_train, pcf_train, r = r) ### compute dissimilarity matrix for the training data
    Dtest <- dissim_matrix(pcf_income, pcf_train, r = r) ### compute dissimilarity matrix training vs testing data
  
    classint <- class_knn_lcv(Classes = classes, Dtrain = Dtrain, Dtest = Dtest)$Predicted_classnumber ### predict the labels 
  
    outint11 <- sum(classint[1:Ntest] != incomeclasses[1:Ntest]) ### return the total number of misclassified patterns
    outint12 <- sum(classint[Mhalf:Mwhole] != incomeclasses[Mhalf:Mwhole])
    
    return(c(outint11, outint12))
} 
  
  indices_rep <- comm.chunk(NRep, form = "vector")
  resuls_ind <- do.call(rbind, lapply(indices_rep, doTheJob))
    
  results <- do.call(rbind, allgather(resuls_ind))
  
  finalize()
  








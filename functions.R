#############################################################
#' Required packages
#############################################################
if(!require(Icens)){
  if(!requireNamespace("BiocManager",quietly=TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("Icens")
  require(Icens)
}
library(VGAM) # kendall.tau
library(MLEcens)
library(interval)
#############################################################
#' General notation
#' @n : the number of observations
#' @A : the 0/1 (n x n) matrix of which the i-th row represents
#'      the possible ranks of the i-the observation
#' 
#' Some functions
#' 
#' @is_possible_rank   : Check whether a rank vector is in the
#'                       restricted permutation space
#' @is_possible_rank_c : c version of @is_possible_rank
########################################################
#' input  : @rank = a rank vector with length @n
#'          @A    = a restricted permutation space
#' output : logical value
########################################################
#' 
#' @permute : Permute MCMC samples
########################################################
#' input  : @df_rank = @M X @n matrix of which each row is a rank vector
#'          @perm    = a permutation
#' output : @M X @n permuted matrix
########################################################
#' 
#' @init_rank_c : Initial point from a restricted permutation space
########################################################
#' input  : @A = a restricted permutation space
#' output : a rank vector with length @n
########################################################
#' 
#' @sampler   : Generate MCMC samples on the restricted permutation space
#' @sampler_c : c version of @sampler
########################################################
#' input  : @A       = a restricted permutation space
#'          @B       = the size of MCMC sample
#'          @burnin  = the number of samples to discard as burnin
#'          @thining = the length of thinning interval
#' output : @B X @n matrix
########################################################
#' 
#' @RP   : Compute RP statistic and generate its approximated null distribution
#' @RP_c : c version of @RP
########################################################
#' input  : @Ax      = a restricted permutation space of X
#'          @Ay      = a restricted permutation space of Y
#'          @perm    = logical. If @perm=TRUE, then approximate 
#'          the null distribution using permutation.
#'          @P       = the number of permutations
#'          @B       = the size of MCMC sample
#'          @burnin  = the number of samples to discard as burnin
#'          @thining = the length of thinning interval
#'          @seed    = seed number
#'          @mes     = message is printed if @mes=TRUE
#' output : the RP test statistic
#'          If @perm=TRUE, then return the approximated null distribution.
########################################################
#' 
#' @calc_tau : Calculate Kendall's tau.
########################################################
#' input  : @X = a rank vector of X
#'          @Y = a rank vector of Y
#' output : Kendall's tau value.
########################################################
#' 
#' @calc_tau_vec : Calculate Kendall's tau.
########################################################
#' input  : @X = a @n X @B matrix of which each row is a rank vector of X
#'          @Y = a @n X @B matrix of which each row is a rank vector of Y
#' output : a vector of Kendall's tau values of length @B
########################################################
#'
#' @data.gen : Generate simulation data set.
########################################################
#' input  : @n          = size of data set
#'          @lam        = rate parameter of the exponential distribution
#'          @tau        = Kendall's tau
#'          @family     = copula model
#'          @seed       = seed number
#'          @cl         = left-censoring parameter
#'          @cr         = right-censoring parameter
#'          @lthr       = left threshold
#'          @rthr       = right threshold
#'          @all_cens=F, @rcens=F : interval-censored data
#'          @all_cens=F, @rcens=T : right-censored data
#'          @all_cens=T, @rcens=F : censored data
#'          @all_cens=T, @rcens=T : wrong input
#' output : @X          = @n X 5 matrix with
#'           1st, 2nd col : left and right endpoint of the observed interval
#'           3rd, 4th col : 0 or 1.
#'                          (1,1) : interval-censored or not censored
#'                          (1,0) : right-censored
#'                          (0,1) : left-censored
#'           5th col : L or R or I. censoring type
#'          @Y          = @n X 5 matrix
#'          @cens_ratio = a vector of length 2 indicating 
#'                        the censoring rate of X and Y
#'          @Xsamp      = a vector of length @n with complete X
#'          @Ysamp      = a vector of length @n with complete Y
#'          @s_tau      = Kendall's tau of @Xsamp and @Ysamp
#'          @LR_X       = @n X 2 matrix of censoring times for X
#'          @LR_Y       = @n X 2 matrix of censoring times for Y
########################################################
#'
#' @define_df : return a data.frame, one output of @data.gen
########################################################
#' input  : @X    = a @n X 2 matrix of censored data
#'          @rthr = right threshold
#'          @lthr = left threshold
#' output : @X    = @n X 5 matrix with
#'           1st, 2nd col : left and right endpoint of the observed interval
#'           3rd, 4th col : 0 or 1.
#'                          (1,1) : interval-censored or not censored
#'                          (1,0) : right-censored
#'                          (0,1) : left-censored
#'           5th col : L or R or I. censoring type
########################################################
#'
#' @define_adj_mat : return the restricted permutation space
#'                   given an output of @define_df
########################################################
#' input  : @xdf = an output of @define_df
#' output : the restricted permutation space
########################################################
is_possible_rank <- function(rank,A){
  n <- length(rank)
  s <- sum(sapply(1:n,function(j){
    return(A[j,rank[j]] == 1)
  }))
  return( (s==n) & (sum(rank)==n*(n+1)/2) )
}
is_possible_rank_c <- function(rank, A){
  n <- length(rank)
  dyn.load(paste0('permDCEN',.Platform$dynlib.ext))
  out <- .C('is_possible_rank',n=as.integer(n),rank=as.double(rank),
            A=as.double(A),value=as.integer(1))
  dyn.unload(paste0('permDCEN',.Platform$dynlib.ext))
  return(as.logical(out$value))
}

permute <- function(df_rank,perm){
  n <- ncol(df_rank)
  M <- nrow(df_rank)
  df_rank_perm <- c(df_rank)
  df_rank_perm[df_rank_perm==perm[1]] <- 0
  for(i in 2:n){
    df_rank_perm[df_rank_perm==perm[i]] <- perm[i-1]
  }
  df_rank_perm[df_rank_perm==0] <- perm[n]
  return(matrix(df_rank_perm,nrow=M,ncol=n))
}

init_rank_c <- function(A){
  n <- nrow(A)
  dyn.load(paste0('permDCEN',.Platform$dynlib.ext))
  out <- .C('init_rank',n = as.integer(n),
            A=as.double(A),value=as.double(rep(0,n)),flag=as.integer(0))
  dyn.unload(paste0('permDCEN',.Platform$dynlib.ext))
  return(list(rank=out$value,flag=as.logical(out$flag)))
}

sampler <- function(A,B=1000,burnin=1000,thining=100){
  n <- nrow(A)
  initial_rank <- apply(A,1,which.max)
  for(i in 1:(n-1)){
    A_rowsum <- apply(A[,i:n],1,sum)
    A_rowsum[initial_rank!=i] <- n+1
    initial_rank[initial_rank==i] <- i+1
    initial_rank[which.min(A_rowsum)] <- i
  }
  stopifnot(is_possible_rank(initial_rank,A))
  Rank_sample <- matrix(0,nrow=B,ncol=n)
  rank_update <- initial_rank
  BB <- thining*B + burnin
  for(b in 1:BB){
    rank_temp <- rank_update
    cd <- sample(1:n,2)
    rank_temp[cd] <- rank_temp[rev(cd)]
    if(is_possible_rank(rank_temp,A)){
      rank_update <- rank_temp
    }
    bb <- (b-burnin)%/%thining
    bres <- (b-burnin)%%thining
    if((b>burnin)&(bres==0)){
      Rank_sample[bb,] <- rank_update
    }
  }
  return(Rank_sample)
}
sampler_c <- function(A,B=1000,burnin=1000,thining=100){
  n <- nrow(A)
  dyn.load(paste0('permDCEN',.Platform$dynlib.ext))
  out <- .C('sampler', n = as.integer(n),
            A= as.double(A),B=as.integer(B),
            burnin=as.integer(burnin),thining=as.integer(thining),
            value=as.double(rep(0,n*B)))
  dyn.unload(paste0('permDCEN',.Platform$dynlib.ext))
  return((matrix(out$value,n,B)))
}

RP <- function(Ax,Ay,perm=TRUE,P=5000,B=1000,
               burnin=1000,thining=100,seed=1234,mes=T){
  set.seed(seed)
  n <- nrow(Ax)
  if(mes) cat("Sampling from Sx\n")
  MCMCsampleX <- sampler(Ax,B,burnin,thining)
  if(mes) cat("Sampling from Sy\n")
  MCMCsampleY <- sampler(Ay,B,burnin,thining)
  if(mes) cat("Compute RP\n")
  atop <- mean(sapply(1:B,function(b){
    return(kendall.tau(MCMCsampleX[b,],MCMCsampleY[b,]))
  }))
  if(perm==TRUE)
  {
    if(mes) cat("Compute RP.perm\n")
    atop.perm <- sapply(1:P,function(p){
      perm <- sample(1:n,n)
      MCMCsampleY.perm <- permute(MCMCsampleY,perm)
      ap <- mean(sapply(1:B,function(b){
        return(kendall.tau(MCMCsampleX[b,],MCMCsampleY.perm[b,]))
      }))
      if((p/P*100)%%10==0){
        cat(p/P*100,"% ",sep="")
      }
      return(ap)
    })
    if(mes) cat("\n")
    return(list(atop=atop,atop.perm=atop.perm))
  } else
  {
    return(list(atop=atop))
  }
}
RP_c <- function(Ax,Ay,perm=TRUE,P=1000,B=1000,
                 burnin=1000,thining=100,seed=1234,mes=T){
  set.seed(seed)
  n <- nrow(Ax)
  if(mes) cat("Sampling from Sx\n")
  MCMCsampleX <- sampler_c(Ax,B,burnin,thining)
  if(mes) cat("Sampling from Sy\n")
  MCMCsampleY <- sampler_c(Ay,B,burnin,thining)
  if(mes) cat("Compute RP\n")
  tau_vec <- calc_tau_vec(MCMCsampleX,MCMCsampleY)
  atop_m <- mean(tau_vec)
  if(perm==TRUE){
    if(mes) cat("Compute RP.perm\n")
    atop.perm <- sapply(1:P,function(p){
      perm <- sample(1:n,n)
      MCMCsampleY.perm <- MCMCsampleY[perm,]
      tau_vec_perm <- calc_tau_vec(MCMCsampleX,MCMCsampleY.perm)
      ap_m <- mean(tau_vec_perm)
      if((p/P*100)%%10==0){
        cat(p/P*100,"% ",sep="")
      }
      return(ap_m)
    })
    if(mes) cat("\n")
    return(list(atop=atop_m,atop.perm=atop.perm,tau_vec=tau_vec))
  }else{
    return(list(atop=atop_m))
  }
}

calc_tau <- function(X,Y){
  n <- as.integer(length(X))
  X <- as.double(X)
  Y <- as.double(Y)
  tau <- as.double(0)
  dyn.load(paste0('permDCEN',.Platform$dynlib.ext))
  res <- .C("k_tau",n,X,Y,tau=tau)
  dyn.unload(paste0('permDCEN',.Platform$dynlib.ext))
  tau <- res$tau
  return(tau)
}
calc_tau_vec <- function(X,Y){
  n <- as.integer(nrow(X))
  B <- as.integer(ncol(X))
  X <- as.double(X)
  Y <- as.double(Y)
  tau <- as.double(rep(0,B))
  dyn.load(paste0('permDCEN',.Platform$dynlib.ext))
  res <- .C("calc_tau_vec",n,B,X,Y,tau=tau)
  dyn.unload(paste0('permDCEN',.Platform$dynlib.ext))
  tau <- res$tau
  return(tau)
}

data.gen <- function(n=50,lam=0.1,tau=0,family="Clayton",
                     seed=1234,cl=3.6,cr=6.0,
                     lthr=0,rthr=100,all_cens=FALSE,rcens=FALSE){
  if(!require(copula)){
    install.packages("copula")
    library(copula)
  }
  set.seed(seed)
  if(tau==0){
    X = rexp(n,rate=lam)
    Y = rexp(n,rate=lam)
  }else{
    # copula parameter
    alpha <- iTau(archmCopula(family), tau)
    # define the copula to sample
    cop <- archmCopula(family, param = alpha, dim = 2)
    # sample
    U <- rCopula(n, cop)
    X <- log(1-U[,1])/(-lam)
    Y <- log(1-U[,2])/(-lam)
  }
  s_tau <- cor(X,Y,method='kendall')
  Lx <- runif(n,0,cl)
  Qx <- runif(n,0,cr)
  Rx <- Lx + Qx
  Ly <- runif(n,0,cl)
  Qy <- runif(n,0,cr)
  Ry <- Ly + Qy
  if(!all_cens & !rcens){
    x_cen <- sum(X > Lx & X <= Rx)/n
    y_cen <- sum(Y > Ly & Y <= Ry)/n
    Xl <- X
    Xr <- X
    Yl <- Y
    Yr <- Y
    ind_x <- X>Lx & X < Rx
    ind_y <- Y>Ly & Y < Ry
    Xl[ind_x] <- Lx[ind_x]
    Xr[ind_x] <- Rx[ind_x]
    Yl[ind_y] <- Ly[ind_y]
    Yr[ind_y] <- Ry[ind_y]
  }else if(rcens){
    x_cen <- sum(X > Rx)/n
    y_cen <- sum(Y > Ry)/n
    Xl <- X
    Xr <- X
    Yl <- Y
    Yr <- Y
    ind_x <- X > Rx
    ind_y <- Y > Ry
    Xl[ind_x] <- Rx[ind_x]
    Xr[ind_x] <- rthr
    Yl[ind_y] <- Ry[ind_y]
    Yr[ind_y] <- rthr
  }else if(all_cens){
    x_cen <- 1
    y_cen <- 1
    Xl <- X
    Xr <- X
    Yl <- Y
    Yr <- Y
    ind_x <- X<= Lx
    ind_y <- Y<= Ly
    Xl[ind_x] <- lthr
    Xr[ind_x] <- Lx[ind_x]
    Yl[ind_y] <- lthr
    Yr[ind_y] <- Ly[ind_y]
    ind_x <- X > Rx
    ind_y <- Y > Ry
    Xl[ind_x] <- Rx[ind_x]
    Xr[ind_x] <- rthr
    Yl[ind_y] <- Ry[ind_y]
    Yr[ind_y] <- rthr
    ind_x <- X>Lx & X <= Rx
    ind_y <- Y>Ly & Y <= Ry
    Xl[ind_x] <- Lx[ind_x]
    Xr[ind_x] <- Rx[ind_x]
    Yl[ind_y] <- Ly[ind_y]
    Yr[ind_y] <- Ry[ind_y]
  }
  Xmat <- cbind(Xl,Xr)
  Ymat <- cbind(Yl,Yr)
  xdf <- define_df(Xmat,rthr=rthr,lthr=lthr)
  ydf <- define_df(Ymat,rthr=rthr,lthr=lthr)
  return(list(X=xdf,Y=ydf,cens_ratio=c(x_cen,y_cen),
              Xsamp=X,Ysamp=Y,s_tau=s_tau,
              LR_X=cbind(Lx,Rx),LR_Y=cbind(Ly,Ry)))
}

define_df <- function(X,rthr=100,lthr=-100){
  n <- nrow(X)
  mat <- matrix(0,n,4)
  cens <- rep("I",n)
  cens[X[,1]==lthr] <- "L"
  cens[X[,2]==rthr] <- "R"
  cens[X[,1]==X[,2]] <- "I"
  mat[,1:2] <- X
  mat[cens=="I",3:4] <- 1
  mat[cens=="L",4] <- 1
  mat[cens=="R",3] <- 1
  mat[cens=="I",3:4] <- 1
  x_df <- data.frame(mat=mat,cens=cens)
  return(x_df)
}

define_adj_mat <- function(xdf){
  n <- nrow(xdf)
  Adj_X <- matrix(0,n,n)
  m1 <- sum(xdf[,5]=="L")
  for(i in 1:n){
    if(xdf[i,5]=="L"){
      left_x <- 1
      right_x <- sum(xdf[xdf[,5]=="R",1]<xdf[i,2]) +
        sum(xdf[xdf[,5]=="I" & xdf[,1]!=xdf[,2],1]<xdf[i,2]) +
        sum(xdf[xdf[,5]=="I" & xdf[,1]==xdf[,2],1]<=xdf[i,2]) + m1
      Adj_X[i,left_x:right_x] <- 1
    }else if(xdf[i,5]=="I"){
      if(xdf[i,1]==xdf[i,2]){
        left_x <- sum(xdf[xdf[,5]!="R",2]<xdf[i,1]) + 1 
      }else{
        left_x <- sum(xdf[xdf[,5]!="R",2]<=xdf[i,1]) + 1
      } 
      right_x <- n - sum(xdf[xdf[,5]=="R",1]>=xdf[i,2]) -
        sum(xdf[xdf[,5]=="I" & xdf[,1]!=xdf[,2],1]>=xdf[i,2]) -
        sum(xdf[xdf[,5]=="I" & xdf[,1]==xdf[,2],1]>xdf[i,2])
      Adj_X[i,left_x:right_x] <- 1
    }else if(xdf[i,5]=="R"){
      left_x <- sum(xdf[xdf[,5]!="R",2]<=xdf[i,1]) + 1
      Adj_X[i,left_x:n] <- 1
    }
  }
  return(Adj_X)
}
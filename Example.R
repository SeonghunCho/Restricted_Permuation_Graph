source("functions.R")
#############################################################
#' Toy example
#############################################################
## setting
n <- 50 # sample size
d <- 2
lam <- 0.1
family <- "Clayton" # copula family
tau <- 0.2
cl <- 3.9
cr <- 6.0
## data generation
dat <- data.gen(n=n,d=d,lam=lam,tau=tau,family=family,seed=1234,cl=cl,cr=cr)
xdf <- dat$X
ydf <- dat$Y
Ax <- define_adj_mat(xdf)
Ay <- define_adj_mat(ydf)
rankX <- init_rank_c(Ax)$rank
rankY <- init_rank_c(Ay)$rank

tic <- Sys.time()
for(i in 1:1000)
  calc_tau(rankX,rankY)
toc <- Sys.time()
toc-tic

tic <- Sys.time()
for(i in 1:1000)
  kendall.tau(rankX,rankY)
toc <- Sys.time()
toc-tic

rankXmat <- matrix(rep(rankX,1000),n,1000)
rankYmat <- matrix(rep(rankY,1000),n,1000)
tic <- Sys.time()
tt <- calc_tau_vec(rankXmat,rankYmat)
toc <- Sys.time()
toc-tic

is_possible_rank_c(1:n,Ax)
is_possible_rank_c(rankX,Ax)

mc_sample <- sampler_c(Ax,B=10,burnin=1000,thining=100)

tic <- Sys.time()
res1 <- RP(Ax,Ay,B=100,burnin=1000,thining=100,perm=F)
toc <- Sys.time()
toc-tic

tic <- Sys.time()
res2 <- RP_c(Ax,Ay,B=100,burnin=1000,thining=100,perm=F)
toc <- Sys.time()
toc-tic
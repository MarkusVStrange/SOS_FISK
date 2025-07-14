set.seed(123)
library(RTMB)
library(Matrix)
setwd("/Users/mavast/Documents/Quantitative modelling/module 9/")

pol <- list()
pol$x <- c(rep(1984,11),rep(2022,11),1984:2022,1984:2022)
pol$y <- c(1:11,1:11,rep(1,39),rep(11,39))

# define grid
DY <- .25
ygr <- seq(2, 10, by=DY)
ycen <- ygr[-1]-0.5*DY
DX <- 1
xgr <- seq(min(cohorts[-c(1:3,41)])-0.5*DX, max(cohorts[-c(1:3,41)])+0.5*DX, by=DX)
xcen <- xgr[-1]-0.5*DX

# setup cell table 
CT <- matrix(rep(NA,((length(xgr)-1)*(length(ygr)-1))), ncol=(length(xgr)-1))
grd <- expand.grid(idx=1:ncol(CT), idy=1:nrow(CT))

# number the cells within polygon
idc <- which(sp::point.in.polygon(xcen[grd$idx],ycen[grd$idy],pol$x, pol$y)==1)
grd <- grd[idc,]
CT[cbind(grd$idy,grd$idx)] <- 1:nrow(grd)

plot(pol$x, pol$y, type="n",xlab="cohort",ylab="age")
abline(v=xgr, col="gray")
abline(h=ygr, col="gray")
#xyc<-cbind(xcen[grd$idx],ycen[grd$idy], 1:nrow(grd))
#text(xyc[,1], xyc[,2], labels=xyc[,3], cex=1, col="darkblue", font=2)

# find neighbors 
AT <- array(c(CT,cbind(CT[,-1],NA), cbind(NA,CT[,-ncol(CT)]), rbind(CT[-1,],NA), rbind(NA,CT[-nrow(CT),])),dim=c(nrow(CT),ncol(CT),5))
nextTo <- apply(AT,3,function(x)x)
nextTo <- nextTo[!is.na(nextTo[,1]),]
nextTo <- nextTo[order(nextTo[,1]),]
D <- rowSums(!is.na(nextTo[,-1]), na.rm=TRUE)

# Setup Q0 and I 
Q0 <- spMatrix(length(D), length(D))
diag(Q0) <- D
dummy <- apply(nextTo, 1, function(x){xx<-na.omit(x[-1]); cbind(rep(x[1],length(xx)),xx)})
nn <- do.call(rbind,dummy)
Q0[nn] <- -1

I<-spMatrix(nrow(Q0), nrow(Q0))
diag(I) <- 1 

# Identify which cell each observation is in  
dat <- list(x=rep(cohorts[-c(1:3,41)],each=nrow(mat)),y=rep(ages[-(14:16)],ncol(mat)),sd=as.vector(mat))

xcut <- cut(dat$x,breaks=xgr)
ycut <- cut(dat$y,breaks=ygr)
dat$cellidx <- CT[cbind(as.integer(ycut), as.integer(xcut))]
dat$inside <- !is.na(dat$cellidx)
dat$Q0 <- Q0
dat$I <- I

# parameters, likelihood and all that ...
par<-list()
par$mu <- ycen*0
par$logDelta <- 0  
par$logSigma <- 0  
par$logLambda <- rep(0, nrow(dat$Q0))
par$logSdObs =0

jnll <- function(par){
  getAll(par, dat)
  na_idx <- !is.na(sd)
  sdObs <- exp(logSdObs)
  delta <- exp(logDelta)
  sigma <- exp(logSigma)
  mu <- rep(mu,length(xcen))
  Q <- Q0 + delta*I
  ret <- -dgmrf(logLambda-mu, Q=Q, log=TRUE, scale=sigma)

  ret <- ret-sum(dnorm(sd[na_idx],mean=exp(logLambda[cellidx][na_idx]),sd=sdObs,log=TRUE))
  ret
}

obj<-MakeADFun(jnll, par, random="logLambda")
fit <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)
sdr
pl <- as.list(sdr,"Est")
plsd <- as.list(sdr,"Std")

cc <- viridis::mako(500)
lambda <- exp(pl$logLambda)
fitmat <- t(matrix(lambda[CT], nrow=nrow(CT), ncol=ncol(CT)))
fields::image.plot(xcen, ycen, fitmat, col=cc, xlab="", ylab="")
points(dat$x, dat$y, cex=ifelse(dat$count==0, 1, dat$count/max(dat$count)*5),
       pch=ifelse(dat$count==0,1,16), col=ifelse(dat$count==0,"red","blue"))

# Final model. ADDITIVE effects with refinding efficiency TRUE
#####

par <- list(logitP_cod = rep(0,length(unique(dat$length.class[dat$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dat$length.class[dat$species=="skrubbe"]))),
            logitTag.diff = 0,logitTime.diff = rep(0,7),
            logitRefind.eff = 0,logitloc.diff = rep(0,2))


nll <- function(par){
  getAll(par, dat)
  
  
  found <- OBS(found)
  n <- OBS(n)
  
  c_sizes <- unique(length.class[species=="torsk"])
  f_sizes <- unique(length.class[species=="skrubbe"])
  
  jnll <- -sum(dbinom(feed.exp$found,feed.exp$n,plogis(c(logitRefind.eff,logitRefind.eff+logitTag.diff)),log=TRUE))
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- (0:1)[unique(PIT)==PIT[i]]
    time.idx <- (0:7)[unique(time)==time[i]]
    loc.idx <- (0:2)[unique(release.loca)==release.loca[i]]
    
    
    diff.PIT <- 0
    diff.time <- 0
    diff.loc <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    if(time.idx>0) diff.time <- logitTime.diff[time.idx]
    if(loc.idx>0) diff.loc <- logitloc.diff[loc.idx]
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.time+diff.loc)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.time+diff.loc)
    }
  }
  jnll <- jnll -sum(dbinom(found,n,p,log=TRUE))
  
  codMarch22 <- log(plogis(logitP_cod)/plogis(logitRefind.eff))
  codApril22 <- log(plogis(logitP_cod+logitTime.diff[1])/plogis(logitRefind.eff))
  codMay22 <- log(plogis(logitP_cod+logitTime.diff[2])/plogis(logitRefind.eff))
  codJune22 <- log(plogis(logitP_cod+logitTime.diff[3])/plogis(logitRefind.eff))
  codMarch24 <- log(plogis(logitP_cod+logitTime.diff[4])/plogis(logitRefind.eff))
  codApril24 <- log(plogis(logitP_cod+logitTime.diff[5])/plogis(logitRefind.eff))
  codMay24 <- log(plogis(logitP_cod+logitTime.diff[6])/plogis(logitRefind.eff))
  codJune24 <- log(plogis(logitP_cod+logitTime.diff[7])/plogis(logitRefind.eff))
  
  flounderMarch22 <- log(plogis(logitP_flounder)/plogis(logitRefind.eff))
  flounderApril22 <- log(plogis(logitP_flounder+logitTime.diff[1])/plogis(logitRefind.eff))
  flounderMay22 <- log(plogis(logitP_flounder+logitTime.diff[2])/plogis(logitRefind.eff))
  flounderJune22 <- log(plogis(logitP_flounder+logitTime.diff[3])/plogis(logitRefind.eff))
  flounderMarch24 <- log(plogis(logitP_flounder+logitTime.diff[4])/plogis(logitRefind.eff))
  flounderApril24 <- log(plogis(logitP_flounder+logitTime.diff[5])/plogis(logitRefind.eff))
  flounderMay24 <- log(plogis(logitP_flounder+logitTime.diff[6])/plogis(logitRefind.eff))
  flounderJune24 <- log(plogis(logitP_flounder+logitTime.diff[7])/plogis(logitRefind.eff))
  
  cod22 <- (codMarch22+codApril22+codMay22+codJune22)/4
  cod24 <- (codMarch24+codApril24+codMay24+codJune24)/4
  flounder22 <- (flounderMarch22+flounderApril22+flounderMay22+flounderJune22)/4
  flounder24 <- (flounderMarch24+flounderApril24+flounderMay24+flounderJune24)/4
  
  cod_total <- (codMarch22+codApril22+codMay22+codJune22+codMarch24+codApril24+codMay24+codJune24)/8
  flounder_total <- (flounderMarch22+flounderApril22+flounderMay22+flounderJune22+flounderMarch24+flounderApril24+flounderMay24+flounderJune24)/8
  
  ADREPORT(codMarch22)
  ADREPORT(codApril22)
  ADREPORT(codMay22)
  ADREPORT(codJune22)
  ADREPORT(codMarch24)
  ADREPORT(codApril24)
  ADREPORT(codMay24)
  ADREPORT(codJune24)
  
  ADREPORT(flounderMarch22)
  ADREPORT(flounderApril22)
  ADREPORT(flounderMay22)
  ADREPORT(flounderJune22)
  ADREPORT(flounderMarch24)
  ADREPORT(flounderApril24)
  ADREPORT(flounderMay24)
  ADREPORT(flounderJune24)
  
  ADREPORT(cod22)
  ADREPORT(cod24)
  ADREPORT(flounder22)
  ADREPORT(flounder24)
  
  ADREPORT(cod_total)
  ADREPORT(flounder_total)
  
  ADREPORT(p)
  
  
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

M4 <- opt
# 1. Number of estimated parameters
k <- length(opt$par)
# 2. Log-likelihood is the negative of obj$fn at optimum
logL <- -obj$fn(opt$par)
# 3. Calculate AIC
AIC4 <- 2 * k - 2 * logL
print(AIC4)

res <- oneStepPredict(obj,method="oneStepGeneric",
                      discrete=TRUE)
res$n <- dat$n
res4 <- res
qqnorm(res$residual)
qqline(res$residual)
#####
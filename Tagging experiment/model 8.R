# Random release time - same release-SD pr. year
#####
release_date <- sort(unique(dat$release))
n_release <- length(release_date)
dat$time_df <- data.frame(release_date,
                          year=as.numeric(as.character(str_sub(release_date,1,4))),
                          month=as.numeric(as.character(str_sub(release_date,6,7))),
                          day=as.numeric(as.character(str_sub(release_date,9,10))))

#dat$release.loca[dat$release.loca!="Kalvø"] <- "other"

par <- list(logitP_cod = rep(0,length(unique(dat$length.class[dat$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dat$length.class[dat$species=="skrubbe"]))),
            logitTag.diff = 0,logitYear.diff=0,
            logitRelease.diff = rep(0,n_release),logReleaseSD = 0,
            logitRefind.eff = 0,logitloc.diff = rep(0,1))


nll <- function(par){
  getAll(par, dat)
  
  releaseSD <- exp(logReleaseSD)
  found <- OBS(found)
  n <- OBS(n)
  
  c_sizes <- unique(length.class[species=="torsk"])
  f_sizes <- unique(length.class[species=="skrubbe"])
  release_date <- sort(unique(release))
  n_time <- length(release_date)
  
  jnll <- -sum(dbinom(feed.exp$found,feed.exp$n,plogis(c(logitRefind.eff,logitRefind.eff+logitTag.diff)),log=TRUE))
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- (0:1)[unique(PIT)==PIT[i]]
    year.idx <- (0:1)[unique(year)==year[i]]
    time.idx <- (1:n_time)[release_date==release[i]]
    loc.idx <- (0:1)[unique(release.loca)==release.loca[i]]
    
    
    
    diff.PIT <- 0
    diff.year <- 0
    diff.time <- 0
    diff.loc <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    if(year.idx>0) diff.year <- logitYear.diff[year.idx]
    diff.time <- logitRelease.diff[time.idx]
    if(loc.idx>0) diff.loc <- logitloc.diff[loc.idx]
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.year+diff.time+diff.loc)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.year+diff.time+diff.loc)
    }
  }
  # index for which sd to use for random effect, different by year
  #sd.idx <- ifelse(str_sub(unique(release),1,4)=="2022",1,2)
  jnll <- jnll-sum(dnorm(logitRelease.diff,mean=0,sd=releaseSD,log=TRUE))
  jnll <- jnll -sum(dbinom(found,n,p,log=TRUE))
  
  cod22 <- log(plogis(logitP_cod)/plogis(logitRefind.eff))
  cod24 <- log(plogis(logitP_cod+logitYear.diff)/plogis(logitRefind.eff))
  
  flounder22 <- log(plogis(logitP_flounder)/plogis(logitRefind.eff))
  flounder24 <- log(plogis(logitP_flounder+logitYear.diff)/plogis(logitRefind.eff))
  
  cod_total <- (cod22+cod24)/2
  flounder_total <- (flounder22+flounder24)/2
  
  idx.22 <- which(str_sub(release_date,1,4)=="2022")
  idx.24 <- which(str_sub(release_date,1,4)=="2024")
  cod22_1 <- log(plogis(logitP_cod[1]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  cod22_2 <- log(plogis(logitP_cod[2]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  cod22_3 <- log(plogis(logitP_cod[3]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  cod22_4 <- log(plogis(logitP_cod[4]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  cod22_5 <- log(plogis(logitP_cod[5]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  cod_rel22 <- c(cod22_1+cod22_2+cod22_3+cod22_4+cod22_5)/5
  
  flounder22_1 <- log(plogis(logitP_flounder[1]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  flounder22_2 <- log(plogis(logitP_flounder[2]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  flounder22_3 <- log(plogis(logitP_flounder[3]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  flounder22_4 <- log(plogis(logitP_flounder[4]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  flounder22_5 <- log(plogis(logitP_flounder[5]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  flounder_rel22 <- c(flounder22_1+flounder22_2+flounder22_3+flounder22_4+flounder22_5)/5
  
  cod24_1 <- log(plogis(logitP_cod[1]+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  cod24_2 <- log(plogis(logitP_cod[2]+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  cod24_3 <- log(plogis(logitP_cod[3]+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  cod24_4 <- log(plogis(logitP_cod[4]+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  cod24_5 <- log(plogis(logitP_cod[5]+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  cod_rel24 <- c(cod24_1+cod24_2+cod24_3+cod24_4+cod24_5)/5
  
  flounder24_1 <- log(plogis(logitP_flounder[1]+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  flounder24_2 <- log(plogis(logitP_flounder[2]+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  flounder24_3 <- log(plogis(logitP_flounder[3]+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  flounder24_4 <- log(plogis(logitP_flounder[4]+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  flounder24_5 <- log(plogis(logitP_flounder[5]+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  flounder_rel24 <- c(flounder24_1+flounder24_2+flounder24_3+flounder24_4+flounder24_5)/5
  
  ADREPORT(cod_rel22)
  ADREPORT(cod_rel24)
  ADREPORT(flounder_rel22)
  ADREPORT(flounder_rel24)
  
  ADREPORT(cod22)
  ADREPORT(cod24)
  ADREPORT(flounder22)
  ADREPORT(flounder24)
  ADREPORT(cod_total)
  ADREPORT(flounder_total)
  
  
  release_effect <- exp(logitRelease.diff)
  ADREPORT(release_effect)
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE, random=c("logitRelease.diff"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

M8 <- opt
# 1. Number of estimated parameters
k <- length(opt$par)
# 2. Log-likelihood is the negative of obj$fn at optimum
logL <- -obj$fn(opt$par)
# 3. Calculate AIC
AIC8 <- 2 * k - 2 * logL
print(AIC8)

res <- oneStepPredict(obj,method="oneStepGeneric",
                      discrete=TRUE)
res$n <- dat$n
res8 <- res
par(mfrow=c(1,2))
qqnorm(res$residual)
qqline(res$residual)

rel <- as.list(sdr,"Est",report=TRUE)$release_effect
plot(1:28,rel)
lines(c(12.5,12.5),c(-5,5))
mean(rel[1:12])
mean(rel[13:28])
#####


# numbers for paper
# Three release locations
cat("sandvig = ",c(-0.50041094 -1.96*0.26525343,-0.50041094 +1.96*0.26525343))
cat("SB Hoved = ",c(0.03320152 -1.96*0.19057540,0.03320152 +1.96*0.19057540))

# "Other" release locations
cat("other = ",c(-0.1499922 -1.96*0.16849299,-0.1499922 +1.96*0.16849299))

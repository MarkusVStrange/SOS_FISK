library("readxl")
library(RTMB)
library(dplyr)
library(LaplacesDemon)
library(ggplot2)
library(egg)
library(stringr)
#################
# read and prepare data
#################
#####
tagged24 <- as.data.frame(read_excel("tagging 2024.xlsx", sheet = "Mærkede"))
found24 <- as.data.frame(read_excel("tagging 2024.xlsx", sheet = "Fundne"))

tagged22 <- as.data.frame(read_excel("tagging 2022.xlsx", sheet = "Mærkede"))
found22 <- as.data.frame(read_excel("tagging 2022.xlsx", sheet = "Fundne"))

taggedFeed <- as.data.frame(read_excel("feeding experiment.xlsx", sheet = "Mærkede"))
foundFeed <- as.data.frame(read_excel("feeding experiment.xlsx", sheet = "Fundne"))

df <- data.frame(PIT_nr = c(tagged24$PIT_nummer,tagged22$PIT_nummer,taggedFeed$PIT_nummer),
                 PIT = c(tagged24$PIT,tagged22$PIT,taggedFeed$PIT),
                 release.loca = c(tagged24$Udsætningslokalitet,tagged22$Udsætningslokalitet,taggedFeed$Udsætningslokalitet),
                 release = c(tagged24$Udsætningsdato,tagged22$Udsætningsdato,taggedFeed$Udsætningsdato),
                 species = c(tagged24$Art,tagged22$Art,taggedFeed$Art),
                 length = c(tagged24$Totallængde,tagged22$Totallængde,taggedFeed$Totallængde),
                 year = c(rep(2024,length(tagged24$PIT)),rep(2022,length(tagged22$PIT)),rep('feeding',length(taggedFeed$PIT))))

df$found <- 0
df$colony <- NA

for(i in 1:length(df$PIT)){
  year <- df$year[i]
  PIT <- df$PIT_nr[i]
  
  if(year=='2022'){
    if(PIT %in% found22$PIT_nummer){
      df$found[i] <- 1
      df$colony[i] <- found22$Skannet[found22$PIT_nummer==PIT]
    }
  }
  if(year=='2024'){
    if(PIT %in% found24$PIT_nummer){
      df$found[i] <- 1
      df$colony[i] <- found24$Skannet[found24$PIT_nummer==PIT]
    }
  }
  if(year=='feeding'){
    if(PIT %in% foundFeed$PIT_nummer){
      df$found[i] <- 1
      df$colony[i] <- foundFeed$Skannet[foundFeed$PIT_nummer==PIT]
    }
  }
}
df$species[df$species=="Torsk"] <- "torsk"
df$species[df$species!="torsk"] <- "skrubbe"

tagging.exp <- df %>% filter(year!='feeding')
feeding.exp <- df %>% filter(year=='feeding')

tagging.exp$month <- factor(format(tagging.exp$release, "%B"),levels = c("marts","april","maj","juni"))

# create size classes
tagging.exp$length.class[tagging.exp$species=="torsk"] <- 
  floor((tagging.exp$length[tagging.exp$species=="torsk"]+10)/40)*40-10
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- 
  floor((tagging.exp$length[tagging.exp$species=="skrubbe"]-10)/30)*30+10
table(tagging.exp$length.class[tagging.exp$species=="torsk"])
table(tagging.exp$length.class[tagging.exp$species=="skrubbe"])

tagging.exp$length.class[which(tagging.exp$length.class>220 & tagging.exp$species=="skrubbe")] <- 220
tagging.exp$length.class[which(tagging.exp$length.class>270 & tagging.exp$species=="torsk")] <- 270

tagging.exp$release.loca[tagging.exp$release.loca!="Kalvø"] <- "other"

dexp <- aggregate(found~PIT+release+species+length.class+release.loca+month+year,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+release+species+length.class+release.loca+month+year,data=tagging.exp,FUN = length)$found
dexp$time <- paste(dexp$year,dexp$month)

dFeed <- aggregate(found~PIT,data=feeding.exp,FUN = sum)
dFeed$n <- aggregate(found~PIT,data=feeding.exp,FUN = length)$found

dat <- as.list(dexp)
dat$feed.exp <- dFeed
rm(list=setdiff(ls(),c('dat')))
#####

#################
# Run models
#################
# Model 0: p = p0
#####

par <- list(logitP_cod = rep(0,length(unique(dat$length.class[dat$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dat$length.class[dat$species=="skrubbe"]))),
            logitRefind.eff = 0)


nll <- function(par){
  getAll(par, dat)
  
  c_sizes <- unique(length.class[species=="torsk"])
  f_sizes <- unique(length.class[species=="skrubbe"])
  
  jnll <- -sum(dbinom(feed.exp$found,feed.exp$n,plogis(c(logitRefind.eff,logitRefind.eff)),log=TRUE))
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]])
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]])
    }
  }
  
  cod <- log(plogis(logitP_cod)/plogis(logitRefind.eff))
  
  
  flounder <- log(plogis(logitP_flounder)/plogis(logitRefind.eff))
  
  
  ADREPORT(cod)
  ADREPORT(flounder)
  
  found <- OBS(found)
  n <- OBS(n)
  
  jnll <- jnll -sum(dbinom(found,n,p,log=TRUE))
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr
res <- oneStepPredict(obj,method="oneStepGeneric",
                      discrete=TRUE)
M0 <- opt

# 1. Number of estimated parameters
k <- length(opt$par)
# 2. Log-likelihood is the negative of obj$fn at optimum
logL <- -obj$fn(opt$par)
# 3. Calculate AIC
AIC0 <- 2 * k - 2 * logL
print(AIC0)

res <- oneStepPredict(obj,method="oneStepGeneric",
                      discrete=TRUE)
res$n <- dat$n

res0 <- res
qqnorm(res$residual)
qqline(res$residual)
#####
# Model 1: p = p0 + PIT
#####
par <- list(logitP_cod = rep(0,length(unique(dat$length.class[dat$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dat$length.class[dat$species=="skrubbe"]))),
            logitTag.diff = 0,
            logitRefind.eff = 0)


nll <- function(par){
  getAll(par, dat)
  
  c_sizes <- unique(length.class[species=="torsk"])
  f_sizes <- unique(length.class[species=="skrubbe"])
  
  jnll <- -sum(dbinom(feed.exp$found,feed.exp$n,plogis(c(logitRefind.eff,logitRefind.eff+logitTag.diff)),log=TRUE))
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- (0:1)[unique(PIT)==PIT[i]]
    
    
    diff.PIT <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT)
    }
  }
  
  cod <- log(plogis(logitP_cod)/plogis(logitRefind.eff))
  
  
  flounder <- log(plogis(logitP_flounder)/plogis(logitRefind.eff))
  
  
  ADREPORT(cod)
  ADREPORT(flounder)
  
  found <- OBS(found)
  n <- OBS(n)
  
  jnll <- jnll -sum(dbinom(found,n,p,log=TRUE))
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

M1 <- opt
# 1. Number of estimated parameters
k <- length(opt$par)
# 2. Log-likelihood is the negative of obj$fn at optimum
logL <- -obj$fn(opt$par)
# 3. Calculate AIC
AIC1 <- 2 * k - 2 * logL
print(AIC1)

res <- oneStepPredict(obj,method="oneStepGeneric",
                      discrete=TRUE)
res$n <- dat$n

res1 <- res
qqnorm(res$residual)
qqline(res$residual)

#####
# Model 2: p = p0 + PIT + year
#####
par <- list(logitP_cod = rep(0,length(unique(dat$length.class[dat$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dat$length.class[dat$species=="skrubbe"]))),
            logitTag.diff = 0,logitYear.diff = 0,
            logitRefind.eff = 0)


nll <- function(par){
  getAll(par, dat)
  
  c_sizes <- unique(length.class[species=="torsk"])
  f_sizes <- unique(length.class[species=="skrubbe"])
  
  jnll <- -sum(dbinom(feed.exp$found,feed.exp$n,plogis(c(logitRefind.eff,logitRefind.eff+logitTag.diff)),log=TRUE))
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- (0:1)[unique(PIT)==PIT[i]]
    year.idx <- (0:1)[unique(year)==year[i]]
    
    
    diff.PIT <- 0
    diff.year <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    if(year.idx>0) diff.year <- logitYear.diff[year.idx]
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.year)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.year)
    }
  }
  
  cod22 <- log(plogis(logitP_cod)/plogis(logitRefind.eff))
  cod24 <- log(plogis(logitP_cod+logitYear.diff[1])/plogis(logitRefind.eff))
  
  
  flounder22 <- log(plogis(logitP_flounder)/plogis(logitRefind.eff))
  flounder24 <- log(plogis(logitP_flounder+logitYear.diff[1])/plogis(logitRefind.eff))
  
  
  cod_total <- (cod22+cod24)/2
  flounder_total <- (flounder22+flounder24)/2
  
  ADREPORT(cod22)
  ADREPORT(cod24)
  ADREPORT(flounder22)
  ADREPORT(flounder24)
  
  ADREPORT(cod_total)
  ADREPORT(flounder_total)
  
  found <- OBS(found)
  n <- OBS(n)
  
  jnll <- jnll -sum(dbinom(found,n,p,log=TRUE))
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

M2 <- opt

res <- oneStepPredict(obj,method="oneStepGeneric",
                      discrete=TRUE)
res$n <- dat$n

res2 <- res
qqnorm(res$residual)
qqline(res$residual)

# 1. Number of estimated parameters
k <- length(opt$par)
# 2. Log-likelihood is the negative of obj$fn at optimum
logL <- -obj$fn(opt$par)
# 3. Calculate AIC
AIC2 <- 2 * k - 2 * logL
print(AIC2)
#####
# Model 3: p = p0 + PIT + year_month
#####

par <- list(logitP_cod = rep(0,length(unique(dat$length.class[dat$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dat$length.class[dat$species=="skrubbe"]))),
            logitTag.diff = 0,logitTime.diff = rep(0,7),
            logitRefind.eff = 0)


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
    
    
    diff.PIT <- 0
    diff.time <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    if(time.idx>0) diff.time <- logitTime.diff[time.idx]
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.time)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.time)
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

M3 <- opt
# 1. Number of estimated parameters
k <- length(opt$par)
# 2. Log-likelihood is the negative of obj$fn at optimum
logL <- -obj$fn(opt$par)
# 3. Calculate AIC
AIC3 <- 2 * k - 2 * logL
print(AIC3)
print(paste("delta AIC between model with and without release location is +",round((584.1139/AIC3-1)*100),"%",sep=""))

res <- oneStepPredict(obj,method="oneStepGeneric",
                      discrete=TRUE)
res$n <- dat$n
res3 <- res
qqnorm(res$residual)
qqline(res$residual)
#####
# Model 4: p = p0 + PIT + year_month + release_location
#####

par <- list(logitP_cod = rep(0,length(unique(dat$length.class[dat$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dat$length.class[dat$species=="skrubbe"]))),
            logitTag.diff = 0,logitTime.diff = rep(0,7),
            logitRefind.eff = 0,logitloc.diff = rep(0,1))


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
    loc.idx <- (0:1)[unique(release.loca)==release.loca[i]]
    
    
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
# Model 5: p = p0 + PIT + Random(year_month)
#####

par <- list(logitP_cod = rep(0,length(unique(dat$length.class[dat$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dat$length.class[dat$species=="skrubbe"]))),
            logitTag.diff = 0,logitTime.diff = rep(0,8),logTimeSD = 0,
            logitRefind.eff = 0)


nll <- function(par){
  getAll(par, dat)
  
  timeSD <- exp(logTimeSD)
  found <- OBS(found)
  n <- OBS(n)
  
  c_sizes <- unique(length.class[species=="torsk"])
  f_sizes <- unique(length.class[species=="skrubbe"])
  
  jnll <- -sum(dbinom(feed.exp$found,feed.exp$n,plogis(c(logitRefind.eff,logitRefind.eff+logitTag.diff)),log=TRUE))
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- (0:1)[unique(PIT)==PIT[i]]
    time.idx <- (1:8)[unique(time)==time[i]]
    
    
    diff.PIT <- 0
    diff.time <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    diff.time <- logitTime.diff[time.idx]
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.time)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.time)
    }
  }
  jnll <- jnll-sum(dnorm(logitTime.diff,mean=0,sd=timeSD,log=TRUE))
  jnll <- jnll -sum(dbinom(found,n,p,log=TRUE))
  
  codMarch22 <- log(plogis(logitP_cod+logitTime.diff[1])/plogis(logitRefind.eff))
  codApril22 <- log(plogis(logitP_cod+logitTime.diff[2])/plogis(logitRefind.eff))
  codMay22 <- log(plogis(logitP_cod+logitTime.diff[3])/plogis(logitRefind.eff))
  codJune22 <- log(plogis(logitP_cod+logitTime.diff[4])/plogis(logitRefind.eff))
  codMarch24 <- log(plogis(logitP_cod+logitTime.diff[5])/plogis(logitRefind.eff))
  codApril24 <- log(plogis(logitP_cod+logitTime.diff[6])/plogis(logitRefind.eff))
  codMay24 <- log(plogis(logitP_cod+logitTime.diff[7])/plogis(logitRefind.eff))
  codJune24 <- log(plogis(logitP_cod+logitTime.diff[8])/plogis(logitRefind.eff))
  
  flounderMarch22 <- log(plogis(logitP_flounder+logitTime.diff[1])/plogis(logitRefind.eff))
  flounderApril22 <- log(plogis(logitP_flounder+logitTime.diff[2])/plogis(logitRefind.eff))
  flounderMay22 <- log(plogis(logitP_flounder+logitTime.diff[3])/plogis(logitRefind.eff))
  flounderJune22 <- log(plogis(logitP_flounder+logitTime.diff[4])/plogis(logitRefind.eff))
  flounderMarch24 <- log(plogis(logitP_flounder+logitTime.diff[5])/plogis(logitRefind.eff))
  flounderApril24 <- log(plogis(logitP_flounder+logitTime.diff[6])/plogis(logitRefind.eff))
  flounderMay24 <- log(plogis(logitP_flounder+logitTime.diff[7])/plogis(logitRefind.eff))
  flounderJune24 <- log(plogis(logitP_flounder+logitTime.diff[8])/plogis(logitRefind.eff))
  
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
obj <- MakeADFun(nll, par, silent=TRUE, random=c("logitTime.diff"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

M5 <- opt
# 1. Number of estimated parameters
k <- length(opt$par)
# 2. Log-likelihood is the negative of obj$fn at optimum
logL <- -obj$fn(opt$par)
# 3. Calculate AIC
AIC5 <- 2 * k - 2 * logL
print(AIC5)
print(paste("delta AIC between model with and without release location is +",round((584.1139/AIC5-1)*100),"%",sep=""))

res <- oneStepPredict(obj,method="oneStepGeneric",
                      discrete=TRUE)
res$n <- dat$n
res5 <- res
#qqnorm(res$residual)
#qqline(res$residual)
#####
###### Model 6: p = p0 + PIT + Random(release_day)
#####
release_date <- sort(unique(dat$release))
n_release <- length(release_date)
dat$time_df <- data.frame(release_date,
                          year=as.numeric(as.character(str_sub(release_date,1,4))),
                          month=as.numeric(as.character(str_sub(release_date,6,7))),
                          day=as.numeric(as.character(str_sub(release_date,9,10))))

par <- list(logitP_cod = rep(0,length(unique(dat$length.class[dat$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dat$length.class[dat$species=="skrubbe"]))),
            logitTag.diff = 0,logitTime.diff = rep(0,n_release),logTimeSD = 0,
            logitRefind.eff = 0)


nll <- function(par){
  getAll(par, dat)
  
  timeSD <- exp(logTimeSD)
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
    time.idx <- (1:n_time)[release_date==release[i]]
    
    
    
    diff.PIT <- 0
    diff.time <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    diff.time <- logitTime.diff[time.idx]
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.time)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.time)
    }
  }
  # index for which sd to use for random effect, different by year
  #sd.idx <- ifelse(str_sub(unique(release),1,4)=="2022",1,2)
  jnll <- jnll-sum(dnorm(logitTime.diff,mean=0,sd=timeSD,log=TRUE))
  jnll <- jnll -sum(dbinom(found,n,p,log=TRUE))
  
  release_effect <- exp(logitTime.diff)
  ADREPORT(release_effect)
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE, random=c("logitTime.diff"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

M6 <- opt
# 1. Number of estimated parameters
k <- length(opt$par)
# 2. Log-likelihood is the negative of obj$fn at optimum
logL <- -obj$fn(opt$par)
# 3. Calculate AIC
AIC6 <- 2 * k - 2 * logL
print(AIC6)

res <- oneStepPredict(obj,method="oneStepGeneric",
                      discrete=TRUE)
res$n <- dat$n
res6 <- res
qqnorm(res$residual)
qqline(res$residual)

rel <- as.list(sdr,"Est",report=TRUE)$release_effect
plot(1:28,rel)
lines(c(12.5,12.5),c(-5,5))
#####
# Model 7: p = p0 + PIT + year + Random(release_day)
#####
release_date <- sort(unique(dat$release))
n_release <- length(release_date)
dat$time_df <- data.frame(release_date,
                          year=as.numeric(as.character(str_sub(release_date,1,4))),
                          month=as.numeric(as.character(str_sub(release_date,6,7))),
                          day=as.numeric(as.character(str_sub(release_date,9,10))))

par <- list(logitP_cod = rep(0,length(unique(dat$length.class[dat$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dat$length.class[dat$species=="skrubbe"]))),
            logitTag.diff = 0,logitYear.diff=0,
            logitRelease.diff = rep(0,n_release),logReleaseSD = 0,
            logitRefind.eff = 0)


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
    
    
    
    diff.PIT <- 0
    diff.year <- 0
    diff.time <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    if(year.idx>0) diff.year <- logitYear.diff[year.idx]
    diff.time <- logitRelease.diff[time.idx]
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.year+diff.time)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.year+diff.time)
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

M7 <- opt
# 1. Number of estimated parameters
k <- length(opt$par)
# 2. Log-likelihood is the negative of obj$fn at optimum
logL <- -obj$fn(opt$par)
# 3. Calculate AIC
AIC7 <- 2 * k - 2 * logL
print(AIC7)

res <- oneStepPredict(obj,method="oneStepGeneric",
                      discrete=TRUE)
res$n <- dat$n
res7 <- res
par(mfrow=c(1,2))
qqnorm(res$residual)
qqline(res$residual)

rel <- as.list(sdr,"Est",report=TRUE)$release_effect
plot(1:28,rel)
lines(c(12.5,12.5),c(-5,5))
mean(rel[1:12])
mean(rel[13:28])


#####
# Model 8: p = p0 + PIT + year + Random(release_day) + release_location
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
res8 <- resobj
par(mfrow=c(1,2))
qqnorm(res$residual)
qqline(res$residual)

rel <- as.list(sdr,"Est",report=TRUE)$release_effect
plot(1:28,rel)
lines(c(12.5,12.5),c(-5,5))
mean(rel[1:12])
mean(rel[13:28])
#####

###################
# Compare models
###################
G_M1.M0 <- 2*(M0$obj-M1$obj)
cat('Model 1 better than model 0; p-value =',round(pchisq(G_M1.M0,df=length(M1$par)-length(M0$par),lower.tail = FALSE),5),'\n') 

G_M2.M1 <- 2*(M1$obj-M2$obj)
cat('Model 2 better than model 1; p-value =',round(pchisq(G_M2.M1,df=length(M2$par)-length(M1$par),lower.tail = FALSE),5),'\n') 

G_M3.M2 <- 2*(M2$obj-M3$obj)
cat('Model 3 better than model 2; p-value =',round(pchisq(G_M3.M2,df=length(M3$par)-length(M2$par),lower.tail = FALSE),5),'\n') 

G_M4.M3 <- 2*(M3$obj-M4$obj)
cat('Model 4 not better than model 3; p-value =',round(pchisq(G_M4.M3,df=length(M4$par)-length(M3$par),lower.tail = FALSE),5),'\n') 

# Random effect models
G_M5.M3 <- 2*(M5$obj-M3$obj) # not sure if Model 5 is a sub-model of model 3
cat('Model 3 better than model 5; p-value =',round(pchisq(G_M5.M3,df=length(M3$par)-length(M5$par),lower.tail = FALSE),5),'\n') 

G_M7.M6 <- 2*(M6$obj-M7$obj)
cat('Model 7 not better than model 6; p-value =',round(pchisq(G_M7.M6,df=length(M7$par)-length(M6$par),lower.tail = FALSE),5),'\n') 

G_M8.M7 <- 2*(M7$obj-M8$obj)
cat('Model 8 not better than model 7; p-value =',round(pchisq(G_M8.M7,df=length(M8$par)-length(M7$par),lower.tail = FALSE),5),'\n') 

# AICc
######
l0=M0[["objective"]]
k0=length(M0[["par"]])
AICc0 = 2*k0 + 2*l0 + 2*k0*(k0+1)/(length(dat$PIT)-k0-1)
l1=M1[["objective"]]
k1=length(M1[["par"]])
AICc1 = 2*k1 + 2*l1 + 2*k1*(k1+1)/(length(dat$PIT)-k1-1)
l2=M2[["objective"]]
k2=length(M2[["par"]])
AICc2 = 2*k2 + 2*l2 + 2*k2*(k2+1)/(length(dat$PIT)-k2-1)
l3=M3[["objective"]]
k3=length(M3[["par"]])
AICc3 = 2*k3 + 2*l3 + 2*k3*(k3+1)/(length(dat$PIT)-k3-1)
l4=M4[["objective"]]
k4=length(M4[["par"]])
AICc4 = 2*k4 + 2*l4 + 2*k4*(k4+1)/(length(dat$PIT)-k4-1)
l5=M5[["objective"]]
k5=length(M5[["par"]])
AICc5 = 2*k5 + 2*l5 + 2*k5*(k5+1)/(length(dat$PIT)-k5-1)
l6=M6[["objective"]]
k6=length(M6[["par"]])
AICc6 = 2*k6 + 2*l6 + 2*k6*(k6+1)/(length(dat$PIT)-k6-1)
l7=M7[["objective"]]
k7=length(M7[["par"]])
AICc7 = 2*k7 + 2*l7 + 2*k7*(k7+1)/(length(dat$PIT)-k7-1)
l8=M8[["objective"]]
k8=length(M8[["par"]])
AICc8= 2*k8 + 2*l8 + 2*k8*(k8+1)/(length(dat$PIT)-k8-1)

cat('Model 0 AICc =',round(AICc0,1),'\n') 
cat('Model 1 AICc =',round(AICc1,1),'\n') 
cat('Model 2 AICc =',round(AICc2,1),'\n') 
cat('Model 3 AICc =',round(AICc3,1),'\n') 
cat('Model 4 AICc =',round(AICc4,1),'\n') 
cat('Model 5 AICc =',round(AICc5,1),'\n') 
cat('Model 6 AICc =',round(AICc6,1),'\n') 
cat('Model 7 AICc =',round(AICc7,1),'\n')
cat('Model 8 AICc =',round(AICc8,1),'\n')

#####

# Plot model diagnostics
#####
par(mfrow=c(3,3))
qqnorm(res0$residual,main=paste("model 0: p = ",round(shapiro.test(res0$residual)$p.value,3)))
qqline(res0$residual)
qqnorm(res1$residual,main=paste("model 1: p = ",round(shapiro.test(res1$residual)$p.value,3)))
qqline(res1$residual)
qqnorm(res2$residual,main=paste("model 2: p = ",round(shapiro.test(res2$residual)$p.value,3)))
qqline(res2$residual)
qqnorm(res3$residual,main=paste("model 3: p = ",round(shapiro.test(res3$residual)$p.value,3)))
qqline(res3$residual)
qqnorm(res4$residual,main=paste("model 4: p = ",round(shapiro.test(res4$residual)$p.value,3)))
qqline(res4$residual)
qqnorm(res5$residual,main=paste("model 5: p = ",round(shapiro.test(res5$residual)$p.value,3)))
qqline(res5$residual)
qqnorm(res6$residual,main=paste("model 6: p = ",round(shapiro.test(res6$residual)$p.value,3)))
qqline(res6$residual)
qqnorm(res7$residual,main=paste("model 7: p = ",round(shapiro.test(res7$residual)$p.value,3)))
qqline(res7$residual)
qqnorm(res8$residual,main=paste("model 8: p = ",round(shapiro.test(res8$residual)$p.value,3)))
qqline(res8$residual)
#####
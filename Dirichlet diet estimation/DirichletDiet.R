library(RTMB)
library(tidyverse)
library(LaplacesDemon)

data_wd <- paste(dirname(dirname(getwd())),"/SOS data/",sep="") # data working directory

# Prepare seal data - month
#####
seal.diet <- read.table(paste(data_wd,"gSeal_diet_prop.csv",sep=""),header=TRUE,sep=';')
seal.diet <- seal.diet %>% filter(Site!="Utklippan")
names(seal.diet)[names(seal.diet) %in% c("Year","Month","Species")] <- c("year","month","species")
prey <- c("cod","flatfish")
seal.diet$species[!(seal.diet$species %in% prey)] <- "other"

sdn1 <- aggregate(n~ymp+month+year,data=seal.diet,FUN = mean)
sd <- aggregate(FW.cor~month+year+species,data=seal.diet,FUN = sum)
sd <- left_join(sd,aggregate(n~month+year,data=sdn1,FUN = sum))
sd$ym <- paste(sd$year,sd$month)

full.diet <- as.character((as.data.frame(table(sd$ym)) %>% filter(Freq==3))$Var1)
sd <- sd %>% filter(ym %in% full.diet)

sealFood <- read.table(paste(data_wd,"sealFood_samplings.csv",sep=""),header=TRUE,sep=';')
sealFood <- sealFood %>% filter(species!="herring")
sealFood$species[!(sealFood$species %in% prey)] <- "flatfish"

sealFood$year <- as.numeric(substr(sealFood$sampling,1,4))
sealFood$month <- substr(sealFood$sampling,6,8)

sf <- aggregate(B_index~month+year+species,data=sealFood,FUN = sum)
sf$sampling <- paste(sf$year,sf$month)
sf <- sf %>% filter(sampling %in% sd$ym)
#meanS.B <- aggregate(biomass~species,data=sf,FUN=mean)$biomass
meanS.B <- c(12270972.32,20092.71) # means from the food time series in "prepareDirichletDiet"
S.idx <- match(sf$species,prey)
sf$Bnorm <- (sf$B_index)/meanS.B[S.idx]

rm(list=setdiff(ls(),c('sd','sf','data_wd')))
names(sd)[names(sd)=="FW.cor"] <- "B"
#####
# seal configuration
par <- list(logK=c(0,0,0))
dat <- list(diet=sd,food=sf)

# run diet model
#####
f <- function(par){
  getAll(par,dat)
  k <- exp(logK)
  #r <- exp(logR)
  prey <- unique(diet$species)
  
  ddirichlet <- function(x, alpha, log = TRUE) {
    logB <- sum(lgamma(alpha)) - lgamma(sum(alpha))
    ret <- sum((alpha - 1) * log(x)) - logB
    if(log){
      ret
    } else {
      exp(ret)
    }
  }
  
  samplings <- unique(diet$ym)
  p_cod <- rep(0,length(samplings))
  p_flat <- rep(0,length(samplings))
  p_other <- rep(0,length(samplings))
  ret <- 0
  
  for(i in 1:length(samplings)){
    f.idx <- food$sampling==samplings[i]
    d.idx <- diet$ym==samplings[i]
    n_prey <- length(prey)
    diet.i <- rep(0,n_prey)
    diet.sum <- sum(diet$B[d.idx])
    
    for(j in 1:n_prey){
      diet.i[j] <- diet$B[which(diet$ym==samplings[i] & diet$species==prey[j])]/diet.sum
      
    }
    
    food.i <- log(c(c(food$Bnorm[f.idx])*100,100))
    p <- (food.i^k)/sum(food.i^k)
    ret <- ret -ddirichlet(x=p,alpha=diet.i*diet$n[d.idx][1])
    
    p_cod[i] <- p[1]
    p_flat[i] <- p[2]
    p_other[i] <- p[3]
  }
  #logitPcod <- log(p_cod/(1-p_cod))
  #logitPflat <- log(p_flat/(1-p_flat))
  #logitPother <- log(p_other/(1-p_other))
  ADREPORT(p_cod)
  ADREPORT(p_flat)
  ADREPORT(p_other)
  
  ret
}
obj <- MakeADFun(f,par,silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

pl <- as.list(sdr,"Est",report = TRUE)
plsd <- as.list(sdr,"Std",report=TRUE)


par(mfrow=c(1,3))
for(i in 1:3){
  plot(1:length(pl[[i]]),(pl[[i]]), xlab = "Sampling index" , ylab = "Diet proportion",ylim=c(0,1),type='l',lwd=2)
  lines(1:length(pl[[i]]),(pl[[i]] -2 * plsd[[i]]), lty = "dotted" , lwd =1.5  )
  lines(1:length(pl[[i]]),(pl[[i]] +2 * plsd[[i]]), lty = "dotted" , lwd =1.5 )    
  
}
(pl[[1]])+(pl[[2]])+(pl[[3]]) # test for conservation

#seals
par(mfrow=c(1,2))
x <- log((10:5000)/10)
kx <- exp(0.3775862) #cod
ky <- exp(0.2061222) # flatfish
kz <- exp(-21.8424387) # other
f <- (x^kx)/(x^kx+log(100)^ky+log(100)^kz)
plot(x, f,type = 'l',lwd=2,xlab="log-scale cod biomass index",ylim=c(0,1),
     main="Cod with constant flatfish")
lines(log(c(100,100)),c(-2,2),type = 'l',lty='dashed')
text(3,0.9,"mean cod biomass")

f <- (x^ky)/(x^ky+log(100)^kx+log(100)^kz)
plot(x, f,type = 'l',lwd=2,xlab="log-scale cod biomass index",ylim=c(0,1),
     main="Flatfish with constant cod")
lines(log(c(100,100)),c(-2,2),type = 'l',lty='dashed')
text(3,0.9,"mean flatfish biomass")
#####

sd$month <- factor(sd$month,levels = c("Jan","Mar","Apr","May","Jun","Aug","Oct","Nov"))
plot(sf$Bnorm[sf$species=="cod"],sd$B[sd$species=="cod"]/sd$n[sd$species=="cod"])
plot(sf$Bnorm[sf$species=="flatfish"],sd$B[sd$species=="flatfish"]/sd$n[sd$species=="flatfish"])

plot(scale(sd$B[sd$species=="flatfish"]),
     scale(sf$Bnorm[sf$species=="flatfish"]),xlab="Diet",ylab = "Food")
lines(c(0,0),c(-5,5),lty="dashed")
lines(c(-5,5),c(0,0),lty="dashed")

plot(scale(sd$B[sd$species=="cod"]),
     scale(sf$Bnorm[sf$species=="cod"]),xlab="Diet",ylab = "Food")
lines(c(0,0),c(-5,5),lty="dashed")
lines(c(-5,5),c(0,0),lty="dashed")

cor(scale(sd$B[sd$species=="flatfish"]),
    scale(sf$Bnorm[sf$species=="flatfish"]))

cor(scale(sd$B[sd$species=="cod"]),
    scale(sf$Bnorm[sf$species=="cod"]))



# Prepare cormorant data - month
#####
cormFood <- read.table(paste(data_wd,"cormorantFood_samplings.csv",sep=""),header=TRUE,sep=';')
cormFood <- cormFood %>% filter(species!="herring")
prey <- c("cod","flatfish")
cormFood$species[!(cormFood$species %in% prey)] <- "flatfish"
cormFood <- aggregate(B_index~sampling+species,data=cormFood,FUN = sum)
corm.diet <-  read.table(paste(data_wd,"corm_diet_prop.csv",sep=""),header=TRUE,sep=';')
corm.diet$species[!(corm.diet$species %in% prey)] <- "other"
corm.diet$month <- factor(corm.diet$month,levels = month <- c("Feb","Mar","Apr","May","Jun","Jul","Aug",
                                                              "Sep","Oct") )

cdn1 <- aggregate(n~dmyl+month+year,data=corm.diet,FUN = mean)
cd <- aggregate(B~species+month+year,data=corm.diet,FUN = sum)
cd <- left_join(cd,aggregate(n~month+year,data=cdn1,FUN = sum))
cd$ym <- paste(cd$year,cd$month)


cormFood$year <- cd$year[match(cormFood$sampling,cd$ym)]
cormFood$month <- cd$month[match(cormFood$sampling,cd$ym)]
cormFood$month <- factor(cormFood$month,levels = month <- c("Feb","Mar","Apr","May","Jun","Jul","Aug",
                                                            "Sep","Oct") )
cormFood$samp_count <- 1

cf <- aggregate(B_index~sampling+year+month+species,data=cormFood,FUN = sum)
cf$n_samp <- aggregate(samp_count~sampling+year+month+species,data=cormFood,FUN = sum)$samp_count

meanC.B <- aggregate(B_index~species,data=cf,FUN=mean)$B_index  
meanC.B <- c(12592402.04,49774.32)

C.idx <- match(cf$species,prey)
cf$Bnorm <- (cf$B_index)/meanC.B[C.idx]

rm(list=setdiff(ls(),c('cd','cf','data_wd')))
#####
# cormorant configuration
par <- list(logK=c(0,0,0))
dat <- list(diet=cd,food=cf)

# run diet model
#####
f <- function(par){
  getAll(par,dat)
  k <- exp(logK)
  #r <- exp(logR)
  prey <- unique(diet$species)
  
  ddirichlet <- function(x, alpha, log = TRUE) {
    logB <- sum(lgamma(alpha)) - lgamma(sum(alpha))
    ret <- sum((alpha - 1) * log(x)) - logB
    if(log){
      ret
    } else {
      exp(ret)
    }
  }
  
  samplings <- unique(diet$ym)
  p_cod <- rep(0,length(samplings))
  p_flat <- rep(0,length(samplings))
  p_other <- rep(0,length(samplings))
  ret <- 0
  
  cod.d <- rep(0,length(samplings))
  cod.f <- rep(0,length(samplings))
  for(i in 1:length(samplings)){
    f.idx <- food$sampling==samplings[i]
    d.idx <- diet$ym==samplings[i]
    n_prey <- length(prey)
    diet.i <- rep(0,n_prey)
    diet.sum <- sum(diet$B[d.idx])
    
    for(j in 1:n_prey){
      diet.i[j] <- diet$B[which(diet$ym==samplings[i] & diet$species==prey[j])]/diet.sum
      
    }
    
    food.i <- log(c(c(food$Bnorm[f.idx])*100,100))
    p <- (food.i^k)/sum(food.i^k)
    ret <- ret -ddirichlet(x=p,alpha=diet.i*diet$n[d.idx][1])
    
    p_cod[i] <- p[1]
    p_flat[i] <- p[2]
    p_other[i] <- p[3]
    
    cod.d[i] <- diet.i[1]
    cod.f[i] <- food.i[1]/sum(food.i)
    print(food.i)
  }
  #logitPcod <- log(p_cod/(1-p_cod))
  #logitPflat <- log(p_flat/(1-p_flat))
  #logitPother <- log(p_other/(1-p_other))
  ADREPORT(p_cod)
  ADREPORT(p_flat)
  ADREPORT(p_other)
  
  ret
}
obj <- MakeADFun(f,par,silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

pl <- as.list(sdr,"Est",report = TRUE)
plsd <- as.list(sdr,"Std",report=TRUE)


par(mfrow=c(1,3))
for(i in 1:3){
  plot(1:length(pl[[i]]),(pl[[i]]), xlab = "Sampling index" , ylab = "Diet proportion",ylim=c(0,1),type='l',lwd=2)
  lines(1:length(pl[[i]]),(pl[[i]] -2 * plsd[[i]]), lty = "dotted" , lwd =1.5  )
  lines(1:length(pl[[i]]),(pl[[i]] +2 * plsd[[i]]), lty = "dotted" , lwd =1.5 )    
  
}
(pl[[1]])+(pl[[2]])+(pl[[3]]) # test for conservation

#cormorants
par(mfrow=c(1,2))
x <- log((10:5000)/10)
kx <- exp(1.144524) #cod
ky <- exp(1.540565) # flatfish
kz <- exp(1.313913) # other
f <- (x^kx)/(x^kx+log(100)^ky+log(100)^kz)
plot(x, f,type = 'l',lwd=2,xlab="log-scale cod biomass index",ylim=c(0,1),
     main="Cod with constant flatfish")
lines(log(c(100,100)),c(-2,2),type = 'l',lty='dashed')
text(3,0.9,"mean cod biomass")

f <- (x^ky)/(x^ky+log(100)^kx+log(100)^kz)
plot(x, f,type = 'l',lwd=2,xlab="log-scale cod biomass index",ylim=c(0,1),
     main="Flatfish with constant cod")
lines(log(c(100,100)),c(-2,2),type = 'l',lty='dashed')
text(3,0.9,"mean flatfish biomass")
#####


plot(scale(log(cd$B[cd$species=="flatfish"])),
     scale(log(cf$Bnorm[cf$species=="flatfish"])),xlab="Diet",ylab = "Food")
lines(c(0,0),c(-5,5),lty="dashed")
lines(c(-5,5),c(0,0),lty="dashed")

plot(scale(log(cd$B[cd$species=="cod"])),
     scale(log(cf$Bnorm[cf$species=="cod"])),xlab="Diet",ylab = "Food")
lines(c(0,0),c(-5,5),lty="dashed")
lines(c(-5,5),c(0,0),lty="dashed")

cor(scale(log(cd$B[cd$species=="flatfish"])),
    scale(log(cf$Bnorm[cf$species=="flatfish"])))

cor(scale(log(cd$B[cd$species=="cod"])),
    scale(log(cf$Bnorm[cf$species=="cod"])))


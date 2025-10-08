####################
library(LaplacesDemon)
library(surveyIndex)
library(RTMB)

# working directory to coefficient-file
wd_coef <- "C:/Users/mavast/Documents/GitHub/SOS data/"

dAll <- readRDS(paste(wd_coef,"DatrasExchange.R",sep=""))
WB <- st_read(paste(wd_coef,"shapefiles/ICES_areas.shp",sep=""))
plot(WB %>% filter(ICES_SUB %in% 22:24))
WB <- as(WB,"Spatial")

dAll <- addSpatialData(dAll,WB)

# Cod ALK QI
#####
d<-subset(dAll, Species=="Gadus morhua",ICES_SUB %in% 22:24)
dAll<-NULL; gc(); ## lose dAll because it takes up a lot of memory
d<-addSpectrum(d,by=1)

ca <- d[[1]]
hh <- d[[2]]
hl <- d[[3]]

hl <- aggregate(HLNoAtLngt~haul.id+LngtCm,data=hl,FUN=sum)

ca$lengthID <- paste(ca$haul.id,ca$LngtCm)
hl$lengthID <- paste(hl$haul.id,hl$LngtCm)

ca_cod <- left_join(ca,hl %>% select(lengthID,HLNoAtLngt))
ca_cod$cohort <- as.numeric(as.character(ca_cod$Year))-ca_cod$Age
#ca_cod <- ca_cod %>% filter(Age>0 & !is.na(HLNoAtLngt))
ca_cod <- ca_cod %>% filter(!is.na(HLNoAtLngt))
hh$Jday <- yday(make_date(year=hh$Year,month=hh$Month,day=hh$Day))/yday(make_date(year=hh$Year,month=12,day=31))
ca_cod <- left_join(ca_cod,hh %>% select(haul.id,Jday))
ca_cod <- ca_cod %>% filter(!(is.na(Jday) | is.na(Age)))
ca_cod <- ca_cod %>% filter(!(Age<1 & Quarter==1))
hist(ca_cod$LngtCm[ca_cod$Age==0])
ca_cod <- ca_cod %>% filter(!(Age<1 & Country=="DE" & Year %in% c(2015,2018)))
hist(ca_cod$LngtCm[ca_cod$Age==0])

cohorts <- sort(unique(ca_cod$cohort))[-c(1:6,(length(unique(ca_cod$cohort))-1):length(unique(ca_cod$cohort)))] # cohort with sufficient observations to fit
ca_cod <- ca_cod %>% filter(cohort %in% cohorts)
ca_cod$logLength <- log(ca_cod$LngtCm)
ca_cod$aq <- paste(ca_cod$Age,ca_cod$Quarter)


par <- list(logL1_1=log(20),logL1Rand=rep(log(20),length(cohorts)-1),logLSD=log(10),
            logL2_1=log(40),logL2Rand=rep(log(40),length(cohorts)-1),
            logL3_1=log(50),logL3Rand=rep(log(50),length(cohorts)-1),
            logSD=log(10))

dat <- list(ca=ca_cod %>% select(haul.id,Age,cohort,logLength,HLNoAtLngt,Jday,Quarter)
            )#%>% filter(Age>0))
#rm(list=setdiff(ls(),c('par','dat')))
# run diet model
#####
f <- function(par){
  getAll(par,dat)
  L1 <- exp(c(logL1_1,logL1Rand))
  L2 <- exp(c(logL2_1,logL2Rand))
  L3 <- exp(c(logL3_1,logL3Rand))
  #L2<- L1+exp(c(logLdiff2_1,logdiff2Rand))
  #L3<- L2+exp(c(logLdiff2_1,logdiff2Rand))-exp(c(neglogLdiff3_1,neglogLdiff3Rand))
  r <- (L3-L2)/(L2-L1)
  t1 <- 1
  t3 <- 6
  t2 <- (t1+t3)/2
  LSD <-exp(logLSD)
  #L2SD <-exp(logL2SD)
  #L3SD <-exp(logL3SD)
  
  SD <- exp(logSD)
  
  cohorts <- sort(unique(ca$cohort))
  ret <- 0
  for(i in 1:length(cohorts)){
    cohort.idx <- which(ca$cohort==cohorts[i])
    ages.i <- sort(unique(ca$Age[cohort.idx]+ca$Jday[cohort.idx]))
    if(i>1){
      ret <- ret-dnorm(log(L1[i]),mean=log(L1[i-1]),sd=LSD,log=TRUE)
      ret <- ret-dnorm(log(L2[i]),mean=log(L2[i-1]),sd=LSD,log=TRUE)
      ret <- ret-dnorm(log(L3[i]),mean=log(L3[i-1]),sd=LSD,log=TRUE)
      
    }
    #plot(ca$Age[cohort.idx]+ca$Jday[cohort.idx],ca$logLength[cohort.idx],main=cohorts[i],col="white")
    for(j in 1:length(ages.i)){
      age.idx <- which((ca$Age[cohort.idx]+ca$Jday[cohort.idx])==ages.i[j])
      t <- ages.i[j]
      L.i <- L1[i]+(L3[i]-L1[i])*(1-r[i]^(2*(t-t1)/(t3-t1)))*(1-r[i]^2)^(-1)
      #vbgr.sd <- calculate_var()
      #sd.fit <- vbgr(ages.i[j],k[i],Linf[i],t_hatch)
      if(sum(ca$HLNoAtLngt[cohort.idx][age.idx]>1)){
        ret <- ret-sum(dnorm(ca$logLength[cohort.idx][age.idx],mean=log(L.i),sd=SD,log=TRUE)*ca$HLNoAtLngt[cohort.idx][age.idx]/sum(log(ca$HLNoAtLngt[cohort.idx][age.idx])))
      }#  points(ages.i[j],log(vbgr.fit),cex=3,col="darkred",pch=19)
      # points(ca$Age[cohort.idx][age.idx]+ca$Jday[cohort.idx][age.idx],ca$logLength[cohort.idx][age.idx],cex=1,col="black",pch=1)
    }
    #lines((0:100)/10,log(vbgr((0:100)/10,0.1603331,112.2677,0.2917808)),lwd=3,col="darkred")
    #readline() # press enter to continue
  }
  x <- (0:1000)/100
  c <- which(cohorts==1997)
  logfit <- log(L1[c]+(L3[c]-L1[c])*(1-r[c]^(2*(x-t1)/(t3-t1)))*(1-r[c]^2)^(-1))
  ADREPORT(logfit)
  ret
}
#obj <- MakeADFun(f,par,random=c('logL1Rand','logdiff2Rand','neglogLdiff3Rand'))
obj <- MakeADFun(f,par,random=c('logL1Rand','logL2Rand','logL3Rand'))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

SDest <- exp(as.list(sdr,"Est")$logSD)
pl <- as.list(sdr,"Est",report=TRUE)$logfit
plsd <- as.list(sdr,"Std",report=TRUE)$logfit

cohortX <- ca_cod %>% filter(cohort==1993)

plot((0:1000)/100,exp(pl) , lwd =2 ,type='l',ylim=c(0,100),col="black")
lines((0:1000)/100,exp(pl -sqrt(SDest^2)), lty = "dotted" , lwd =1.5,col="blue")
lines((0:1000)/100,exp(pl +sqrt(SDest^2)), lty = "dotted" , lwd =1.5,col="blue")  
#lines((0:1000)/100,exp(pl - SD), lty = "dotted" , lwd =1.5,col="red")
#lines((0:1000)/100,exp(pl +SD), lty = "dotted" , lwd =1.5,col="red")  
#lines((0:1000)/100,exp(pl -sqrt(SD^2+SDest^2)), lty = "dotted" , lwd =1.5,col="orange")
#lines((0:1000)/100,exp(pl +sqrt(SD^2+SDest^2)), lty = "dotted" , lwd =1.5,col="orange")  

#lines((0:1000)/100,exp(pl -2 * plsd), lty = "dotted" , lwd =1.5  ,col="red")
#lines((0:1000)/100,exp(pl +2 * plsd), lty = "dotted" , lwd =1.5  ,col="red") 
#lines((0:1000)/100,exp(pl -2 * exp(-1.71186124)), lty = "dotted" , lwd =1.5 ,col="orange" )
#lines((0:1000)/100,exp(pl +2 * exp(-1.71186124)), lty = "dotted" , lwd =1.5 ,col="orange") 
points(cohortX$Age+cohortX$Jday,cohortX$LngtCm)

hist(rnorm(699,mean=exp(pl[81]),sd=2*exp(sqrt(plsd^2+SDest^2)[81])))

colors <- colorRampPalette(c("red", "blue"))(length(cohorts))
params <- as.list(sdr,"Est")

for(i in 1:length(cohorts)){
  L1 <- exp(c(params$logL1_1,params$logL1Rand))[i]
  L2 <- exp(c(params$logL2_1,params$logL2Rand))[i]
  L3 <- exp(c(params$logL3_1,params$logL3Rand))[i]
  r <- (L3-L2)/(L2-L1)
  t1 <- 1
  t3 <- 6
  t2 <- (t1+t3)/2
  t <- (0:1000)/100
  L <- L1+(L3-L1)*(1-r^(2*(t-t1)/(t3-t1)))*(1-r^2)^(-1)
  if(i==1){
    plot(t,L, lwd =1.5 ,type='l',ylim=c(0,120),col=colors[i])
  }
  if(i>1){
    lines(t,L, lwd =1.5 ,col=colors[i])
  }
  #cohortX <- ca_cod %>% filter(cohort==cohorts[i])
  #points(cohortX$Age+cohortX$Jday,cohortX$LngtCm)
  #readline() # press enter to continue
}

for(i in 1:length(cohorts)){
  L1 <- exp(c(params$logL1_1,params$logL1Rand))[i]
  L2 <- exp(c(params$logL2_1,params$logL2Rand))[i]
  L3 <- exp(c(params$logL3_1,params$logL3Rand))[i]
  r <- (L3-L2)/(L2-L1)
  t1 <- 1
  t3 <- 6
  t2 <- (t1+t3)/2
  t <- (0:1000)/100
  L <- L1+(L3-L1)*(1-r^(2*(t-t1)/(t3-t1)))*(1-r^2)^(-1)
  if(i==1){
    plot(t,L, lwd =1.5 ,type='l',ylim=c(0,20),col=colors[i],xlim=c(0,1))
  }
  if(i>1){
    lines(t,L, lwd =1.5 ,col=colors[i])
  }
}



################
# Estimated hatching time
################
par <- list(t1=0,logLSD2=log(10),logLSD3=log(10),
            logL2_1=log(40),logL2Rand=rep(log(40),length(cohorts)-1),
            logL3_1=log(50),logL3Rand=rep(log(50),length(cohorts)-1),
            logSD=log(10),logSDcons=0)

dat <- list(ca=ca_cod %>% select(haul.id,Age,cohort,logLength,HLNoAtLngt,Jday,Quarter)
)#%>% filter(Age>0))
#rm(list=setdiff(ls(),c('par','dat')))
# run diet model
#####
f <- function(par){
  getAll(par,dat)
  L1 <- 0#exp(c(logL1_1,logL1Rand))
  L2 <- exp(c(logL2_1,logL2Rand))
  L3 <- exp(c(logL3_1,logL3Rand))
  #L2<- L1+exp(c(logLdiff2_1,logdiff2Rand))
  #L3<- L2+exp(c(logLdiff2_1,logdiff2Rand))-exp(c(neglogLdiff3_1,neglogLdiff3Rand))
  r <- (L3-L2)/(L2-L1)
  #t1 <- 1
  t1<- plogis(t1)*0.4958904+0.04109589
  t3 <- 6
  t2 <- (t1+t3)/2
  LSD2 <-exp(logLSD2)
  LSD3 <-exp(logLSD3)
  #t1SD <- exp(logSDt1)
  #L2SD <-exp(logL2SD)
  #L3SD <-exp(logL3SD)
  
  SD <- exp(logSD)
  SDcons <- exp(logSDcons)
  
  cohorts <- sort(unique(ca$cohort))
  ret <- 0
  for(i in 1:length(cohorts)){
    cohort.idx <- which(ca$cohort==cohorts[i])
    ages.i <- sort(unique(ca$Age[cohort.idx]+ca$Jday[cohort.idx]))
    if(i>1){
      #ret <- ret-dnorm(t1[i],mean=t1[i-1],sd=t1SD,log=TRUE)
      ret <- ret-dnorm(log(L2[i]),mean=log(L2[i-1]),sd=LSD2,log=TRUE)
      ret <- ret-dnorm(log(L3[i]),mean=log(L3[i-1]),sd=LSD3,log=TRUE)
      
    }
    #plot(ca$Age[cohort.idx]+ca$Jday[cohort.idx],ca$logLength[cohort.idx],main=cohorts[i],col="white")
    for(j in 1:length(ages.i)){
      age.idx <- which((ca$Age[cohort.idx]+ca$Jday[cohort.idx])==ages.i[j])
      t <- ages.i[j]
      L.i <- L1+(L3[i]-L1)*(1-r[i]^(2*(t-t1)/(t3-t1)))*(1-r[i]^2)^(-1)
      #vbgr.sd <- calculate_var()
      #sd.fit <- vbgr(ages.i[j],k[i],Linf[i],t_hatch)
      if(sum(ca$HLNoAtLngt[cohort.idx][age.idx]>1)){
        ret <- ret-sum(dnorm(ca$logLength[cohort.idx][age.idx],
                             mean=log(L.i),sd=SD+(-log(L.i) + log(exp(log(L.i)) + SDcons)),
                             log=TRUE)*ca$HLNoAtLngt[cohort.idx][age.idx]/
                         sum(log(ca$HLNoAtLngt[cohort.idx][age.idx])))
      }#  points(ages.i[j],log(vbgr.fit),cex=3,col="darkred",pch=19)
      # points(ca$Age[cohort.idx][age.idx]+ca$Jday[cohort.idx][age.idx],ca$logLength[cohort.idx][age.idx],cex=1,col="black",pch=1)
    }
    #lines((0:100)/10,log(vbgr((0:100)/10,0.1603331,112.2677,0.2917808)),lwd=3,col="darkred")
    #readline() # press enter to continue
  }
  x <- (0:1000)/100
  c <- which(cohorts==2016)
  logfit <- log(L1+(L3[c]-L1)*(1-r[c]^(2*(x-t1)/(t3-t1)))*(1-r[c]^2)^(-1))
  ADREPORT(logfit)
  ret
}
#obj <- MakeADFun(f,par,random=c('logL1Rand','logdiff2Rand','neglogLdiff3Rand'))
obj <- MakeADFun(f,par,random=c('logL2Rand','logL3Rand'),
                map = list(t1 = factor(NA)))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

SDest <- exp(as.list(sdr,"Est")$logSD)
SDcons <- exp(as.list(sdr,"Est")$logSDcons)
pl <- as.list(sdr,"Est",report=TRUE)$logfit
plsd <- as.list(sdr,"Std",report=TRUE)$logfit


plot(exp(pl), type = 'l')
lines(exp(pl - plsd*2))
lines(exp(pl + plsd*2))
cohortX <- ca_cod %>% filter(cohort==2016)

SD <- -pl + log(exp(pl) + SDcons)


plot((0:1000)/100,exp(pl) , lwd =2 ,type='l',ylim=c(0,100),col="black")
lines((0:1000)/100,exp(pl -sqrt(SDest^2)), lty = "dotted" , lwd =1.5,col="blue")
lines((0:1000)/100,exp(pl +sqrt(SDest^2)), lty = "dotted" , lwd =1.5,col="blue")  
lines((0:1000)/100,exp(pl - SD), lty = "dotted" , lwd =1.5,col="red")
lines((0:1000)/100,exp(pl +SD), lty = "dotted" , lwd =1.5,col="red")  
lines((0:1000)/100,exp(pl -sqrt(SD^2+SDest^2)), lty = "dotted" , lwd =1.5,col="orange")
lines((0:1000)/100,exp(pl +sqrt(SD^2+SDest^2)), lty = "dotted" , lwd =1.5,col="orange")  
points(cohortX$Age+cohortX$Jday,cohortX$LngtCm)

hist(rnorm(699,mean=exp(pl[100]),sd=2*exp(sqrt(SD^2+SDest^2)[100])))

colors <- colorRampPalette(c("red", "blue"))(length(cohorts))
params <- as.list(sdr,"Est")

for(i in 1:length(cohorts)){
  L1 <- 0
  L2 <- exp(c(params$logL2_1,params$logL2Rand))[i]
  L3 <- exp(c(params$logL3_1,params$logL3Rand))[i]
  r <- (L3-L2)/(L2-L1)
  t1<- plogis(params$t1)*0.4958904+0.04109589
  t3 <- 6
  t2 <- (t1+t3)/2
  t <- (0:1000)/100
  L <- L1+(L3-L1)*(1-r^(2*(t-t1)/(t3-t1)))*(1-r^2)^(-1)
  if(i==1){
    plot(t,L, lwd =1.5 ,type='l',ylim=c(0,120),col=colors[i])
  }
  if(i>1){
    lines(t,L, lwd =1.5 ,col=colors[i])
  }
  #cohortX <- ca_cod %>% filter(cohort==cohorts[i])
  #points(cohortX$Age+cohortX$Jday,cohortX$LngtCm)
  #readline() # press enter to continue
}

for(i in 1:length(cohorts)){
  L1 <- 0
  L2 <- exp(c(params$logL2_1,params$logL2Rand))[i]
  L3 <- exp(c(params$logL3_1,params$logL3Rand))[i]
  r <- (L3-L2)/(L2-L1)
  t1<- plogis(params$t1)*0.4958904+0.04109589
  t3 <- 6
  t2 <- (t1+t3)/2
  t <- (0:1000)/100
  L <- L1+(L3-L1)*(1-r^(2*(t-t1)/(t3-t1)))*(1-r^2)^(-1)
  if(i==1){
    plot(t,L, lwd =1.5 ,type='l',ylim=c(0,20),col=colors[i],xlim=c(0,1))
  }
  if(i>1){
    lines(t,L, lwd =1.5 ,col=colors[i])
  }
}

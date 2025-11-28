####################
library(LaplacesDemon)
library(surveyIndex)
library(RTMB)
library(sf)
library(lubridate)
library(tidyverse)
# prepare data
#####
# working directory to coefficient-file
wd_coef <- "C:/Users/mavast/Documents/GitHub/SOS data/"

dAll <- readRDS(paste(wd_coef,"DatrasExchange.R",sep=""))
WB <- st_read(paste(wd_coef,"shapefiles/ICES_areas.shp",sep=""))
plot(WB %>% filter(ICES_SUB %in% 22:24))
WB <- as(WB,"Spatial")

dAll <- addSpatialData(dAll,WB)


d<-subset(dAll, Species=="Gadus morhua",ICES_SUB %in% 22:24)
dAll<-NULL; gc(); ## lose dAll because it takes up a lot of memory
d<-addSpectrum(d,by=1)
options(warn = -1)
ca <- d[['CA']]
hh <- d[['HH']]
hl <- d[['HL']]

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
ca_cod$yq <- paste(ca_cod$Year,ca_cod$Quarter)
t <- aggregate(Jday~yq,data=ca_cod,FUN=mean)


ca_cod$qday <- t$Jday[match(ca_cod$yq,t$yq)]


wd_larvae <- "C:/Users/mavast/Downloads/EggsAndLarvaeData_0911281322"
d1 <- read.table(paste(wd_larvae,"EggsAndLarvae_EH_0911281322.csv",sep="/"),header=TRUE,sep=',')
d2 <- read.table(paste(wd_larvae,"EggsAndLarvae_EM_0911281322.csv",sep="/"),header=TRUE,sep=',')
d3 <- read.table(paste(wd_larvae,"EggsnAndLarvae_Historical_0911281322.csv",sep="/"),header=TRUE,sep=',', fill = TRUE)
x <- as.numeric(unique(substr(d1$HaulID,1,4)))

cod <- d3 %>% filter(Species=="Gadus morhua" & Stage=="LV" & !(is.na(Length)))

#rm(list=setdiff(ls(),c('par','dat')))
#####

# model growth functions
#####
#plot(yday(cod$DateTime),cod$Length)

larvae <- data.frame(haul.id="larvae", Age=yday(cod$DateTime)/365,cohort=1989,
                     logLength=log(cod$Length/10),HLNoAtLngt=cod$NoCounted,qday=0,Quarter=NA) # NB larvae lengths not on log scale


dat <- list(ca=rbind(larvae %>% filter(!is.na(logLength)),ca_cod %>% select(haul.id,Age,cohort,logLength,HLNoAtLngt,qday,Quarter))
            %>% filter(cohort==2000))
dat$ca <- dat$ca %>% filter(Age<=10) # Remove an outlier

plot(dat$ca$Age+dat$ca$qday,exp(dat$ca$logLength),xlim=c(0,13))
rm(list=setdiff(ls(),c('dat','larvae','cod','ca_cod','params')))
#%>% filter(Age>0))
#rm(list=setdiff(ls(),c('par','dat')))
par <- list()
par$logL1 <- log(15)
par$logL2 <- log(35)
par$logL3 <- log(50)
par$logSD_m <- 0 # deviation from mean fit
par$logSD_c <- 0 # year class SD
par$logSDspawn <- 0
par$tHatch <- 0.1


f <- function(par){
  getAll(par,dat)
  L1 <- exp(logL1)
  L2 <- exp(logL2)
  L3 <- exp(logL3)

  r <- (L3-L2)/(L2-L1)

  t1 <- 1
  t3 <- 6
  #t2 <- (t1+t3)/2
  SD_m <- exp(logSD_m)
  SD_c <- exp(logSD_c)
  SDspawn <- exp(logSDspawn)
  #tHatch <- 0
  
  VbgrF <- expression(L1+(L3-L1)*(1-((L3-L2)/(L2-L1))^(2*(t-t1)/(t3-t1)))*(1-((L3-L2)/(L2-L1))^2)^(-1))
  logVbgrF <- expression(log(L1+(L3-L1)*(1-((L3-L2)/(L2-L1))^(2*(t-t1)/(t3-t1)))*(1-((L3-L2)/(L2-L1))^2)^(-1)))
  difflogVBGR <- D(logVbgrF,"t")
  diffVBGR <- D(logVbgrF,"t")
  
  cohorts <- sort(unique(ca$cohort))
  ret <- 0
  for(i in 1:length(cohorts)){
    cohort.idx <- which(ca$cohort==cohorts[i])
    times.i <- sort(unique(ca$Age[cohort.idx]+ca$qday[cohort.idx]))
    
    for(j in 1:length(times.i)){
      age.idx <- which((ca$Age[cohort.idx]+ca$qday[cohort.idx])==times.i[j])
      if(length(age.idx)<2) break
      age.j <- times.i[j]-tHatch
      mean_size <- sum(ca$logLength[cohort.idx][age.idx]*ca$HLNoAtLngt[cohort.idx][age.idx])/
        sum(ca$HLNoAtLngt[cohort.idx][age.idx])
      
      weight <- ca$HLNoAtLngt[cohort.idx][age.idx]/(sum((ca$HLNoAtLngt[cohort.idx][age.idx])))

      if(ca$haul.id[cohort.idx][age.idx][1]!="larvae"){
        dL_dt <- eval(difflogVBGR,list(t=age.j,L1=L1,L2=L2,L3=L3,t1=t1,t3=t3))
        logL.pred <- eval(logVbgrF,list(t=age.j,L1=L1,L2=L2,L3=L3,t1=t1,t3=t3))
        
        LSD <- sqrt(SD_c^2+(SDspawn*dL_dt)^2)
        
        ret <- ret-sum(dnorm(mean_size,mean=logL.pred,sd=SD_m,log=TRUE)) # fit mean length
        
        ret <- ret-sum(dnorm(ca$logLength[cohort.idx][age.idx],
                             mean=mean_size,sd=LSD,
                             log=TRUE)*weight)
      }
      
      if(ca$haul.id[cohort.idx][age.idx][1]=="larvae"){
        dL_dt <- eval(diffVBGR,list(t=age.j,L1=L1,L2=L2,L3=L3,t1=t1,t3=t3))
        L.pred <- eval(VbgrF,list(t=age.j,L1=L1,L2=L2,L3=L3,t1=t1,t3=t3))
        
        LSD <- SDspawn*dL_dt
        
        ret <- ret-sum(dnorm(exp(ca$logLength[cohort.idx][age.idx]),
                             mean=L.pred,sd=LSD,
                             log=TRUE)*weight)
      }
    }
    logL.Hatch <- eval(logVbgrF,list(t=tHatch,L1=L1,L2=L2,L3=L3,t1=t1,t3=t3))
    ret <- ret-dnorm(logL.Hatch,log(0.4),SD_m)
  }
  x <- (0:1000)/100
  c <- which(cohorts==1989)
  logfit <- log(eval(VbgrF,list(t=x,L1=L1,L2=L2[c],L3=L3[c],t1=t1[c],t3=t3)))
  ADREPORT(logfit)
  #dLdSpawn <- eval(diffVBGR,list(t=x,L1=L1,L2=L2,L3=L3,t1=t1,t3=t3))
  #LogLSD <- sqrt(SD_c^2+(SDspawn*dLdSpawn)^2)
  #ADREPORT(LogLSD)
  
  ret
}
obj <- MakeADFun(f,par)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr
M1 <- opt

pl <- as.list(sdr,"Est",report=TRUE)$logfit


coefs <- exp(summary(sdr)[1:6,1])
# define gowth- and differential equations to calculate total SD
logVbgrF <- expression(log(L1+(L3-L1)*(1-((L3-L2)/(L2-L1))^(2*(t-t1)/(t3-t1)))*(1-((L3-L2)/(L2-L1))^2)^(-1)))
diffVBGR <- D(logVbgrF,"t")
dL_dt <- eval(diffVBGR,list(t=(0:1000)/100,
                            L1=coefs[1],L2=coefs[2],L3=coefs[3],t1=1,t3=6))

totSD <- sqrt(coefs[5]^2+(coefs[6]*dL_dt)^2)
d <- aggregate(logLength~Age+qday,data=dat$ca,FUN=mean)
d.sd <- aggregate(logLength~Age+qday,data=dat$ca,FUN=sd)
dL_dt.data <- eval(diffVBGR,list(t=d$Age+d$qday,
                            L1=coefs[1],L2=coefs[2],L3=coefs[3],t1=1,t3=6))

totSD.data <- sqrt(coefs[5]^2+(coefs[6]*dL_dt.data)^2)

par(mfrow=c(1,1))
plot((0:1000)/100,exp(pl) , lwd =3 ,type='l',
     ylim=c(0,100),col="black",ylab="Length [cm]",xlab="Age [years]",
     main="2000 cohort")
lines((0:1000)/100,exp(pl -totSD), lty = "dotted" , lwd =2.5,col="orange")
lines((0:1000)/100,exp(pl +totSD), lty = "dotted" , lwd =2.5,col="orange")  
points(dat$ca$Age+dat$ca$qday,exp(dat$ca$logLength))
points(d$Age+d$qday,exp(d$logLength),pch=19,cex=2)
points(d$Age+d$qday,exp(d$logLength-totSD.data),pch=19,cex=1,col="orange")
points(d$Age+d$qday,exp(d$logLength+totSD.data),pch=19,cex=1,col="orange")
points(d$Age+d$qday,exp(d$logLength-d.sd$logLength),pch=19,cex=1,col="red")
points(d$Age+d$qday,exp(d$logLength+d.sd$logLength),pch=19,cex=1,col="red")

legend(7,40,legend=c("data mean","data SD","modelled SD"),pch=19,col=c("black","red","orange"))

# Year 1
plot((0:1000)/100,exp(pl) , lwd =3 ,type='l',
     ylim=c(0,25),col="black",ylab="Length [cm]",xlab="Age [years]",
     main="2000 cohort",xlim=c(0,1))
lines((0:1000)/100,exp(pl -totSD), lty = "dotted" , lwd =2.5,col="orange")
lines((0:1000)/100,exp(pl +totSD), lty = "dotted" , lwd =2.5,col="orange")  
points(dat$ca$Age+dat$ca$qday,exp(dat$ca$logLength))
points(d$Age+d$qday,exp(d$logLength),pch=19,cex=2)
points(d$Age+d$qday,exp(d$logLength-totSD.data),pch=19,cex=1,col="orange")
points(d$Age+d$qday,exp(d$logLength+totSD.data),pch=19,cex=1,col="orange")
points(d$Age+d$qday,exp(d$logLength-d.sd$logLength),pch=19,cex=1,col="red")
points(d$Age+d$qday,exp(d$logLength+d.sd$logLength),pch=19,cex=1,col="red")

legend(0,15,legend=c("data mean","data SD","modelled SD"),pch=19,col=c("black","red","orange"))
#####


########################
#Conseptual simulations
#######################
# cohort spread
#####
vbgr <- function(t,L1=exp(2.737346),L2=exp(3.863767),L3=exp(4.256434),t1=1,t3=6) {
  # exponent terms
  
  t2 <- (t1+t3)/2
  
  r <- (L3-L2)/(L2-L1)
  # main formula: length zero at age 0
  L <- L1+(L3-L1)*(1-r^(2*(t-t1)/(t3-t1)))*(1-r^2)^(-1)
  
  return(L)
}
logSD <- 0.0#exp(-0.86121831)



x <- (0:(3650))/365

L <- matrix(rep(NA,length(x)*3000),ncol=3000)
par(mfrow=c(1,1))
vbgrM <- vbgr(x)
plot(x,log(vbgrM),
     type='l',lwd=2,ylab="length",xlab="age")
lines(c(-10,100),c(0,0),lty="dashed")
lines(c(-10,100),log(c(exp(2.737346),exp(2.737346))),lty="dashed")
lines(c(-10,100),log(c(exp(3.863767),exp(3.863767))),lty="dashed")
lines(c(-10,100),log(c(exp(4.256434),exp(4.256434))),lty="dashed")

for(i in 1:3000){
  t <- x+rnorm(1,mean=0,sd=31/365)
  SD <- rnorm(1,0,logSD)
  L[,i] <- (vbgr(t))+SD
  #lines(x,(L[,i]),type='l',lwd=2,ylab="cm",xlab="age")
}
#m <- apply((L), 1, mean,na.rm = TRUE)
#L[which(L<0)] <- NA
SDs <- apply((L), 1, sd,na.rm = TRUE)
par(mfrow=c(1,1))

plot(x,SDs,type='l',lwd=2)
logVbgrF <- expression((L1+(L3-L1)*(1-((L3-L2)/(L2-L1))^(2*(t-t1)/(t3-t1)))*(1-((L3-L2)/(L2-L1))^2)^(-1)))
diffVBGR <- D(logVbgrF,"t")
dL_dt <- eval(diffVBGR,list(t=x,L1=exp(2.737346),L2=exp(3.863767),L3=exp(4.256434),t1=1,t3=6))

lines(x,sqrt((dL_dt*(31/365))^2+logSD^2),lwd=2,col="red")
#####       

# seasonality
#####
x <- (1:100)/100
A <- 0
phi <- 0.5

s <- A*sin(2*pi*(x-phi))
plot(x,s)
plot((0:1000)/100,exp(pl)+s, lwd =3 ,type='l',
     ylim=c(0,25),col="black",ylab="Length [cm]",xlab="Age [years]",
     main="2000 cohort",xlim=c(0,1))
plot((0:1000)/100,exp(pl)+s, lwd =3 ,type='l',
     ylim=c(0,100),col="black",ylab="Length [cm]",xlab="Age [years]",
     main="2000 cohort")
#####
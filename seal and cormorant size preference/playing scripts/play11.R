library(RTMB)
library(tidyverse)
library(sf)
library(DATRAS)
# prepare data
#####
# working directory to coefficient-file
wd_coef <- "C:/Users/mavast/Documents/GitHub/SOS data/"

#dAll <- readRDS(paste(wd_coef,"DatrasExchange.R",sep=""))
dAll <- readRDS("datras.R")
WB <- st_read(paste(wd_coef,"shapefiles/ICES_areas.shp",sep=""))
#plot(WB %>% filter(ICES_SUB %in% 22:24))
WB <- as(WB,"Spatial")

dAll <- addSpatialData(dAll,WB)

# cod
d<-subset(dAll, Species=="Gadus morhua",ICES_SUB %in% 22:24)
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
#hist(ca_cod$LngtCm[ca_cod$Age==0])
ca_cod <- ca_cod %>% filter(!(Age<1 & Country=="DE" & Year %in% c(2015,2018)))
#hist(ca_cod$LngtCm[ca_cod$Age==0])

cohorts <- sort(unique(ca_cod$cohort))[-c(1:6)] # cohort with sufficient observations to fit
ca_cod <- ca_cod %>% filter(cohort %in% cohorts)
ca_cod$logLength <- log(ca_cod$LngtCm)
ca_cod$aq <- paste(ca_cod$Age,ca_cod$Quarter)
ca_cod$yq <- paste(ca_cod$Year,ca_cod$Quarter)
t <- aggregate(Jday~yq,data=ca_cod,FUN=mean)
ca_cod$qday <- t$Jday[match(ca_cod$yq,t$yq)]
rm(list=setdiff(ls(),c('ca_cod','dAll')))

#####
cohortX<- 1990:2022
ca_cod$t <- ca_cod$Age+ca_cod$qday
ca_cod$Length <- exp(ca_cod$logLength)

df <- ca_cod %>% select(Length,t,cohort,HLNoAtLngt,haul.id,Year)

ca_cod$n <- 1
df.n <- aggregate(n~haul.id+Length,data=ca_cod,FUN=sum)

df <- left_join(df,df.n,)
df$mul <- df$HLNoAtLngt/df$n


dat <- as.list(df %>% filter(cohort %in% cohortX))

n_cohorts <- length(unique(cohortX))
par <- list()
par$logL1.1 <- log(15)
par$logL2.1 <- log(35)
par$logL3.1 <- log(50)
par$logL1.rand <- rep(log(15),n_cohorts-1)
par$logL2.rand <- rep(log(35),n_cohorts-1)
par$logL3.rand <- rep(log(50),n_cohorts-1)

par$logSD_m <- 0 # deviation from mean fit
par$logSD_c <- 0 # year class SD
par$logSDage <- 0
par$logRWsd <- 0
#par$tHatch <- 0.1

f <- function(par){
  getAll(par,dat)
  L1 <- exp(c(logL1.1,logL1.rand))
  L2 <- exp(c(logL2.1,logL2.rand))
  L3 <- exp(c(logL3.1,logL3.rand))
  t_hatch <- 61/365
  
  r <- (L3-L2)/(L2-L1)
  
  t1 <- 1-t_hatch
  t3 <- 5-t_hatch
  #t2 <- (t1+t3)/2
  SD_m <- exp(logSD_m)
  SD_c <- exp(logSD_c)
  SD_age <- exp(logSDage)
  SD_RW <- exp(logRWsd)
  
  VbgrF <- expression(L1+(L3-L1)*(1-((L3-L2)/(L2-L1))^(2*(age-t1)/(t3-t1)))*(1-((L3-L2)/(L2-L1))^2)^(-1))
  diffVBGR <- D(VbgrF,"age")
  
  cohorts <- sort(unique(cohort))
  
  df_mean <- matrix(NA,ncol=4)
  ret <- 0
  
  for(i in 1:length(cohorts)){
    cohort.idx <- which(cohort==cohorts[i])
    times.i <- sort(unique(t[cohort.idx]))
    if(i>1){
      ret <- ret-dnorm(log(L1[i]),mean=log(L1[i-1]),sd=SD_RW,log=TRUE)
      ret <- ret-dnorm(log(L2[i]),mean=log(L2[i-1]),sd=SD_RW,log=TRUE)
      ret <- ret-dnorm(log(L3[i]),mean=log(L3[i-1]),sd=SD_RW,log=TRUE)
    }
    
    for(j in 1:length(times.i)){
      age.idx <- which((t[cohort.idx])==times.i[j])
      #mean_idx <- which(t_mean==times.i[j] & cohort_mean==cohorts[i])
      age.j <- times.i[j]-t_hatch
      #      mean_size <- L_mean[mean_idx]
      #weighted mean
      mean_size <- sum(Length[cohort.idx][age.idx]*mul[cohort.idx][age.idx])/sum(mul[cohort.idx][age.idx])
      #weighted std. dev.
      mean_sd <- sqrt( sum(mul[cohort.idx][age.idx]*(Length[cohort.idx][age.idx]-mean_size)^2)/sum(mul[cohort.idx][age.idx]))
      
      sd_fac <- ifelse(length(age.idx)>1,1,0)
      
      weight <- mul[cohort.idx][age.idx]/sqrt(sum((mul[cohort.idx][age.idx])))*sd_fac
      
      weight_m <- sum(mul[cohort.idx][age.idx])/(sum(mul[cohort.idx][age.idx]))# or /sqrt(sum(mul[cohort.idx][age.idx]))
      
      dL_dt <- eval(diffVBGR,list(age=age.j,L1=L1[i],L2=L2[i],L3=L3[i],
                                  t1=t1,t3=t3))
      
      L.pred <- eval(VbgrF,list(age=age.j,L1=L1[i],L2=L2[i],L3=L3[i],
                                t1=t1,t3=t3))
      
      LSD <- sqrt((SD_age*dL_dt)^2+(L.pred*SD_c)^2)
      
      ret <- ret-dnorm(mean_size,mean=L.pred,sd=SD_m*L.pred,log=TRUE)*weight_m # fit mean length
      if(length(age.idx)<2) next
      
      ret <- ret-sum(dnorm(Length[cohort.idx][age.idx],
                           mean=mean_size,sd=LSD,
                           log=TRUE)*weight)
      
      df_mean <- rbind(df_mean,c(round(cohorts[i]),times.i[j],mean_size,mean_sd))
    }
    L.Hatch <- eval(VbgrF,list(age=t_hatch,L1=L1[i],L2=L2[i],L3=L3[i],
                               t1=t1,t3=t3))
    ret <- ret-dnorm(L.Hatch,0.4,0.01,log=TRUE)#*10^4 play with sd here if the model doesn't converge
    
  }
  x <- (0:1000)/100
  c <- which(cohorts==cohortX)
  # evaluate for x
  dL_dt <- eval(diffVBGR,list(age=x,L1=L1[c],L2=L2[c],L3=L3[c],
                              t1=t1,t3=t3))
  fit <- eval(VbgrF,list(age=x,L1=L1[c],L2=L2[c],L3=L3[c],
                         t1=t1,t3=t3))
  tot_sd <- sqrt((SD_age*dL_dt)^2+(fit*SD_c)^2)
  
  #evaluate for data
  df_mean <- df_mean[-1,]
  
  ADREPORT(fit)
  ADREPORT(tot_sd)
  ADREPORT(df_mean)
  
  ret
}
obj <- MakeADFun(f,par,silent=TRUE,random=c("logL1.rand","logL2.rand","logL3.rand"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

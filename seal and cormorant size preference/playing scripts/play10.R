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

# flounder
d<-subset(dAll, Species=="Platichthys flesus",ICES_SUB %in% 22:24)
d<-addSpectrum(d,by=1)
ca <- d[['CA']]
hh <- d[['HH']]
hl <- d[['HL']]

hl <- aggregate(HLNoAtLngt~haul.id+LngtCm,data=hl,FUN=sum)

ca$lengthID <- paste(ca$haul.id,ca$LngtCm)
hl$lengthID <- paste(hl$haul.id,hl$LngtCm)
ca_flounder <- left_join(ca,hl %>% select(lengthID,HLNoAtLngt))
ca_flounder$cohort <- as.numeric(as.character(ca_flounder$Year))-ca_flounder$Age
#ca_flounder <- ca_flounder %>% filter(Age>0 & !is.na(HLNoAtLngt))
ca_flounder <- ca_flounder %>% filter(!is.na(HLNoAtLngt))
hh$Jday <- yday(make_date(year=hh$Year,month=hh$Month,day=hh$Day))/yday(make_date(year=hh$Year,month=12,day=31))
ca_flounder <- left_join(ca_flounder,hh %>% select(haul.id,Jday))
ca_flounder <- ca_flounder %>% filter(!(is.na(Jday) | is.na(Age)))
ca_flounder <- ca_flounder %>% filter(!(Age<1 & Quarter==1))
hist(ca_flounder$LngtCm[ca_flounder$Age==0])
plot(ca_flounder$Age+ca_flounder$Jday,ca_flounder$LngtCm)

ca_flounder <- ca_flounder %>% filter(Age<50)
table(ca_flounder$cohort)
cohorts <- sort(unique(ca_flounder$cohort))[-c(1:4,(length(unique(ca_flounder$cohort))-3):length(unique(ca_flounder$cohort)))] # cohort with sufficient observations to fit
ca_flounder <- ca_flounder %>% filter(cohort %in% cohorts)
ca_flounder$logLength <- log(ca_flounder$LngtCm)
ca_flounder$aq <- paste(ca_flounder$Age,ca_flounder$Quarter)
ca_flounder$yq <- paste(ca_flounder$Year,ca_flounder$Quarter)
t <- aggregate(Jday~yq,data=ca_flounder,FUN=mean)
ca_flounder$qday <- t$Jday[match(ca_flounder$yq,t$yq)]
rm(list=setdiff(ls(),c('ca_cod','ca_flounder','dAll')))

#####

# cod growth functions
#####

ca_cod <- ca_cod %>% filter(Age<11) # Remove an outliers

ca_cod$t <- ca_cod$Age+ca_cod$qday
ca_cod$Length <- exp(ca_cod$logLength)


df <- ca_cod %>% select(Length,t,cohort,HLNoAtLngt,haul.id,Year)
cohortX <- min(ca_cod$cohort):(max(ca_cod$cohort)-2)

ca_cod$n <- 1
df.n <- aggregate(n~haul.id+Length,data=ca_cod,FUN=sum)

df <- left_join(df,df.n,)
df$mul <- df$HLNoAtLngt/df$n


dat <- as.list(df %>% filter(cohort %in% cohortX))
dat$cohortX <- 2020
#plot(dat$t,dat$Length)

rm(list=setdiff(ls(),c('dat','ca_cod','cod')))
n_cohorts <- length(unique(dat$cohort))
#%>% filter(Age>0))
#rm(list=setdiff(ls(),c('par','dat')))
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
  
  t1 <- 1
  t3 <- 5
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

fit <- as.list(sdr,"Est",report=TRUE)$fit
tot_sd <- as.list(sdr,"Est",report=TRUE)$tot_sd
dat_means <- as.list(sdr,"Est",report=TRUE)$df_mean
cohortX <- dat$cohortX
m.idx <- which(dat_means[,1]==cohortX)
coefs <- exp(summary(sdr)[1:7,1])
VbgrF <- expression(L1+(L3-L1)*(1-((L3-L2)/(L2-L1))^(2*(age-t1)/(t3-t1)))*(1-((L3-L2)/(L2-L1))^2)^(-1))
diffVBGR <- D(VbgrF,"age")
dL_dt <- eval(diffVBGR,list(age=dat_means[m.idx,2],L1=coefs[1],L2=coefs[2],L3=coefs[3],
                            t1=1,t3=6))
dat_sd <- sqrt((coefs[6]*dL_dt)^2+(dat_means[m.idx,3]*coefs[5])^2)



par(mfrow=c(1,1))
plot((0:1000)/100,(fit) , lwd =3 ,type='l',
     ylim=c(0,120),col="black",ylab="Length [cm]",xlab="Age [years]",
     main=paste(cohortX,"cohort"))
lines((0:1000)/100,(fit -tot_sd), lty = "dotted" , lwd =2.5,col="orange")
lines((0:1000)/100,(fit +tot_sd), lty = "dotted" , lwd =2.5,col="orange")  
points(dat$t[dat$cohort==cohortX],dat$Length[dat$cohort==cohortX])
points(dat_means[m.idx,2],dat_means[m.idx,3],pch=19,cex=2)
points(dat_means[m.idx,2],dat_means[m.idx,3]-dat_sd,pch=19,cex=1,col="orange")
points(dat_means[m.idx,2],dat_means[m.idx,3]+dat_sd,pch=19,cex=1,col="orange")
points(dat_means[m.idx,2],dat_means[m.idx,3]-dat_means[m.idx,4],pch=19,cex=1,col="red")
points(dat_means[m.idx,2],dat_means[m.idx,3]+dat_means[m.idx,4],pch=19,cex=1,col="red")
lines(c(-10,100),c(0,0),lty="dashed")
lines(c( 61/365, 61/365),c(-10,200),lty="dashed")

legend(7,40,legend=c("data mean","data SD","modelled SD"),pch=19,col=c("black","red","orange"))
#####

plot((-100:1000)/100,(fit) , lwd =3 ,type='l',xlim=c(-1,2),
     ylim=c(0,40),col="black",ylab="Length [cm]",xlab="Age [years]",
     main=paste(cohortX,"cohort"))
lines((-100:1000)/100,(fit -tot_sd), lty = "dotted" , lwd =2.5,col="orange")
lines((-100:1000)/100,(fit +tot_sd), lty = "dotted" , lwd =2.5,col="orange")  
points(dat$t[dat$cohort==cohortX],dat$Length[dat$cohort==cohortX])
points(dat_means[m.idx,2],dat_means[m.idx,3],pch=19,cex=2)
points(dat_means[m.idx,2],dat_means[m.idx,3]-dat_sd,pch=19,cex=1,col="orange")
points(dat_means[m.idx,2],dat_means[m.idx,3]+dat_sd,pch=19,cex=1,col="orange")
points(dat_means[m.idx,2],dat_means[m.idx,3]-dat_means[m.idx,4],pch=19,cex=1,col="red")
points(dat_means[m.idx,2],dat_means[m.idx,3]+dat_means[m.idx,4],pch=19,cex=1,col="red")
lines(c(-10,100),c(0,0),lty="dashed")
lines(c( 61/365, 61/365),c(-10,200),lty="dashed")

legend(1.5,10,legend=c("data mean","data SD","modelled SD"),pch=19,col=c("black","red","orange"))


for(i in 1:10){
  
}

#####

# flounder growth functions
#####

ca_flounder$t <- ca_flounder$Age+ca_flounder$qday
ca_flounder$Length <- exp(ca_flounder$logLength)
plot(ca_flounder$t,ca_flounder$Length)

df <- ca_flounder %>% select(Length,t,cohort,HLNoAtLngt,haul.id,Year)
ca_flounder$n <- 1
df.n <- aggregate(n~haul.id+Length,data=ca_flounder,FUN=sum)

df <- left_join(df,df.n,)
df$mul <- df$HLNoAtLngt/df$n

cohortX <- 2000:2001
dat <- as.list(df %>% filter(cohort %in% cohortX) %>% select(-c(HLNoAtLngt,haul.id,n)))
dat$cohortX <- cohortX
#plot(dat$t,dat$Length)

rm(list=setdiff(ls(),c('dat','ca_flounder','cod')))
n_cohorts <- length(unique(dat$cohort))
#%>% filter(Age>0))
#rm(list=setdiff(ls(),c('par','dat')))
par <- list()
par$logL1.1 <- log(10)
par$logL2.1 <- log(30)
par$logL3.1 <- log(40)
par$logL1.rand <- rep(log(10),n_cohorts-1)
par$logL2.rand <- rep(log(30),n_cohorts-1)
par$logL3.rand <- rep(log(40),n_cohorts-1)

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
  t_hatch <- ((75+197)/2)/365
  
  r <- (L3-L2)/(L2-L1)
  
  t1 <- 1
  t3 <- 6
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
    ret <- ret-dnorm(L.Hatch,0.6,0.01,log=TRUE)#*10^4 play with sd here if the model doesn't converge
    
  }
  x <- (0:2000)/100
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
obj <- MakeADFun(f,par,silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

fit <- as.list(sdr,"Est",report=TRUE)$fit
tot_sd <- as.list(sdr,"Est",report=TRUE)$tot_sd
dat_means <- as.list(sdr,"Est",report=TRUE)$df_mean
cohortX <- 1997
m.idx <- which(dat_means[,1]==cohortX)
coefs <- exp(summary(sdr)[1:7,1])

VbgrF <- expression(L1+(L3-L1)*(1-((L3-L2)/(L2-L1))^(2*(t-t1)/(t3-t1)))*(1-((L3-L2)/(L2-L1))^2)^(-1))
diffVBGR <- D(VbgrF,"t")
dL_dt <- eval(diffVBGR,list(t=dat_means[m.idx,2],L1=coefs[1],L2=coefs[2],L3=coefs[3],
                            t1=1,t3=6))
dat_sd <- sqrt((coefs[6]*dL_dt)^2+(dat_means[m.idx,3]*coefs[5])^2)

par(mfrow=c(1,1))
plot((-100:2000)/100,(fit) , lwd =3 ,type='l',
     ylim=c(0,120),col="black",ylab="Length [cm]",xlab="Age [years]",
     main=paste(cohortX,"cohort"))
lines((-100:2000)/100,(fit -tot_sd), lty = "dotted" , lwd =2.5,col="orange")
lines((-100:2000)/100,(fit +tot_sd), lty = "dotted" , lwd =2.5,col="orange")  
points(dat$t[dat$cohort==cohortX],dat$Length[dat$cohort==cohortX])
points(dat_means[m.idx,2],dat_means[m.idx,3],pch=19,cex=2)
points(dat_means[m.idx,2],dat_means[m.idx,3]-dat_sd,pch=19,cex=1,col="orange")
points(dat_means[m.idx,2],dat_means[m.idx,3]+dat_sd,pch=19,cex=1,col="orange")
points(dat_means[m.idx,2],dat_means[m.idx,3]-dat_means[m.idx,4],pch=19,cex=1,col="red")
points(dat_means[m.idx,2],dat_means[m.idx,3]+dat_means[m.idx,4],pch=19,cex=1,col="red")
lines(c(-10,100),c(0,0),lty="dashed")
lines(c(((75+197)/2)/365,((75+197)/2)/365),c(-10,200),lty="dashed")

legend(7,40,legend=c("data mean","data SD","modelled SD"),pch=19,col=c("black","red","orange"))


plot((-100:2000)/100,(fit) , lwd =3 ,type='l',xlim=c(-1,2),
     ylim=c(0,40),col="black",ylab="Length [cm]",xlab="Age [years]",
     main=paste(cohortX,"cohort"))
lines((-100:2000)/100,(fit -tot_sd), lty = "dotted" , lwd =2.5,col="orange")
lines((-100:2000)/100,(fit +tot_sd), lty = "dotted" , lwd =2.5,col="orange")  
points(dat$t[dat$cohort==cohortX],dat$Length[dat$cohort==cohortX])
points(dat_means[m.idx,2],dat_means[m.idx,3],pch=19,cex=2)
points(dat_means[m.idx,2],dat_means[m.idx,3]-dat_sd,pch=19,cex=1,col="orange")
points(dat_means[m.idx,2],dat_means[m.idx,3]+dat_sd,pch=19,cex=1,col="orange")
points(dat_means[m.idx,2],dat_means[m.idx,3]-dat_means[m.idx,4],pch=19,cex=1,col="red")
points(dat_means[m.idx,2],dat_means[m.idx,3]+dat_means[m.idx,4],pch=19,cex=1,col="red")
lines(c(-10,100),c(0,0),lty="dashed")
lines(c(((75+197)/2)/365,((75+197)/2)/365),c(-10,200),lty="dashed")

legend(1.5,10,legend=c("data mean","data SD","modelled SD"),pch=19,col=c("black","red","orange"))


#####



########################
#Conseptual simulations - cod
#######################
#####
vbgr <- function(t,L1=exp(2.737346),L2=exp(3.863767),L3=exp(4.256434),t1=1,t3=6) {
  # exponent terms
  
  t2 <- (t1+t3)/2
  
  r <- (L3-L2)/(L2-L1)
  # main formula: length zero at age 0
  L <- L1+(L3-L1)*(1-r^(2*(t-t1)/(t3-t1)))*(1-r^2)^(-1)
  
  return(L)
}
logSD <- 0.1#exp(-0.86121831)
t_sd <- 31/365



x <- (0:(3650*2))/365

L <- matrix(rep(NA,length(x)*10000),ncol=10000)
par(mfrow=c(1,1))
vbgrM <- vbgr(x)
plot(x,(vbgrM),
     type='l',lwd=2,ylab="length",xlab="age",ylim=c(0,200))
lines(c(-10,100),c(0,0),lty="dashed")
lines(c(-10,100),(c(exp(2.737346),exp(2.737346))),lty="dashed")
lines(c(-10,100),(c(exp(3.863767),exp(3.863767))),lty="dashed")
lines(c(-10,100),(c(exp(4.256434),exp(4.256434))),lty="dashed")
si <- rep(0,1000)
for(i in 1:10000){
  t <- x+rnorm(1,mean=0,sd=t_sd)
  SD <- rnorm(1,1,logSD) # this is exp(SD_log), where SD_log is the standard deviation of length on log-scale
  L[,i] <- (vbgr(t))*SD
  #lines(x,(L[,i]),type='l',lwd=2,ylab="cm",xlab="age")
  si[i] <- SD
}
#m <- apply((L), 1, mean,na.rm = TRUE)
#L[which(L<0)] <- NA
SDs <- apply((L), 1, sd,na.rm = TRUE)
par(mfrow=c(1,1))

plot(x,SDs,type='l',lwd=2)
logVbgrF <- expression((L1+(L3-L1)*(1-((L3-L2)/(L2-L1))^(2*(t-t1)/(t3-t1)))*(1-((L3-L2)/(L2-L1))^2)^(-1)))
diffVBGR <- D(logVbgrF,"t")
dL_dt <- eval(diffVBGR,list(t=x,L1=exp(2.737346),L2=exp(3.863767),L3=exp(4.256434),t1=1,t3=6))

lines(x,sqrt((dL_dt*(t_sd))^2+(vbgrM*logSD)^2),lwd=2,col="red")


#####  

########################
#Conseptual simulations - cod
#######################
#####
vbgr <- function(t,L1=11,L2=30,L3=40,t1=1,t3=10) {
  # exponent terms
  
  t2 <- (t1+t3)/2
  
  r <- (L3-L2)/(L2-L1)
  # main formula: length zero at age 0
  L <- L1+(L3-L1)*(1-r^(2*(t-t1)/(t3-t1)))*(1-r^2)^(-1)
  
  return(L)
}
plot(dat$t_mean,dat$L_mean,ylim=c(0,60),xlim=c(0,20))
lines((0:200)/10,vbgr((0:200)/10),lwd=2,col="red")

#####  
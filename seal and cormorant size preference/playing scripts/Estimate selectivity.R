library(LaplacesDemon)
library(surveyIndex)
library(RTMB)

# working directory to coefficient-file
wd_coef <- "C:/Users/mavast/Documents/GitHub/SOS data/"
coefs <- readRDS(paste(wd_coef,"coefficients.R",sep=""))

# Function to calculate fish length as a function of age.
vbgr <- function(age, k,Linf) {
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf
  
  return(L_t)
} 

last_index <- function(v,i){ # v is a vector and i is how many steps from the last observation you wish
  idx <- (length(v)-i+1):(length(v))
  idx
}

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
ca_cod <- ca_cod %>% filter(Age<9 & !is.na(HLNoAtLngt))
#ca_cod <- ca_cod %>% filter(Age>=1 & Age<9 & !is.na(HLNoAtLngt))
#ca_cod <- ca_cod %>% filter(!(Age==1 & Quarter==1))
hh$Jday <- yday(make_date(year=hh$Year,month=hh$Month,day=hh$Day))/yday(make_date(year=hh$Year,month=12,day=31))
ca_cod <- left_join(ca_cod,hh %>% select(haul.id,Jday))
ca_cod <- ca_cod %>% filter(!is.na(Jday))
ca_cod <- ca_cod %>% filter(!is.na(Age))
ca_cod <- ca_cod %>% filter(!(Age<1 & LngtCm>30))


table(ca_cod$cohort)
cohorts <- sort(unique(ca_cod$cohort))[-c(1:6,38:40)] # cohort with sufficient observations to fit

ca_cod$logLength <- log(ca_cod$LngtCm)

# parameters
par <- list(logK=rep(log(0.1603331),length(cohorts)),
            logLinf=rep(log(112.2677),length(cohorts)),
            #transT_hatch=rep(0.02210022,length(cohorts)),
            transT_hatch = 0,logitW=0,logitU=0,
            logSD=0)

dat <- list(ca=ca_cod %>% select(haul.id,Age,cohort,logLength,HLNoAtLngt,Jday,Quarter) %>% filter(cohort %in% cohorts))
#rm(list=setdiff(ls(),c('par','dat')))
# run diet model
#####
f <- function(par){
  getAll(par,dat)
  k <- exp(logK)
  Linf <- exp(logLinf)
  t_hatch <- plogis(transT_hatch)*0.4958904+0.04109589 # confine the hatching time between Jan. 15th (Julian day=0.0411) and July 15th (Julian day=0.5369) (Fiskeatlas)
  #t_hatch <- 0.2917808 
  SD <- exp(logSD)
  u <- plogis(logitU)
  w <- plogis(logitW)
  #u <- 0.23 #exp(logU)
  #w <- 0.91 #exp(logW)
  
  
  # Function to calculate fish length as a function of t, Age (years) + time of year.
  vbgr <- function(t, k,Linf,t_hatch) {
    age <- t-t_hatch
    X_t <- t-floor(t) # time of year [0,1]
    phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
    phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching
    L_hatch <- 0.4
    
    L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf
    
    return(L_t)
  }
  
  
  cohorts <- sort(unique(ca$cohort))
  ret <- 0
  for(i in 1:length(cohorts)){
    cohort.idx <- which(ca$cohort==cohorts[i])
    ages.i <- sort(unique(ca$Age[cohort.idx]+ca$Jday[cohort.idx]))
    ages.i <- ages.i[ages.i>0.5]
    #ages.i <- sort(unique(ca$Age[cohort.idx]+as.numeric(as.character(ca$Quarter[cohort.idx]))*0.25-0.25/2))
    
    #plot(ca$Age[cohort.idx]+ca$Jday[cohort.idx],ca$logLength[cohort.idx],main=cohorts[i],col="white")
    for(j in 1:length(ages.i)){
      age.idx <- which((ca$Age[cohort.idx]+ca$Jday[cohort.idx])==ages.i[j])
      vbgr.fit <- vbgr(ages.i[j],k[i],Linf[i],t_hatch)
      #vbgr.sd <- calculate_var()
      #sd.fit <- vbgr(ages.i[j],k[i],Linf[i],t_hatch)
      ret <- ret-sum(dnorm(ca$logLength[cohort.idx][age.idx],mean=log(vbgr.fit),sd=SD,log=TRUE)*ca$HLNoAtLngt[cohort.idx][age.idx]/(log(sum(ca$HLNoAtLngt[cohort.idx][age.idx]))+0.1))                             
      #  points(ages.i[j],log(vbgr.fit),cex=3,col="darkred",pch=19)
      # points(ca$Age[cohort.idx][age.idx]+ca$Jday[cohort.idx][age.idx],ca$logLength[cohort.idx][age.idx],cex=1,col="black",pch=1)
    }
    #lines((0:100)/10,log(vbgr((0:100)/10,0.1603331,112.2677,0.2917808)),lwd=3,col="darkred")
    #readline() # press enter to continue
  }
  cohort.idx <- 2016-min(ca$cohort)+1
  logVbgr1997 <- log(vbgr((0:100)/10,k[cohort.idx],Linf[cohort.idx],t_hatch))
  ADREPORT(logVbgr1997)
  ret
}
obj <- MakeADFun(f,par,silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

SDest <- exp(as.list(sdr,"Est")$logSD)
pl <- as.list(sdr,"Est",report=TRUE)$logVbgr1997
plsd <- as.list(sdr,"Std",report=TRUE)$logVbgr1997

cohortX <- ca_cod %>% filter(cohort==2016)

plot((0:100)/10,exp(pl),type='l',lwd=3,xlab="Age",ylab="cm",xlim=c(0,2),ylim=c(0,40))
lines(c(0.85,0.85),c(-5,200))
lines(c(-10,20),c(0,0))
lines(c(-10,20),c(15,15))
points(cohortX$Age+cohortX$Jday,cohortX$LngtCm)


plot((0:100)/10,exp(pl) , lwd =2 ,type='l',ylim=c(0,100))
lines((0:100)/10,exp(pl -2 * sqrt(plsd^2+SDest^2)), lty = "dotted" , lwd =1.5  )
lines((0:100)/10,exp(pl +2 * sqrt(plsd^2+SDest^2)), lty = "dotted" , lwd =1.5 )  
lines((0:100)/10,exp(pl -2 * plsd), lty = "dotted" , lwd =1.5  ,col="red")
lines((0:100)/10,exp(pl +2 * plsd), lty = "dotted" , lwd =1.5  ,col="red") 
lines((0:100)/10,exp(pl -2 * exp(-1.48400125  )), lty = "dotted" , lwd =1.5 ,col="orange" )
lines((0:100)/10,exp(pl +2 * exp(-1.48400125  )), lty = "dotted" , lwd =1.5 ,col="orange") 
points(cohortX$Age+cohortX$Jday,cohortX$LngtCm)

plot((0:100)/10,exp(pl),type='l',lwd=3,xlab="Age",ylab="cm",xlim=c(0,2),ylim=c(0,40))
lines(c(0.85,0.85),c(-5,200))
lines(c(-10,20),c(0,0))
lines(c(-10,20),c(15,15))
points(ca_cod$Age+ca_cod$Jday,ca_cod$LngtCm)
#points(df_ag$Age+df_ag$Jday,df_ag$LngtCm,pch=19,cex=1,col="orange")


####################
#selectivity analysis
####################
x <- as.data.frame(summary(sdr))
x <- x[1:(length(x[,1])-length(pl)),]
x$cohort <- c(rep(cohorts,2),rep(NA,4))

ca_cod <- left_join(ca,hl %>% select(lengthID,HLNoAtLngt))
ca_cod$cohort <- as.numeric(as.character(ca_cod$Year))-ca_cod$Age
ca_cod <- ca_cod %>% filter(Age<2 & !is.na(HLNoAtLngt))
ca_cod <- ca_cod %>% filter(!(Age==1 & Quarter==4))
hh$Jday <- yday(make_date(year=hh$Year,month=hh$Month,day=hh$Day))/yday(make_date(year=hh$Year,month=12,day=31))
ca_cod <- left_join(ca_cod,hh %>% select(haul.id,Jday))
ca_cod <- ca_cod %>% filter(!is.na(Jday))
ca_cod <- ca_cod %>% filter(!is.na(Age))
ca_cod <- ca_cod %>% filter(!(Age<1 & LngtCm>30))




ca_cod$logLength <- log(ca_cod$LngtCm)


rm(list=setdiff(ls(),c('x','ca_cod','cohorts')))


# parameters
par <- list(logL50=log(9.01),
            logSR=log(1.61))
dat <- list(ca=ca_cod %>% select(haul.id,Age,cohort,logLength,HLNoAtLngt,Jday,Quarter) %>% 
              filter(cohort %in% cohorts),
            parameters=x)
#rm(list=setdiff(ls(),c('par','dat')))
# run diet model
#####
f <- function(par){
  getAll(par,dat)
  x_cohor <- unique(x$cohort)[!is.na(unique(x$cohort))]
  k <- exp(x$Estimate[1:length(x_cohor)])
  Linf <- exp(x$Estimate[(length(x_cohor)+1):(2*length(x_cohor))])
  t_hatch <- plogis(x[length(x[,1])-3,1])*0.4958904+0.04109589
  SD <- exp(x[length(x[,1]),1])
  L50 <- exp(logL50)
  SR <- exp(logSR)*5
  u <- plogis(x[length(x[,1])-1,1]) 
  w <- plogis(x[length(x[,1])-2,1])
  
  
  
  sel <- function(x){
    1/(1+3^(-2/SR*(x-L50)))
  }
  
  # Function to calculate fish length as a function of t, Age (years) + time of year.
  vbgr <- function(t, k,Linf,t_hatch) {
    age <- t-t_hatch
    X_t <- t-floor(t) # time of year [0,1]
    phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
    phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching
    L_hatch <- 0.4
    
    L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf
    
    return(L_t)
  }
  
  
  cohorts <- sort(unique(ca$cohort))
  
  ret <- 0
  for(i in 1:length(cohorts)){
    cohort.idx <- which(ca$cohort==cohorts[i])
    coef.idx <- which(x_cohor==cohorts[i])
    ages.i <- sort(unique(ca$Age[cohort.idx]+ca$Jday[cohort.idx]))

    Q_ages <- c(mean(ages.i[ages.i<1]),mean(ages.i[ages.i>1]))
    Q_ages <- Q_ages[!is.na(Q_ages)]
    for(j in 1:length(Q_ages)){
      age.idx <- which(ca$Age[cohort.idx]==floor(Q_ages[j]))
      vbgr.fit <- vbgr(Q_ages[j],k[coef.idx],Linf[coef.idx],t_hatch)

      ret <- ret-sum(dnorm(ca$logLength[cohort.idx][age.idx],
                       mean=log(vbgr.fit),sd=SD,log=TRUE)*
                       ca$HLNoAtLngt[cohort.idx][age.idx]*
                       sel(exp(ca$logLength[cohort.idx][age.idx]))/
                       sum(ca$HLNoAtLngt[cohort.idx][age.idx]))
            
      
    }
    ss <- ca %>% filter(cohort %in% cohorts[i] & Age==1)
    r <- round(rnorm(sum(ss$HLNoAtLngt)*1.2,mean = vbgr.fit,sd = 7)) # simulated data
    dr <- data.frame(l = (0:300)/10,
                     Freq=dnorm(log((0:300)/10),mean=log(vbgr.fit),sd=SD)/
                       max(dnorm(log((0:300)/10),mean=log(vbgr.fit),sd=SD))*
                       150) #make simulation a df
    
    ggplot(ss, aes(x = exp(logLength), y = HLNoAtLngt)) + # plot the observed and simulated length-frequencies together
      geom_col(fill = "steelblue") +geom_vline(xintercept =vbgr.fit, linewidth = 2, linetype = "dashed")+
      geom_line(data=dr,aes(x=l,y=Freq),linewidth=2)+
      labs(x = "Length Class", y = "Frequency") +
      theme_minimal()
    #readline() # press enter to continue
  }
  #cohort.idx <- 1997-min(ca$cohort)+1
  #logVbgr1997 <- log(vbgr((0:100)/10,k[cohort.idx],Linf[cohort.idx],t_hatch))
  #ADREPORT(logVbgr1997)
  ret
}
obj <- MakeADFun(f,par,silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr













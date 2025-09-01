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
ca_cod <- ca_cod %>% filter(Age>1 & !is.na(HLNoAtLngt))
hh$Jday <- yday(make_date(year=hh$Year,month=hh$Month,day=hh$Day))/yday(make_date(year=hh$Year,month=12,day=31))
ca_cod <- left_join(ca_cod,hh %>% select(haul.id,Jday))
ca_cod <- ca_cod %>% filter(!is.na(Jday))


cohorts <- sort(unique(ca_cod$cohort))[-c(1:6,39:40)] # cohort with sufficient observations to fit

ca_cod$logLength <- log(ca_cod$LngtCm)

# parameters
par <- list(logK=rep(log(0.1603331),length(cohorts)),
            logLinf=rep(log(112.2677),length(cohorts)),
            #transT_hatch=rep(0.02210022,length(cohorts)),
            transT_hatch = 0,
            logSD=0,logitW=0,logitU=0)

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
    #plot(ca$Age[cohort.idx]+ca$Jday[cohort.idx],ca$logLength[cohort.idx],main=cohorts[i],col="white")
    for(j in 1:length(ages.i)){
      age.idx <- which((ca$Age[cohort.idx]+ca$Jday[cohort.idx])==ages.i[j])
      vbgr.fit <- vbgr(ages.i[j],k,Linf,t_hatch)+err[i]
      
      ret <- RET+dnorm(err[i], 0, sd =SD, log = TRUE)
      #vbgr.sd <- calculate_var()
      #sd.fit <- vbgr(ages.i[j],k[i],Linf[i],t_hatch)
      ret <- ret-sum(dnorm(ca$logLength[cohort.idx][age.idx],mean=log(vbgr.fit),sd=SD_samp,log=TRUE)*ca$HLNoAtLngt[cohort.idx][age.idx])
    #  points(ages.i[j],log(vbgr.fit),cex=3,col="darkred",pch=19)
     # points(ca$Age[cohort.idx][age.idx]+ca$Jday[cohort.idx][age.idx],ca$logLength[cohort.idx][age.idx],cex=1,col="black",pch=1)
      
    }
    #lines((0:100)/10,log(vbgr((0:100)/10,0.1603331,112.2677,0.2917808)),lwd=3,col="darkred")
    #readline() # press enter to continue
  }
  cohort.idx <- 1997-min(ca$cohort)+1
  logVbgr1997 <- log(vbgr((0:100)/10,k[cohort.idx],Linf[cohort.idx],t_hatch))
  ADREPORT(logVbgr1997)
  ret
}
obj <- MakeADFun(f,par,silent=TRUE, random = c('err'))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

SDest <- exp(as.list(sdr,"Est")$logSD)
pl <- as.list(sdr,"Est",report=TRUE)$logVbgr1997
plsd <- as.list(sdr,"Std",report=TRUE)$logVbgr1997



cohortX <- ca_cod %>% filter(cohort==1997)
n_hauls <- length(unique(cohortX$haul.id))
haul_mean <- numeric(n_hauls)
haul_sd <- numeric(n_hauls)
haul_t <- numeric(n_hauls) 
for(i in 1:n_hauls){
  haul.i <- unique(cohortX$haul.id)[i]
  dada <- cohortX %>% filter(haul.id==haul.i)

  haul_mean[i] <- sum(dada$logLength*dada$HLNoAtLngt)/sum(dada$HLNoAtLngt)
  haul_sd[i] <- sqrt(sum(dada$HLNoAtLngt * (dada$logLength - haul_mean[i])^2) / sum(dada$HLNoAtLngt))
  haul_t[i] <-dada$Age[1]+dada$Jday[1]
  if(length(unique(dada$logLength))==1) print(i)
}


plot((0:100)/10,exp(pl) , lwd =2 ,type='l',ylim=c(0,100))
lines((0:100)/10,exp(pl -2 * sqrt(plsd^2+SDest^2)), lty = "dotted" , lwd =1.5  )
lines((0:100)/10,exp(pl +2 * sqrt(plsd^2+SDest^2)), lty = "dotted" , lwd =1.5 )  
lines((0:100)/10,exp(pl -2 * plsd), lty = "dotted" , lwd =1.5  ,col="red")
lines((0:100)/10,exp(pl +2 * plsd), lty = "dotted" , lwd =1.5  ,col="red") 
lines((0:100)/10,exp(pl -2 * exp(-1.71186124)), lty = "dotted" , lwd =1.5 ,col="orange" )
lines((0:100)/10,exp(pl +2 * exp(-1.71186124)), lty = "dotted" , lwd =1.5 ,col="orange") 
points(cohortX$Age+cohortX$Jday,cohortX$LngtCm)


for(i in 1:n_hauls){
  print(i)
  idx <- which(cohortX$haul.id==unique(cohortX$haul.id)[i])
  plot((0:100)/10,exp(pl) , lwd =2 ,type='l',ylim=c(0,100))
  lines((0:100)/10,exp(pl -sqrt(plsd^2+SDest^2)), lty = "dotted" , lwd =1.5  )
  lines((0:100)/10,exp(pl +sqrt(plsd^2+SDest^2)), lty = "dotted" , lwd =1.5 )  
  #lines((0:100)/10,exp(pl -2 * plsd), lty = "dotted" , lwd =1.5  ,col="red")
  #lines((0:100)/10,exp(pl +2 * plsd), lty = "dotted" , lwd =1.5  ,col="red") 
  #lines((0:100)/10,exp(pl -2 * exp(-1.71186124)), lty = "dotted" , lwd =1.5 ,col="orange" )
  #lines((0:100)/10,exp(pl +2 * exp(-1.71186124)), lty = "dotted" , lwd =1.5 ,col="orange") 
  points(cohortX$Age[idx]+cohortX$Jday[idx],cohortX$LngtCm[idx])
  points(haul_t[i],exp(haul_mean[i]),pch=19,cex=2,col="darkred")
  points(haul_t[i],exp(haul_mean[i]+haul_sd[i]),pch=19,cex=1.5,col="orange")
  points(haul_t[i],exp(haul_mean[i]-haul_sd[i]),pch=19,cex=1.5,col="orange")
  
  readline() # press enter to continue
  
}


length(unique(cohortX$haul.id))

###########
# cod
##########
# define growth parameters
L_hatch <- coefs$L_hatch_cod # length at hatching, cm (Pepin et al, 1997)
u <- coefs$u # Amplitude of seasonality (same as cod)
w <- coefs$w # Timing of peak growth
t_hatch <- coefs$t_hatch_cod # time of hatching, Julian day / 365 - spawning from Jan. to May (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching

cod <- df_all %>% filter(species=="cod" & metric=="mean") # cod

before <- 1:3
after <- 37:39
c_cohort <- sort(unique(cod$cohort))[-c(before,after)]
c_coef <- data.frame(k=rep(0,length(c_cohort)),Linf=rep(0,length(c_cohort)),
                     cohort=c_cohort,species="cod")

par <- list(logK = 0,logitTransLinf=0)
t.c <- (0:1000)/100
colors <- colorRampPalette(c("red", "blue"))(length(c_cohort))

for (i in 1:length(c_cohort)){
  dada <- cod %>% filter(cohort ==c_cohort[i] & !is.na(value))
  
  obj_fn <- function(par) {
    k <- exp(par[["logK"]])
    Linf <- plogis(par[["logitTransLinf"]])*150
    pred <- vbgr(dada$age, k,Linf)
    sum((dada$value - pred)^2)
  }
  opt <- nlminb(par, objective = obj_fn)
  c_coef[i,1:2] <- c(exp(opt$par[1]),plogis(opt$par[2])*150)
  
  if(i==1){
    plot(t.c, vbgr(t.c,c_coef[i,1],c_coef[i,2]),type = 'l',lwd=2,col=colors[i],ylim=c(0,120),xlim=c(0,8),
         ylab = "length [cm]",xlab = "Fish age [years]",main="Cod cohorts 1989-2017")
  }
  #main=paste("cod",i+1988)
  if(i>1){
    lines(t.c, vbgr(t.c,c_coef[i,1],c_coef[i,2]),type = 'l',lwd=2,col=colors[i])
  }
  #points(dada$age,dada$value)
  #readline() # press enter to continue
}
coef_before <- data.frame(k=rep(mean(c_coef$k[1:3]),length(before)),
                          Linf=rep(mean(c_coef$Linf[1:3]),length(before)),
                          cohort=min(c_coef$cohort)-rev(before),
                          species=c_coef$species[1])

coef_after <- data.frame(k=rep(mean(c_coef$k[last_index(c_coef$k,3)]),length(after)),
                         Linf=rep(mean(c_coef$Linf[last_index(c_coef$Linf,3)]),length(after)),
                         cohort=max(c_coef$cohort)+1:length(after),
                         species=c_coef$species[1])
c_coef <- rbind(coef_before,c_coef,coef_after)
coefs$cod_varying_growth <- c_coef
saveRDS(coefs,paste(wd_coef,"coefficients.R",sep=""))
#######

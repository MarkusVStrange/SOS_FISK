# cod growth functions
#####

ca_cod <- ca_cod %>% filter(Age<=10) # Remove an outliers

ca_cod$t <- ca_cod$Age+ca_cod$qday
ca_cod$Length <- exp(ca_cod$logLength)

df_m <- data.frame(t_mean = rep(ca_cod$t,times=round(ca_cod$HLNoAtLngt)),
                   L_mean = rep(ca_cod$Length,times=round(ca_cod$HLNoAtLngt)),
                   cohort_mean = rep(ca_cod$cohort,times=round(ca_cod$HLNoAtLngt)),
                   n_mean = 1)
cohortX <- 1988:2021
df_mean <- aggregate(L_mean~t_mean+cohort_mean,data=df_m %>% filter(cohort_mean %in% cohortX),FUN=mean)
df_mean$n_mean <- aggregate(n_mean~t_mean+cohort_mean,data=df_m %>% filter(cohort_mean %in% cohortX),FUN=sum)$n_mean
df_mean$sd_mean <- aggregate(L_mean~t_mean+cohort_mean,data=df_m %>% filter(cohort_mean %in% cohortX),FUN=sd)$L_mean
#dat <- list(ca=rbind(larvae %>% filter(!is.na(logLength)),ca_cod %>% select(haul.id,Age,cohort,logLength,HLNoAtLngt,qday,Quarter)))#%>% filter(cohort %in% 2000:2010)))
dat <- as.list(rbind(ca_cod %>% select(haul.id,Age,cohort,Length,HLNoAtLngt,qday,t)) %>%
                 filter(cohort %in% cohortX))
dat$t_mean <- df_mean$t_mean
dat$cohort_mean <- df_mean$cohort_mean
dat$L_mean <- df_mean$L_mean
dat$n_mean <- df_mean$n_mean
dat$sd_mean <- df_mean$sd_mean

#plot(dat$t,dat$Length)
rm(list=setdiff(ls(),c('dat','ca_cod','cod')))
n_cohorts <- length(unique(dat$cohort[dat$haul.id!="larvae"]))
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
par$logSDage.1 <- 0
par$logSDage.rand <- rep(0,n_cohorts-1)
par$logRWsd <- 0
par$logRWAgesd <- 0
#par$tHatch <- 0.1

f <- function(par){
  getAll(par,dat)
  L1 <- exp(c(logL1.1,logL1.rand))
  L2 <- exp(c(logL2.1,logL2.rand))
  L3 <- exp(c(logL3.1,logL3.rand))
  t_hatch <- 61/365
  
  r <- (L3-L2)/(L2-L1)
  
  t1 <- 1
  t3 <- 6
  #t2 <- (t1+t3)/2
  SD_m <- exp(logSD_m)
  SD_c <- exp(logSD_c)
  SD_age <- exp(c(logSDage.1,logSDage.rand))
  SD_RW <- exp(logRWsd)
  SD_RWage <- exp(logRWAgesd)
  
  VbgrF <- expression(L1+(L3-L1)*(1-((L3-L2)/(L2-L1))^(2*(t-t1)/(t3-t1)))*(1-((L3-L2)/(L2-L1))^2)^(-1))
  diffVBGR <- D(VbgrF,"t")
  
  cohorts <- sort(unique(cohort))
  ret <- 0
  
  for(i in 1:length(cohorts)){
    cohort.idx <- which(cohort==cohorts[i] | haul.id=="larvae")
    times.i <- sort(unique(t[cohort.idx]))
    if(i>1){
      ret <- ret-dnorm(log(L1[i]),mean=log(L1[i-1]),sd=SD_RW,log=TRUE)
      ret <- ret-dnorm(log(L2[i]),mean=log(L2[i-1]),sd=SD_RW,log=TRUE)
      ret <- ret-dnorm(log(L3[i]),mean=log(L3[i-1]),sd=SD_RW,log=TRUE)
      ret <- ret-dnorm(log(SD_age[i]),mean=log(SD_age[i-1]),sd=SD_RWage,log=TRUE)
    }
    
    for(j in 1:length(times.i)){
      age.idx <- which((t[cohort.idx])==times.i[j])
      mean_idx <- which(t_mean==times.i[j] & cohort_mean==cohorts[i])
      
      mean_size <- L_mean[mean_idx]
      
      sd_fac <- ifelse(length(age.idx)>1,1,0)
      
      weight <- HLNoAtLngt[cohort.idx][age.idx]/sqrt(sum((HLNoAtLngt[cohort.idx][age.idx])))*sd_fac
      weight_m <- n_mean[mean_idx]/(sum(n_mean[mean_idx]))
      
      dL_dt <- eval(diffVBGR,list(t=times.i[j],L1=L1[i],L2=L2[i],L3=L3[i],
                                  t1=t1,t3=t3))
      
      
      L.pred <- eval(VbgrF,list(t=times.i[j],L1=L1[i],L2=L2[i],L3=L3[i],
                                t1=t1,t3=t3))
      
      LSD <- sqrt((SD_c*dL_dt)^2+(L.pred*SD_age[i])^2)
      
      ret <- ret-sum(dnorm(Length[cohort.idx][age.idx],
                           mean=mean_size,sd=LSD,
                           log=TRUE)*weight)
      ret <- ret-dnorm(mean_size,mean=L.pred,sd=SD_m*L.pred,log=TRUE)*weight_m # fit mean length
    }
    L.Hatch <- eval(VbgrF,list(t=t_hatch,L1=L1[i],L2=L2[i],L3=L3[i],
                               t1=t1,t3=t3))
    ret <- ret-dnorm(L.Hatch,0.4,0.01,log=TRUE)#*10^4 play with sd here if the model doesn't converge
    
  }
  cohortX <- 2009
  x <- (-100:1000)/100
  c <- which(cohorts==cohortX)
  # evaluate for x
  dL_dt <- eval(diffVBGR,list(t=x,L1=L1[c],L2=L2[c],L3=L3[c],
                              t1=t1,t3=t3))
  fit <- eval(VbgrF,list(t=x,L1=L1[c],L2=L2[c],L3=L3[c],
                         t1=t1,t3=t3))
  tot_sd <- sqrt((SD_c*dL_dt)^2+(fit*SD_age[c])^2)
  ADREPORT(fit)
  ADREPORT(tot_sd)
  #evaluate for data
  times <- t_mean[cohort_mean==cohortX]
  
  dL_dt <- eval(diffVBGR,list(t=times,L1=L1[c],L2=L2[c],L3=L3[c],
                              t1=t1,t3=t3))
  means <- matrix(c(times,L_mean[cohort_mean==cohortX],sd_mean[cohort_mean==cohortX]),ncol=3)
  dat_sd <- sqrt((SD_c*dL_dt)^2+(means[,2]*SD_age[c])^2)
  ADREPORT(means)
  ADREPORT(dat_sd)
  ADREPORT(L1)
  
  ret
}
obj <- MakeADFun(f,par,silent=TRUE,random=c("logL1.rand","logL2.rand","logL3.rand","logSDage.rand"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

M1 <- opt

k <- length(opt$par)
# 2. Log-likelihood is the negative of obj$fn at optimum
logL <- -obj$fn(opt$par)
# 3. Calculate AIC
AIC1 <- 2 * k - 2 * logL
AIC1

fit <- as.list(sdr,"Est",report=TRUE)$fit
tot_sd <- as.list(sdr,"Est",report=TRUE)$tot_sd
dat_means <- as.list(sdr,"Est",report=TRUE)$means
dat_sd <- as.list(sdr,"Est",report=TRUE)$dat_sd
cohortX <- 2009

par(mfrow=c(1,1))
plot((-100:1000)/100,(fit) , lwd =3 ,type='l',
     ylim=c(0,120),col="black",ylab="Length [cm]",xlab="Age [years]",
     main=paste(cohortX,"cohort"))
lines((-100:1000)/100,(fit -tot_sd), lty = "dotted" , lwd =2.5,col="orange")
lines((-100:1000)/100,(fit +tot_sd), lty = "dotted" , lwd =2.5,col="orange")  
points(dat$Age[dat$cohort==cohortX]+dat$qday[dat$cohort==cohortX],dat$Length[dat$cohort==cohortX])
points(dat_means[,1],dat_means[,2],pch=19,cex=2)
points(dat_means[,1],dat_means[,2]-dat_sd,pch=19,cex=1,col="orange")
points(dat_means[,1],dat_means[,2]+dat_sd,pch=19,cex=1,col="orange")
points(dat_means[,1],dat_means[,2]-dat_means[,3],pch=19,cex=1,col="red")
points(dat_means[,1],dat_means[,2]+dat_means[,3],pch=19,cex=1,col="red")
lines(c(-10,100),c(0,0),lty="dashed")

legend(7,40,legend=c("data mean","data SD","modelled SD"),pch=19,col=c("black","red","orange"))


plot((-100:1000)/100,(fit) , lwd =3 ,type='l',xlim=c(-1,2),
     ylim=c(0,40),col="black",ylab="Length [cm]",xlab="Age [years]",
     main=paste(cohortX,"cohort"))
lines((-100:1000)/100,(fit -tot_sd), lty = "dotted" , lwd =2.5,col="orange")
lines((-100:1000)/100,(fit +tot_sd), lty = "dotted" , lwd =2.5,col="orange")  
points(dat$Age+dat$qday,(dat$Length))
points(dat_means[,1],dat_means[,2],pch=19,cex=2)
points(dat_means[,1],dat_means[,2]-dat_sd,pch=19,cex=1,col="orange")
points(dat_means[,1],dat_means[,2]+dat_sd,pch=19,cex=1,col="orange")
points(dat_means[,1],dat_means[,2]-dat_means[,3],pch=19,cex=1,col="red")
points(dat_means[,1],dat_means[,2]+dat_means[,3],pch=19,cex=1,col="red")

legend(1.5,10,legend=c("data mean","data SD","modelled SD"),pch=19,col=c("black","red","orange"))
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
t_sd <- 0#31/365



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


plot(sqrt((dL_dt*(t_sd))^2+(vbgrM*logSD)^2)/SDs*100-100,lwd=1.5)
lines(live/3000,lwd=1.5,col="red")
#####  



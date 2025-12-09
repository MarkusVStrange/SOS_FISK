# model growth functions
#####

ca_cod <- ca_cod %>% filter(Age<=10) # Remove an outliers

ca_cod$t <- ca_cod$Age+ca_cod$qday
ca_cod$Length <- exp(ca_cod$logLength)

df_m <- data.frame(t_mean = rep(ca_cod$t,times=round(ca_cod$HLNoAtLngt)),
                   L_mean = rep(ca_cod$Length,times=round(ca_cod$HLNoAtLngt)),
                   cohort_mean = rep(ca_cod$cohort,times=round(ca_cod$HLNoAtLngt)),
                   n_mean = 1)
cohortX <- 2013
df_mean <- aggregate(L_mean~t_mean+cohort_mean,data=df_m %>% filter(cohort_mean==cohortX),FUN=mean)
df_mean$n_mean <- aggregate(n_mean~t_mean+cohort_mean,data=df_m %>% filter(cohort_mean==cohortX),FUN=sum)$n_mean
df_mean$sd_mean <- aggregate(L_mean~t_mean+cohort_mean,data=df_m %>% filter(cohort_mean==cohortX),FUN=sd)$L_mean
#dat <- list(ca=rbind(larvae %>% filter(!is.na(logLength)),ca_cod %>% select(haul.id,Age,cohort,logLength,HLNoAtLngt,qday,Quarter)))#%>% filter(cohort %in% 2000:2010)))
dat <- as.list(rbind(ca_cod %>% select(haul.id,Age,cohort,Length,HLNoAtLngt,qday,t)) %>%
                 filter(cohort==cohortX))
dat$t_mean <- df_mean$t_mean
dat$cohort_mean <- df_mean$cohort_mean
dat$L_mean <- df_mean$L_mean
dat$n_mean <- df_mean$n_mean
dat$sd_mean <- df_mean$sd_mean

plot(dat$t,dat$Length)


rm(list=setdiff(ls(),c('dat','ca_cod','larvae','cod')))
n_cohorts <- length(unique(dat$cohort[dat$haul.id!="larvae"]))
#%>% filter(Age>0))
#rm(list=setdiff(ls(),c('par','dat')))
par <- list()
par$logK.1 <- log(0.3)
par$logLinf.1 <- log(70)
#par$t0.1 <- 0
#par$logK.rand <- rep(log(0.3),n_cohorts-1) # not including 1st year
#par$logLinf.rand <- rep(log(70),n_cohorts-1)

par$logSD_m <- -1 # deviation from mean fit
#par$logSD_k <- -5 # k SD
par$logSD_Linf <- 2 # k SD
par$logSD_t0 <- -2 # Age SD
#par$tHatch <- 0.1 # time of hatching
#logRWSD <- -1 # Random walk SD
#par$U <- 0
#par$logitW <- 0

f <- function(par){
  getAll(par,dat)
  k <- c(exp(logK.1))#,exp(logK.rand))
  Linf <- c(exp(logLinf.1))#,exp(logLinf.rand))
  #t0 <- c(t0.1)
  
  SD_m <- exp(logSD_m)
  SD_k <- 0#exp(logSD_k)
  SD_Linf <- exp(logSD_Linf)
  SD_t0 <- exp(logSD_t0)
  u <- 0.23 #plogis(logitU)#0.23 # Amplitude of seasonality
  w <- 0.91 #plogis(logitW) Timing of peak growth
  t0 <- 60/365 #((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to May (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
  
  VbgrF <- expression(Linf*(1-exp(-k*(t-t0+phi_t))))
  diffVBGR_k <- D(VbgrF,"k")
  diffVBGR_Linf <- D(VbgrF,"Linf")
  diffVBGR_age <- D(VbgrF,"t0")
  
  cohorts <- sort(unique(cohort[haul.id!="larvae"]))
  
  ret <- 0
  for(i in 1:length(cohorts)){
    cohort.idx <- which(cohort==cohorts[i] | haul.id=="larvae")
    times.i <- sort(unique(t[cohort.idx]))
    
    #if(i>1){
    #  ret <- ret-dnorm(log(L1[i]),mean=log(L1[i-1]),sd=RWSD,log=TRUE)
    #  ret <- ret-dnorm(log(L2[i]),mean=log(L2[i-1]),sd=RWSD,log=TRUE)
    #  ret <- ret-dnorm(log(L3[i]),mean=log(L3[i-1]),sd=RWSD,log=TRUE)
    #}
    
    for(j in 1:length(times.i)){
      age.idx <- which((t[cohort.idx])==times.i[j])
      mean_idx <- which(t_mean==times.i[j] & cohort_mean==cohorts[i])
      #if(length(age.idx)<2) next
      age.j <- times.i[j]-t0
      phi_t.j <- u*sin(2*pi*(times.i[j]-w))/(2*pi) # seasonal variability in growth
      mean_size <- L_mean[mean_idx]
      
      sd_fac <- ifelse(length(age.idx)>1,1,0)
      
      weight <- HLNoAtLngt[cohort.idx][age.idx]/sqrt(sum((HLNoAtLngt[cohort.idx][age.idx])))*sd_fac
      weight_m <- n_mean[mean_idx]/sqrt(sum(n_mean[mean_idx]))
      
      dL_dk <- eval(diffVBGR_k,list(t=times.i[j],t0=t0[i],k=k[i],Linf=Linf[i],
                                    phi_t =phi_t.j))
      dL_dLinf <- eval(diffVBGR_Linf,list(t=times.i[j],t0=t0[i],k=k[i],Linf=Linf[i],
                                          phi_t =phi_t.j))
      dL_dt0 <- eval(diffVBGR_age,list(t=times.i[j],t0=t0[i],k=k[i],Linf=Linf[i],
                                        phi_t =phi_t.j))
      
      L.pred <- eval(VbgrF,list(t=times.i[j],t0=t0[i],k=k[i],Linf=Linf[i],
                                phi_t =phi_t.j))
      
      LSD <- sqrt((SD_k*dL_dk)^2+(SD_Linf*dL_dLinf)^2+(SD_t0*dL_dt0)^2)
      
      ret <- ret-sum(dnorm(Length[cohort.idx][age.idx],
                           mean=mean_size,sd=LSD,
                           log=TRUE)*weight)
      ret <- ret-dnorm(mean_size,mean=L.pred,sd=SD_m,log=TRUE)*weight_m # fit mean length
    }
  }
  cohortX <- 2013
  x <- (-100:1000)/100
  c <- which(cohorts==cohortX)
  phi_t <- u*sin(2*pi*(x-w))/(2*pi) # seasonal variability in growth
  # evaluate for x
  dL_dk <- eval(diffVBGR_k,list(t=x,t0=t0[c],k=k[c],Linf=Linf[c],phi_t =phi_t))
  dL_dLinf <- eval(diffVBGR_Linf,list(t=x,t0=t0[c],k=k[c],Linf=Linf[c],phi_t =phi_t))
  dL_dt0 <- eval(diffVBGR_age,list(t=x,t0=t0[c],k=k[c],Linf=Linf[c],phi_t =phi_t))
  fit <- (eval(VbgrF,list(t=x,t0=t0[c],k=k[c],Linf=Linf[c],phi_t =phi_t)))
  tot_sd <- sqrt((dL_dLinf*SD_Linf)^2+(dL_dt0*SD_t0)^2+(dL_dk*SD_k))
  ADREPORT(fit)
  ADREPORT(tot_sd)
  #evaluate for data
  times <- t_mean[cohort_mean==cohortX]
  phi_t <- u*sin(2*pi*(times-w))/(2*pi) # seasonal variability in growth
  
  dL_dk <- eval(diffVBGR_k,list(t=times,t0=t0[c],k=k[c],Linf=Linf[c],phi_t =phi_t))
  dL_dLinf <- eval(diffVBGR_Linf,list(t=times,t0=t0[c],k=k[c],Linf=Linf[c],phi_t =phi_t))
  dL_dt0 <- eval(diffVBGR_age,list(t=times,t0=t0[c],k=k[c],Linf=Linf[c],phi_t =phi_t))
  means <- matrix(c(times,L_mean[cohort_mean==cohortX],sd_mean[cohort_mean==cohortX]),ncol=3)
  dat_sd <- sqrt((dL_dLinf*SD_Linf)^2+(dL_dt0*SD_t0)^2+(dL_dk*SD_k))
  ADREPORT(means)
  ADREPORT(dat_sd)
  
  ret
}
obj <- MakeADFun(f,par,silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr


fit <- as.list(sdr,"Est",report=TRUE)$fit
tot_sd <- as.list(sdr,"Est",report=TRUE)$tot_sd
dat_means <- as.list(sdr,"Est",report=TRUE)$means
dat_sd <- as.list(sdr,"Est",report=TRUE)$dat_sd

par(mfrow=c(1,1))
plot((-100:1000)/100,(fit) , lwd =3 ,type='l',
     ylim=c(0,120),col="black",ylab="Length [cm]",xlab="Age [years]",
     main="2002 cohort")
lines((-100:1000)/100,(fit -tot_sd), lty = "dotted" , lwd =2.5,col="orange")
lines((-100:1000)/100,(fit +tot_sd), lty = "dotted" , lwd =2.5,col="orange")  
points(dat$Age+dat$qday,(dat$Length))
points(dat_means[,1],dat_means[,2],pch=19,cex=2)
points(dat_means[,1],dat_means[,2]-dat_sd,pch=19,cex=1,col="orange")
points(dat_means[,1],dat_means[,2]+dat_sd,pch=19,cex=1,col="orange")
points(dat_means[,1],dat_means[,2]-dat_means[,3],pch=19,cex=1,col="red")
points(dat_means[,1],dat_means[,2]+dat_means[,3],pch=19,cex=1,col="red")

legend(7,40,legend=c("data mean","data SD","modelled SD"),pch=19,col=c("black","red","orange"))


plot((-100:1000)/100,(fit) , lwd =3 ,type='l',xlim=c(-1,2),
     ylim=c(0,40),col="black",ylab="Length [cm]",xlab="Age [years]",
     main="2002 cohort")
lines((-100:1000)/100,(fit -tot_sd), lty = "dotted" , lwd =2.5,col="orange")
lines((-100:1000)/100,(fit +tot_sd), lty = "dotted" , lwd =2.5,col="orange")  
points(dat$Age+dat$qday,(dat$Length))
points(dat_means[,1],dat_means[,2],pch=19,cex=2)
points(dat_means[,1],dat_means[,2]-dat_sd,pch=19,cex=1,col="orange")
points(dat_means[,1],dat_means[,2]+dat_sd,pch=19,cex=1,col="orange")
points(dat_means[,1],dat_means[,2]-dat_means[,3],pch=19,cex=1,col="red")
points(dat_means[,1],dat_means[,2]+dat_means[,3],pch=19,cex=1,col="red")

legend(1.5,10,legend=c("data mean","data SD","modelled SD"),pch=19,col=c("black","red","orange"))
lines(c(1,1),c(-10,100),lty="dashed")
fit[197]


########################
#Conseptual simulations
#######################
# cohort spread
#####
vbgr <- function(t,k=k.i,Linf=Linf.i,t0=t0.i) {
  # exponent terms
  u <- 0.23 # Amplitude of seasonality
  w <- 0.91 # Timing of peak growth
  #t0 <-  # time of hatching, Julian day / 365 - spawning from Jan. to May (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
  phi_t <- u*sin(2*pi*(t-w))/(2*pi) # seasonal variability in growth
  

  L <- Linf*(1-exp(-k*(t-t0+phi_t)))

  return(L)
}
k.i <- exp(-1.342149)
Linf.i <- exp(4.427219)
t0.i <- ((16+197)/2)/365
k_sd <- 0.0
Linf_sd <- exp(2.429177 )
t0_sd <- exp(-1.412180 )

x <- (0:(3650*2))/365
L <- matrix(rep(NA,length(x)*10000),ncol=10000)
par(mfrow=c(1,1))
vbgrM <- vbgr(x)
plot(x,(vbgrM),
     type='l',lwd=2,ylab="length",xlab="age",ylim=c(0,120))
lines(c(-10,20),c(0,0),lty="dashed")

for(i in 1:10000){
  t0.id <- ((16+197)/2)/365+rnorm(1,mean=0,sd=t0_sd)
  k.id <- rnorm(1,mean=k.i,sd=k_sd)
  Linf.id <- rnorm(1,mean=Linf.i,sd=Linf_sd)
  L[,i] <- (vbgr(x,Linf=Linf.id,k=k.id,t0=t0.id))
  #lines(x,(L[,i]),type='l',lwd=2,ylab="cm",xlab="age")
}
#m <- apply((L), 1, mean,na.rm = TRUE)
#L[which(L<0)] <- NA
SDs <- apply((L), 1, sd,na.rm = TRUE)
par(mfrow=c(1,1))

plot(x,SDs,type='l',lwd=2)
VbgrF <- expression(Linf*(1-exp(-k*(t-t0+u*sin(2*pi*(t-w))/(2*pi)))))
diffVBGR.t0 <- D(VbgrF,"t0")
diffVBGR.k <- D(VbgrF,"k")
diffVBGR.Linf <- D(VbgrF,"Linf")

dL_dt0 <- eval(diffVBGR.t0,list(t=x,t0=t0.i,k=k.i,Linf=Linf.i,
                                  u =u, w=w,X_t=x))
dL_dk <- eval(diffVBGR.k,list(t=x,t0=t0.i,k=k.i,Linf=Linf.i,
                              u =u, w=w,X_t=x))
dL_dLinf <- eval(diffVBGR.Linf,list(t=x,t0=t0.i,k=k.i,Linf=Linf.i,
                                    u =u, w=w,X_t=x))
u <- 0.23 # Amplitude of seasonality
w <- 0.91 # Timing of peak growth
lines(x,sqrt((dL_dt0*t0_sd)^2+(dL_dk*k_sd)^2+(dL_dLinf*Linf_sd)^2),lwd=2,col="red")
lines(c(0,0),c(-10,10),lty="dashed")
lines(c(1,1),c(-10,10),lty="dashed")


plot(x,sqrt((dL_dt0*t0_sd)^2+(dL_dk*k_sd)^2+(dL_dLinf*Linf_sd)^2)/SDs*100-100,type='l',lwd=2,
     xlab="mean age", ylab="relative estimation error [%]")
##### 

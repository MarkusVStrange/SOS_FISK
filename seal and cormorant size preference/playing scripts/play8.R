# model growth functions
#####

ca_cod <- ca_cod %>% filter(Age<=10) # Remove an outliers

ca_cod$t <- ca_cod$Age+ca_cod$qday
ca_cod$Length <- exp(ca_cod$logLength)

df_m <- data.frame(t_mean = rep(ca_cod$t,times=round(ca_cod$HLNoAtLngt)),
                      L_mean = rep(ca_cod$Length,times=round(ca_cod$HLNoAtLngt)),
                      cohort_mean = rep(ca_cod$cohort,times=round(ca_cod$HLNoAtLngt)),
                      n_mean = 1)

df_mean <- aggregate(L_mean~t_mean+cohort_mean,data=df_m %>% filter(cohort_mean==2002),FUN=mean)
df_mean$n_mean <- aggregate(n_mean~t_mean+cohort_mean,data=df_m %>% filter(cohort_mean==2002),FUN=sum)$n_mean
#dat <- list(ca=rbind(larvae %>% filter(!is.na(logLength)),ca_cod %>% select(haul.id,Age,cohort,logLength,HLNoAtLngt,qday,Quarter)))#%>% filter(cohort %in% 2000:2010)))
dat <- as.list(rbind(ca_cod %>% select(haul.id,Age,cohort,Length,HLNoAtLngt,qday,t)) %>%
              filter(cohort==2002))
dat$t_mean <- df_mean$t_mean
dat$cohort_mean <- df_mean$cohort_mean
dat$L_mean <- df_mean$L_mean
dat$n_mean <- df_mean$n_mean
plot(dat$t,dat$Length)


rm(list=setdiff(ls(),c('dat','ca_cod','larvae','cod')))
n_cohorts <- length(unique(dat$cohort[dat$haul.id!="larvae"]))
#%>% filter(Age>0))
#rm(list=setdiff(ls(),c('par','dat')))
par <- list()
par$logK.1 <- log(0.3)
par$logLinf.1 <- log(70)
#par$logK.rand <- rep(log(0.3),n_cohorts-1) # not including 1st year
#par$logLinf.rand <- rep(log(70),n_cohorts-1)

par$logSD_m <- -1 # deviation from mean fit
#par$logSD_k <- -5 # k SD
par$logSD_Linf <- 2 # k SD
par$logSDage <- -2 # Age SD
#par$tHatch <- 0.1 # time of hatching
#logRWSD <- -1 # Random walk SD
#par$U <- 0
par$logitW <- 0

f <- function(par){
  getAll(par,dat)
  k <- c(exp(logK.1))#,exp(logK.rand))
  Linf <- c(exp(logLinf.1))#,exp(logLinf.rand))
  
  SD_m <- exp(logSD_m)
  SD_k <- 0#exp(logSD_k)
  SD_Linf <- exp(logSD_Linf)
  SDage <- exp(logSDage)
  u <- 0.23 #plogis(logitU)#0.23 # Amplitude of seasonality
  w <- plogis(logitW)#0.91 # Timing of peak growth
  t_hatch <- 60/365 #((16+197)/2)/365 # time of hatching, Julian day / 365 - spawning from Jan. to May (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning

  VbgrF <- expression((-Linf)*exp(-k*(phi_t+age))+Linf)
  diffVBGR_k <- D(VbgrF,"k")
  diffVBGR_Linf <- D(VbgrF,"Linf")
  diffVBGR_age <- D(VbgrF,"age")
  
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
      age.j <- times.i[j]-t_hatch
      phi_t.j <- u*sin(2*pi*(times.i[j]-w))/(2*pi) # seasonal variability in growth
      mean_size <- L_mean[mean_idx]
      
      sd_fac <- ifelse(length(age.idx)>1,1,0)
      
      weight <- HLNoAtLngt[cohort.idx][age.idx]/sqrt(sum((HLNoAtLngt[cohort.idx][age.idx])))*sd_fac
      weight_m <- n_mean[mean_idx]/(sum(n_mean[mean_idx]))
      
      dL_dk <- eval(diffVBGR_k,list(age=age.j,k=k[i],Linf=Linf[i],
                                          phi_t =phi_t.j))
      dL_dLinf <- eval(diffVBGR_Linf,list(age=age.j,k=k[i],Linf=Linf[i],
                                          phi_t =phi_t.j))
      dL_dage <- eval(diffVBGR_age,list(age=age.j,k=k[i],Linf=Linf[i],
                                       phi_t =phi_t.j))
      
      L.pred <- eval(VbgrF,list(age=age.j,k=k[i],Linf=Linf[i],
                                        phi_t =phi_t.j))
        
      LSD <- sqrt((SD_k*dL_dk)^2+(SD_Linf*dL_dLinf)^2+(SDage*dL_dage)^2)
        
      ret <- ret-sum(dnorm(Length[cohort.idx][age.idx],
                             mean=mean_size,sd=LSD,
                           log=TRUE)*weight)
      ret <- ret-dnorm(mean_size,mean=L.pred,sd=SD_m,log=TRUE)*weight_m # fit mean length
    }
  }
  x <- (0:1000)/100
  c <- which(cohorts==2002)
  phi_t.j <- u*sin(2*pi*(x-w))/(2*pi) # seasonal variability in growth

  logfit <- (eval(VbgrF,list(age=x-t_hatch,k=k[i],Linf=Linf[i],
                                phi_t =phi_t.j)))
  ADREPORT(logfit)
  ret
}
obj <- MakeADFun(f,par,silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

pl <- as.list(sdr,"Est",report=TRUE)$logfit
coefs <- exp(summary(sdr)[1:7,1])
# define gowth- and differential equations to calculate total SD
logVbgrF <- expression(((-Linf)*exp(-k*(u*sin(2*pi*(X_t-w))/(2*pi)+age))+Linf))
diffVBGR_age <- D(logVbgrF,"age")
diffVBGR_Linf <- D(logVbgrF,"Linf")
u <- 0.23#plogis(log(coefs[6])) # Amplitude of seasonality
w <- 0.91#plogis(log(coefs[7]))  # Timing of peak growth
t_hatch=60/365#((16+197)/2)/365
dL_dage <- eval(diffVBGR_age,list(age=(0:1000)/100-t_hatch,k=coefs[1],Linf=coefs[2],
                            u =u, w=w,X_t=(0:1000)/100))
dL_dLinf <- eval(diffVBGR_Linf,list(age=(0:1000)/100-t_hatch,k=coefs[1],Linf=coefs[2],
                                  u =u, w=w,X_t=(0:1000)/100))

totSD <- sqrt((coefs[4]*dL_dLinf)^2+(coefs[5]*dL_dage)^2)
df <- dat
df <- data.frame(t=rep(df$Age+df$qday,times=df$HLNoAtLngt),
                 Length=rep(df$Length,times=df$HLNoAtLngt))
d <- aggregate(Length~t,data = df,FUN=mean)
d$sd <- aggregate(Length~t,data = df,FUN=sd)$Length


dL_dage.data <- eval(diffVBGR_age,list(age=d$t-t_hatch,k=coefs[1],Linf=coefs[2],
                                 u =u, w=w,X_t=d$t))
dL_dLinf.data <- eval(diffVBGR_Linf,list(age=d$t-t_hatch,k=coefs[1],Linf=coefs[2],
                                        u =u, w=w,X_t=d$t))

totSD.data <- sqrt((coefs[4]*dL_dLinf.data)^2+(coefs[5]*dL_dage.data)^2)

par(mfrow=c(1,1))
plot((0:1000)/100,(pl) , lwd =3 ,type='l',
     ylim=c(0,120),col="black",ylab="Length [cm]",xlab="Age [years]",
     main="2000 cohort")
lines((0:1000)/100,(pl -totSD), lty = "dotted" , lwd =2.5,col="orange")
lines((0:1000)/100,(pl +totSD), lty = "dotted" , lwd =2.5,col="orange")  
points(dat$Age+dat$qday,(dat$Length))
points(d$t,(d$Length),pch=19,cex=2)
points(d$t,(d$Length-totSD.data),pch=19,cex=1,col="orange")
points(d$t,(d$Length+totSD.data),pch=19,cex=1,col="orange")
points(d$t,(d$Length-d$sd),pch=19,cex=1,col="red")
points(d$t,(d$Length+d$sd),pch=19,cex=1,col="red")

legend(7,40,legend=c("data mean","data SD","modelled SD"),pch=19,col=c("black","red","orange"))


par(mfrow=c(1,1))
plot((0:1000)/100,(pl) , lwd =3 ,type='l',
     ylim=c(0,40),xlim=c(0,2),col="black",ylab="Length [cm]",xlab="Age [years]",
     main="2000 cohort")
lines((0:1000)/100,(pl -totSD), lty = "dotted" , lwd =2.5,col="orange")
lines((0:1000)/100,(pl +totSD), lty = "dotted" , lwd =2.5,col="orange")  
points(dat$Age+dat$qday,(dat$Length))
points(d$t,(d$Length),pch=19,cex=2)
points(d$t,(d$Length-totSD.data),pch=19,cex=1,col="orange")
points(d$t,(d$Length+totSD.data),pch=19,cex=1,col="orange")
points(d$t,(d$Length-d$sd),pch=19,cex=1,col="red")
points(d$t,(d$Length+d$sd),pch=19,cex=1,col="red")

legend(1.5,10,legend=c("data mean","data SD","modelled SD"),pch=19,col=c("black","red","orange"))
lines(c(1,1),c(-10,100),lty="dashed")
pl[100]


########################
#Conseptual simulations
#######################
# cohort spread
#####
vbgr <- function(t,k=k.i,Linf=Linf.i,t_hatch=((16+197)/2)/365) {
  # exponent terms
  u <- 0.23 # Amplitude of seasonality
  w <- 0.91 # Timing of peak growth
  #t_hatch <-  # time of hatching, Julian day / 365 - spawning from Jan. to May (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
  phi_t <- u*sin(2*pi*(t-w))/(2*pi) # seasonal variability in growth
  
  age <- t-t_hatch
  
  L <- (-Linf)*exp(-k*(phi_t+age))+Linf
  
  return(L)
}
k.i <- exp(-1.342149)
Linf.i <- exp(4.427219)
k_sd <- 0.0
Linf_sd <- exp(2.429177 )
t_hatch_sd <- exp(-1.412180 )

x <- (0:(3650*2))/365
L <- matrix(rep(NA,length(x)*10000),ncol=10000)
par(mfrow=c(1,1))
vbgrM <- vbgr(x)
plot(x,(vbgrM),
     type='l',lwd=2,ylab="length",xlab="age",ylim=c(0,120))
lines(c(-10,20),c(0,0),lty="dashed")

for(i in 1:10000){
  t_hatch.id <- ((16+197)/2)/365+rnorm(1,mean=0,sd=t_hatch_sd)
  k.id <- rnorm(1,mean=k.i,sd=k_sd)
  Linf.id <- rnorm(1,mean=Linf.i,sd=Linf_sd)
  L[,i] <- (vbgr(x,Linf=Linf.id,k=k.id,t_hatch=t_hatch.id))
  #lines(x,(L[,i]),type='l',lwd=2,ylab="cm",xlab="age")
}
#m <- apply((L), 1, mean,na.rm = TRUE)
#L[which(L<0)] <- NA
SDs <- apply((L), 1, sd,na.rm = TRUE)
par(mfrow=c(1,1))

plot(x,SDs,type='l',lwd=2)
VbgrF <- expression(((-Linf)*exp(-k*(u*sin(2*pi*(X_t-w))/(2*pi)+age))+Linf))
diffVBGR.age <- D(VbgrF,"age")
diffVBGR.k <- D(VbgrF,"k")
diffVBGR.Linf <- D(VbgrF,"Linf")

dL_dage <- eval(diffVBGR.age,list(age=x-t_hatch,k=k.i,Linf=Linf.i,
                                  u =u, w=w,X_t=x))
dL_dk <- eval(diffVBGR.k,list(age=x-t_hatch,k=k.i,Linf=Linf.i,
                              u =u, w=w,X_t=x))
dL_dLinf <- eval(diffVBGR.Linf,list(age=x-t_hatch,k=k.i,Linf=Linf.i,
                                    u =u, w=w,X_t=x))
u <- 0.23 # Amplitude of seasonality
w <- 0.91 # Timing of peak growth
t_hatch=((16+197)/2)/365
lines(x,sqrt((dL_dage*t_hatch_sd)^2+(dL_dk*k_sd)^2+(dL_dLinf*Linf_sd)^2),lwd=2,col="red")
lines(c(0,0),c(-10,10),lty="dashed")
lines(c(1,1),c(-10,10),lty="dashed")


plot(x-t_hatch,sqrt((dL_dage*t_hatch_sd)^2+(dL_dk*k_sd)^2+(dL_dLinf*Linf_sd)^2)/SDs*100-100,type='l',lwd=2,
     xlab="mean age", ylab="relative estimation error [%]")
##### 

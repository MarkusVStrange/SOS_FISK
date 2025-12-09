############
# read data - read in an sdr-object (mainVBGR.R)
############
d <- rbind(readRDS("N_Cod1.RData"),readRDS("N_Cod4.RData"))
ca_cod <- d %>% filter(N!=0)
ca_cod$cohort <- ca_cod$Year-ca_cod$Age
ca_cod$Age[ca_cod$Age>7] <- 7
ca_cod$t <- ca_cod$Age+ca_cod$qday



names(ca_cod)[names(ca_cod)=="LngtCm"] <- "Length"

df <- ca_cod %>% select(Length,t,cohort,Quarter,N,Year,Species,Age)

d1 <- aggregate(N~t+cohort+Quarter+Age+Year+Species,data=df %>% filter(t>0.5 & Age<7),FUN=sum)
d2 <- aggregate(N~t+Quarter+Age+Year+Species,data=df %>% filter(t>0.5 & Age==7),FUN=sum)
d2$cohort <- d2$Year-d2$Age
d <- rbind(d1,d2)
# some plots
par(mfrow=c(6,6))
for(i in 1:36){
  dada <- d %>% filter(Quarter==4 & cohort==1987+i)
  #plot(dada$t,log(dada$N),pch=19,ylab="log(N)",xlab="Age",
  #     main=paste("Q",dada$Quarter[1],"cohort",unique(dada$cohort)))
  
}


for(i in 1:36){
  dada <- d %>% filter(Quarter==1 & cohort==1987+i)
  #plot(dada$t,log(dada$N),pch=19,ylab="log(N)",xlab="Age",
  #     main=paste("Q",dada$Quarter[1],"cohort",unique(dada$cohort)))
}

for(i in 1:36){
  dada <- d %>% filter(Quarter==4 & cohort==1987+i)
  plot(dada$t,log(dada$N),pch=19,ylab="log(N)",xlab="Age",
       main=paste("cohort",unique(dada$cohort)))
  dada <- d %>% filter(Quarter==1 & cohort==1987+i)
  points(dada$t,log(dada$N),pch=19,col="red")
}


#cohortX <- 1988:2022
#coh.idx <- which(cohortX %in% 1988:2022)
dat <- as.list(d %>% filter(t>0.5))

# assumption: everything above age 1.5 is fully selected

rm(list=setdiff(ls(),c('dat','ca_cod','sdr','coh.idx')))

df_L <- data.frame(L1 = as.list(sdr,"Est",report=TRUE)$L1,
                   L2 = as.list(sdr,"Est",report=TRUE)$L2,
                   L3 = as.list(sdr,"Est",report=TRUE)$L3,
                   cohort=1988:2022)

L.i <- rep(0,length(df_L$L1))
VbgrF <- expression(L1+(L3-L1)*(1-((L3-L2)/(L2-L1))^(2*(age-t1)/(t3-t1)))*(1-((L3-L2)/(L2-L1))^2)^(-1))

for(i in 1:length(L.i)){
  L.i[i] <- eval(VbgrF,list(age=8-61/365,L1=df_L$L1[i],L2=df_L$L2[i],L3=df_L$L3[i],
                            t1=1,t3=5))
}

LW <- read.table("C:/Users/mavast/Documents/GitHub/SOS data/length-weight.csv",header=TRUE,sep=';')

dat$Lmax <- mean(L.i)
dat$LW <- LW %>% filter(species==unique(dat$Species))
species <- unique(dat$Species)
dat$t_hatches <- c(cod=61/365,flounder=((75+197)/2)/365,plaice=0.03972603,dab=0.4561644)
dat$L_hatches <-c(cod=0.4,flounder=mean(c(0.697,0.662)),plaice=mean(c(0.697,0.662)),dab=mean(c(0.697,0.662)))
dat$t1s <- c(cod=1,flounder=2,plaice=2,dab=2)
dat$t3s <- c(cod=5,flounder=10,plaice=10,dab=10)
dat$Ls <- df_L
dat$obsID <- paste(dat$Year,dat$Quarter,dat$Age)
dat$mo <- c(0,0.06,0.6,0.84,0.86,0.9,0.94,1)
############

############
# parameters and likelihood function
############
par <- list()
par$logQdiff <- 0
par$logQ1_4 <- c(-2,-1,-1) # catchability for Age 0 and 1. Age 2+ log-catchability is 0
par$logM.rand <- rep(-5,length(unique(dat$Year)))
par$logF.rand <- rep(-4,length(unique(dat$Year)))
par$logN <- matrix(5, nrow=length(0:max(dat$Age)), ncol=length(unique(dat$Year))+1)

par$logSD_F <- 0 
par$logSD_M <- 0
par$logSD_N <- c(0)
par$logSDobs <- 0
par$logSD_R <- 0

f <- function(par){
  getAll(par,dat)
  L1 <- Ls$L1
  L2 <- Ls$L2
  L3 <- Ls$L3
  t_hatch <- as.numeric(t_hatches[unique(Species)])
  L_hatch <- as.numeric(L_hatches[unique(Species)])
  N <- OBS(N)
  
  r <- (L3-L2)/(L2-L1)
  
  t1 <- as.numeric(t1s[unique(Species)])
  t3 <- as.numeric(t3s[unique(Species)])
  #t2 <- (t1+t3)/2
  #M_slope <- -exp(rep(logM.1,length(Year)))
  M_slope <- -plogis(logM.rand)
  #F_scale <- exp(rep(logF.1,length(Year)))
  F_scale <- exp(logF.rand)
  SD_M <- exp(logSD_M)
  SD_F <- exp(logSD_F)
  SD_N <- exp(logSD_N)
  SD_obs <- exp(logSDobs)
  SD_R <- exp(logSD_R)
  
  #qDiff <- exp(logQdiff) # use Q1 as reference
  
  VbgrF <- expression(L1+(L3-L1)*(1-((L3-L2)/(L2-L1))^(2*(age-t1)/(t3-t1)))*(1-((L3-L2)/(L2-L1))^2)^(-1))
  
  cohorts <- sort(unique(Ls$cohort))
  
  ret <- 0
  
  na <- max(Age)-(min(Age))+1
  ny <- max(Year)-min(Year)+1+1
  
  
  # N
  for(y in 2:ny){
    ret <- ret + -dnorm(logN[1,y], logN[1,y-1], SD_R, log=TRUE)
    if(y>2){
      ret <- ret-dnorm(logM.rand[y-1],mean=logM.rand[y-2],sd=SD_M,log=TRUE)
      ret <- ret-dnorm(logF.rand[y-1],mean=logF.rand[y-2],sd=SD_F,log=TRUE)
    }
    
    for(a in 2:na){ # for ages 1-7+
      age <- a-1
      if((sort(unique(Year))[y]-age) %in% cohorts){
        coh.idx <- which(cohorts %in% (sort(unique(Year))[y]-age))
      }
      
      days <- (1:365)/365
      
      L <- ifelse(days<t_hatch,0,eval(VbgrF,list(age=days+age-1-t_hatch,L1=L1[coh.idx],L2=L2[coh.idx],L3=L3[coh.idx],
                                                 t1=t1,t3=t3)))
      
      M.i <-  ifelse(L>Lmax | L<L_hatch,0,M_slope[y-1]*L-M_slope[y-1]*Lmax)
      
      Z <- sum((M.i+L*F_scale[y-1])/365)
      
      if(a<na){
        pred <- logN[a-1,y-1]-Z    
        ret <- ret + -dnorm(logN[a,y], pred, SD_N, log=TRUE)  
      }else{
        pred <- log(exp(logN[a-1,y-1]-Z)+exp(logN[a-1,y-1]-Z))
        ret <- ret + -dnorm(logN[a,y], pred, SD_N, log=TRUE) 
        
      }
    }
  }
  
  
  # Obs 
  
  logPred <- numeric(length(obsID))
  
  for(i in 1:length(obsID)){
    t.i <- t[i]
    y.idx <- Year[i]-min(Year)+1
    
    if((Year[i]-Age[i]) %in% cohorts){
      coh.idx <- which(cohorts %in% (Year[i]-Age[i]))
    }
    
    days <- seq((1/365),t.i-Age[i],length=round((t.i-Age[i])*365))
    
    
    L <- eval(VbgrF,list(age=days+Age[i]-t_hatch,L1=L1[coh.idx],L2=L2[coh.idx],L3=L3[coh.idx],
                         t1=t1,t3=t3))
    
    M.i <-  ifelse(L>Lmax | L<L_hatch,0,M_slope[y.idx]*L-M_slope[y.idx]*Lmax)
    
    Z <- sum((M.i+L*F_scale[y.idx])/length(days))
    
    logQ <- ifelse(Quarter[i]==1,0,logQdiff)
    
    if(t.i<2){
      idx <- which(c("0 4","1 1","1 4","2 1") %in% paste(Age[i],Quarter[i]))
      logQ <- logQ+logQ1_4[idx]
    }
    #    logQ <- ifelse(t.i>1.5,logQ,logQ+logQ1_4[Age[i]+1])
    logPred[i] <- logQ-Z+logN[Age[i]+1,y.idx]
  }
  ret <- ret-sum(dnorm(log(N),mean=logPred,sd=SD_obs,log=TRUE))
  
  
  # calculate TSB and M'ish and F'ish
  fish_weight <- matrix(NA,ncol=na,nrow=ny-1) # not including 2026
  fish_length <- matrix(NA,ncol=na,nrow=ny-1)
  fish_weight_J1 <- matrix(NA,ncol=na,nrow=ny-1) #fish weight on January 1st
  for(y in 1:(ny-1)){
    for(a in 1:na){ 
      age <- a-1
      if((sort(unique(Year))[y]-age) %in% cohorts){
        coh.idx <- which(cohorts %in% (sort(unique(Year))[y]-age))
      }
      fish_length[y,a] <- eval(VbgrF,list(age=age+0.5-t_hatch,L1=L1[coh.idx],L2=L2[coh.idx],L3=L3[coh.idx],
                                          t1=t1,t3=t3))
      L_J1 <- ifelse(age==0,0,eval(VbgrF,list(age=age-t_hatch,L1=L1[coh.idx],L2=L2[coh.idx],L3=L3[coh.idx],
                                              t1=t1,t3=t3)))
      
      if(sort(unique(Year))[y] %in% LW$year){
        lw.idx <- which(LW$year %in% (sort(unique(Year))[y]))
      }
      fish_weight[y,a] <- LW$a[lw.idx]*fish_length[y,a]^LW$b[lw.idx]
      fish_weight_J1[y,a] <- LW$a[lw.idx]*L_J1^LW$b[lw.idx]
    }
  }
  
  N <- exp(logN)
  TSB <- colSums(N[,-ny]*t(fish_weight_J1))
  SSB <- colSums(mo*N[,-ny]*t(fish_weight_J1)) # For comparison with stock assessment
  k <- 100 # for numerical stability in  AD
  S <- 1.0 / (1.0 + exp(k * (fish_length - Lmax)))   # ≈1 when L < Lmax, ≈0 when L > Lmax
  M.ish <- t(M_slope * (fish_length - Lmax) * S)
  F.ish <- t(fish_length*F_scale)
  Z.ish <- M.ish+F.ish
  logR <- logN[2,] # Age 1 for comparability
  
  
  
  ADREPORT(logN)
  ADREPORT(TSB)
  ADREPORT(M.ish)
  ADREPORT(F.ish)
  ADREPORT(Z.ish)
  ADREPORT(logPred)
  ADREPORT(logR)
  ADREPORT(SSB)
  ret
}
obj <- MakeADFun(f,par,random=c("logN","logM.rand","logF.rand"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdrep <- sdreport(obj)
sdrep
###############


df_pred <- data.frame(logpred_N=(as.list(sdrep,"Est",report=TRUE)$logPred),
                      t=dat$t,cohort=dat$cohort,
                      sd=as.list(sdrep,"Std",report=TRUE)$logPred)

# plot on log-scale. Read d from top of page
ggplot(data=df_pred,aes(x=t,y=logpred_N)) +  
  geom_line(linewidth=1)+
  facet_wrap(~cohort,scales="free_y")+geom_point(data=d,aes(x=t,y=log(N)),color="red",fill="red",size=1.2)+
  geom_ribbon(aes(ymin = logpred_N-1.96*sd, ymax = logpred_N+1.96*sd), alpha = 0.2,linetype=0) +
  ggtitle("WB cod Surba")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    legend.title = element_blank())


logN <- melt(as.list(sdrep,"Est",report=TRUE)$logN)
names(logN) <- c("Age","Year","logN")
logN$Age <- (0:7)[logN$Age]
logN$Year <- (min(dat$Year):(max(dat$Year)+1))[logN$Year]
logN$cohort <- logN$Year-logN$Age
Nsd <- melt(as.list(sdrep,"Std",report=TRUE)$logN)
logN$sd <- Nsd$value


ggplot(data=logN,aes(x=Year,y=log10(exp(logN)),col=factor(Age))) +  
  geom_line(linewidth=1)+facet_wrap(~Age)+
  geom_ribbon(aes(ymin = log10(exp(logN-1.96*sd)), ymax = log10(exp(logN+1.96*sd))), alpha = 0.2) +
  ggtitle("WB cod Surba")+ylab("N (log10-scale)")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    legend.title = element_blank(),
                    legend.position = "none")



TSB <- data.frame(TSB=as.list(sdrep,"Est",report=TRUE)$TSB,
                  year=sort(unique(dat$Year)),
                  sd=as.list(sdrep,"Std",report=TRUE)$TSB )
ggplot(data=TSB,aes(x=year,y=TSB/mean(TSB))) +  
  geom_line(linewidth=1)+
  geom_ribbon(aes(ymin = (TSB-1.96*sd)/mean(TSB), ymax = (TSB+1.96*sd)/mean(TSB)), alpha = 0.2) +
  ggtitle("WB cod Surba")+ylab("Relative TSB")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    legend.title = element_blank())

SSB <- data.frame(SSB=as.list(sdrep,"Est",report=TRUE)$SSB,
                  year=sort(unique(dat$Year)),
                  sd=as.list(sdrep,"Std",report=TRUE)$SSB )
ggplot(data=SSB,aes(x=year,y=SSB/mean(SSB))) +  
  geom_line(linewidth=1)+
  geom_ribbon(aes(ymin = (SSB-1.96*sd)/mean(SSB), ymax = (SSB+1.96*sd)/mean(SSB)), alpha = 0.2) +
  ggtitle("WB cod Surba")+ylab("Relative SSB")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    legend.title = element_blank())

logR <- data.frame(logR=as.list(sdrep,"Est",report=TRUE)$logR,
                  year=c(sort(unique(dat$Year)),2026),
                  sd=as.list(sdrep,"Std",report=TRUE)$logR )

ggplot(data=logR,aes(x=year,y=(logR)/mean((logR)))) +  
  geom_line(linewidth=1)+
  geom_ribbon(aes(ymin = (logR-1.96*sd)/mean((logR)), ymax = (logR+1.96*sd)/mean((logR))), alpha = 0.2) +
  ggtitle("WB cod Surba")+ylab("Relative Recruitment (log-scale)")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    legend.title = element_blank())




Z <- melt(as.list(sdrep,"Est",report=TRUE)$Z.ish)
names(Z) <- c("Age","Year","Z")
Z$Age <- (0:7)[Z$Age]
Z$Year <- (min(dat$Year):(max(dat$Year)+1))[Z$Year]
Z$cohort <- Z$Year-Z$Age
Zsd <- melt(as.list(sdrep,"Std",report=TRUE)$Z.ish)
Z$sd <- Zsd$value


ggplot(data=Z,aes(x=Year,y=Z,col=factor(Age))) +  
  geom_line(linewidth=1)+facet_wrap(~Age,scales = "free_y")+
  geom_ribbon(aes(ymin = Z-1.96*sd, ymax = Z+1.96*sd), alpha = 0.2,linetype=0) +
  ggtitle("WB cod Surba")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    legend.title = element_blank(),
                    legend.position = "none")


M <- melt(as.list(sdrep,"Est",report=TRUE)$M.ish)
names(M) <- c("Age","Year","M")
M$Age <- (0:7)[M$Age]
M$Year <- (min(dat$Year):(max(dat$Year)+1))[M$Year]
M$cohort <- M$Year-M$Age
Msd <- melt(as.list(sdrep,"Std",report=TRUE)$M.ish)
M$sd <- Msd$value

ggplot(data=M,aes(x=Year,y=M,col=factor(Age))) +  
  geom_line(linewidth=1)+facet_wrap(~Age)+
  geom_ribbon(aes(ymin = M-1.96*sd, ymax = M+1.96*sd), alpha = 0.2,linetype=0) +
  ggtitle("WB cod Surba")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    legend.title = element_blank())


F <- melt(as.list(sdrep,"Est",report=TRUE)$F.ish)
names(F) <- c("Age","Year","F")
F$Age <- (0:7)[F$Age]
F$Year <- (min(dat$Year):(max(dat$Year)+1))[F$Year]
F$cohort <- F$Year-F$Age
Fsd <- melt(as.list(sdrep,"Std",report=TRUE)$F.ish)
F$sd <- Fsd$value

ggplot(data=F,aes(x=Year,y=F,col=factor(Age))) +  
  geom_line(linewidth=1)+facet_wrap(~Age)+
  geom_ribbon(aes(ymin = F-1.96*sd, ymax = F+1.96*sd), alpha = 0.2,linetype=0) +
  ggtitle("WB cod Surba")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    legend.title = element_blank())




Q1 <- d %>% filter(Quarter==1) %>% select(N,Age,Year)
ggplot(data=Q1,aes(x=Year,y=(N),col=factor(Age))) +  
  geom_line(linewidth=1)+facet_wrap(~Age,scales = "free_y")+
  ggtitle("Survey Q1")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    legend.title = element_blank())

Q4 <- d %>% filter(Quarter==4) %>% select(N,Age,Year)
ggplot(data=Q4,aes(x=Year,y=(N),col=factor(Age))) +  
  geom_line(linewidth=1)+facet_wrap(~Age,scales = "free_y")+
  ggtitle("Survey Q4")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    legend.title = element_blank())


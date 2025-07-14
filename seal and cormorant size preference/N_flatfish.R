source('prepare DATRAS.R')
rm(list=setdiff(ls(),c('hl_N')))

#hl_N <- read.table('hl_N.CB.csv',header=TRUE,sep=';')
#quarter.day <- data.frame(Quarter=1:4,jday=c(0.1676808,0.2902867,0.6282315,0.8768811))
#hl_N$jday <- quarter.day$jday[hl_N$Quarter]
#hl_N$Species <- recode(hl_N$Species,"Clupea harengus"="herring",
#                       "Clupea harengus"="cod","Limanda limanda"="dab",
#                       "Platichthys flesus"="flounder","Pleuronectes platessa"="plaice")
#hl_N$CPUE <- hl_N$N # for convenience
# prepare data
#####

hl_N$cohort <- hl_N$Year-hl_N$Age
hl_N$yd <- hl_N$Age+hl_N$jday
names(hl_N)[names(hl_N)=="species"] <- "Species"
#hl_N$N_age <- hl_N$N
# conseptual plot
#cohort.y <- hl_N %>% filter(cohort==1990)
#cohort.y_ag <- aggregate(N_age~Age+Quarter+Species+yd,data=cohort.y,FUN = sum)
#cohort.0_15 <- cohort.y_ag %>% filter(Age<=15)
#ggplot(data=cohort.0_15,aes(x=yd,y=log(N_age),color = Species))+
#  geom_point(pch=19,cex=2)+
#  labs(x='Age [years]',y=~paste('CPUE [#  ',h^-1,']'))+
#  theme(axis.title.x = element_text(size = 12),
#        axis.title.y = element_text(size = 12),
#        axis.text = element_text(size = 10),
#        legend.position =  'none')
# estimate mortality, Z (M+F), 1985-2014
hl_N.ag <-  aggregate(N_age~Age+Quarter+Species+cohort+yd+Year,data=hl_N %>% filter(Species %in% c("flounder","plaice","dab")),FUN = sum)
haulDura <- aggregate(HaulDur~haulID+Quarter+Year,data=hl_N,FUN = mean)
hd <- aggregate(HaulDur~Quarter+Year,data=haulDura,FUN = sum)
hl_N.ag <- hl_N.ag %>% left_join(hd)
hl_N.ag$CPUE <- hl_N.ag$N_age/hl_N.ag$HaulDur*60 # CPUE in ind. pr. hour
rm(list = c('haulDura','hd'))
sort(table(hl_N.ag$yd))
hl_N.ag <- hl_N.ag[round((round(hl_N.ag$yd,2)-floor(hl_N.ag$yd)),2)!=0.29,]
sort(table(hl_N.ag$yd))

## THis IS ALL DONE IN ALK_SurveyIndex


#hl_N.ag <-  aggregate(CPUE~Age+Quarter+Species+cohort+yd+Year,data=hl_N %>% filter(Species %in% c("flounder","plaice","dab")),FUN = sum)

# another conceptual plot
cohort.y <- hl_N.ag %>% filter(cohort==1991)
cohort.y_ag <- aggregate(CPUE~Age+Quarter+Species+yd,data=cohort.y,FUN = sum)
cohort.0_15 <- cohort.y_ag %>% filter(Age<10)
ggplot(data=cohort.0_15,aes(x=yd,y=log(CPUE),color = Species))+
  geom_point(pch=19,cex=2)+ylim(-10,5)+
  facet_wrap(~Species)+geom_point(data = cohort.0_15 %>% filter(Age==8),aes(x=yd,y=log(CPUE)),color = "black")+
  labs(x='Age [years]',y=~paste('CPUE [#  ',h^-1,']'))+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position =  'none')

ggplot(data=cohort.0_15,aes(x=yd,y=log(CPUE),color = Species))+
  geom_point(pch=19,cex=2)+ylim(-10,5)+
  facet_wrap(~Species)+geom_point(data = cohort.0_15 %>% filter(Age==8),aes(x=yd,y=log(CPUE)),color = "black")+
  labs(x='Age [years]',y=~paste('CPUE [#  ',h^-1,']'))+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position =  'none')
#####

# estimate survey index for cohorts with sufficient information (i.e, before 2021)
#####
cohort <- 1979:2020 # cohorts to estimate
d <- hl_N.ag %>% filter(Quarter %in% c(1,4) & yd>=3 & Age<15)
N_est <- data.frame(cohort=NA,Age=NA,Species="O",N=NA)

us.ages <- sort(unique(hl_N.ag$yd[which(hl_N.ag$yd<3)])) # ages for undersampling investigation
us <-  array(NA, dim = c(length(cohort), length(us.ages), 3))

flat.info <- array(NA, dim = c(length(cohort), 4, 3)) 

#d <- d[-which(d$cohort==2017 & d$Species=="flounder" & d$CPUE<1),] # remove an outlier
# dim: cohort, statistic, Species. Statistics include: 1=slope (M), 2=intercept (R), 3=adj. R^2, 4=n (available ages)

colors <- colorRampPalette(c("red", "blue"))(length(cohort))
for(i in 1:length(cohort)){
  cohort.i <- cohort[i]
  dada <- d %>% filter(cohort==cohort.i)
  fit.flounder <- lm(log(CPUE)~yd,data=dada %>% filter(Species=="flounder" & Age<=(min(Age)+7) & CPUE>0))
  fit.plaice <- lm(log(CPUE)~yd,data=dada %>% filter(Species=="plaice" & Age<=(min(Age)+7) & CPUE>0))
  fit.dab <- lm(log(CPUE)~yd,data=dada %>% filter(Species=="dab" & Age<=(min(Age)+7) & CPUE>0))
  flat.info[i,,1] <- c(coef(fit.flounder)[2],coef(fit.flounder)[1],summary(fit.flounder)$adj.r.squared,length(unique(dada$Age[dada$Species=="flounder"])))
  flat.info[i,,2] <- c(coef(fit.plaice)[2],coef(fit.plaice)[1],summary(fit.plaice)$adj.r.squared,length(unique(dada$Age[dada$Species=="plaice"])))
  flat.info[i,,3] <- c(coef(fit.dab)[2],coef(fit.dab)[1],summary(fit.dab)$adj.r.squared,length(unique(dada$Age[dada$Species=="dab"])))
  
  ages <- seq(0,20*365,by=1)/365
  n_ages <- length(ages)
  cohorts <- rep(cohort.i,n_ages)
  
  flounder <- data.frame(cohort=cohorts,Age=ages,Species=rep("flounder",n_ages),
                         N=exp(predict(fit.flounder,newdata = data.frame(yd=ages))))
  plaice <- data.frame(cohort=cohorts,Age=ages,Species=rep("plaice",n_ages),
                       N=exp(predict(fit.plaice,newdata = data.frame(yd=ages))))
  dab <- data.frame(cohort=cohorts,Age=ages,Species=rep("dab",n_ages),
                    N=exp(predict(fit.dab,newdata = data.frame(yd=ages))))
  N_flat <- rbind(flounder,plaice,dab)
  N_est <- rbind(N_est,flounder,plaice,dab)
  
  usd <- hl_N.ag %>% filter(cohort==cohort.i & yd %in% us.ages)
  
  
  if("flounder" %in% usd$Species){
    f <- usd %>% filter(Species=="flounder")
    f.idx <- us.ages %in% f$yd
    est.f <- flounder$N[round(ages*365) %in% round(us.ages*365)]
    us[i,f.idx,1] <- f$CPUE/est.f[f.idx]
  }
  if("plaice" %in% usd$Species){
    p <- usd %>% filter(Species=="plaice")
    p.idx <- us.ages %in% p$yd
    est.p <- plaice$N[round(ages*365) %in% round(us.ages*365)]
    us[i,p.idx,2] <- p$CPUE/est.p[p.idx]
  }
  if("dab" %in% usd$Species){
    da <- usd %>% filter(Species=="dab")
    d.idx <- us.ages %in% da$yd
    est.d <- dab$N[round(ages*365) %in% round(us.ages*365)]
    us[i,d.idx,3] <- da$CPUE/est.d[d.idx]
  }
  if(i==1){
    #plot(flounder$Age,log(flounder$N),type = 'l',lwd=2,col=colors[1],ylim=c(-10,5),
    #     ylab = "log(CPUE)",xlab = "Fish age [years]",main="flounder")
  }
  if(i>1){
    #lines(flounder$Age,log(flounder$N),type = 'l',lwd=2,col=colors[i])
  }
}
N_est <- N_est[-1,]
#before1989 <- d %>% filter(cohort<1989) %>% select(cohort,yd,Species,CPUE)
#names(before1989)[c(2,4)] <- c("Age","N")
#N_est <- rbind(N_est,before1989)
N_est$year <- floor(N_est$cohort+N_est$Age)
N_est <- N_est %>% filter(year<=2024)
N_est <- N_est %>% filter(year>=1991)
N_est$jd <- N_est$Age-floor(N_est$Age)



# for checking
ggplot(data=dada,aes(x=yd,y=log(CPUE),color = Species))+
  geom_point(pch=19,cex=2)+xlim(0,15)+
  facet_wrap(~Species)+geom_line(data=N_flat,aes(x=Age,y=log(N),color = Species))+
  labs(x='Age [years]',y=~paste('CPUE [#  ',h^-1,']'))+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position =  'none')

# plot flatfish statistrics
par(mfrow=c(1,3))
# M
plot(cohort,-flat.info[,1,1],ylim=c(0,1.3),main="Flounder",ylab="M",type='l',lwd=2)
plot(cohort,-flat.info[,1,2],ylim=c(0,1.3),main="Plaice",ylab="M",type='l',lwd=2)
plot(cohort,-flat.info[,1,3],ylim=c(0,1.3),main="Dab",ylab="M",type='l',lwd=2)
# R
plot(cohort,flat.info[,2,1],ylim=c(3,10),main="Flounder",ylab="R-index",type='l',lwd=2)
plot(cohort,flat.info[,2,2],ylim=c(3,10),main="Plaice",ylab="R-index",type='l',lwd=2)
plot(cohort,flat.info[,2,3],ylim=c(3,10),main="Dab",ylab="R-index",type='l',lwd=2)
# adj. R^2
plot(cohort,flat.info[,3,1],ylim=c(0,1),main="Flounder",ylab="Adj. R^2",type='l',lwd=2)
plot(cohort,flat.info[,3,2],ylim=c(0,1),main="Plaice",ylab="Adj. R^2",type='l',lwd=2)
plot(cohort,flat.info[,3,3],ylim=c(0,1),main="Dab",ylab="Adj. R^2",type='l',lwd=2)
# R
plot(-flat.info[,1,1],flat.info[,2,1],xlim=c(0,1.3),ylim=c(3,10),main="Flounder",ylab="R-index",xlab="M")
plot(-flat.info[,1,2],flat.info[,2,2],xlim=c(0,1.3),ylim=c(3,10),main="Plaice",ylab="R-index",xlab="M")
plot(-flat.info[,1,3],flat.info[,2,3],xlim=c(0,1.3),ylim=c(3,10),main="Dab",ylab="R-index",xlab="M")



#####

# estimate undersampling
#####
#us[us>2] <- NA
undersampling <- data.frame(Species = rep(c("flounder", "plaice","dab"),each=length(us.ages)),
                            age = rep(us.ages,3),us=c(colMeans(us[,,1],na.rm = TRUE),
                                                      colMeans(us[,,2],na.rm = TRUE),
                                                      colMeans(us[,,3],na.rm = TRUE)))

undersampling <- data.frame(Species = rep(c("flounder", "plaice","dab"),each=length(us.ages)),
                            age = rep(us.ages,3),us=c(apply(us[,,1],2,median,na.rm = TRUE),
                                                      apply(us[,,2],2,median,na.rm = TRUE),
                                                      apply(us[,,3],2,median,na.rm = TRUE)))
#plot undersampling
ggplot(data=undersampling,aes(x=age,y=us,color = Species))+
  geom_point(pch=19,cex=2)+ylim(0,1.5)+
  facet_wrap(~Species)+
  labs(x='Age [years]',y=~paste('CPUE [#  ',h^-1,']'))+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position =  'none')

func <- function(x,A,B,k){
  -A * (exp(-k * x))+B
}
par(mfrow=c(1,3))
x <- (0:400)/100
d.dab <- undersampling %>% filter(Species=="dab")
dab.us <-  nls(log(us) ~ -A*exp(-k*age)+B,
               data=d.dab, start = list(A=12,B=0, k=1))
plot(undersampling$age[undersampling$Species=="dab"],
     log(undersampling$us[undersampling$Species=="dab"]))
lines(x,func(x,coef(dab.us)[1],coef(dab.us)[2],coef(dab.us)[3]),type='l',lwd=2)


d.flounder <- undersampling %>% filter(Species=="flounder")
flounder.us <-  nls(log(us) ~ -A*exp(-k*age)+B,
                    data=d.flounder, start = list(A=12,B=0, k=1))
plot(undersampling$age[undersampling$Species=="flounder"],
     log(undersampling$us[undersampling$Species=="flounder"]))
lines(x,func(x,coef(flounder.us)[1],coef(flounder.us)[2],coef(flounder.us)[3]),type='l',lwd=2)

d.plaice <- undersampling %>% filter(Species=="plaice")
plaice.us <-    nls(log(us) ~ -A*exp(-k*age)+B,
                    data=d.plaice, start = list(A=12,B=0, k=1))
plot(d.plaice$age,
     log(d.plaice$us))
lines(x,func(x,coef(plaice.us)[1],coef(plaice.us)[2],coef(plaice.us)[3]),type='l',lwd=2)


par(mfrow=c(1,1))
plot(x,func(x,coef(dab.us)[1],coef(dab.us)[2],coef(dab.us)[3]),type='l',lwd=2,col="red",xlim=c(1,3),ylim=c(-10,0))
lines(x,func(x,coef(flounder.us)[1],coef(flounder.us)[2],coef(flounder.us)[3]),type='l',lwd=2,col="orange")
lines(x,func(x,coef(plaice.us)[1],coef(plaice.us)[2],coef(plaice.us)[3]),type='l',lwd=2,col="yellow3")
points(undersampling$age[undersampling$Species=="dab"],
       log(undersampling$us[undersampling$Species=="dab"]),col="red")
points(undersampling$age[undersampling$Species=="flounder"],
       log(undersampling$us[undersampling$Species=="flounder"]),col="orange")
points(d.plaice$age,
       log(d.plaice$us),col="yellow3")


#####

# estimate estimate survey index for cohorts 2021 to 24 through undersampling functions
#####
du <- hl_N.ag %>% filter(Quarter %in% c(1,4) & cohort>2020)
N_est.recent <- data.frame(cohort=NA,Age=NA,Species="O",N=NA)
years <- 2021:2024
par(mfrow=c(2,2))
for(i in 1:length(years)){
  dada <- du %>% filter(cohort==years[i] & Species %in% c("flounder","plaice","dab"))
  x <- (0:4000)/1000
  d.f <- du %>% filter(cohort==years[i] & Species=="flounder")
  d.f$CPUEcor <- c(d.f$CPUE[d.f$yd<3]/exp(func(d.f$yd[d.f$yd<3],coef(flounder.us)[1],coef(flounder.us)[2],coef(flounder.us)[3])),
                   d.f$CPUE[d.f$yd>3])
  fit.flounder <- lm(log(CPUEcor)~yd,data=d.f)
  #plot(d.f$yd,log(d.f$CPUEcor),ylim=c(-1,10),main="Flounder",xlab="Age [years] + Julian day",ylab="Estimated log(CPUE)")
  #lines(x,predict(fit.flounder,newdata = data.frame(yd=x)))
  d.p <- du %>% filter(cohort==years[i] & Species=="plaice")
  d.p$CPUEcor <- c(d.p$CPUE[d.p$yd<3]/exp(func(d.p$yd[d.p$yd<3],coef(plaice.us)[1],coef(plaice.us)[2],coef(plaice.us)[3])),
                   d.p$CPUE[d.p$yd>3])
  fit.plaice <- lm(log(CPUEcor)~yd,data=d.p)
  plot(d.p$yd,log(d.p$CPUEcor),ylim=c(-1,10),main="Plaice",xlab="Age [years] + Julian day",ylab="Estimated log(CPUE)")
  lines(x,predict(fit.plaice,newdata = data.frame(yd=x)))
  d.d <- du %>% filter(cohort==years[i] & Species=="dab")
  d.d$CPUEcor <- c(d.d$CPUE[d.d$yd<3]/exp(func(d.d$yd[d.d$yd<3],coef(dab.us)[1],coef(dab.us)[2],coef(dab.us)[3])),
                   d.d$CPUE[d.d$yd>3])
  fit.dab <- lm(log(CPUEcor)~yd,data=d.d)
  #plot(d.d$yd,log(d.d$CPUEcor),ylim=c(-1,10),main="Dab",xlab="Age [years] + Julian day",ylab="Estimated log(CPUE)")
  #lines(x,predict(fit.dab,newdata = data.frame(yd=x)))
  
  ages <- seq(0,(max(years)-years[i]+1)*365,by=1)/365
  n_ages <- length(ages)
  cohorts <- rep(years[i],n_ages)
  
  flounder <- data.frame(cohort=cohorts,Age=ages,Species=rep("flounder",n_ages),
                         N=exp(predict(fit.flounder,newdata = data.frame(yd=ages))))
  plaice <- data.frame(cohort=cohorts,Age=ages,Species=rep("plaice",n_ages),
                       N=exp(predict(fit.plaice,newdata = data.frame(yd=ages))))
  dab <- data.frame(cohort=cohorts,Age=ages,Species=rep("dab",n_ages),
                    N=exp(predict(fit.dab,newdata = data.frame(yd=ages))))
  N_est.recent <- rbind(N_est.recent,flounder,plaice,dab)
  
}
N_est.recent <- N_est.recent[-1,]
N_est.recent$year <- floor(N_est.recent$cohort+N_est.recent$Age)
N_est.recent <- N_est.recent %>% filter(year<=2024)
N_est.recent$jd <- N_est.recent$Age-floor(N_est.recent$Age)
#####

# combine survey index estimates and plot N and TSB
#####
N_est.fin <- rbind(N_est,N_est.recent)
write.table(N_est.fin,"data/N_est.flatfish.csv",sep = ';')
N_est.fin <- read.table("data/N_est.flatfish.csv",sep = ';',header=TRUE)
rm(list=setdiff(ls(),c('N_est.fin','hl_N')))

N_jan1 <- N_est.fin %>% filter(Age %in% 0:30)
TN_est <- aggregate(N~Species+year,data = N_jan1 %>% filter(Age>=3),FUN=sum)



t <- aggregate(CPUE~Species+Year+Quarter,data=hl_N.ag %>% filter(Age>=3),FUN=sum)
day_q1 <- unique(round(hl_N.ag$yd*365)-round(floor(hl_N.ag$yd)*365))[1]
N_est.fin$Jday <- round(N_est.fin$jd*365)
N_Q1 <- N_est.fin %>% filter(Jday==day_q1)
T1N_est <- aggregate(N~Species+year,data = N_Q1 %>% filter(Age>=3),FUN=sum)
# plot N
ggplot(data=TN_est,aes(x=year,y=N,color = Species))+
  geom_point(pch=19,cex=2)+ylim(0,10)+
  facet_wrap(~Species)+geom_line(data=t %>% filter(Quarter==1),aes(x=Year,y=CPUE,color = Species))+
  labs(x='Age [years]',y=~paste('CPUE [#  ',h^-1,']'))+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position =  'none')

ggplot(data=T1N_est,aes(x=year,y=N,color = Species))+
  geom_point(pch=19,cex=2)+
  facet_wrap(~Species)+geom_line(data=t %>% filter(Quarter==1),aes(x=Year,y=CPUE,color = Species))+
  labs(x='Age [years]',y=~paste('CPUE [#  ',h^-1,']'))+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position =  'none')


ggplot(data=t,aes(x=Year,y=CPUE,color = Species))+
  geom_line(pch=19,cex=2)+
  facet_wrap(~Species)+
  labs(x='Age [years]',y=~paste('CPUE [#  ',h^-1,']'))+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position =  'none')

t13 <- hl_N.ag %>% filter(Year==2013 & Age>=3 & Quarter==1)
N13 <- N_Q1 %>% filter(year==2013 & Age>=3)
#####

N_est.fin <- read.table("data/N_est.flatfish.csv",sep = ';',header=TRUE)
LW <- read.table("data/length-weight.csv",sep = ';',header=TRUE)
exp_fit <- function(L,a,b){
  a*L^b
}
# Define Flatfish functions and read data
#####
# define flounder growth function
vbgrFlounder <- function(age) {
  Linf <- 34.0031653 # estimated in "fit flatfish growths from DATRAS" 
  k <- 0.4537312  # estimated in "fit flatfish growths from DATRAS" 
  L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
  u <- 0.23 # Amplitude of seasonality (same as cod)
  w <- 0.91 # Timing of peak growth
  t_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching
  
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf
  
  return(L_t)
} #Original

# define plaice growth function
vbgrPlaice <- function(age) {
  Linf <- 41.637260  # estimated in "fit flatfish growths from DATRAS" 
  k <- 0.240703   # estimated in "fit flatfish growths from DATRAS" 
  L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
  u <- 0.23 # Amplitude of seasonality (same as cod)
  w <- 0.91 # Timing of peak growth
  t_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching
  
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf
  
  return(L_t)
} #Original

# define dab growth function
vbgrDab <- function(age) {
  Linf <- 26.5093927   # estimated in "fit flatfish growths from DATRAS" 
  k <- 0.4846326    # estimated in "fit flatfish growths from DATRAS" 
  L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
  u <- 0.23 # Amplitude of seasonality (same as cod)
  w <- 0.91 # Timing of peak growth
  t_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching
  
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf
  
  return(L_t)
} #Original

#####

Flounder_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Plaice_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
Dab_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 

flounder <- N_jan1 %>% filter(Species=="flounder" & Age>0)
flounder$L <- vbgrFlounder(flounder$Age-Flounder_hatch )
LW.f <- LW %>% filter(species=="flounder")
f.idx <- match(flounder$year,as.numeric(LW.f$year))
flounder$B <- exp_fit(flounder$L,LW.f$a[f.idx],LW.f$b[f.idx])*flounder$N

plaice <- N_jan1 %>% filter(Species=="plaice" & Age>0)
plaice$L <- vbgrPlaice(plaice$Age-Plaice_hatch )
LW.p <- LW %>% filter(species=="plaice")
p.idx <- match(plaice$year,as.numeric(LW.p$year))
plaice$B <- exp_fit(plaice$L,LW.p$a[p.idx],LW.p$b[p.idx])*plaice$N

dab <- N_jan1 %>% filter(Species=="dab" & Age>0)
dab$L <- vbgrDab(dab$Age-Dab_hatch )
LW.d <- LW %>% filter(species=="dab")
d.idx <- match(dab$year,as.numeric(LW.d$year))
dab$B <- exp_fit(dab$L,LW.d$a[d.idx],LW.d$b[d.idx])*dab$N

B <- aggregate(B~Species+year,data=rbind(flounder,plaice,dab),FUN=sum)

# plot N
ggplot(data=B,aes(x=year,y=B,color = Species))+
  geom_point(pch=19,cex=2)+
  facet_wrap(~Species)+
  labs(x='Age [years]',y=~paste('CPUE [#  ',h^-1,']'))+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position =  'none')

plot((0:500)/10,exp_fit((0:500)/10,0.017196246,2.862030))
lines((0:500)/10,exp_fit((0:500)/10,0.012088384 ,3.018525 ))

exp_fit(10,0.017196246,2.862030)




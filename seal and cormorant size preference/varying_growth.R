library(LaplacesDemon)
source('sd estimation.R')

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

###########
# herring
##########
# define growth parameters
L_hatch <- coefs$L_hatch_herring # length at hatching, cm (Pepin et al, 1997)
u <- coefs$u # Amplitude of seasonality (same as herring)
w <- coefs$w # Timing of peak growth
t_hatch <- coefs$t_hatch_herring # time of hatching, Julian day / 365 - spawning from Jan. to May (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching

herring <- df_all %>% filter(species=="herring" & metric=="mean") # herring

before <- 1:3
after <- 39:40
h_cohort <- sort(unique(herring$cohort))[-c(before,after)]
h_coef <- data.frame(k=rep(0,length(h_cohort)),Linf=rep(0,length(h_cohort)),
                     cohort=h_cohort,species="herring")

par <- list(logK = 0,logLinf=0)
t.c <- (0:1000)/100
colors <- colorRampPalette(c("red", "blue"))(length(h_cohort))

for (i in 1:length(h_cohort)){
  dada <- herring %>% filter(cohort ==h_cohort[i] & !is.na(value))
  
  obj_fn <- function(par) {
    k <- exp(par[["logK"]])
    Linf <- exp(par[["logLinf"]])
    pred <- vbgr(dada$age, k,Linf)
    sum((dada$value - pred)^2)
  }
  opt <- nlminb(par, objective = obj_fn)
  h_coef[i,1:2] <- exp(opt$par)
  
  if(i==1){
    plot(t.c, vbgr(t.c,exp(opt$par[1]),exp(opt$par[2])),type = 'l',lwd=2,col=colors[i],ylim=c(0,40),xlim=c(0,8),
         ylab = "length [cm]",xlab = "Fish age [years]",main="Herring cohorts 1986-2020")
  }
  if(i>1){
    lines(t.c, vbgr(t.c,exp(opt$par[1]),exp(opt$par[2])),type = 'l',lwd=2,col=colors[i])
  }
  #points(dada$age,dada$value)
  #readline() # press enter to continue
}
coef_before <- data.frame(k=rep(mean(h_coef$k[1:3]),length(before)),
                          Linf=rep(mean(h_coef$Linf[1:3]),length(before)),
                          cohort=min(h_coef$cohort)-rev(before),
                          species=h_coef$species[1])

coef_after <- data.frame(k=rep(mean(h_coef$k[last_index(h_coef$k,3)]),length(after)),
                         Linf=rep(mean(h_coef$Linf[last_index(h_coef$Linf,3)]),length(after)),
                         cohort=max(h_coef$cohort)+1:length(after),
                         species=h_coef$species[1])
h_coef <- rbind(coef_before,h_coef,coef_after)


coefs$herring_varying_growth <- h_coef
saveRDS(coefs,paste(wd_coef,"coefficients.R",sep=""))
#######

###########
# flounder
##########
# define growth parameters
L_hatch <- coefs$L_hatch_flounder # mean from Kennedy et al., 2007
u <- coefs$u # Amplitude of seasonality (same as cod)
w <- coefs$w # Timing of peak growth
t_hatch <- coefs$t_hatch_flounder# time of hatching, Julian day / 365 - from fiskeatlas 
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching

flounder <- df_all %>% filter(species=="flounder" & metric=="mean") # flounder

before <- 1:3
after <- 41:42
f_cohort <- sort(unique(flounder$cohort))[-c(before,after)]
f_coef <- data.frame(k=rep(0,length(f_cohort)),Linf=rep(0,length(f_cohort)),
                     cohort=f_cohort,species="flounder")

par <- list(logK = 0,logLinf=0)
t.c <- (0:1500)/100
colors <- colorRampPalette(c("red", "blue"))(length(f_cohort))

for (i in 1:length(f_cohort)){
  dada <- flounder %>% filter(cohort ==f_cohort[i] & !is.na(value))
  
  obj_fn <- function(par) {
    k <- exp(par[["logK"]])
    Linf <- exp(par[["logLinf"]])
    pred <- vbgr(dada$age, k,Linf)
    sum((dada$value - pred)^2)
  }
  opt <- nlminb(par, objective = obj_fn)
  f_coef[i,1:2] <- exp(opt$par)
  
  if(i==1){
    plot(t.c, vbgr(t.c,exp(opt$par[1]),exp(opt$par[2])),type = 'l',lwd=2,col=colors[i],ylim=c(0,50),xlim=c(0,8),
              ylab = "length [cm]",xlab = "Fish age [years]",main="Flounder cohorts 1984-2020")
  }
  if(i>1){
    lines(t.c, vbgr(t.c,exp(opt$par[1]),exp(opt$par[2])),type = 'l',lwd=2,col=colors[i])
  }
  #points(dada$age,dada$value)
  #readline() # press enter to continue
}
coef_before <- data.frame(k=rep(mean(f_coef$k[1:3]),length(before)),
                          Linf=rep(mean(f_coef$Linf[1:3]),length(before)),
                          cohort=min(f_coef$cohort)-rev(before),
                          species=f_coef$species[1])

coef_after <- data.frame(k=rep(mean(f_coef$k[last_index(f_coef$k,3)]),length(after)),
                         Linf=rep(mean(f_coef$Linf[last_index(f_coef$Linf,3)]),length(after)),
                         cohort=max(f_coef$cohort)+1:length(after),
                         species=f_coef$species[1])
f_coef <- rbind(coef_before,f_coef,coef_after)

coefs$flounder_varying_growth <- f_coef
saveRDS(coefs,paste(wd_coef,"coefficients.R",sep=""))
#######

###########
# plaice
##########
# define growth parameters
L_hatch <- coefs$L_hatch_plaice # mean from Kennedy et al., 2007
u <- coefs$u # Amplitude of seasonality (same as cod)
w <- coefs$w # Timing of peak growth
t_hatch <- coefs$t_hatch_plaice # time of hatching, Julian day / 365 - from fiskeatlas 
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching

plaice <- df_all %>% filter(species=="plaice" & metric=="mean") # plaice

before <- 1:4
after <- 41:42
p_cohort <- sort(unique(plaice$cohort))[-c(before,after)]
p_coef <- data.frame(k=rep(0,length(p_cohort)),Linf=rep(0,length(p_cohort)),
                     cohort=p_cohort,species="plaice")

par <- list(logK = 0,logLinf=0)
t.c <- (0:1500)/100
colors <- colorRampPalette(c("red", "blue"))(length(p_cohort))

for (i in 1:length(p_cohort)){
  dada <- plaice %>% filter(cohort ==p_cohort[i] & !is.na(value))
  dada <- dada %>% filter(!(age>10 & value<28)) # remove some outliers in recent cohorts
  
  obj_fn <- function(par) {
    k <- exp(par[["logK"]])
    Linf <- exp(par[["logLinf"]])
    pred <- vbgr(dada$age, k,Linf)
    sum((dada$value - pred)^2)
  }
  opt <- nlminb(par, objective = obj_fn)
  p_coef[i,1:2] <- exp(opt$par)
  if(i==1){
    plot(t.c, vbgr(t.c,exp(opt$par[1]),exp(opt$par[2])),type = 'l',lwd=2,col=colors[i],ylim=c(0,50),xlim=c(0,8),
         ylab = "length [cm]",xlab = "Fish age [years]",main="Plaice cohorts 1985-2020")
  }
  if(i>1){
    lines(t.c, vbgr(t.c,exp(opt$par[1]),exp(opt$par[2])),type = 'l',lwd=2,col=colors[i])
  }
  #points(dada$age,dada$value)
  #readline() # press enter to continue
}
coef_before <- data.frame(k=rep(mean(p_coef$k[1:3]),length(before)),
                          Linf=rep(mean(p_coef$Linf[1:3]),length(before)),
                          cohort=min(p_coef$cohort)-rev(before),
                          species=p_coef$species[1])

coef_after <- data.frame(k=rep(mean(p_coef$k[last_index(p_coef$k,3)]),length(after)),
                         Linf=rep(mean(p_coef$Linf[last_index(p_coef$Linf,3)]),length(after)),
                         cohort=max(p_coef$cohort)+1:length(after),
                         species=p_coef$species[1])
p_coef <- rbind(coef_before,p_coef,coef_after)

coefs$plaice_varying_growth <- p_coef
saveRDS(coefs,paste(wd_coef,"coefficients.R",sep=""))
#######

###########
# dab
##########
# define growth parameters
L_hatch <- coefs$L_hatch_dab # mean from Kennedy et al., 2007
u <- coefs$u # Amplitude of seasonality (same as cod)
w <- coefs$w # Timing of peak growth
t_hatch <- coefs$t_hatch_dab # time of hatching, Julian day / 365 - from fiskeatlas 
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching

dab <- df_all %>% filter(species=="dab" & metric=="mean") # dab

before <- 1:4
after <- 38:39
d_cohort <- sort(unique(dab$cohort))[-c(before,after)]
d_coef <- data.frame(k=rep(0,length(d_cohort)),Linf=rep(0,length(d_cohort)),
                     cohort=d_cohort,species="dab")
par <- list(logK = 0,logLinf=0)
t.c <- (0:1500)/100
colors <- colorRampPalette(c("red", "blue"))(length(d_cohort))

for (i in 1:length(d_cohort)){
  dada <- dab %>% filter(cohort ==d_cohort[i] & !is.na(value))
  
  obj_fn <- function(par) {
    k <- exp(par[["logK"]])
    Linf <- exp(par[["logLinf"]])
    pred <- vbgr(dada$age, k,Linf)
    sum((dada$value - pred)^2)
  }
  opt <- nlminb(par, objective = obj_fn)
  d_coef[i,1:2] <- exp(opt$par)
  
  if(i==1){
    plot(t.c, vbgr(t.c,exp(opt$par[1]),exp(opt$par[2])),type = 'l',lwd=2,col=colors[i],ylim=c(0,50),xlim=c(0,8),
         ylab = "length [cm]",xlab = "Fish age [years]",main="Dab cohorts 1988-2020")
  }
  if(i>1){
    lines(t.c, vbgr(t.c,exp(opt$par[1]),exp(opt$par[2])),type = 'l',lwd=2,col=colors[i])
  }
  #points(dada$age,dada$value)
  #readline() # press enter to continue
}
coef_before <- data.frame(k=rep(mean(d_coef$k[1:3]),length(before)),
                          Linf=rep(mean(d_coef$Linf[1:3]),length(before)),
                          cohort=min(d_coef$cohort)-rev(before),
                          species=d_coef$species[1])

coef_after <- data.frame(k=rep(mean(d_coef$k[last_index(d_coef$k,3)]),length(after)),
                         Linf=rep(mean(d_coef$Linf[last_index(d_coef$Linf,3)]),length(after)),
                         cohort=max(d_coef$cohort)+1:length(after),
                         species=d_coef$species[1])
d_coef <- rbind(coef_before,d_coef,coef_after)

coefs$dab_varying_growth <- d_coef
saveRDS(coefs,paste(wd_coef,"coefficients.R",sep=""))
#######


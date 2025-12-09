source('sd estimation.R')
rm(list=setdiff(ls(),c('df_all')))
# working directory to coefficient-file
wd_coef <- "C:/Users/mavast/Documents/GitHub/SOS data/"
coefs <- readRDS(paste(wd_coef,"coefficients.R",sep=""))


# SD function and simulation test
#####
# Function to calculate fish length SD as a function of age.
calculate_sd <- function(age, k_sd,age_sd,Linf_sd) {
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  # Partial derivatives
  dL_dLinf <- -exp(-k * (phi_t - phi_hatch +age)) + 1
  dL_dk <- (L_hatch - Linf) * (phi_hatch - phi_t -age ) * exp(-k * (phi_t +age - phi_hatch))
  dL_dage <- -(L_hatch - Linf)*k*exp(-k*(phi_t + age - phi_hatch))
  
  # Variance
  var_L <- (dL_dLinf * Linf_sd)^2 + (dL_dk * k_sd)^2 + (dL_dage * age_sd)^2
  
  # Standard deviation
  sd_L <- sqrt(var_L)
  
  return(sd_L)
}
# simulation test
#n <- 100000
#sd.0 <- matrix(rep(0,length(t.c)*n),nrow = n)
#t.i <- (t.c+t_hatch)
#X_t <- t.i-floor(t.i) # time of year [0,1]
#phi_t.i <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth

#for (i in 1:n){
#  k.i <- rnorm(1,mean=k,sd=exp(opt$par[1]))
#  Linf.i <- rnorm(1,mean=Linf,sd=exp(opt$par[3]))
#  fish.age <- t.c+rnorm(1,0,exp(opt$par[2]))
#  L <- (L_hatch-Linf.i)*exp(-k.i*(phi_t.i+fish.age-phi_hatch))+Linf.i
#lines(t,L,col='red')
#  sd.0[i,] <- L

#}
#lines(t.c,apply(sd.0,2,sd),lwd=2,col="red",lty = 'dashed')
#####


###########
# cod 
##########
cod <- df_all %>% filter(species=="cod" & metric=="sd")
# define growth parameters
L_hatch <- coefs$L_hatch_cod  # length at hatching, cm (Pepin et al, 1997)
u <- coefs$u # Amplitude of seasonality
w <- coefs$w # Timing of peak growth
t_hatch <- coefs$t_hatch_cod # time of hatching, Julian day / 365 - spawning from Jan. to July (a bit arbitrary from Hüssy et al., 2011), average set to peak spawning

#simulate larvae distribution by assuming 95 % CI between Jan. 15th and July 15th 
r <- round(rnorm(1000000,mean=t_hatch*365,sd=(197-t_hatch*365)/1.96))
idx <- r<=197 # index below the 97.5% quantile
t <- sort(unique(r[idx]))/365 # larvae hatching days
age <- 197/365-t # larvae age on May 15th
phi_t <- u*sin(2*pi*(t-w))/(2*pi) # seasonal variability in growth
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching

cohorts <- sort(unique(cod$cohort))
n_cohorts <- length(cohorts)
cod_sd_coef <- data.frame(k_sd=rep(0,n_cohorts),
                          A_sd=0,Linf_sd=0, cohort=0)

for(i in 1:n_cohorts){
  cohort.i <- cohorts[i]
  k <- coefs$cod_varying_growth$k[coefs$cod_varying_growth$cohort==cohort.i]
  Linf <- coefs$cod_varying_growth$Linf[coefs$cod_varying_growth$cohort==cohort.i]
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf # larvae length on July 15th
  df.hatch <- data.frame(age = age,length=L_t,n=as.matrix(table(r[idx]/365))[,1]) # combine in dataframe
  mean.weighted <- sum(df.hatch$length*df.hatch$n)/sum(df.hatch$n) # calculate mean larvae length on July 15th
  sd.weighted <- sqrt(sum(df.hatch$n*(df.hatch$length-mean.weighted)^2)/sum(df.hatch$n)) # calculate sd for larvae length on July 15th
  #sd(rep(df.hatch$length,times=df.hatch$n)) # test sd is calculated correctly
  #sd.weighted
  
  # Now estimate growth parameter uncertainties
  cod.sd <- cod$value[cod$cohort==cohorts[i]] # cod sd at age
  cod.sd <- c(sd.weighted,cod.sd) # add the calculated sd after spawning
  cod.age <- cod$age[cod$cohort==cohorts[i]]-t_hatch # get ages
  cod.age <- c(197/365-t_hatch,cod.age)  # add age of the calculated sd after spawning
  
  na.idx <- is.na(cod.sd)
  cod.sd <- cod.sd[!na.idx]
  cod.age <- cod.age[!na.idx]
  #plot(cod.age,cod.sd,ylim=c(0,20),xlab="Age [years]",ylab="SD [cm]]",main="Cod")
  #points(cod.age[1],cod.sd[1],pch=19,col="red")

  par <- list(logK_sd = -4,logA_sd=-2,logLinf_sd=3)
  
  obj_fn <- function(par) {
    k_sd <- exp(par[["logK_sd"]])
    age_sd <- exp(par[["logA_sd"]])
    Linf_sd <- exp(par[["logLinf_sd"]])
    pred <- calculate_sd(cod.age, k_sd, age_sd,Linf_sd)
    sum((cod.sd - pred)^2)
  }
  opt <- nlminb(par, objective = obj_fn)
  cod_sd_coef[i,] <- c(exp(opt$par),cohorts[i])
  #
  t.c <- (0:1500)/100
  plot(t.c, calculate_sd(t.c,exp(opt$par[1]),exp(opt$par[2]),exp(opt$par[3])),
       type = "l", col = "black", lwd = 2,
       ylab = 'Standard deviation',xlab='Age [years]',xlim=c(0,8.5),ylim=c(0,25),main=paste("cod -",cohorts[i]))
  points(cod.age,cod.sd)
  readline() # press enter to continue
}
coefs$cod_varying_sd <- cod_sd_coef
saveRDS(coefs,paste(wd_coef,"coefficients.R",sep=""))
###########

###########
# herring 
##########
herring <- df_all %>% filter(species=="herring" & metric=="sd")
# define growth parameters
L_hatch <- coefs$L_hatch_herring # figure 2 in Bauer et al., 2014
u <- coefs$u # Amplitude of seasonality
w <- coefs$w # Timing of peak growth
t_hatch <- coefs$t_hatch_herring # From figure 4 in Polte et al., 2021

#simulate larvae distribution
hatch_dur <- mean(c(33,37,43,44,46,48,49,52,54,56,57,63,65,66,68,69))/365 # From figure 4B in Polte et al., 2021
r <- round(rnorm(1000000,mean=t_hatch*365,sd=365*((hatch_dur/2)/1.96)))/365
#hist(r)
#quantile(r,c(0.025,0.975))
idx <- r<=t_hatch+hatch_dur/2 # index below the 97.5% quantile
t <- sort(unique(r[idx])) # larvae hatching days

age <- t_hatch+hatch_dur/2-t # larvae age on May 15th
phi_t <- u*sin(2*pi*(t-w))/(2*pi) # seasonal variability in growth
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching


cohorts <- sort(unique(herring$cohort))
n_cohorts <- length(cohorts)
herring_sd_coef <- data.frame(k_sd=rep(0,n_cohorts),
                              A_sd=0,Linf_sd=0, cohort=0)

for(i in 1:n_cohorts){
  cohort.i <- cohorts[i]
  k <- coefs$herring_varying_growth$k[coefs$herring_varying_growth$cohort==cohort.i]
  Linf <- coefs$herring_varying_growth$Linf[coefs$herring_varying_growth$cohort==cohort.i]
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf # larvae length on July 15th
  df.hatch <- data.frame(age = age,length=L_t,n=as.data.frame(table(r[idx]/365))$Freq) # combine in dataframe
  mean.weighted <- sum(df.hatch$length*df.hatch$n)/sum(df.hatch$n) # calculate mean larvae length on July 15th
  sd.weighted <- sqrt(sum(df.hatch$n*(df.hatch$length-mean.weighted)^2)/sum(df.hatch$n)) # calculate sd for larvae length on July 15th
  #sd(rep(df.hatch$length,times=df.hatch$n)) # test sd is calculated correctly
  #sd.weighted
  
  # Now estimate growth parameter uncertainties
  herring.sd <- herring$value[herring$cohort==cohorts[i]] # herring sd at age
  herring.sd <- c(sd.weighted,herring.sd) # add the calculated sd after spawning
  herring.age <- herring$age[herring$cohort==cohorts[i]]-t_hatch # get ages
  herring.age <- c(hatch_dur/2,herring.age)  # add age of the calculated sd after spawning
  
  na.idx <- is.na(herring.sd)
  herring.sd <- herring.sd[!na.idx]
  herring.age <- herring.age[!na.idx]
  #plot(herring.age,herring.sd,ylim=c(0,20),xlab="Age [years]",ylab="SD [cm]]",main="herring")
  #points(herring.age[1],herring.sd[1],pch=19,col="red")
  
  par <- list(logK_sd = -4,logA_sd=-2,logLinf_sd=3)
  
  obj_fn <- function(par) {
    k_sd <- exp(par[["logK_sd"]])
    age_sd <- exp(par[["logA_sd"]])
    Linf_sd <- exp(par[["logLinf_sd"]])
    pred <- calculate_sd(herring.age, k_sd, age_sd,Linf_sd)
    sum((herring.sd - pred)^2)
  }
  opt <- nlminb(par, objective = obj_fn)
  herring_sd_coef[i,] <- c(exp(opt$par),cohorts[i])
  #
  t.c <- (0:1500)/100
  plot(t.c, calculate_sd(t.c,exp(opt$par[1]),exp(opt$par[2]),exp(opt$par[3])),
       type = "l", col = "black", lwd = 2,
       ylab = 'Standard deviation',xlab='Age [years]',xlim=c(0,12),ylim=c(0,5),main=paste("herring -",cohorts[i]))
  points(herring.age,herring.sd)
  readline() # press enter to continue
}
herring_sd_coef[2,1:3] <- (herring_sd_coef[1,1:3]+herring_sd_coef[3,1:3])/2
coefs$herring_varying_sd <- herring_sd_coef
saveRDS(coefs,paste(wd_coef,"coefficients.R",sep=""))
###########

###########
# flounder 
##########
flounder <- df_all %>% filter(species=="flounder" & metric=="sd")
# define growth parameters
L_hatch <- coefs$L_hatch_flounder # figure 2 in Bauer et al., 2014
u <- coefs$u # Amplitude of seasonality
w <- coefs$w # Timing of peak growth
t_hatch <- coefs$t_hatch_flounder # From figure 4 in Polte et al., 2021

# simulate larvae distribution by assuming 95 % CI between Mar. 15th and Jul. 15th 
r <- round(rnorm(1000000,mean=t_hatch*365,sd=(197-t_hatch*365)/1.96))/365
#hist(r)
#quantile(r,c(0.025,0.975))*365 # test
idx <- r<= 197/365# index below the 97.5% quantile
t <- sort(unique(r[idx])) # larvae hatching days

age <- 197/365-t # larvae age on May 15th
phi_t <- u*sin(2*pi*(t-w))/(2*pi) # seasonal variability in growth
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching


cohorts <- sort(unique(flounder$cohort))
n_cohorts <- length(cohorts)
flounder_sd_coef <- data.frame(k_sd=rep(0,n_cohorts),
                               A_sd=0,Linf_sd=0, cohort=0)

for(i in 1:n_cohorts){
  cohort.i <- cohorts[i]
  k <- coefs$flounder_varying_growth$k[coefs$flounder_varying_growth$cohort==cohort.i]
  Linf <- coefs$flounder_varying_growth$Linf[coefs$flounder_varying_growth$cohort==cohort.i]
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf # larvae length on July 15th
  df.hatch <- data.frame(age = age,length=L_t,n=as.data.frame(table(r[idx]/365))$Freq) # combine in dataframe
  mean.weighted <- sum(df.hatch$length*df.hatch$n)/sum(df.hatch$n) # calculate mean larvae length on July 15th
  sd.weighted <- sqrt(sum(df.hatch$n*(df.hatch$length-mean.weighted)^2)/sum(df.hatch$n)) # calculate sd for larvae length on July 15th
  #sd(rep(df.hatch$length,times=df.hatch$n)) # test sd is calculated correctly
  #sd.weighted
  
  # Now estimate growth parameter uncertainties
  flounder.sd <- flounder$value[flounder$cohort==cohorts[i]] # flounder sd at age
  flounder.sd <- c(sd.weighted,flounder.sd) # add the calculated sd after spawning
  flounder.age <- flounder$age[flounder$cohort==cohorts[i]]-t_hatch # get ages
  flounder.age <- c(197/365-t_hatch,flounder.age)  # add age of the calculated sd after spawning
  
  na.idx <- is.na(flounder.sd)
  flounder.sd <- flounder.sd[!na.idx]
  flounder.age <- flounder.age[!na.idx]
  #plot(flounder.age,flounder.sd,ylim=c(0,20),xlab="Age [years]",ylab="SD [cm]]",main="flounder")
  #points(flounder.age[1],flounder.sd[1],pch=19,col="red")
  
  par <- list(logK_sd = -4,logA_sd=-2,logLinf_sd=3)
  
  obj_fn <- function(par) {
    k_sd <- exp(par[["logK_sd"]])
    age_sd <- exp(par[["logA_sd"]])
    Linf_sd <- exp(par[["logLinf_sd"]])
    pred <- calculate_sd(flounder.age, k_sd, age_sd,Linf_sd)
    sum((flounder.sd - pred)^2)
  }
  opt <- nlminb(par, objective = obj_fn)
  flounder_sd_coef[i,] <- c(exp(opt$par),cohorts[i])
  #
  t.c <- (0:1500)/100
  plot(t.c, calculate_sd(t.c,exp(opt$par[1]),exp(opt$par[2]),exp(opt$par[3])),
       type = "l", col = "black", lwd = 2,
       ylab = 'Standard deviation',xlab='Age [years]',xlim=c(0,15),ylim=c(0,10),main=paste("flounder -",cohorts[i]))
  points(flounder.age,flounder.sd)
  #readline() # press enter to continue
}
flounder_sd_coef[40:42,1:3] <- as.numeric(rep(flounder_sd_coef[39,1:3],each=3))
coefs$flounder_varying_sd <- flounder_sd_coef
saveRDS(coefs,paste(wd_coef,"coefficients.R",sep=""))
###########

###########
# plaice 
##########
plaice <- df_all %>% filter(species=="plaice" & metric=="sd")
# define growth parameters
L_hatch <- coefs$L_hatch_plaice # figure 2 in Bauer et al., 2014
u <- coefs$u # Amplitude of seasonality
w <- coefs$w # Timing of peak growth
t_hatch <- coefs$t_hatch_plaice # From figure 4 in Polte et al., 2021

# simulate larvae distribution by assuming 95 % CI between Nov. 15th and Mar. 15th 
r <- round(rnorm(1000000,mean=t_hatch*365,sd=(75-t_hatch*365)/1.96))/365
#hist(r)
#quantile(r,c(0.025,0.975))*365 # test
idx <- r<= 75/365# index below the 97.5% quantile
t <- sort(unique(r[idx])) # larvae hatching days

age <- 75/365-t # larvae age on Mar 15th
phi_t <- u*sin(2*pi*(t-w))/(2*pi) # seasonal variability in growth
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching

cohorts <- sort(unique(plaice$cohort))
n_cohorts <- length(cohorts)
plaice_sd_coef <- data.frame(k_sd=rep(0,n_cohorts),
                             A_sd=0,Linf_sd=0, cohort=0)
for(i in 1:n_cohorts){
  cohort.i <- cohorts[i]
  k <- coefs$plaice_varying_growth$k[coefs$plaice_varying_growth$cohort==cohort.i]
  Linf <- coefs$plaice_varying_growth$Linf[coefs$plaice_varying_growth$cohort==cohort.i]
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf # larvae length on July 15th
  df.hatch <- data.frame(age = age,length=L_t,n=as.data.frame(table(r[idx]/365))$Freq) # combine in dataframe
  mean.weighted <- sum(df.hatch$length*df.hatch$n)/sum(df.hatch$n) # calculate mean larvae length on July 15th
  sd.weighted <- sqrt(sum(df.hatch$n*(df.hatch$length-mean.weighted)^2)/sum(df.hatch$n)) # calculate sd for larvae length on July 15th
  #sd(rep(df.hatch$length,times=df.hatch$n)) # test sd is calculated correctly
  #sd.weighted
  
  # Now estimate growth parameter uncertainties
  plaice.sd <- plaice$value[plaice$cohort==cohorts[i]] # plaice sd at age
  plaice.sd <- c(sd.weighted,plaice.sd) # add the calculated sd after spawning
  plaice.age <- plaice$age[plaice$cohort==cohorts[i]]-t_hatch # get ages
  plaice.age <- c(75/365-t_hatch,plaice.age)  # add age of the calculated sd after spawning
  
  na.idx <- is.na(plaice.sd)
  plaice.sd <- plaice.sd[!na.idx]
  plaice.age <- plaice.age[!na.idx]
  #plot(plaice.age,plaice.sd,ylim=c(0,20),xlab="Age [years]",ylab="SD [cm]]",main="plaice")
  #points(plaice.age[1],plaice.sd[1],pch=19,col="red")
  
  par <- list(logK_sd = -4,logA_sd=-2,logLinf_sd=3)
  
  obj_fn <- function(par) {
    k_sd <- exp(par[["logK_sd"]])
    age_sd <- exp(par[["logA_sd"]])
    Linf_sd <- exp(par[["logLinf_sd"]])
    pred <- calculate_sd(plaice.age, k_sd, age_sd,Linf_sd)
    sum((plaice.sd - pred)^2)
  }
  opt <- nlminb(par, objective = obj_fn)
  plaice_sd_coef[i,] <- c(exp(opt$par),cohorts[i])
  #
  t.c <- (0:1500)/100
  plot(t.c, calculate_sd(t.c,exp(opt$par[1]),exp(opt$par[2]),exp(opt$par[3])),
       type = "l", col = "black", lwd = 2,
       ylab = 'Standard deviation',xlab='Age [years]',xlim=c(0,15),ylim=c(0,10),main=paste("plaice -",cohorts[i]))
  points(plaice.age,plaice.sd)
  #readline() # press enter to continue
}
plaice_sd_coef[42,1:3] <- plaice_sd_coef[41,1:3]
coefs$plaice_varying_sd <- plaice_sd_coef
saveRDS(coefs,paste(wd_coef,"coefficients.R",sep=""))
###########

###########
# dab 
##########
dab <- df_all %>% filter(species=="dab" & metric=="sd")
# define growth parameters
L_hatch <- coefs$L_hatch_dab # figure 2 in Bauer et al., 2014
u <- coefs$u # Amplitude of seasonality
w <- coefs$w # Timing of peak growth
t_hatch <- coefs$t_hatch_dab # From figure 4 in Polte et al., 2021

# simulate larvae distribution by assuming 95 % CI between Apr. 15th and Aug. 15th 
r <- round(rnorm(1000000,mean=t_hatch*365,sd=(228-t_hatch*365)/1.96))/365
#hist(r)
#quantile(r,c(0.025,0.975))*365 # test
idx <- r<= 228/365# index below the 97.5% quantile
t <- sort(unique(r[idx])) # larvae hatching days

age <- 228/365-t # larvae age on Mar 15th
phi_t <- u*sin(2*pi*(t-w))/(2*pi) # seasonal variability in growth
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching


cohorts <- sort(unique(dab$cohort))
n_cohorts <- length(cohorts)
dab_sd_coef <- data.frame(k_sd=rep(0,n_cohorts),
                          A_sd=0,Linf_sd=0, cohort=0)
for(i in 1:n_cohorts){
  cohort.i <- cohorts[i]
  k <- coefs$dab_varying_growth$k[coefs$dab_varying_growth$cohort==cohort.i]
  Linf <- coefs$dab_varying_growth$Linf[coefs$dab_varying_growth$cohort==cohort.i]
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf 
  df.hatch <- data.frame(age = age,length=L_t,n=as.data.frame(table(r[idx]/365))$Freq) # combine in dataframe
  mean.weighted <- sum(df.hatch$length*df.hatch$n)/sum(df.hatch$n) 
  sd.weighted <- sqrt(sum(df.hatch$n*(df.hatch$length-mean.weighted)^2)/sum(df.hatch$n)) 
  #sd(rep(df.hatch$length,times=df.hatch$n)) # test sd is calculated correctly
  #sd.weighted
  
  # Now estimate growth parameter uncertainties
  dab.sd <- dab$value[dab$cohort==cohorts[i]] # dab sd at age
  dab.sd <- c(sd.weighted,dab.sd) # add the calculated sd after spawning
  dab.age <- dab$age[dab$cohort==cohorts[i]]-t_hatch # get ages
  dab.age <- c(228/365-t_hatch,dab.age)  # add age of the calculated sd after spawning
  
  na.idx <- is.na(dab.sd)
  dab.sd <- dab.sd[!na.idx]
  dab.age <- dab.age[!na.idx]
  #plot(dab.age,dab.sd,ylim=c(0,20),xlab="Age [years]",ylab="SD [cm]]",main="dab")
  #points(dab.age[1],dab.sd[1],pch=19,col="red")
  
  par <- list(logK_sd = -4,logA_sd=-2,logLinf_sd=3)
  
  obj_fn <- function(par) {
    k_sd <- exp(par[["logK_sd"]])
    age_sd <- exp(par[["logA_sd"]])
    Linf_sd <- exp(par[["logLinf_sd"]])
    pred <- calculate_sd(dab.age, k_sd, age_sd,Linf_sd)
    sum((dab.sd - pred)^2)
  }
  opt <- nlminb(par, objective = obj_fn)
  dab_sd_coef[i,] <- c(exp(opt$par),cohorts[i])
  #
  t.c <- (0:1500)/100
  plot(t.c, calculate_sd(t.c,exp(opt$par[1]),exp(opt$par[2]),exp(opt$par[3])),
       type = "l", col = "black", lwd = 2,
       ylab = 'Standard deviation',xlab='Age [years]',xlim=c(0,8),ylim=c(0,10),main=paste("dab -",cohorts[i]))
  points(dab.age,dab.sd)
  readline() # press enter to continue
}
dab_sd_coef[1:2,1:3] <- rep(dab_sd_coef[3,1:3],2)
coefs$dab_varying_sd <- dab_sd_coef
saveRDS(coefs,paste(wd_coef,"coefficients.R",sep=""))
###########


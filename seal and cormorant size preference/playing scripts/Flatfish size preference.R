setwd("/Users/mavast/Desktop/Markus SOS")
library(sicegar)
library(ggplot2)
library(lubridate)
library(minpack.lm)

##############################
# flatfish size preference based on surveys
##############################
# Growth model, cohort spread by length, and individuals by age and time, functions: vbgr.fixed(), vbgr.sd(), N_age
#####
vbgr.dab <- function(t){ # von Bertalanffy growth rate. Estimated from DATRAS survey
  Linf <- 29.91456 # cm, estimated from survey
  L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
  k <- 0.50193 # vbgr growth rate, estimated from survey
  u <- 0.22 # Amplitude of seasonality (same as cod)
  w <- 0.90 # Timing of peak growth
  
  t_hatch <- 167/365 # time of hatching, Julian day / 365 - from fiskeatlas 
  X_t <- t-floor(t) # time of year [0,1]
  fish.age <- t-t_hatch
  
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+fish.age-phi_hatch))+Linf # vbgr estimated length
  return(L_t)
}

vbgr.flounder <- function(t){ # von Bertalanffy growth rate. Estimated from DATRAS survey
  Linf <- 38.57021 # cm, estimated from survey
  L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
  k <- 0.35311 # vbgr growth rate, estimated from survey
  u <- 0.22 # Amplitude of seasonality (same as cod)
  w <- 0.90 # Timing of peak growth
  
  t_hatch <- 75/365 # time of hatching, Julian day / 365 - from fiskeatlas 
  X_t <- t-floor(t) # time of year [0,1]
  fish.age <- t-t_hatch
  
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+fish.age-phi_hatch))+Linf # vbgr estimated length
  return(L_t)
}

vbgr.plaice <- function(t){ # von Bertalanffy growth rate. Estimated from DATRAS survey
  Linf <- 45.47677 # cm, estimated from survey
  L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
  k <- 0.22814 # vbgr growth rate, estimated from survey
  u <- 0.22 # Amplitude of seasonality (same as cod)
  w <- 0.90 # Timing of peak growth
  
  t_hatch <- 15/365 # time of hatching, Julian day / 365 - from fiskeatlas 
  X_t <- t-floor(t) # time of year [0,1]
  fish.age <- t-t_hatch
  
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+fish.age-phi_hatch))+Linf # vbgr estimated length
  return(L_t)
}


dab.sd <- function(t,n){ # von Bertalanffy growth rate. Parameters from McQueen et al., 2019
  # t in fraction year, Julian day/365
  Linf <- rnorm(n,mean=154.56,sd=16.75) # Assymptotic length - with stochasticity
  #Linf <- 154.56 # cm
  L_hatch <- 0.4 # length at hatching, cm (Pepin et al, 1997)
  t_hatch <- 32/365
  k <- rnorm(n,mean=0.11,sd=0.01508671) # vbgr growth rate
  u <- 0.22 # Amplitude of seasonality
  w <- 0.90 # Timing of peak growth
  #t_hatch <- 61/365 # time of hatching, Julian day / 365 - march 1st (van Deurs and Hüssy, 2011)
  #NB t_hatch is now more or less arbitrary
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth
  sd <- rep(0,length(t))
  for (i in 1:length(sd)){
    fish.age <- t[i]-t_hatch+rnorm(n,mean=0,sd=0.1) # time of hatching with stochasticity, Julian day / 365 - mean = march 1st
    L_t <- (L_hatch-Linf)*exp(-k*(phi_t[i]+fish.age-phi_hatch))+Linf # vbgr estimated length
    sd[i] <- sd(L_t)
  }
  return(sd)
}



N_age <- function(M_age,age,yrs){ # Number of individuals by age
  M_age$Z <- M_age$M_var+M_age$Fest
  Ns <- matrix(rep(0,(length(age)+1)*length(yrs)*365),ncol = length(age)+1)
  Ns[,length(age)+1] <- rep(yrs,each=365)
  x <- seq(1/365,1,length=365)
  for(i in 1:length(yrs)){
    for(j in 1:length(age)){
      di <- M_age %>% filter(ages==age[j] & years==yrs[i])
      dip1 <- M_age %>% filter(ages==(age[j]+1) & years==(yrs[i]+1))
      
      Ns[(365*(i-1)+1):(365*i),j] <- di$N*exp(-x*di$Z)
    }
  }
  return(Ns)
}
#####
source('prepare DATRAS.R')
source('sd estimation.R')
rm(list=setdiff(ls(),c('df')))

###########
# flounder
##########
# define growth parameters
L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
u <- 0.23 # Amplitude of seasonality (same as cod)
w <- 0.91 # Timing of peak growth
t_hatch <- ((75+197)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching

flounder.mean <- df$mean[df$species=="flounder"] # flounder mean length at age
flounder.age <- df$yd[df$species=="flounder"]-t_hatch # get ages
plot(flounder.age,flounder.mean)
flounder.age <- flounder.age[-(25:26)]
flounder.mean <- flounder.mean[-(25:26)]

# Function to calculate sd of L as a function of age. Based on partial derivatives of L (delta method)
vbgr <- function(age, k,Linf) {
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf
  
  return(L_t)
} #Original

par <- list(logK = 0,logLinf=0)

obj_fn <- function(par) {
  k <- exp(par[["logK"]])
  Linf <- exp(par[["logLinf"]])
  pred <- vbgr(flounder.age, k,Linf)
  sum((flounder.mean - pred)^2)
}
opt <- nlminb(par, objective = obj_fn)

exp(opt$par)
par(mfrow=c(1,1))
t.c <- (0:1500)/100
plot(t.c, vbgr(t.c,exp(opt$par[1]),exp(opt$par[2])),
     type = "l", col = "black", lwd = 2,
     ylab = 'Length [cm]',xlab='Age [years]',xlim=c(0,12),ylim=c(0,40),main="Flounder")
points(flounder.age,flounder.mean)
#######

###########
# plaice
##########
# define growth parameters
L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
u <- 0.23 # Amplitude of seasonality (same as cod)
w <- 0.91 # Timing of peak growth
t_hatch <- ((-46+75)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching

plaice.mean <- df$mean[df$species=="plaice"] # plaice mean length at age
plaice.age <- df$yd[df$species=="plaice"]-t_hatch # get ages
plot(plaice.age,plaice.mean)


# Function to calculate sd of L as a function of age. Based on partial derivatives of L (delta method)
vbgr <- function(age, k,Linf) {
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf
  
  return(L_t)
} #Original

par <- list(logK = 0,logLinf=0)

obj_fn <- function(par) {
  k <- exp(par[["logK"]])
  Linf <- exp(par[["logLinf"]])
  pred <- vbgr(plaice.age, k,Linf)
  sum((plaice.mean - pred)^2)
}
opt <- nlminb(par, objective = obj_fn)

exp(opt$par)
par(mfrow=c(1,1))
t.c <- (0:1500)/100
plot(t.c, vbgr(t.c,exp(opt$par[1]),exp(opt$par[2])),
     type = "l", col = "black", lwd = 2,
     ylab = 'Length [cm]',xlab='Age [years]',xlim=c(0,12),ylim=c(0,40),main="plaice")
points(plaice.age,plaice.mean)
#######

###########
# dab
##########
# define growth parameters
L_hatch <- mean(c(0.697,0.662)) # mean from Kennedy et al., 2007
u <- 0.23 # Amplitude of seasonality (same as cod)
w <- 0.91 # Timing of peak growth
t_hatch <- ((105+228)/2)/365 # time of hatching, Julian day / 365 - from fiskeatlas 
phi_hatch <- u*sin(2*pi*(t_hatch-w))/(2*pi) # seasonal variability in growth at hatching

dab.mean <- df$mean[df$species=="dab"] # dab mean length at age
dab.age <- df$yd[df$species=="dab"]-t_hatch # get ages
plot(dab.age,dab.mean)


# Function to calculate sd of L as a function of age. Based on partial derivatives of L (delta method)
vbgr <- function(age, k,Linf) {
  t <- age+t_hatch
  X_t <- t-floor(t) # time of year [0,1]
  phi_t <- u*sin(2*pi*(X_t-w))/(2*pi) # seasonal variability in growth
  
  L_t <- (L_hatch-Linf)*exp(-k*(phi_t+age-phi_hatch))+Linf
  
  return(L_t)
} #Original

par <- list(logK = 0,logLinf=0)

obj_fn <- function(par) {
  k <- exp(par[["logK"]])
  Linf <- exp(par[["logLinf"]])
  pred <- vbgr(dab.age, k,Linf)
  sum((dab.mean - pred)^2)
}
opt <- nlminb(par, objective = obj_fn)

exp(opt$par)
par(mfrow=c(1,1))
t.c <- (0:1500)/100
plot(t.c, vbgr(t.c,exp(opt$par[1]),exp(opt$par[2])),
     type = "l", col = "black", lwd = 2,
     ylab = 'Length [cm]',xlab='Age [years]',xlim=c(0,11),ylim=c(0,40),main="dab")
points(dab.age,dab.mean)
#######

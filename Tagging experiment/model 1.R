# Final model. ADDITIVE effects with refinding efficiency TRUE
#####
par <- list(logitP_cod = rep(0,length(unique(dat$length.class[dat$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dat$length.class[dat$species=="skrubbe"]))),
            logitTag.diff = 0,
            logitRefind.eff = 0)


nll <- function(par){
  getAll(par, dat)
  
  c_sizes <- unique(length.class[species=="torsk"])
  f_sizes <- unique(length.class[species=="skrubbe"])
  
  jnll <- -sum(dbinom(feed.exp$found,feed.exp$n,plogis(c(logitRefind.eff,logitRefind.eff+logitTag.diff)),log=TRUE))
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- (0:1)[unique(PIT)==PIT[i]]

    
    diff.PIT <- 0

    if(PIT.idx>0) diff.PIT <- logitTag.diff

    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT)
    }
  }
  
  cod <- log(plogis(logitP_cod)/plogis(logitRefind.eff))

  
  flounder <- log(plogis(logitP_flounder)/plogis(logitRefind.eff))


  ADREPORT(cod)
  ADREPORT(flounder)

  found <- OBS(found)
  n <- OBS(n)
  
  jnll <- jnll -sum(dbinom(found,n,p,log=TRUE))
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

M1 <- opt
# 1. Number of estimated parameters
k <- length(opt$par)
# 2. Log-likelihood is the negative of obj$fn at optimum
logL <- -obj$fn(opt$par)
# 3. Calculate AIC
AIC1 <- 2 * k - 2 * logL
print(AIC1)

res <- oneStepPredict(obj,method="oneStepGeneric",
                      discrete=TRUE)
res$n <- dat$n

res1 <- res
qqnorm(res$residual)
qqline(res$residual)

#####

#Looks at model outputs
#####

V_i <- res$n * res$mean * (1 - res$mean / res$n)
dispersion <- sum((res$observation - res$mean)^2 / V_i) / (length(res$observation) - length(obj$par))
paste("dispersion statitic =",round(dispersion,2)," i.e., underdispersed")
est <- as.list(sdr, "Est",report=TRUE)

sd <- as.list(sdr, "Std",report=TRUE)


df_year <- data.frame(species = rep(c("cod","flounder"),each=5*2),
                      length_class=c(rep(factor(c("11-15","16-20","21-25","26-30",">30"),levels=c("11-15","16-20","21-25","26-30",">30")),2),
                                     rep(factor(c("10-12","13-15","16-18","19-21",">21"),levels=c("10-12","13-15","16-18","19-21",">21")),2)),
                      year = rep(rep(c(2022,2024),each=5),2),
                      est = c(est$cod22,est$cod24,est$flounder22,est$flounder24),
                      sd = c(sd$cod22,sd$cod24,sd$flounder22,sd$flounder24))

df_total <- data.frame(species = rep(c("cod","flounder"),each=5),
                       length_class=c(factor(c("11-15","16-20","21-25","26-30",">30"),levels=c("11-15","16-20","21-25","26-30",">30")),
                                      factor(c("10-12","13-15","16-18","19-21",">21"),levels=c("10-12","13-15","16-18","19-21",">21"))),
                       est = c(est$cod_total,est$flounder_total),
                       sd = c(sd$cod_total,sd$flounder_total))

#####
dexp$fish <- paste(dexp$species,dexp$length.class)
dexp$year <- factor(dexp$year)


model <- glmmTMB(cbind(found, n - found) ~ fish+PIT+time,
                 family = binomial(link = "logit"),
                 data = dexp)
summary(model)

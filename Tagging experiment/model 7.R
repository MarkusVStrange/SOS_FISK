# Random release time - same release-SD pr. year
#####
release_date <- sort(unique(dat$release))
n_release <- length(release_date)
dat$time_df <- data.frame(release_date,
                          year=as.numeric(as.character(str_sub(release_date,1,4))),
                          month=as.numeric(as.character(str_sub(release_date,6,7))),
                          day=as.numeric(as.character(str_sub(release_date,9,10))))

par <- list(logitP_cod = rep(0,length(unique(dat$length.class[dat$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dat$length.class[dat$species=="skrubbe"]))),
            logitTag.diff = 0,logitYear.diff=0,
            logitRelease.diff = rep(0,n_release),logReleaseSD = 0,
            logitRefind.eff = 0)


nll <- function(par){
  getAll(par, dat)
  
  releaseSD <- exp(logReleaseSD)
  found <- OBS(found)
  n <- OBS(n)
  
  c_sizes <- unique(length.class[species=="torsk"])
  f_sizes <- unique(length.class[species=="skrubbe"])
  release_date <- sort(unique(release))
  n_time <- length(release_date)
  
  jnll <- -sum(dbinom(feed.exp$found,feed.exp$n,plogis(c(logitRefind.eff,logitRefind.eff+logitTag.diff)),log=TRUE))
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- (0:1)[unique(PIT)==PIT[i]]
    year.idx <- (0:1)[unique(year)==year[i]]
    time.idx <- (1:n_time)[release_date==release[i]]
    
    
    
    diff.PIT <- 0
    diff.year <- 0
    diff.time <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    if(year.idx>0) diff.year <- logitYear.diff[year.idx]
    diff.time <- logitRelease.diff[time.idx]
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.year+diff.time)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.year+diff.time)
    }
  }
  # index for which sd to use for random effect, different by year
  #sd.idx <- ifelse(str_sub(unique(release),1,4)=="2022",1,2)
  jnll <- jnll-sum(dnorm(logitRelease.diff,mean=0,sd=releaseSD,log=TRUE))
  jnll <- jnll -sum(dbinom(found,n,p,log=TRUE))
  
  cod22 <- log(plogis(logitP_cod)/plogis(logitRefind.eff))
  cod24 <- log(plogis(logitP_cod+logitYear.diff)/plogis(logitRefind.eff))
  
  flounder22 <- log(plogis(logitP_flounder)/plogis(logitRefind.eff))
  flounder24 <- log(plogis(logitP_flounder+logitYear.diff)/plogis(logitRefind.eff))
  
  cod_total <- (cod22+cod24)/2
  flounder_total <- (flounder22+flounder24)/2
  
  idx.22 <- which(str_sub(release_date,1,4)=="2022")
  idx.24 <- which(str_sub(release_date,1,4)=="2024")
  cod22_1 <- log(plogis(logitP_cod[1]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  cod22_2 <- log(plogis(logitP_cod[2]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  cod22_3 <- log(plogis(logitP_cod[3]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  cod22_4 <- log(plogis(logitP_cod[4]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  cod22_5 <- log(plogis(logitP_cod[5]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  cod_rel22 <- c(cod22_1+cod22_2+cod22_3+cod22_4+cod22_5)/5
  
  flounder22_1 <- log(plogis(logitP_flounder[1]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  flounder22_2 <- log(plogis(logitP_flounder[2]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  flounder22_3 <- log(plogis(logitP_flounder[3]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  flounder22_4 <- log(plogis(logitP_flounder[4]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  flounder22_5 <- log(plogis(logitP_flounder[5]+logitRelease.diff[idx.22])/plogis(logitRefind.eff))
  flounder_rel22 <- c(flounder22_1+flounder22_2+flounder22_3+flounder22_4+flounder22_5)/5
  
  cod24_1 <- log(plogis(logitP_cod[1]+logitYear.diff+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  cod24_2 <- log(plogis(logitP_cod[2]+logitYear.diff+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  cod24_3 <- log(plogis(logitP_cod[3]+logitYear.diff+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  cod24_4 <- log(plogis(logitP_cod[4]+logitYear.diff+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  cod24_5 <- log(plogis(logitP_cod[5]+logitYear.diff+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  cod_rel24 <- c(cod24_1+cod24_2+cod24_3+cod24_4+cod24_5)/5
  
  flounder24_1 <- log(plogis(logitP_flounder[1]+logitYear.diff+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  flounder24_2 <- log(plogis(logitP_flounder[2]+logitYear.diff+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  flounder24_3 <- log(plogis(logitP_flounder[3]+logitYear.diff+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  flounder24_4 <- log(plogis(logitP_flounder[4]+logitYear.diff+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  flounder24_5 <- log(plogis(logitP_flounder[5]+logitYear.diff+logitRelease.diff[idx.24])/plogis(logitRefind.eff))
  flounder_rel24 <- c(flounder24_1+flounder24_2+flounder24_3+flounder24_4+flounder24_5)/5
  
  ADREPORT(cod_rel22)
  ADREPORT(cod_rel24)
  ADREPORT(flounder_rel22)
  ADREPORT(flounder_rel24)
  
  ADREPORT(cod22)
  ADREPORT(cod24)
  ADREPORT(flounder22)
  ADREPORT(flounder24)
  ADREPORT(cod_total)
  ADREPORT(flounder_total)
  
  
  release_effect <- exp(logitRelease.diff)
  ADREPORT(release_effect)
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE, random=c("logitRelease.diff"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

M7 <- opt
# 1. Number of estimated parameters
k <- length(opt$par)
# 2. Log-likelihood is the negative of obj$fn at optimum
logL <- -obj$fn(opt$par)
# 3. Calculate AIC
AIC7 <- 2 * k - 2 * logL
print(AIC7)


res <- oneStepPredict(obj,method="oneStepGeneric",
                      discrete=TRUE)

res$n <- dat$n
res7 <- res
#par(mfrow=c(1,2))
qqnorm(res$residual)
qqline(res$residual)

rel <- as.list(sdr,"Est",report=TRUE)$release_effect
plot(1:28,rel)
lines(c(12.5,12.5),c(-5,5))
mean(rel[1:12])
mean(rel[13:28])


#####

#Prepare results
#####
est <- as.list(sdr, "Est",report=TRUE)

sd <- as.list(sdr, "Std",report=TRUE)

df_year <- data.frame(species = rep(c("cod","flounder"),each=5*2),
                      length_class=c(rep(factor(c("[11,15)","[15,19)","[19-23)","[23,27)","27+"),
                                                levels=c("[11,15)","[15,19)","[19-23)","[23,27)","27+")),2),
                                     rep(factor(c("[10-13)","[13-16)","[16-19)","[19-22)","22+"),
                                                levels=c("[10-13)","[13-16)","[16-19)","[19-22)","22+")),2)),
                      year = rep(rep(c(2022,2024),each=5),2),
                      est = c(est$cod22,est$cod24,est$flounder22,est$flounder24),
                      sd = c(sd$cod22,sd$cod24,sd$flounder22,sd$flounder24))

df_total <- data.frame(species = rep(c("cod","flounder"),each=5),
                       length_class=c(factor(c("[11,15)","[15,19)","[19-23)","[23,27)","27+"),
                                             levels=c("[11,15)","[15,19)","[19-23)","[23,27)","27+")),
                                      factor(c("[10-13)","[13-16)","[16-19)","[19-22)","22+"),
                                             levels=c("[10-13)","[13-16)","[16-19)","[19-22)","22+"))),
                       est = c(est$cod_total,est$flounder_total),
                       sd = c(sd$cod_total,sd$flounder_total))

dev.off()
par(mfrow=c(2,2))

plot(exp(est$cod_rel22),type='l',lwd=4,ylim=c(0,1.5), 
     ylab="Average predation",xlab="",main="Cod 2022")
lines(exp(est$cod_rel22+2*sd$cod_rel22),lty="dashed",lwd=2)
lines(exp(est$cod_rel22-2*sd$cod_rel22),lty="dashed",lwd=2)

plot(exp(est$cod_rel24),type='l',lwd=4,ylim=c(0,1.5),xlab=""
     , ylab="",main="Cod 2024")
lines(exp(est$cod_rel24+2*sd$cod_rel24),lty="dashed",lwd=2)
lines(exp(est$cod_rel24-2*sd$cod_rel24),lty="dashed",lwd=2)

plot(exp(est$flounder_rel22),type='l',lwd=4,ylim=c(0,1.5),xlab="Release event",
     ylab="Average predation",main="Flounder 2022")
lines(exp(est$flounder_rel22+2*sd$flounder_rel22),lty="dashed",lwd=2)
lines(exp(est$flounder_rel22-2*sd$flounder_rel22),lty="dashed",lwd=2)

plot(exp(est$flounder_rel24),type='l',lwd=4,ylim=c(0,1.5),xlab="Release event",
     ylab="",main="Flounder 2024")
lines(exp(est$flounder_rel24+2*sd$flounder_rel24),lty="dashed",lwd=2)
lines(exp(est$flounder_rel24-2*sd$flounder_rel24),lty="dashed",lwd=2)

png(filename="release.png",width = 480,height=480)
#####

#plot years
#####
pc22 <- ggplot(df_year %>% filter(species=="cod" & year==2022),aes(x = length_class, y = exp(est))) +
  geom_bar(stat = "identity", position="dodge", width = 0.5,fill="steelblue3") +
  geom_errorbar(aes(ymin = exp(est-sd), ymax = exp(est+sd)), 
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Cod")+
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.5)+
  labs(y = "predation proportion", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1.3),expand=c(0,0))+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "n",
                    legend.title = element_blank(),
                    plot.title = element_text(hjust = 0.5))+
  annotate("text", x = -Inf, y = Inf, label = "A", 
           hjust = -0.5, vjust = 1.5, size = 4.5)

pc24 <- ggplot(df_year %>% filter(species=="cod" & year==2024),aes(x = length_class, y = exp(est))) +
  geom_bar(stat = "identity", position="dodge", width = 0.5,fill="steelblue3") +
  geom_errorbar(aes(ymin = exp(est-sd), ymax = exp(est+sd)), 
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.5)+
  labs(y = "predation proportion", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1.3),expand=c(0,0))+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "n",
                    legend.title = element_blank(),
                    plot.title = element_text(hjust = 0.5))+
  annotate("text", x = -Inf, y = Inf, label = "B", 
           hjust = -0.5, vjust = 1.5, size = 4.5)

pf22 <- ggplot(df_year %>% filter(species=="flounder" & year==2022),aes(x = length_class, y = exp(est))) +
  geom_bar(stat = "identity", position="dodge", width = 0.5,fill="steelblue3") +
  geom_errorbar(aes(ymin = exp(est-sd), ymax = exp(est+sd)), 
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Flounder")+
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.5)+
  labs(y = "predation proportion", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1.3),expand=c(0,0))+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "n",
                    legend.title = element_blank(),
                    plot.title = element_text(hjust = 0.5))+
  annotate("text", x = -Inf, y = Inf, label = "C", 
           hjust = -0.5, vjust = 1.5, size = 4.5)

pf24 <- ggplot(df_year %>% filter(species=="flounder" & year==2024),aes(x = length_class, y = exp(est))) +
  geom_bar(stat = "identity", position="dodge", width = 0.5,fill="steelblue3") +
  geom_errorbar(aes(ymin = exp(est-sd), ymax = exp(est+sd)), 
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.5)+
  labs(y = "predation proportion", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1.3),expand=c(0,0))+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "n",
                    legend.title = element_blank(),
                    plot.title = element_text(hjust = 0.5))+
  annotate("text", x = -Inf, y = Inf, label = "D", 
           hjust = -0.5, vjust = 1.5, size = 4.5)




p <- ggarrange(pc22,pf22,pc24,pf24,ncol=2)
ggsave("predation_years.png",p,height=8,width=10,dpi=600)
#####

#plot total
#####
pt1 <- ggplot(df_total %>% filter(species=="cod"),aes(x = length_class, y = exp(est))) +
  geom_bar(stat = "identity", position="dodge", width = 0.5,fill="steelblue3") +
  geom_errorbar(aes(ymin = exp(est-sd), ymax = exp(est+sd)), 
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Cod")+
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.5)+
  labs(y = "predation proportion", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1.3),expand=c(0,0))+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "n",
                    legend.title = element_blank(),
                    plot.title = element_text(hjust = 0.5))+
  annotate("text", x = -Inf, y = Inf, label = "A", 
           hjust = -0.5, vjust = 1.5, size = 4.5)


pt2 <- ggplot(df_total %>% filter(species=="flounder"),aes(x = length_class, y = exp(est))) +
  geom_bar(stat = "identity", position="dodge", width = 0.5,fill="steelblue3") +
  geom_errorbar(aes(ymin = exp(est-sd), ymax = exp(est+sd)), 
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Flounder")+
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.5)+
  labs(y = "predation proportion", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1.3),expand=c(0,0))+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "n",
                    legend.title = element_blank(),
                    plot.title = element_text(hjust = 0.5),
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.line.y = element_blank())+
  annotate("text", x = -Inf, y = Inf, label = "B", 
           hjust = -0.5, vjust = 1.5, size = 4.5)


p <- ggarrange(pt1,pt2,ncol=2)
ggsave("predation_total.png",p,height=5,width=10,dpi = 600)
#####


# numbers for paper
#####

# av. cod in 22 and 24
c(exp(mean(df_year$est[1:5])),exp(mean(df_year$est[6:10])))
# av. flounder in 22 and 24
c(exp(mean(df_year$est[11:15])),exp(mean(df_year$est[16:20])))

# av. cod year difference
mean(exp(df_year$est[1:5])/exp(df_year$est[6:10]))
# av. flounder year difference
mean(exp(df_year$est[11:15])/exp(df_year$est[16:20]))

mean(exp(c(df_year$est[1:5],df_year$est[11:15]))/exp(c(df_year$est[6:10],df_year$est[16:20])))

# cod 2022 predation range
cat("cod 22 predation range: ",round(exp(c(min(df_year$est[1:5]),max(df_year$est[1:5])))*100,1),"%")
# cod 2024 predation range
cat("cod 24 predation range: ",round(exp(c(min(df_year$est[6:10]),max(df_year$est[6:10])))*100,1),"%")
# cod 2022 predation range
cat("flounder 22 predation range: ",round(exp(c(min(df_year$est[11:15]),max(df_year$est[11:15])))*100,1),"%")
# cod 2024 predation range
cat("cod 24 predation range: ",round(exp(c(min(df_year$est[16:20]),max(df_year$est[16:20])))*100,1),"%")

# Total experiment
#cod
cat("cod predation: ",round(exp(df_total$est[1:5])*100,1),"%")
#flounder
cat("flounder predation: ",round(exp(df_total$est[6:10])*100,1),"%")




# PIT odds
exp(sdr$par.fixed[11])-1

# Year odds
exp(sdr$par.fixed[12])-1

# Chance of detection, CoD
#CoD14
plogis(sdr$par.fixed[14])

odds14 <- plogis(sdr$par.fixed[14])/(1-plogis(sdr$par.fixed[14]))
odds23 <- odds14*exp(sdr$par.fixed[11])

#CoD23
odds23/(1+odds23)

# CI of model parameters
coefs <- summary(sdr)[1:14,]
cat("p_cod11-14 = ",plogis(coefs[1,1]),' [',plogis(c(coefs[1,1]-1.96*coefs[1,2],coefs[1,1]+1.96*coefs[1,2])),']')
cat("p_cod15-18 = ",plogis(coefs[2,1]),' [',plogis(c(coefs[2,1]-1.96*coefs[2,2],coefs[2,1]+1.96*coefs[2,2])),']')
cat("p_cod19-22 = ",plogis(coefs[3,1]),' [',plogis(c(coefs[3,1]-1.96*coefs[3,2],coefs[3,1]+1.96*coefs[3,2])),']')
cat("p_cod23-26 = ",plogis(coefs[4,1]),' [',plogis(c(coefs[4,1]-1.96*coefs[4,2],coefs[4,1]+1.96*coefs[4,2])),']')
cat("p_cod27+ = ",plogis(coefs[5,1]),' [',plogis(c(coefs[5,1]-1.96*coefs[5,2],coefs[5,1]+1.96*coefs[5,2])),']')


cat("p_flo10-12 = ",plogis(coefs[6,1]),' [',plogis(c(coefs[6,1]-1.96*coefs[6,2],coefs[6,1]+1.96*coefs[6,2])),']')
cat("p_flo13-15 = ",plogis(coefs[7,1]),' [',plogis(c(coefs[7,1]-1.96*coefs[7,2],coefs[7,1]+1.96*coefs[7,2])),']')
cat("p_flo16-18 = ",plogis(coefs[8,1]),' [',plogis(c(coefs[8,1]-1.96*coefs[8,2],coefs[8,1]+1.96*coefs[8,2])),']')
cat("p_flo19-21 = ",plogis(coefs[9,1]),' [',plogis(c(coefs[9,1]-1.96*coefs[9,2],coefs[9,1]+1.96*coefs[9,2])),']')
cat("p_flo22+ = ",plogis(coefs[10,1]),' [',plogis(c(coefs[10,1]-1.96*coefs[10,2],coefs[10,1]+1.96*coefs[10,2])),']')

cat("PIT23 = ",coefs[11,1],' [',c(coefs[11,1]-1.96*coefs[11,2],coefs[11,1]+1.96*coefs[11,2]),']')

cat("Year24 = ",coefs[12,1],' [',c(coefs[12,1]-1.96*coefs[12,2],coefs[12,1]+1.96*coefs[12,2]),']')

cat("sigma = ",exp(coefs[13,1]),' [',exp(c(coefs[13,1]-1.96*coefs[13,2],coefs[13,1]+1.96*coefs[13,2])),']')

cat("alpha = ",plogis(coefs[14,1]),' [',plogis(c(coefs[14,1]-1.96*coefs[14,2],coefs[14,1]+1.96*coefs[14,2])),']')

cat("not Kalvø = ",-0.1248690,' [',c(-0.1248690-1.96*0.1717503,-0.1248690+1.96*0.1717503),']')

#####


# data for paper - Get tagging.exp data frame
#######
tagging.exp$n <- 1
# cod
cod <- aggregate(found~length.class,data=tagging.exp %>% filter(species=="torsk"),FUN=sum)
cod$n <- aggregate(n~length.class,data=tagging.exp %>% filter(species=="torsk"),FUN=sum)$n
cod$p <- cod$found/cod$n
cod

#Flounder
flounder <- aggregate(found~length.class,data=tagging.exp %>% filter(species=="skrubbe"),FUN=sum)
flounder$n <- aggregate(n~length.class,data=tagging.exp %>% filter(species=="skrubbe"),FUN=sum)$n
flounder$p <- flounder$found/flounder$n
flounder

#PIT
PIT <- aggregate(found~PIT,data=tagging.exp,FUN=sum)
PIT$n <- aggregate(n~PIT,data=tagging.exp,FUN=sum)$n
PIT$p <- PIT$found/PIT$n
PIT


ggplot(data=tagging.exp %>% filter(species=="torsk"),aes(x=length.class,y=n,fill=PIT))+
  geom_bar(stat = "identity", width = 20)
ggplot(data=flounder %>% filter(species=="skrubbe"),aes(x=length.class,y=n,fill=PIT))+
  geom_bar(stat = "identity", width = 20)
#year
year <- aggregate(found~year,data=tagging.exp,FUN=sum)
year$n <- aggregate(n~year,data=tagging.exp,FUN=sum)$n
year$p <- year$found/year$n
year

#release
release <- aggregate(found~release.loca+year,data=tagging.exp,FUN=sum)
release$n <- aggregate(n~release.loca+year,data=tagging.exp,FUN=sum)$n
release$p <- release$found/release$n
release

#release
release <- aggregate(found~release.loca,data=tagging.exp,FUN=sum)
release$n <- aggregate(n~release.loca,data=tagging.exp,FUN=sum)$n
release$p <- release$found/release$n
release


#feeding
feeding.exp$n <- 1
feeding <- aggregate(found~PIT,data=feeding.exp,FUN=sum)
feeding$n <- aggregate(n~PIT,data=feeding.exp,FUN=sum)$n
feeding$p <- feeding$found/feeding$n
feeding




#######
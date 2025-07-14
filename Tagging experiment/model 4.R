# Final model. ADDITIVE effects with refinding efficiency TRUE
#####

# create size classes
tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350

dexp <- aggregate(found~PIT+species+length.class+month+year,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class+month+year,data=tagging.exp,FUN = length)$found
dexp$time <- paste(dexp$year,dexp$month)

dFeed <- aggregate(found~PIT,data=feeding.exp,FUN = sum)
dFeed$n <- aggregate(found~PIT,data=feeding.exp,FUN = length)$found

dat <- as.list(dexp)
dat$feed.exp <- dFeed

par <- list(logitP_cod = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logitTag.diff = 0,logitTime.diff = rep(0,7),
            logitRefind.eff = 0)


nll <- function(par){
  getAll(par, dat)

  
  found <- OBS(found)
  n <- OBS(n)
  
  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  jnll <- -sum(dbinom(feed.exp$found,feed.exp$n,plogis(c(logitRefind.eff,logitRefind.eff+logitTag.diff)),log=TRUE))
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- (0:1)[unique(PIT)==PIT[i]]
    time.idx <- (0:7)[unique(time)==time[i]]

    
    diff.PIT <- 0
    diff.time <- 0

    if(PIT.idx>0) diff.PIT <- logitTag.diff
    if(time.idx>0) diff.time <- logitTime.diff[time.idx]

    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.time)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.time)
    }
  }
  jnll <- jnll -sum(dbinom(found,n,p,log=TRUE))

  codMarch22 <- log(plogis(logitP_cod)/plogis(logitRefind.eff))
  codApril22 <- log(plogis(logitP_cod+logitTime.diff[1])/plogis(logitRefind.eff))
  codMay22 <- log(plogis(logitP_cod+logitTime.diff[2])/plogis(logitRefind.eff))
  codJune22 <- log(plogis(logitP_cod+logitTime.diff[3])/plogis(logitRefind.eff))
  codMarch24 <- log(plogis(logitP_cod+logitTime.diff[4])/plogis(logitRefind.eff))
  codApril24 <- log(plogis(logitP_cod+logitTime.diff[5])/plogis(logitRefind.eff))
  codMay24 <- log(plogis(logitP_cod+logitTime.diff[6])/plogis(logitRefind.eff))
  codJune24 <- log(plogis(logitP_cod+logitTime.diff[7])/plogis(logitRefind.eff))
  
  flounderMarch22 <- log(plogis(logitP_flounder)/plogis(logitRefind.eff))
  flounderApril22 <- log(plogis(logitP_flounder+logitTime.diff[1])/plogis(logitRefind.eff))
  flounderMay22 <- log(plogis(logitP_flounder+logitTime.diff[2])/plogis(logitRefind.eff))
  flounderJune22 <- log(plogis(logitP_flounder+logitTime.diff[3])/plogis(logitRefind.eff))
  flounderMarch24 <- log(plogis(logitP_flounder+logitTime.diff[4])/plogis(logitRefind.eff))
  flounderApril24 <- log(plogis(logitP_flounder+logitTime.diff[5])/plogis(logitRefind.eff))
  flounderMay24 <- log(plogis(logitP_flounder+logitTime.diff[6])/plogis(logitRefind.eff))
  flounderJune24 <- log(plogis(logitP_flounder+logitTime.diff[7])/plogis(logitRefind.eff))
  
  cod22 <- (codMarch22+codApril22+codMay22+codJune22)/4
  cod24 <- (codMarch24+codApril24+codMay24+codJune24)/4
  flounder22 <- (flounderMarch22+flounderApril22+flounderMay22+flounderJune22)/4
  flounder24 <- (flounderMarch24+flounderApril24+flounderMay24+flounderJune24)/4
  
  cod_total <- (codMarch22+codApril22+codMay22+codJune22+codMarch24+codApril24+codMay24+codJune24)/8
  flounder_total <- (flounderMarch22+flounderApril22+flounderMay22+flounderJune22+flounderMarch24+flounderApril24+flounderMay24+flounderJune24)/8
  
  ADREPORT(codMarch22)
  ADREPORT(codApril22)
  ADREPORT(codMay22)
  ADREPORT(codJune22)
  ADREPORT(codMarch24)
  ADREPORT(codApril24)
  ADREPORT(codMay24)
  ADREPORT(codJune24)
  
  ADREPORT(flounderMarch22)
  ADREPORT(flounderApril22)
  ADREPORT(flounderMay22)
  ADREPORT(flounderJune22)
  ADREPORT(flounderMarch24)
  ADREPORT(flounderApril24)
  ADREPORT(flounderMay24)
  ADREPORT(flounderJune24)
  
  ADREPORT(cod22)
  ADREPORT(cod24)
  ADREPORT(flounder22)
  ADREPORT(flounder24)
  
  ADREPORT(cod_total)
  ADREPORT(flounder_total)
  
  ADREPORT(p)

  
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr
# 1. Number of estimated parameters
k <- length(opt$par)
# 2. Log-likelihood is the negative of obj$fn at optimum
logL <- -obj$fn(opt$par)
# 3. Calculate AIC
AIC <- 2 * k - 2 * logL
print(AIC)
563.988

res <- oneStepPredict(obj)
res$n <- dexp$n
qqnorm(res$residual)
qqline(res$residual)
#V_i <- res$n * res$mean * (1 - res$mean / res$n)
#dispersion <- sum((res$observation - res$mean)^2 / V_i) / (length(res$observation) - length(obj$par))
#paste("dispersion statitic =",round(dispersion,2)," i.e., underdispersed")

est <- as.list(sdr, "Est",report=TRUE)

sd <- as.list(sdr, "Std",report=TRUE)


df <- data.frame(species = rep(c("cod","flounder"),each=5*4*2),
                 month=rep(rep(c("March", "April","May","June"),each=5),4),
                 length_class=c(rep(factor(c("11-15","16-20","21-25","26-30",">30"),levels=c("11-15","16-20","21-25","26-30",">30")),8),
                                rep(factor(c("10-12","13-15","16-18","19-21",">21"),levels=c("10-12","13-15","16-18","19-21",">21")),8)),
                 year = rep(rep(c(2022,2024),each=20),2),
                 est = c(est[[1]],est[[2]],est[[3]],est[[4]],est[[5]],est[[6]],est[[7]],est[[8]],
                         est[[9]],est[[10]],est[[11]],est[[12]],est[[13]],est[[14]],est[[15]],est[[16]]),
                 sd = c(sd[[1]],sd[[2]],sd[[3]],sd[[4]],sd[[5]],sd[[6]],sd[[7]],sd[[8]],
                        sd[[9]],sd[[10]],sd[[11]],sd[[12]],sd[[13]],sd[[14]],sd[[15]],sd[[16]]))

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


# OBS tænk over størrelsesklasser of evt. +gruppe
# PIT size confidence interval
est.PIT <- as.list(sdr, "Est")$logitTag.diff
sd.PIT <- as.list(sdr, "Std")$logitTag.diff

z <- est.PIT / sd.PIT
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Pit size CI = [",round(exp(est.PIT-2*sd.PIT),2),",",round(exp(est.PIT+2*sd.PIT),2),"]")

paste("23 mm tags have",round(exp(est.PIT)*100-100,1),"% higher odds of recovery than 14 mm tags")


# time effects
est.time <- as.list(sdr, "Est")$logitTime.diff
sd.time <- as.list(sdr, "Std")$logitTime.diff
paste("April 2022 CI = [",round(exp(est.time[1]-1.96*sd.time[1]),3),",",round(exp(est.time[1]+1.96*sd.time[1]),3),"]")
paste("Release in low season has",round(exp(est.time[1])*100-100,1),"% lower odds of recovery than release in the high season")

paste("May 2022 CI = [",round(exp(est.time[2]-1.96*sd.time[2]),3),",",round(exp(est.time[2]+1.96*sd.time[2]),3),"]")
paste("Release in low season has",round(exp(est.time[2])*100-100,2),"% lower odds of recovery than release in the high season")

paste("June 2022 CI = [",round(exp(est.time[3]-1.96*sd.time[3]),3),",",round(exp(est.time[3]+1.96*sd.time[3]),3),"]")
paste("Release in low season has",round(exp(est.time[3])*100-100,2),"% lower odds of recovery than release in the high season")

paste("March 2024 CI = [",round(exp(est.time[4]-1.96*sd.time[4]),3),",",round(exp(est.time[4]+1.96*sd.time[4]),3),"]")
paste("Release in low season has",round(exp(est.time[4])*100-100,2),"% lower odds of recovery than release in the high season")

paste("April 2024 CI = [",round(exp(est.time[5]-1.96*sd.time[5]),3),",",round(exp(est.time[5]+1.96*sd.time[5]),3),"]")
paste("Release in low season has",round(exp(est.time[5])*100-100,2),"% lower odds of recovery than release in the high season")

paste("May 2024 CI = [",round(exp(est.time[6]-1.96*sd.time[6]),3),",",round(exp(est.time[6]+1.96*sd.time[6]),3),"]")
paste("Release in low season has",round(exp(est.time[6])*100-100,2),"% lower odds of recovery than release in the high season")

paste("June 2024 CI = [",round(exp(est.time[7]-1.96*sd.time[7]),3),",",round(exp(est.time[7]+1.96*sd.time[7]),3),"]")
paste("Release in low season has",round(exp(est.time[7])*100-100,2),"% lower odds of recovery than release in the high season")


#OBS måske undersøg om fodringsforsøg kan relateres til udsætningsmåned.
# Fodring i juni
#####


#plot year and month
#####

df$month <- factor(df$month,levels = c("March", "April", "May",   "June" ))
# plot function and CI for cod

p1 <- ggplot(df %>% filter(species=="cod" & year==2022),aes(x = length_class, y = exp(est),fill=month)) +
  geom_bar(stat = "identity", position="dodge", width = 0.5) +scale_fill_manual(values=c("March"="brown2","April"="orange","May"="green","June"="darkgreen"))+
  geom_errorbar(aes(ymin = exp(est-sd), ymax = exp(est+sd)), 
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Cod 2022")+
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
                    plot.title = element_text(hjust = 0.5))



p2 <- ggplot(df %>% filter(species=="cod" & year==2024),aes(x = length_class, y = exp(est),fill=month)) +
  geom_bar(stat = "identity", position="dodge", width = 0.5) +scale_fill_manual(values=c("March"="brown2","April"="orange","May"="green","June"="darkgreen"))+
  geom_errorbar(aes(ymin = exp(est-sd), ymax = exp(est+sd)),
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Cod 2024")+
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.5)+
  labs(y = "predation proportion", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1.3),expand=c(0,0))+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "n",
                    legend.title = element_blank(),
                    plot.title = element_text(hjust = 0.5))

p3 <- ggplot(df %>% filter(species=="flounder" & year==2022),aes(x = length_class, y = exp(est),fill=month)) +
  geom_bar(stat = "identity", position="dodge", width = 0.5) +scale_fill_manual(values=c("March"="brown2","April"="orange","May"="green","June"="darkgreen"))+
  geom_errorbar(aes(ymin = exp(est-sd), ymax = exp(est+sd)), 
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Flounder 2022")+
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
                    plot.title = element_text(hjust = 0.5))

p4 <- ggplot(df %>% filter(species=="flounder" & year==2024),aes(x = length_class, y = exp(est),fill=month)) +
  geom_bar(stat = "identity", position="dodge", width = 0.5) +scale_fill_manual(values=c("March"="brown2","April"="orange","May"="green","June"="darkgreen"))+
  geom_errorbar(aes(ymin = exp(est-sd), ymax = exp(est+sd)),
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Flounder 2024")+
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
                    legend.position = c(0.87, 0.5),
                    legend.title = element_blank(),
                    plot.title = element_text(hjust = 0.5))

p <- ggarrange(p1,p3,p2,p4,ncol=2)

ggsave("predation.png",p,height=8,width=10)

#####

#plot years
#####
pc22 <- ggplot(df_year %>% filter(species=="cod" & year==2022),aes(x = length_class, y = exp(est))) +
  geom_bar(stat = "identity", position="dodge", width = 0.5,fill="steelblue3") +
  geom_errorbar(aes(ymin = exp(est-sd), ymax = exp(est+sd)), 
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Cod 2022")+
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
                    plot.title = element_text(hjust = 0.5))

pc24 <- ggplot(df_year %>% filter(species=="cod" & year==2024),aes(x = length_class, y = exp(est))) +
  geom_bar(stat = "identity", position="dodge", width = 0.5,fill="steelblue3") +
  geom_errorbar(aes(ymin = exp(est-sd), ymax = exp(est+sd)), 
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Cod 2024")+
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.5)+
  labs(y = "predation proportion", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1.3),expand=c(0,0))+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "n",
                    legend.title = element_blank(),
                    plot.title = element_text(hjust = 0.5))

pf22 <- ggplot(df_year %>% filter(species=="flounder" & year==2022),aes(x = length_class, y = exp(est))) +
  geom_bar(stat = "identity", position="dodge", width = 0.5,fill="steelblue3") +
  geom_errorbar(aes(ymin = exp(est-sd), ymax = exp(est+sd)), 
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Flounder 2022")+
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
                    plot.title = element_text(hjust = 0.5))

pf24 <- ggplot(df_year %>% filter(species=="flounder" & year==2024),aes(x = length_class, y = exp(est))) +
  geom_bar(stat = "identity", position="dodge", width = 0.5,fill="steelblue3") +
  geom_errorbar(aes(ymin = exp(est-sd), ymax = exp(est+sd)), 
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Flounder 2024")+
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
                    plot.title = element_text(hjust = 0.5))




p <- ggarrange(pc22,pf22,pc24,pf24,ncol=2)
p
ggsave("predation_years.png",p,height=8,width=10)
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
                    plot.title = element_text(hjust = 0.5))


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
                    plot.title = element_text(hjust = 0.5))


p <- ggarrange(pt1,pt2,ncol=2)
p
ggsave("predation_total.png",p,height=6,width=10)
#####


# numbers for paper
#####
# av. cod 
exp(c(min(df_total$est[1:5]),max(df_total$est[1:5])))
# av. flounder 
exp(c(min(df_total$est[6:10]),max(df_total$est[6:10])))

# av. cod year difference
mean(exp(df_year$est[1:5])/exp(df_year$est[6:10]))
# av. flounder year difference
mean(exp(df_year$est[11:15])/exp(df_year$est[16:20]))

mean(exp(c(df_year$est[1:5],df_year$est[11:15]))/exp(c(df_year$est[6:10],df_year$est[16:20])))

# PIT odds
exp(sdr$par.fixed[11])
#####
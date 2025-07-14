# Final model. ADDITIVE effects with refinding efficiency TRUE
#####

# create size classes
tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350

dexp <- aggregate(found~PIT+species+length.class+year,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class+year,data=tagging.exp,FUN = length)$found

dFeed <- aggregate(found~PIT,data=feeding.exp,FUN = sum)
dFeed$n <- aggregate(found~PIT,data=feeding.exp,FUN = length)$found

dat <- as.list(dexp)
dat$feed.exp <- dFeed

par <- list(logitP_cod = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logitTag.diff = 0,logitYear.diff = 0,
            logitRefind.eff = 0)


nll <- function(par){
  getAll(par, dat)
  
  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  jnll <- -sum(dbinom(feed.exp$found,feed.exp$n,plogis(c(logitRefind.eff,logitRefind.eff+logitTag.diff)),log=TRUE))
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- (0:1)[unique(PIT)==PIT[i]]
    year.idx <- (0:1)[unique(year)==year[i]]
    
    
    diff.PIT <- 0
    diff.year <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    if(year.idx>0) diff.year <- logitYear.diff[year.idx]
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.year)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.year)
    }
  }
  
  cod22 <- log(plogis(logitP_cod)/plogis(logitRefind.eff))
  cod24 <- log(plogis(logitP_cod+logitYear.diff[1])/plogis(logitRefind.eff))
 
  
  flounder22 <- log(plogis(logitP_flounder)/plogis(logitRefind.eff))
  flounder24 <- log(plogis(logitP_flounder+logitYear.diff[1])/plogis(logitRefind.eff))
  
 
  cod_total <- (cod22+cod24)/2
  flounder_total <- (flounder22+flounder24)/2
  
  ADREPORT(cod22)
  ADREPORT(cod24)
  ADREPORT(flounder22)
  ADREPORT(flounder24)
  
  ADREPORT(cod_total)
  ADREPORT(flounder_total)
  
  found <- OBS(found)
  n <- OBS(n)
  
  jnll <- jnll -sum(dbinom(found,n,p,log=TRUE))
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

res <- oneStepPredict(obj)
res$n <- dexp$n


qqnorm(res$residual)
qqline(res$residual)

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

ggsave("predation_total.png",p,height=6,width=10)
#####

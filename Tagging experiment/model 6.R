d <- aggregate(found~species+month+year+length.class+PIT,data=dexp,FUN=sum) 
d$n <- aggregate(n~species+month+year+length.class+PIT,data=dexp,FUN=sum)$n 

d$p <- d$found/d$n
d$p[d$PIT=="14 mm"] <- plogis(logit(d$p[d$PIT=="14 mm"])+0.363)
d$found_cor <- (d$n*d$p)

df_year <- aggregate(found_cor~species+month+year+length.class,data=d,FUN=sum) 
df_year$n <- aggregate(n~species+month+year+length.class,data=d,FUN=sum)$n 
df_year$p <- df_year$found_cor/df_year$n
names(df_year)[names(df_year)=="length.class"] <- "length_class"

#plot year and month
#####

df_year$month <- factor(df_year$month,levels = c("marts", "april", "maj",   "juni" ))
# plot function and CI for cod

p1 <- ggplot(df_year %>% filter(species=="torsk" & year==2022),aes(x = length_class, y = p,fill=month)) +
  geom_bar(stat = "identity", position="dodge", width = 20) +scale_fill_manual(values=c("marts"="brown2","april"="orange","maj"="green","juni"="darkgreen"))+
  ggtitle("Cod 2022")+
  labs(y = "predation proportion", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1),expand=c(0,0))+
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



p2 <- ggplot(df_year %>% filter(species=="skrubbe" & year==2022),aes(x = length_class, y = p,fill=month)) +
  geom_bar(stat = "identity", position="dodge", width = 20) +scale_fill_manual(values=c("marts"="brown2","april"="orange","maj"="green","juni"="darkgreen"))+
  ggtitle("Flounder 2022")+
  labs(y = "predation proportion", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1),expand=c(0,0))+
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
p3 <- ggplot(df_year %>% filter(species=="torsk" & year==2024),aes(x = length_class, y = p,fill=month)) +
  geom_bar(stat = "identity", position="dodge", width = 20) +scale_fill_manual(values=c("marts"="brown2","april"="orange","maj"="green","juni"="darkgreen"))+
  ggtitle("Cod 2024")+
  labs(y = "predation proportion", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1),expand=c(0,0))+
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
p4 <- ggplot(df_year %>% filter(species=="skrubbe" & year==2024),aes(x = length_class, y = p,fill=month)) +
  geom_bar(stat = "identity", position="dodge", width = 20) +scale_fill_manual(values=c("marts"="brown2","april"="orange","maj"="green","juni"="darkgreen"))+
  ggtitle("Flounder2024")+
  labs(y = "predation proportion", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1),expand=c(0,0))+
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


p <- ggarrange(p1,p3,p2,p4,ncol=2)
#####

#model 6
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
            logitTag.diff = 0,logitTime.diff = rep(0,n_release),logTimeSD = 0,
            logitRefind.eff = 0)


nll <- function(par){
  getAll(par, dat)
  
  timeSD <- exp(logTimeSD)
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
    time.idx <- (1:n_time)[release_date==release[i]]
    
    
    
    diff.PIT <- 0
    diff.time <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    diff.time <- logitTime.diff[time.idx]
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.time)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.time)
    }
  }
  # index for which sd to use for random effect, different by year
  #sd.idx <- ifelse(str_sub(unique(release),1,4)=="2022",1,2)
  jnll <- jnll-sum(dnorm(logitTime.diff,mean=0,sd=timeSD,log=TRUE))
  jnll <- jnll -sum(dbinom(found,n,p,log=TRUE))
  
  release_effect <- exp(logitTime.diff)
  ADREPORT(release_effect)
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE, random=c("logitTime.diff"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

M6 <- opt
# 1. Number of estimated parameters
k <- length(opt$par)
# 2. Log-likelihood is the negative of obj$fn at optimum
logL <- -obj$fn(opt$par)
# 3. Calculate AIC
AIC6 <- 2 * k - 2 * logL
print(AIC6)

res <- oneStepPredict(obj,method="oneStepGeneric",
                      discrete=TRUE)
res$n <- dat$n
res6 <- res
qqnorm(res$residual)
qqline(res$residual)

rel <- as.list(sdr,"Est",report=TRUE)$release_effect
plot(1:28,rel)
lines(c(12.5,12.5),c(-5,5))
#####
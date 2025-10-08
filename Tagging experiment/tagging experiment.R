library("readxl")
library(RTMB)
library(dplyr)
library(LaplacesDemon)
library(ggplot2)
library(egg)

# read and prepare data
#####
tagged24 <- as.data.frame(read_excel("tagging 2024.xlsx", sheet = "Mærkede"))
found24 <- as.data.frame(read_excel("tagging 2024.xlsx", sheet = "Fundne"))

tagged22 <- as.data.frame(read_excel("tagging 2022.xlsx", sheet = "Mærkede"))
found22 <- as.data.frame(read_excel("tagging 2022.xlsx", sheet = "Fundne"))

taggedFeed <- as.data.frame(read_excel("feeding experiment.xlsx", sheet = "Mærkede"))
foundFeed <- as.data.frame(read_excel("feeding experiment.xlsx", sheet = "Fundne"))

df <- data.frame(PIT_nr = c(tagged24$PIT_nummer,tagged22$PIT_nummer,taggedFeed$PIT_nummer),
                 PIT = c(tagged24$PIT,tagged22$PIT,taggedFeed$PIT),
                 release.loca = c(tagged24$Udsætningslokalitet,tagged22$Udsætningslokalitet,taggedFeed$Udsætningslokalitet),
                 release = c(tagged24$Udsætningsdato,tagged22$Udsætningsdato,taggedFeed$Udsætningsdato),
                 species = c(tagged24$Art,tagged22$Art,taggedFeed$Art),
                 length = c(tagged24$Totallængde,tagged22$Totallængde,taggedFeed$Totallængde),
                 year = c(rep(2024,length(tagged24$PIT)),rep(2022,length(tagged22$PIT)),rep('feeding',length(taggedFeed$PIT))))

df$found <- 0
df$colony <- NA

for(i in 1:length(df$PIT)){
  year <- df$year[i]
  PIT <- df$PIT_nr[i]
  
  if(year=='2022'){
    if(PIT %in% found22$PIT_nummer){
      df$found[i] <- 1
      df$colony[i] <- found22$Skannet[found22$PIT_nummer==PIT]
    }
  }
  if(year=='2024'){
    if(PIT %in% found24$PIT_nummer){
      df$found[i] <- 1
      df$colony[i] <- found24$Skannet[found24$PIT_nummer==PIT]
    }
  }
  if(year=='feeding'){
    if(PIT %in% foundFeed$PIT_nummer){
      df$found[i] <- 1
      df$colony[i] <- foundFeed$Skannet[foundFeed$PIT_nummer==PIT]
    }
  }
}
df$species[df$species=="Torsk"] <- "torsk"
df$species[df$species!="torsk"] <- "skrubbe"

tagging.exp <- df %>% filter(year!='feeding')
feeding.exp <- df %>% filter(year=='feeding')

tagging.exp$month <- factor(format(tagging.exp$release, "%B"),levels = c("marts","april","maj","juni"))

# create size classes
tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350

dexp <- aggregate(found~PIT+species+length.class+release.loca+month+year,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class+release.loca+month+year,data=tagging.exp,FUN = length)$found
dexp$time <- paste(dexp$year,dexp$month)

dFeed <- aggregate(found~PIT,data=feeding.exp,FUN = sum)
dFeed$n <- aggregate(found~PIT,data=feeding.exp,FUN = length)$found

dat <- as.list(dexp)
dat$feed.exp <- dFeed
rm(list=setdiff(ls(),c('dat')))
#####

# compare feeding experiments
#####

dFeed <- aggregate(found~PIT+release,data=feeding.exp,FUN = sum)
dFeed$n <- aggregate(found~PIT+release,data=feeding.exp,FUN = length)$found

par <- list(logitp_14=0,logitp_23=0,logitp_14diff=0,logitp_23diff=0)

nll <- function(par){
  getAll(par, dFeed)
  p <- c(0,0,0,0)
  p[1] <- plogis(logitp_14)
  p[2] <- plogis(logitp_23)
  p[3] <- plogis(logitp_14+logitp_14diff)
  p[4] <- plogis(logitp_23+logitp_23diff)
  
  ADREPORT(logitp_14diff)
  ADREPORT(logitp_23diff)
  
  -sum(dbinom(found, n, p, log=TRUE))
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

# CI of logitp_14diff
est <- as.list(sdr, "Est", report=TRUE)$logitp_14diff
sd <- as.list(sdr, "Std", report=TRUE)$logitp_14diff
paste("CI of 14 mm tag difference between experiments")
round(c(est-2*sd,est+2*sd),2) # CI


# CI of logitp_23diff
est <- as.list(sdr, "Est", report=TRUE)$logitp_23diff
sd <- as.list(sdr, "Std", report=TRUE)$logitp_23diff
paste("CI of 23 mm tag difference between experiments")
round(c(est-2*sd,est+2*sd),2) # CI

# No difference between feeding experiments for either tag
#####

# PIT difference
#####

# create size classes
tagging.exp$length.class <- ceiling(tagging.exp$length/50)*50

dexp <- aggregate(found~PIT+species+length.class,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class,data=tagging.exp,FUN = length)$found

dexp$length.class[dexp$length.class<150] <- 150

dat <- as.list(dexp)
#dat <- as.list(dexp %>% filter(species=="torsk" & PIT=="14 mm"))

par <- list(logitP_cod = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logTag.diff = 0)


nll <- function(par){
  getAll(par, dat)
  p_cod <- plogis(logitP_cod)
  p_flounder <- plogis(logitP_flounder)
  tag.diff <- exp(logTag.diff)
  
  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    if(species[i]=='torsk'){
      if(PIT[i]=="14 mm") p[i] <- p_cod[c_sizes==length.class[i]]
      if(PIT[i]=="23 mm") p[i] <- p_cod[c_sizes==length.class[i]]*tag.diff
    }
    if(species[i]=='skrubbe'){
        if(PIT[i]=="14 mm") p[i] <- p_flounder[f_sizes==length.class[i]]
        if(PIT[i]=="23 mm") p[i] <- p_flounder[f_sizes==length.class[i]]*tag.diff
    }
  }
  nll <- -sum(dbinom(found,n,p,log=TRUE))
  nll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

# plot function and CI for cod
est.cod <- as.list(sdr, "Est")$logitP_cod
sd.cod <- as.list(sdr, "Std")$logitP_cod

df.cod <- data.frame(est = plogis(est.cod), 
                     size_class = unique(dexp$length.class[dexp$species=="torsk"]),
                     CI.up = plogis(est.cod  +2*sd.cod),
                     CI.low = plogis(est.cod  -2*sd.cod),
                     n = aggregate(n~length.class,data=dexp %>% filter(species=="torsk"),FUN=sum)$n)

ggplot(df.cod,aes(x = size_class, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 20) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +
  geom_text(aes(x=size_class,y=1.05,label = paste("n =",n)))+
  labs(y = "Mean value ± CI", x = "mm") +
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())

# plot function and CI for flounder
est.flounder <- as.list(sdr, "Est")$logitP_flounder
sd.flounder <- as.list(sdr, "Std")$logitP_flounder

df.flounder <- data.frame(est = plogis(est.flounder), 
                     size_class = unique(dexp$length.class[dexp$species=="skrubbe"]),
                     CI.up = plogis(est.flounder  +2*sd.flounder),
                     CI.low = plogis(est.flounder  -2*sd.flounder),
                     n = aggregate(n~length.class,data=dexp %>% filter(species=="skrubbe"),FUN=sum)$n)

ggplot(df.flounder,aes(x = size_class, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 20) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +
  geom_text(aes(x=size_class,y=1.05,label = paste("n =",n)))+
  labs(y = "Mean value ± CI", x = "mm") +
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())


# PIT size confidence interval

est.PIT <- as.list(sdr, "Est")$logTag.diff
sd.PIT <- as.list(sdr, "Std")$logTag.diff

(CV <- exp(c(est.PIT-2*sd.PIT,est.PIT+2*sd.PIT)))

# 23mm tags are significantly easier to find.
paste("23 mm tags have",round(exp(est.PIT)*100-100,1),"% higher finding probability than 14 mm tags")

#####

# Various variables - 1 cm bins
#####
# create size classes
#tagging.exp$found[which(tagging.exp$colony!="skarvkoloni-Hopsø")] <- 0

tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/10)*10  # 10 
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/10)*10 # 40

tagging.exp$length.class[which(tagging.exp$length.class>240 & tagging.exp$species=="skrubbe")] <- 250  # create >240mm group
tagging.exp$length.class[which(tagging.exp$length.class<130 & tagging.exp$species=="skrubbe")] <- 120 # create <130mm group

tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 310 # create >300mm group
tagging.exp$length.class[which(tagging.exp$length.class<140 & tagging.exp$species=="torsk")] <- 130 # create <140mm group



dexp <- aggregate(found~PIT+species+length.class+month+release.loca+year,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class+month+release.loca+year,data=tagging.exp,FUN = length)$found




dat <- as.list(dexp)
#dat <- as.list(dexp %>% filter(species=="torsk" & PIT=="14 mm"))
par <- list(logitP_cod = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logTag.diff = 0,logYear.diff = 0,logLoca.diff=0,logMonth.diff = c(0,0,0))


nll <- function(par){
  getAll(par, dat)
  p_cod <- plogis(logitP_cod)
  p_flounder <- plogis(logitP_flounder)
  tag.diff <- exp(logTag.diff)
  year.diff <- exp(logYear.diff)
  loca.diff <- exp(logLoca.diff)
  month.diff <- exp(logMonth.diff)
  
  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- rep(0:1)[unique(PIT)==PIT[i]]
    month.idx <- rep(0:3)[unique(month)==month[i]]
    year.idx <- rep(0:1)[unique(year)==year[i]]
    loca.idx <- rep(0:2)[unique(release.loca)==release.loca[i]]
    
    fac.PIT <- 1
    fac.month <- 1
    fac.year <- 1
    fac.loca <- 1
    
    if(PIT.idx>0) fac.PIT <- tag.diff
    if(month.idx>0) fac.month <- month.diff[month.idx]
    if(year.idx>0) fac.year <- year.diff
    if(loca.idx>0) fac.loca <- loca.diff
    
    
    
    if(species[i]=='torsk'){
      p[i] <- p_cod[c_sizes==length.class[i]]*fac.PIT*fac.year*fac.loca*fac.month
    }
    if(species[i]=='skrubbe'){
      p[i] <- p_flounder[f_sizes==length.class[i]]*fac.PIT*fac.year*fac.loca*fac.month
    }
  }
  
  df <- data.frame(p = p, size=length.class,species=species,n=n,
                   p_cor = p*n)
  ag <- aggregate(p_cor~size+species,data=df,FUN=sum)
  ag.n <- aggregate(n~size+species,data=df,FUN=sum)
  ag$w.m <- ag$p_cor/ag.n$n 
  
  weighted.mean_C <- ag$w.m[ag$species=="torsk"]
  weighted.mean_F <- ag$w.m[ag$species=="skrubbe"]
  ADREPORT(weighted.mean_C)
  ADREPORT(weighted.mean_F)

  nll <- -sum(dbinom(found,n,p,log=TRUE))
  nll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

# plot function and CI for cod
est.cod <- as.list(sdr, "Est")$logitP_cod
sd.cod <- as.list(sdr, "Std")$logitP_cod

df.cod <- data.frame(est = plogis(est.cod), 
                     size_class = unique(dexp$length.class[dexp$species=="torsk"]),
                     CI.up = plogis(est.cod  +2*sd.cod),
                     CI.low = plogis(est.cod  -2*sd.cod),
                     n = aggregate(n~length.class,data=dexp %>% filter(species=="torsk"),FUN=sum)$n)

ggplot(df.cod,aes(x = size_class-5, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 8) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +scale_x_continuous(breaks = seq(100, 390, by = 10))+
  geom_text(aes(x=size_class-5,y=1.05,label = paste("n =",n)))+
  labs(y = "Mean value ± CI", x = "mm") + ggtitle("Cod")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())

plotwithci <- function(x, y, sd, col='darkred', alpha=.4, trans=function(x)x, lwd=2, add=FALSE,main="", xlab=deparse(substitute(x)), ylab=deparse(substitute(y))){  
  code <- col2rgb(col)[,1]/255
  tc <- rgb(code[1], code[2], code[3], alpha)  
  xx <- c(x, rev(x))
  yy <- trans(c(y-2*sd, rev(y+2*sd)))
  if(!add)plot(xx, yy, type='n', xlab=xlab, ylab=ylab,ylim=c(0,1),main=main)
  polygon(xx, yy, col=tc, border=NA)
  lines(x, trans(y), col=col, lwd=lwd)
}

plotwithci(df.cod$size_class, est.cod, sd.cod,trans=function(x)plogis(x), ylab="Tag refinding probability", col ="navy",main="Cod")


# plot function and CI for flounder
est.flounder <- as.list(sdr, "Est")$logitP_flounder
sd.flounder <- as.list(sdr, "Std")$logitP_flounder

df.flounder <- data.frame(est = plogis(est.flounder), 
                          size_class = unique(dexp$length.class[dexp$species=="skrubbe"]),
                          CI.up = plogis(est.flounder  +2*sd.flounder),
                          CI.low = plogis(est.flounder  -2*sd.flounder),
                          n = aggregate(n~length.class,data=dexp %>% filter(species=="skrubbe"),FUN=sum)$n)

ggplot(df.flounder,aes(x = size_class-5, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 8) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +scale_x_continuous(breaks = seq(80, 280, by = 10)) +
  geom_text(aes(x=size_class-5,y=1.05,label = paste("n =",n)))+
  labs(y = "Mean value ± CI", x = "mm") +ggtitle("Flounder")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())

plotwithci(df.flounder$size_class, est.flounder, sd.flounder,trans=function(x)plogis(x), ylab="Tag refinding probability", col ="navy",main="Flounder")
# OBS tænk over størrelsesklasser of evt. +gruppe
# PIT size confidence interval
est.PIT <- as.list(sdr, "Est")$logTag.diff
sd.PIT <- as.list(sdr, "Std")$logTag.diff

z <- est.PIT / sd.PIT
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Pit size CI = [",round(exp(est.PIT-2*sd.PIT),2),",",round(exp(est.PIT+2*sd.PIT),2),"]")

paste("23 mm tags have",round(exp(est.PIT)*100-100,1),"% higher finding probability than 14 mm tags")

# Year effect confidence interval
est.year <- as.list(sdr, "Est")$logYear.diff
sd.year <- as.list(sdr, "Std")$logYear.diff

z <- est.year / sd.year
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("year effect CI = [",round(exp(est.year-2*sd.year),2),",",round(exp(est.year+2*sd.year),2),"]")

paste("2024 has",-round(exp(est.year)*100-100,1),"% lower finding probability than 2022")

# location effect confidence interval
est.loca <- as.list(sdr, "Est")$logLoca.diff
sd.loca <- as.list(sdr, "Std")$logLoca.diff

z <- est.loca / sd.loca
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Release at sandvig or Sønderballe Hoved CI = [",round(exp(est.loca[1]-2*sd.loca[1]),2),",",round(exp(est.loca[1]+2*sd.loca[1]),2),"]")

#paste("Release at Sandvig or Sønderballe Hoved has",-round(exp(est.loca[1])*100-100,1),"% lower finding probability than release at Kalvø")
paste("Release at sandvig or Sønderballe Hoved is not significantly different from release at Kalvø")



# month effect confidence interval
est.month <- as.list(sdr, "Est")$logMonth.diff
sd.month <- as.list(sdr, "Std")$logMonth.diff

paste("Release in April CI = [",round(exp(est.month[1]-2*sd.month[1]),2),",",round(exp(est.month[1]+2*sd.month[1]),2),"]")
paste("Release in May CI = [",round(exp(est.month[2]-2*sd.month[2]),2),",",round(exp(est.month[2]+2*sd.month[2]),2),"]")
paste("Release in June CI = [",round(exp(est.month[3]-2*sd.month[3]),2),",",round(exp(est.month[3]+2*sd.month[3]),2),"]")



paste("Release in April has",round(exp(est.month[1])*100-100,1),"% higher finding probability than release in March")
paste("Release in May has",round(exp(est.month[2])*100-100,1),"% higher finding probability than release in March")
#paste("Release in June has",round(exp(est.month[3])*100-100,1),"% higher finding probability than release in March, overlapping CI")
paste("Release in June is not significantly different from release in March")
print("No significant difference between April, May, and June")
#OBS måske undersøg om fodringsforsøg kan relateres til udsætningsmåned.
# Fodring i juni
#####

# Various variables - 5 and 3 cm bins for cod and flounder, respectively. multiplicative effects
#####

# create size classes
#tagging.exp$found[which(tagging.exp$colony!="skarvkoloni-Hopsø")] <- 0

tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350



dexp <- aggregate(found~PIT+species+length.class+month+release.loca+year,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class+month+release.loca+year,data=tagging.exp,FUN = length)$found




dat <- as.list(dexp)
#dat <- as.list(dexp %>% filter(species=="torsk" & PIT=="14 mm"))
par <- list(logitP_cod = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logTag.diff = 0,logYear.diff = 0,logLoca.diff=0,logMonth.diff = c(0,0,0))


nll <- function(par){
  getAll(par, dat)
  p_cod <- plogis(logitP_cod)
  p_flounder <- plogis(logitP_flounder)
  tag.diff <- exp(logTag.diff)
  year.diff <- exp(logYear.diff)
  loca.diff <- exp(logLoca.diff)
  month.diff <- exp(logMonth.diff)
  
  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- rep(0:1)[unique(PIT)==PIT[i]]
    month.idx <- rep(0:3)[unique(month)==month[i]]
    year.idx <- rep(0:1)[unique(year)==year[i]]
    loca.idx <- rep(0:2)[unique(release.loca)==release.loca[i]]
    
    fac.PIT <- 1
    fac.month <- 1
    fac.year <- 1
    fac.loca <- 1
    
    if(PIT.idx>0) fac.PIT <- tag.diff
    if(month.idx>0) fac.month <- month.diff[month.idx]
    if(year.idx>0) fac.year <- year.diff
    if(loca.idx>0) fac.loca <- loca.diff
    
    
    
    if(species[i]=='torsk'){
      p[i] <- p_cod[c_sizes==length.class[i]]*fac.PIT*fac.year*fac.loca*fac.month
    }
    if(species[i]=='skrubbe'){
      p[i] <- p_flounder[f_sizes==length.class[i]]*fac.PIT*fac.year*fac.loca*fac.month
    }
  }
  
  nll <- -sum(dbinom(found,n,p,log=TRUE))
  nll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

# plot function and CI for cod
est.cod <- as.list(sdr, "Est")$logitP_cod
sd.cod <- as.list(sdr, "Std")$logitP_cod

df.cod <- data.frame(est = plogis(est.cod), 
                     size_class = unique(dexp$length.class[dexp$species=="torsk"]),
                     CI.up = plogis(est.cod  +2*sd.cod),
                     CI.low = plogis(est.cod  -2*sd.cod),
                     n = aggregate(n~length.class,data=dexp %>% filter(species=="torsk"),FUN=sum)$n)

ggplot(df.cod,aes(x = size_class-5, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 8) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +scale_x_continuous(breaks = seq(100, 390, by = 10))+
  geom_text(aes(x=size_class-5,y=1.05,label = paste("n =",n)))+
  labs(y = "Mean value ± CI", x = "mm") + ggtitle("Cod")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())



plotwithci <- function(x, y, sd, col='darkred', alpha=.4, trans=function(x)x, lwd=2, add=FALSE,main="", xlab=deparse(substitute(x)), ylab=deparse(substitute(y))){  
  code <- col2rgb(col)[,1]/255
  tc <- rgb(code[1], code[2], code[3], alpha)  
  xx <- c(x, rev(x))
  yy <- trans(c(y-2*sd, rev(y+2*sd)))
  if(!add)plot(xx, yy, type='n', xlab=xlab, ylab=ylab,ylim=c(0,1),main=main)
  polygon(xx, yy, col=tc, border=NA)
  lines(x, trans(y), col=col, lwd=lwd)
}

plotwithci(df.cod$size_class, est.cod, sd.cod,trans=function(x)plogis(x), ylab="Tag refinding probability", col ="navy",main="Cod")


# plot function and CI for flounder
est.flounder <- as.list(sdr, "Est")$logitP_flounder
sd.flounder <- as.list(sdr, "Std")$logitP_flounder

df.flounder <- data.frame(est = plogis(est.flounder), 
                          size_class = unique(dexp$length.class[dexp$species=="skrubbe"]),
                          CI.up = plogis(est.flounder  +2*sd.flounder),
                          CI.low = plogis(est.flounder  -2*sd.flounder),
                          n = aggregate(n~length.class,data=dexp %>% filter(species=="skrubbe"),FUN=sum)$n)

ggplot(df.flounder,aes(x = size_class-5, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 8) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +scale_x_continuous(breaks = seq(80, 280, by = 10)) +
  geom_text(aes(x=size_class-5,y=1.05,label = paste("n =",n)))+
  labs(y = "Mean value ± CI", x = "mm") +ggtitle("Flounder")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())

plotwithci(df.flounder$size_class, est.flounder, sd.flounder,trans=function(x)plogis(x), ylab="Tag refinding probability", col ="navy",main="Flounder")
# OBS tænk over størrelsesklasser of evt. +gruppe
# PIT size confidence interval
est.PIT <- as.list(sdr, "Est")$logTag.diff
sd.PIT <- as.list(sdr, "Std")$logTag.diff

z <- est.PIT / sd.PIT
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Pit size CI = [",round(exp(est.PIT-2*sd.PIT),2),",",round(exp(est.PIT+2*sd.PIT),2),"]")

paste("23 mm tags have",round(exp(est.PIT)*100-100,1),"% higher finding probability than 14 mm tags")

# Year effect confidence interval
est.year <- as.list(sdr, "Est")$logYear.diff
sd.year <- as.list(sdr, "Std")$logYear.diff

z <- est.year / sd.year
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("year effect CI = [",round(exp(est.year-2*sd.year),2),",",round(exp(est.year+2*sd.year),2),"]")

paste("2024 has",-round(exp(est.year)*100-100,1),"% lower finding probability than 2022")

# location effect confidence interval
est.loca <- as.list(sdr, "Est")$logLoca.diff
sd.loca <- as.list(sdr, "Std")$logLoca.diff

z <- est.loca / sd.loca
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Release at sandvig or Sønderballe Hoved CI = [",round(exp(est.loca[1]-2*sd.loca[1]),2),",",round(exp(est.loca[1]+2*sd.loca[1]),2),"]")

#paste("Release at Sandvig or Sønderballe Hoved has",-round(exp(est.loca[1])*100-100,1),"% lower finding probability than release at Kalvø")
paste("Release at sandvig or Sønderballe Hoved is not significantly different from release at Kalvø")



# month effect confidence interval
est.month <- as.list(sdr, "Est")$logMonth.diff
sd.month <- as.list(sdr, "Std")$logMonth.diff

paste("Release in April CI = [",round(exp(est.month[1]-2*sd.month[1]),2),",",round(exp(est.month[1]+2*sd.month[1]),2),"]")
paste("Release in May CI = [",round(exp(est.month[2]-2*sd.month[2]),2),",",round(exp(est.month[2]+2*sd.month[2]),2),"]")
paste("Release in June CI = [",round(exp(est.month[3]-2*sd.month[3]),2),",",round(exp(est.month[3]+2*sd.month[3]),2),"]")



paste("Release in April has",round(exp(est.month[1])*100-100,1),"% higher finding probability than release in March")
paste("Release in May has",round(exp(est.month[2])*100-100,1),"% higher finding probability than release in March")
#paste("Release in June has",round(exp(est.month[3])*100-100,1),"% higher finding probability than release in March, overlapping CI")
paste("Release in June is not significantly different from release in March")
print("No significant difference between April, May, and June")
#OBS måske undersøg om fodringsforsøg kan relateres til udsætningsmåned.
# Fodring i juni
#####

# Estimate predation  2022-24 mean
#####

# create size classes
#tagging.exp$found[which(tagging.exp$colony!="skarvkoloni-Hopsø")] <- 0



tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350


dexp <- aggregate(found~PIT+species+length.class,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class,data=tagging.exp,FUN = length)$found


dFeed <- aggregate(found~PIT,data=feeding.exp,FUN = sum)
dFeed$n <- aggregate(found~PIT,data=feeding.exp,FUN = length)$found

#dexp$length.class[dexp$length.class<150] <- 150

dat <- as.list(dexp)
dat$feed.exp <- dFeed
#dat <- as.list(dexp %>% filter(species=="torsk" & PIT=="14 mm"))
par <- list(logitP_cod = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logTag.diff = 0,logitRecap.eff = 0)


nll <- function(par){
  getAll(par, dat)
  p_cod <- plogis(logitP_cod)
  p_flounder <- plogis(logitP_flounder)
  tag.diff <- exp(logTag.diff)
  recap.eff <- plogis(logitRecap.eff)
  
  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- rep(0:1)[unique(PIT)==PIT[i]]
    
    fac.PIT <- 1
    
    
    if(PIT.idx>0) fac.PIT <- tag.diff
    if(species[i]=='torsk'){
      p[i] <- p_cod[c_sizes==length.class[i]]*fac.PIT*recap.eff
    }
    if(species[i]=='skrubbe'){
      p[i] <- p_flounder[f_sizes==length.class[i]]*fac.PIT*recap.eff
    }
  }
  jnll <- -sum(dbinom(feed.exp$found,feed.exp$n,c(recap.eff,recap.eff*tag.diff),log=TRUE))
  jnll <- jnll-sum(dbinom(found,n,p,log=TRUE))
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

# plot function and CI for cod
est.cod <- as.list(sdr, "Est")$logitP_cod
sd.cod <- as.list(sdr, "Std")$logitP_cod

df.cod <- data.frame(est = plogis(est.cod), 
                     size_class = unique(dexp$length.class[dexp$species=="torsk"]),
                     CI.up = plogis(est.cod  +2*sd.cod),
                     CI.low = plogis(est.cod  -2*sd.cod),
                     n = aggregate(n~length.class,data=dexp %>% filter(species=="torsk"),FUN=sum)$n)

cp <- ggplot(df.cod,aes(x = size_class-25, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 40) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +scale_x_continuous(breaks = seq(100, 400, by = 50))+
  geom_text(aes(x=size_class-25,y=1.05,label = paste("n =",n)))+
  labs(y = "eating probability", x = "fish length [mm]") + ggtitle("Cod")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())
cp
ggsave("cod.png",cp,width=7,height=6)
# plot function and CI for flounder
est.flounder <- as.list(sdr, "Est")$logitP_flounder
sd.flounder <- as.list(sdr, "Std")$logitP_flounder

df.flounder <- data.frame(est = plogis(est.flounder), 
                          size_class = unique(dexp$length.class[dexp$species=="skrubbe"]),
                          CI.up = plogis(est.flounder  +2*sd.flounder),
                          CI.low = plogis(est.flounder  -2*sd.flounder),
                          n = aggregate(n~length.class,data=dexp %>% filter(species=="skrubbe"),FUN=sum)$n)

fp <- ggplot(df.flounder,aes(x = size_class-15, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 25) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +scale_x_continuous(breaks = seq(90, 300, by = 30)) +
  geom_text(aes(x=size_class-15,y=1.05,label = paste("n =",n)))+
  labs(y = "eating probability", x = "fish length [mm]") +ggtitle("Flounder")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())

fp
ggsave("flounder.png",fp,width=7,height=6)
# PIT size confidence interval
est.recap <- as.list(sdr, "Est")$logitRecap.eff
sd.recap <- as.list(sdr, "Std")$logitRecap.eff

paste("Recapture CI = [",round(exp(est.recap-2*sd.recap),2),",",round(exp(est.recap+2*sd.recap),2),"]")

paste("The probability of recovering a 14 mm tag is",round(exp(est.recap)*100,1),"%")

est.PIT <- as.list(sdr, "Est")$logTag.diff
paste("The probability of recovering a 23 mm tag is",round(exp(est.recap)*100*exp(est.PIT),1),"%")
paste("OBS there might be something wrong with the recovery, because quite far from fractions from the feeding experiment")

#####

# Estimate predation 2022-24 separate
#####

# create size classes
#tagging.exp$found[which(tagging.exp$colony!="skarvkoloni-Hopsø")] <- 0



tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350


dexp <- aggregate(found~PIT+species+length.class+year,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class+year,data=tagging.exp,FUN = length)$found


dFeed <- aggregate(found~PIT,data=feeding.exp,FUN = sum)
dFeed$n <- aggregate(found~PIT,data=feeding.exp,FUN = length)$found

#dexp$length.class[dexp$length.class<150] <- 150

dat <- as.list(dexp)
dat$feed.exp <- dFeed
#dat <- as.list(dexp %>% filter(species=="torsk" & PIT=="14 mm"))
par <- list(logitP_cod22 = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder22 = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logitP_cod24 = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder24 = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logTag.diff = 0,logitRecap.eff = 0)


nll <- function(par){
  getAll(par, dat)
  p_cod22 <- plogis(logitP_cod22)
  p_flounder22 <- plogis(logitP_flounder22)
  p_cod24 <- plogis(logitP_cod24)
  p_flounder24 <- plogis(logitP_flounder24)
  tag.diff <- exp(logTag.diff)
  recap.eff <- plogis(logitRecap.eff)
  
  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- rep(0:1)[unique(PIT)==PIT[i]]
    
    fac.PIT <- 1
    
    
    if(PIT.idx>0) fac.PIT <- tag.diff
    
    if(year[i]==2022){
      if(species[i]=='torsk'){
        p[i] <- p_cod22[c_sizes==length.class[i]]*fac.PIT*recap.eff
      }
      if(species[i]=='skrubbe'){
        p[i] <- p_flounder22[f_sizes==length.class[i]]*fac.PIT*recap.eff
      }
    }
    if(year[i]==2024){
      if(species[i]=='torsk'){
        p[i] <- p_cod24[c_sizes==length.class[i]]*fac.PIT*recap.eff
      }
      if(species[i]=='skrubbe'){
        p[i] <- p_flounder24[f_sizes==length.class[i]]*fac.PIT*recap.eff
      }
    }
  }
  jnll <- -sum(dbinom(feed.exp$found,feed.exp$n,c(recap.eff,recap.eff*tag.diff),log=TRUE))
  jnll <- jnll-sum(dbinom(found,n,p,log=TRUE))
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr


# compute n pr. year
n.cod <- aggregate(n~length.class+year,data=dexp %>% filter(species=="torsk"),FUN=sum)
n.flounder <- aggregate(n~length.class+year,data=dexp %>% filter(species=="skrubbe"),FUN=sum)
# plot function and CI for cod 2022
est.cod22 <- as.list(sdr, "Est")$logitP_cod22
sd.cod22 <- as.list(sdr, "Std")$logitP_cod22

df.cod22 <- data.frame(est = plogis(est.cod22), 
                     size_class = unique(dexp$length.class[dexp$species=="torsk"]),
                     CI.up = plogis(est.cod22  +2*sd.cod22),
                     CI.low = plogis(est.cod22  -2*sd.cod22),
                     n = n.cod$n[n.cod$year==2022])

ggplot(df.cod22,aes(x = size_class-25, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 40) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +scale_x_continuous(breaks = seq(100, 400, by = 50))+
  geom_text(aes(x=size_class-25,y=1.05,label = paste("n =",n)))+
  labs(y = "Mean value ± CI", x = "mm") + ggtitle("Cod 2022")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())

# plot function and CI for cod 2024
est.cod24 <- as.list(sdr, "Est")$logitP_cod24
sd.cod24 <- as.list(sdr, "Std")$logitP_cod24

df.cod24 <- data.frame(est = plogis(est.cod24), 
                       size_class = unique(dexp$length.class[dexp$species=="torsk"]),
                       CI.up = plogis(est.cod24  +2*sd.cod24),
                       CI.low = plogis(est.cod24  -2*sd.cod24),
                       n = n.cod$n[n.cod$year==2024])

ggplot(df.cod24,aes(x = size_class-25, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 40) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +scale_x_continuous(breaks = seq(100, 400, by = 50))+
  geom_text(aes(x=size_class-25,y=1.05,label = paste("n =",n)))+
  labs(y = "Mean value ± CI", x = "mm") + ggtitle("Cod 2024")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())

# plot function and CI for flounder
est.flounder22 <- as.list(sdr, "Est")$logitP_flounder22
sd.flounder22 <- as.list(sdr, "Std")$logitP_flounder22

df.flounder22 <- data.frame(est = plogis(est.flounder22), 
                          size_class = unique(dexp$length.class[dexp$species=="skrubbe"]),
                          CI.up = plogis(est.flounder22  +2*sd.flounder22),
                          CI.low = plogis(est.flounder22  -2*sd.flounder22),
                          n = n.flounder$n[n.flounder$year==2022])

ggplot(df.flounder22,aes(x = size_class-15, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 25) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +scale_x_continuous(breaks = seq(90, 300, by = 30)) +
  geom_text(aes(x=size_class-15,y=1.05,label = paste("n =",n)))+
  labs(y = "Mean value ± CI", x = "mm") +ggtitle("Flounder 2022")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())

# plot function and CI for flounder
est.flounder24 <- as.list(sdr, "Est")$logitP_flounder24
sd.flounder24 <- as.list(sdr, "Std")$logitP_flounder24

df.flounder24 <- data.frame(est = plogis(est.flounder24), 
                          size_class = unique(dexp$length.class[dexp$species=="skrubbe"]),
                          CI.up = plogis(est.flounder24  +2*sd.flounder24),
                          CI.low = plogis(est.flounder24  -2*sd.flounder24),
                          n = n.flounder$n[n.flounder$year==2024])

ggplot(df.flounder24,aes(x = size_class-15, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 25) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +scale_x_continuous(breaks = seq(90, 300, by = 30)) +
  geom_text(aes(x=size_class-15,y=1.05,label = paste("n =",n)))+
  labs(y = "Mean value ± CI", x = "mm") +ggtitle("Flounder 2024")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())


# PIT size confidence interval
est.recap <- as.list(sdr, "Est")$logitRecap.eff
sd.recap <- as.list(sdr, "Std")$logitRecap.eff

paste("Recapture CI = [",round(exp(est.recap-2*sd.recap),2),",",round(exp(est.recap+2*sd.recap),2),"]")

paste("The probability of recovering a 14 mm tag is",round(exp(est.recap)*100,1),"%")

est.PIT <- as.list(sdr, "Est")$logTag.diff
paste("The probability of recovering a 23 mm tag is",round(exp(est.recap)*100*exp(est.PIT),1),"%")
paste("OBS there might be something wrong with the recovery, because quite far from fractions from the feeding experiment")

#####

# 2022-24 separate without feeding correction
#####

# create size classes
#tagging.exp$found[which(tagging.exp$colony!="skarvkoloni-Hopsø")] <- 0

tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350


dexp <- aggregate(found~PIT+species+length.class+year,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class+year,data=tagging.exp,FUN = length)$found

dat <- as.list(dexp)
#dat <- as.list(dexp %>% filter(species=="torsk" & PIT=="14 mm"))
par <- list(logitP_cod22 = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder22 = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logitP_cod24 = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder24 = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logTag.diff = 0)


nll <- function(par){
  getAll(par, dat)
  p_cod22 <- plogis(logitP_cod22)
  p_flounder22 <- plogis(logitP_flounder22)
  p_cod24 <- plogis(logitP_cod24)
  p_flounder24 <- plogis(logitP_flounder24)
  tag.diff <- exp(logTag.diff)

  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- rep(0:1)[unique(PIT)==PIT[i]]
    
    fac.PIT <- 1
    
    
    if(PIT.idx>0) fac.PIT <- tag.diff
    
    if(year[i]==2022){
      if(species[i]=='torsk'){
        p[i] <- p_cod22[c_sizes==length.class[i]]*fac.PIT
      }
      if(species[i]=='skrubbe'){
        p[i] <- p_flounder22[f_sizes==length.class[i]]*fac.PIT
      }
    }
    if(year[i]==2024){
      if(species[i]=='torsk'){
        p[i] <- p_cod24[c_sizes==length.class[i]]*fac.PIT
      }
      if(species[i]=='skrubbe'){
        p[i] <- p_flounder24[f_sizes==length.class[i]]*fac.PIT
      }
    }
  }
  jnll <- -sum(dbinom(found,n,p,log=TRUE))
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr


# compute n pr. year
n.cod <- aggregate(n~length.class+year,data=dexp %>% filter(species=="torsk"),FUN=sum)
n.flounder <- aggregate(n~length.class+year,data=dexp %>% filter(species=="skrubbe"),FUN=sum)
# plot function and CI for cod 2022
est.cod22 <- as.list(sdr, "Est")$logitP_cod22
sd.cod22 <- as.list(sdr, "Std")$logitP_cod22

df.cod22 <- data.frame(est = plogis(est.cod22), 
                       size_class = unique(dexp$length.class[dexp$species=="torsk"]),
                       CI.up = plogis(est.cod22  +2*sd.cod22),
                       CI.low = plogis(est.cod22  -2*sd.cod22),
                       n = n.cod$n[n.cod$year==2022])

ggplot(df.cod22,aes(x = size_class-25, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 40) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +scale_x_continuous(breaks = seq(100, 400, by = 50))+
  geom_text(aes(x=size_class-25,y=1.05,label = paste("n =",n)))+
  labs(y = "Mean value ± CI", x = "mm") + ggtitle("Cod 2022")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())


plotwithci <- function(x, y, sd, col='darkred', alpha=.4, trans=function(x)x, lwd=2, add=FALSE,main="", xlab=deparse(substitute(x)), ylab=deparse(substitute(y))){  
  code <- col2rgb(col)[,1]/255
  tc <- rgb(code[1], code[2], code[3], alpha)  
  xx <- c(x, rev(x))
  yy <- trans(c(y-2*sd, rev(y+2*sd)))
  if(!add)plot(xx, yy, type='n', xlab=xlab, ylab=ylab,ylim=c(0,1),main=main)
  polygon(xx, yy, col=tc, border=NA)
  lines(x, trans(y), col=col, lwd=lwd)
}

plotwithci(df.cod22$size_class, est.cod22, sd.cod22,trans=function(x)plogis(x), ylab="Tag refinding probability", col ="navy",main="Cod 2022")

# plot function and CI for cod 2024
est.cod24 <- as.list(sdr, "Est")$logitP_cod24
sd.cod24 <- as.list(sdr, "Std")$logitP_cod24

df.cod24 <- data.frame(est = plogis(est.cod24), 
                       size_class = unique(dexp$length.class[dexp$species=="torsk"]),
                       CI.up = plogis(est.cod24  +2*sd.cod24),
                       CI.low = plogis(est.cod24  -2*sd.cod24),
                       n = n.cod$n[n.cod$year==2024])

ggplot(df.cod24,aes(x = size_class-25, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 40) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +scale_x_continuous(breaks = seq(100, 400, by = 50))+
  geom_text(aes(x=size_class-25,y=1.05,label = paste("n =",n)))+
  labs(y = "Mean value ± CI", x = "mm") + ggtitle("Cod 2024")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())
plotwithci(df.cod24$size_class, est.cod24, sd.cod24,trans=function(x)plogis(x), ylab="Tag refinding probability", col ="darkred",main="Cod 2024")



# plot function and CI for flounder
est.flounder22 <- as.list(sdr, "Est")$logitP_flounder22
sd.flounder22 <- as.list(sdr, "Std")$logitP_flounder22

df.flounder22 <- data.frame(est = plogis(est.flounder22), 
                            size_class = unique(dexp$length.class[dexp$species=="skrubbe"]),
                            CI.up = plogis(est.flounder22  +2*sd.flounder22),
                            CI.low = plogis(est.flounder22  -2*sd.flounder22),
                            n = n.flounder$n[n.flounder$year==2022])

ggplot(df.flounder22,aes(x = size_class-15, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 25) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +scale_x_continuous(breaks = seq(90, 300, by = 30)) +
  geom_text(aes(x=size_class-15,y=1.05,label = paste("n =",n)))+
  labs(y = "Mean value ± CI", x = "mm") +ggtitle("Flounder 2022")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())
plotwithci(df.flounder22$size_class, est.flounder22, sd.flounder22,trans=function(x)plogis(x), ylab="Tag refinding probability", col ="navy",main="Flounder 2022")

# plot function and CI for flounder
est.flounder24 <- as.list(sdr, "Est")$logitP_flounder24
sd.flounder24 <- as.list(sdr, "Std")$logitP_flounder24

df.flounder24 <- data.frame(est = plogis(est.flounder24), 
                            size_class = unique(dexp$length.class[dexp$species=="skrubbe"]),
                            CI.up = plogis(est.flounder24  +2*sd.flounder24),
                            CI.low = plogis(est.flounder24  -2*sd.flounder24),
                            n = n.flounder$n[n.flounder$year==2024])

ggplot(df.flounder24,aes(x = size_class-15, y = est)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 25) +
  geom_errorbar(aes(ymin = CI.low, ymax = CI.up), 
                width = 5, linewidth = 0.8) +scale_x_continuous(breaks = seq(90, 300, by = 30)) +
  geom_text(aes(x=size_class-15,y=1.05,label = paste("n =",n)))+
  labs(y = "Mean value ± CI", x = "mm") +ggtitle("Flounder 2024")+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "top",
                    legend.title = element_blank())

plotwithci(df.flounder24$size_class, est.flounder24, sd.flounder24,trans=function(x)plogis(x), ylab="Tag refinding probability", col ="darkred",main="Flounder 2024")

# PIT size confidence interval
est.recap <- as.list(sdr, "Est")$logitRecap.eff
sd.recap <- as.list(sdr, "Std")$logitRecap.eff

paste("Recapture CI = [",round(exp(est.recap-2*sd.recap),2),",",round(exp(est.recap+2*sd.recap),2),"]")

paste("The probability of recovering a 14 mm tag is",round(exp(est.recap)*100,1),"%")

est.PIT <- as.list(sdr, "Est")$logTag.diff
paste("The probability of recovering a 23 mm tag is",round(exp(est.recap)*100*exp(est.PIT),1),"%")
paste("OBS there might be something wrong with the recovery, because quite far from fractions from the feeding experiment")

#####


# Various variables - 5 and 3 cm bins for cod and flounder, respectively. ADDITIVE effects
#####

# create size classes
#tagging.exp$found[which(tagging.exp$colony!="skarvkoloni-Hopsø")] <- 0

tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350
tagging.exp$release.loca[tagging.exp$release.loca %in% c("Sandvig", "Sønderballe Hoved" )] <- "non-kalvo"


dexp <- aggregate(found~PIT+species+length.class+month+release.loca+year,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class+month+release.loca+year,data=tagging.exp,FUN = length)$found

dat <- as.list(dexp)
#dat <- as.list(dexp %>% filter(species=="torsk" & PIT=="14 mm"))
str(dexp)
par <- list(logitP_cod = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logitTag.diff = 0,logitYear.diff = 0,logitLoca.diff=0,logitMonth.diff = c(0,0,0))


nll <- function(par){
  getAll(par, dat)

  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- rep(0:1)[unique(PIT)==PIT[i]]
    month.idx <- rep(0:3)[unique(month)==month[i]]
    year.idx <- rep(0:1)[unique(year)==year[i]]
    loca.idx <- rep(0:1)[unique(release.loca)==release.loca[i]]
    
    diff.PIT <- 0
    diff.month <- 0
    diff.year <- 0
    diff.loca <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    if(month.idx>0) diff.month <- logitMonth.diff[month.idx]
    if(year.idx>0) diff.year <- logitYear.diff
    if(loca.idx>0) diff.loca <- logitLoca.diff

    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.year+diff.loca+diff.month)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.year+diff.loca+diff.month)
    }
  }
  
  nll <- -sum(dbinom(found,n,p,log=TRUE))
  nll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr


# OBS tænk over størrelsesklasser of evt. +gruppe
# PIT size confidence interval
est.PIT <- as.list(sdr, "Est")$logitTag.diff
sd.PIT <- as.list(sdr, "Std")$logitTag.diff

z <- est.PIT / sd.PIT
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Pit size CI = [",round(exp(est.PIT-2*sd.PIT),2),",",round(exp(est.PIT+2*sd.PIT),2),"]")

paste("23 mm tags have",round(exp(est.PIT)*100-100,1),"% higher odds of recovery than 14 mm tags")

# Year effect confidence interval
est.year <- as.list(sdr, "Est")$logitYear.diff
sd.year <- as.list(sdr, "Std")$logitYear.diff

z <- est.year / sd.year
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("year effect CI = [",round(exp(est.year-2*sd.year),2),",",round(exp(est.year+2*sd.year),2),"]")

paste("2024 has",-round(exp(est.year)*100-100,1),"% lower odds of recovery than 2022")

# location effect confidence interval
est.loca <- as.list(sdr, "Est")$logitLoca.diff
sd.loca <- as.list(sdr, "Std")$logitLoca.diff

z <- est.loca / sd.loca
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Release at sandvig or Sønderballe Hoved CI = [",round(exp(est.loca[1]-2*sd.loca[1]),2),",",round(exp(est.loca[1]+2*sd.loca[1]),2),"]")

#paste("Release at Sandvig or Sønderballe Hoved has",-round(exp(est.loca[1])*100-100,1),"% lower finding probability than release at Kalvø")
paste("Release at sandvig or Sønderballe Hoved is not significantly different from release at Kalvø")



# month effect confidence interval
est.month <- as.list(sdr, "Est")$logitMonth.diff
sd.month <- as.list(sdr, "Std")$logitMonth.diff

paste("Release in April CI = [",round(exp(est.month[1]-1.96*sd.month[1]),3),",",round(exp(est.month[1]+1.96*sd.month[1]),3),"]")
paste("Release in May CI = [",round(exp(est.month[2]-1.96*sd.month[2]),3),",",round(exp(est.month[2]+1.96*sd.month[2]),3),"]")
paste("Release in June CI = [",round(exp(est.month[3]-1.96*sd.month[3]),3),",",round(exp(est.month[3]+1.96*sd.month[3]),3),"]")



paste("Release in April has",round(exp(est.month[1])*100-100,1),"% higher odds of recovery than release in March")
paste("Release in May has",round(exp(est.month[2])*100-100,1),"% higher odds of recovery than release in March")
paste("Release in June has",round(exp(est.month[3])*100-100,1),"% higher finding probability than release in March, overlapping CI")
paste("Release in June is not significantly different from release in March")
print("No significant difference between April, May, and June")
#OBS måske undersøg om fodringsforsøg kan relateres til udsætningsmåned.
# Fodring i juni
#####

# Hele svineriet. ADDITIVE effects without refinding efficiency
#####

# create size classes
#tagging.exp$found[which(tagging.exp$colony!="skarvkoloni-Hopsø")] <- 0

tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350
tagging.exp$release.loca[tagging.exp$release.loca %in% c("Sandvig", "Sønderballe Hoved" )] <- "non-kalvo"

dexp <- aggregate(found~PIT+species+length.class+month+release.loca+year,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class+month+release.loca+year,data=tagging.exp,FUN = length)$found

dat <- as.list(dexp)

par <- list(logitP_cod = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logitTag.diff = 0,logitYear.diff = 0,logitLoca.diff=0,logitMonth.diff = c(0,0,0))


nll <- function(par){
  getAll(par, dat)
  
  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- rep(0:1)[unique(PIT)==PIT[i]]
    month.idx <- rep(0:3)[unique(month)==month[i]]
    year.idx <- rep(0:1)[unique(year)==year[i]]
    loca.idx <- rep(0:1)[unique(release.loca)==release.loca[i]]
    
    diff.PIT <- 0
    diff.month <- 0
    diff.year <- 0
    diff.loca <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    if(month.idx>0) diff.month <- logitMonth.diff[month.idx]
    if(year.idx>0) diff.year <- logitYear.diff
    if(loca.idx>0) diff.loca <- logitLoca.diff
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.year+diff.loca+diff.month)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.year+diff.loca+diff.month)
    }
  }
  
  
  
  
  
  jnll <- -sum(dbinom(found,n,p,log=TRUE))
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr


# OBS tænk over størrelsesklasser of evt. +gruppe
# PIT size confidence interval
est.PIT <- as.list(sdr, "Est")$logitTag.diff
sd.PIT <- as.list(sdr, "Std")$logitTag.diff

z <- est.PIT / sd.PIT
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Pit size CI = [",round(exp(est.PIT-2*sd.PIT),2),",",round(exp(est.PIT+2*sd.PIT),2),"]")

paste("23 mm tags have",round(exp(est.PIT)*100-100,1),"% higher odds of recovery than 14 mm tags")

# Year effect confidence interval
est.year <- as.list(sdr, "Est")$logitYear.diff
sd.year <- as.list(sdr, "Std")$logitYear.diff

z <- est.year / sd.year
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("year effect CI = [",round(exp(est.year-2*sd.year),2),",",round(exp(est.year+2*sd.year),2),"]")

paste("2024 has",-round(exp(est.year)*100-100,1),"% lower odds of recovery than 2022")

# location effect confidence interval
est.loca <- as.list(sdr, "Est")$logitLoca.diff
sd.loca <- as.list(sdr, "Std")$logitLoca.diff

z <- est.loca / sd.loca
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Release at sandvig or Sønderballe Hoved CI = [",round(exp(est.loca[1]-2*sd.loca[1]),2),",",round(exp(est.loca[1]+2*sd.loca[1]),2),"]")

#paste("Release at Sandvig or Sønderballe Hoved has",-round(exp(est.loca[1])*100-100,1),"% lower finding probability than release at Kalvø")
paste("Release at sandvig or Sønderballe Hoved is not significantly different from release at Kalvø")



# month effect confidence interval
est.month <- as.list(sdr, "Est")$logitMonth.diff
sd.month <- as.list(sdr, "Std")$logitMonth.diff

paste("Release in April CI = [",round(exp(est.month[1]-1.96*sd.month[1]),3),",",round(exp(est.month[1]+1.96*sd.month[1]),3),"]")
paste("Release in May CI = [",round(exp(est.month[2]-1.96*sd.month[2]),3),",",round(exp(est.month[2]+1.96*sd.month[2]),3),"]")
paste("Release in June CI = [",round(exp(est.month[3]-1.96*sd.month[3]),3),",",round(exp(est.month[3]+1.96*sd.month[3]),3),"]")



paste("Release in April has",round(exp(est.month[1])*100-100,1),"% higher odds of recovery than release in March")
paste("Release in May has",round(exp(est.month[2])*100-100,1),"% higher odds of recovery than release in March")
paste("Release in June has",round(exp(est.month[3])*100-100,1),"% higher finding probability than release in March, overlapping CI")
paste("Release in June is not significantly different from release in March")
print("No significant difference between April, May, and June")
#OBS måske undersøg om fodringsforsøg kan relateres til udsætningsmåned.
# Fodring i juni
#####

# Without location. ADDITIVE effects without refinding efficiency
#####

# create size classes
#tagging.exp$found[which(tagging.exp$colony!="skarvkoloni-Hopsø")] <- 0

tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350

dexp <- aggregate(found~PIT+species+length.class+month+year,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class+month+year,data=tagging.exp,FUN = length)$found

dat <- as.list(dexp)

par <- list(logitP_cod = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logitTag.diff = 0,logitYear.diff = 0,logitMonth.diff = c(0,0,0))


nll <- function(par){
  getAll(par, dat)
  
  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- rep(0:1)[unique(PIT)==PIT[i]]
    month.idx <- rep(0:3)[unique(month)==month[i]]
    year.idx <- rep(0:1)[unique(year)==year[i]]

    diff.PIT <- 0
    diff.month <- 0
    diff.year <- 0

    if(PIT.idx>0) diff.PIT <- logitTag.diff
    if(month.idx>0) diff.month <- logitMonth.diff[month.idx]
    if(year.idx>0) diff.year <- logitYear.diff

    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.year+diff.month)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.year+diff.month)
    }
  }
  
  
  
  
  
  jnll <- -sum(dbinom(found,n,p,log=TRUE))
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

paste("AIC = ",round(2*length(opt$par)+2*obj$fn(opt$par)))

# OBS tænk over størrelsesklasser of evt. +gruppe
# PIT size confidence interval
est.PIT <- as.list(sdr, "Est")$logitTag.diff
sd.PIT <- as.list(sdr, "Std")$logitTag.diff

z <- est.PIT / sd.PIT
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Pit size CI = [",round(exp(est.PIT-2*sd.PIT),2),",",round(exp(est.PIT+2*sd.PIT),2),"]")

paste("23 mm tags have",round(exp(est.PIT)*100-100,1),"% higher odds of recovery than 14 mm tags")

# Year effect confidence interval
est.year <- as.list(sdr, "Est")$logitYear.diff
sd.year <- as.list(sdr, "Std")$logitYear.diff

z <- est.year / sd.year
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("year effect CI = [",round(exp(est.year-2*sd.year),2),",",round(exp(est.year+2*sd.year),2),"]")

paste("2024 has",-round(exp(est.year)*100-100,1),"% lower odds of recovery than 2022")

# location effect confidence interval
est.loca <- as.list(sdr, "Est")$logitLoca.diff
sd.loca <- as.list(sdr, "Std")$logitLoca.diff

z <- est.loca / sd.loca
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Release at sandvig or Sønderballe Hoved CI = [",round(exp(est.loca[1]-2*sd.loca[1]),2),",",round(exp(est.loca[1]+2*sd.loca[1]),2),"]")

#paste("Release at Sandvig or Sønderballe Hoved has",-round(exp(est.loca[1])*100-100,1),"% lower finding probability than release at Kalvø")
paste("Release at sandvig or Sønderballe Hoved is not significantly different from release at Kalvø")



# month effect confidence interval
est.month <- as.list(sdr, "Est")$logitMonth.diff
sd.month <- as.list(sdr, "Std")$logitMonth.diff

paste("Release in April CI = [",round(exp(est.month[1]-1.96*sd.month[1]),3),",",round(exp(est.month[1]+1.96*sd.month[1]),3),"]")
paste("Release in May CI = [",round(exp(est.month[2]-1.96*sd.month[2]),3),",",round(exp(est.month[2]+1.96*sd.month[2]),3),"]")
paste("Release in June CI = [",round(exp(est.month[3]-1.96*sd.month[3]),3),",",round(exp(est.month[3]+1.96*sd.month[3]),3),"]")



paste("Release in April has",round(exp(est.month[1])*100-100,1),"% higher odds of recovery than release in March")
paste("Release in May has",round(exp(est.month[2])*100-100,1),"% higher odds of recovery than release in March")
paste("Release in June has",round(exp(est.month[3])*100-100,1),"% higher finding probability than release in March, overlapping CI")
paste("Release in June is not significantly different from release in March")
print("No significant difference between April, May, and June")
#OBS måske undersøg om fodringsforsøg kan relateres til udsætningsmåned.
# Fodring i juni
#####

# Final model. ADDITIVE effects without refinding efficiency
#####

# create size classes
#tagging.exp$found[which(tagging.exp$colony!="skarvkoloni-Hopsø")] <- 0

tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350
tagging.exp$release.loca[tagging.exp$release.loca %in% c("Sandvig", "Sønderballe Hoved" )] <- "non-kalvo"

tagging.exp$season[tagging.exp$month %in% c("marts","juni")] <- "low season"
tagging.exp$season[tagging.exp$month %in% c("april","maj")] <- "high season"


dexp <- aggregate(found~PIT+species+length.class+season+year,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class+season+year,data=tagging.exp,FUN = length)$found



dat <- as.list(dexp)

par <- list(logitP_cod = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logitTag.diff = 0,logitYear.diff = 0,logitSeason.diff=0)


nll <- function(par){
  getAll(par, dat)
  
  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- rep(0:1)[unique(PIT)==PIT[i]]
    year.idx <- rep(0:1)[unique(year)==year[i]]
    season.idx <- rep(0:1)[unique(season)==season[i]]
    
    diff.PIT <- 0
    diff.season <- 0
    diff.year <- 0

    if(PIT.idx>0) diff.PIT <- logitTag.diff
    if(season.idx>0) diff.season <- logitSeason.diff[season.idx]
    if(year.idx>0) diff.year <- logitYear.diff

    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.year+diff.season)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.year+diff.season)
    }
  }
  
  found <- OBS(found)
  n <- OBS(n)
  
  jnll <- -sum(dbinom(found,n,p,log=TRUE))
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr


# OBS tænk over størrelsesklasser of evt. +gruppe
# PIT size confidence interval
est.PIT <- as.list(sdr, "Est")$logitTag.diff
sd.PIT <- as.list(sdr, "Std")$logitTag.diff

z <- est.PIT / sd.PIT
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Pit size CI = [",round(exp(est.PIT-2*sd.PIT),2),",",round(exp(est.PIT+2*sd.PIT),2),"]")

paste("23 mm tags have",round(exp(est.PIT)*100-100,1),"% higher odds of recovery than 14 mm tags")

# Year effect confidence interval
est.year <- as.list(sdr, "Est")$logitYear.diff
sd.year <- as.list(sdr, "Std")$logitYear.diff

z <- est.year / sd.year
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("year effect CI = [",round(exp(est.year-2*sd.year),2),",",round(exp(est.year+2*sd.year),2),"]")

paste("2024 has",-round(exp(est.year)*100-100,1),"% lower odds of recovery than 2022")

# location effect confidence interval
est.loca <- as.list(sdr, "Est")$logitLoca.diff
sd.loca <- as.list(sdr, "Std")$logitLoca.diff

z <- est.loca / sd.loca
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Release at sandvig or Sønderballe Hoved CI = [",round(exp(est.loca[1]-2*sd.loca[1]),2),",",round(exp(est.loca[1]+2*sd.loca[1]),2),"]")

#paste("Release at Sandvig or Sønderballe Hoved has",-round(exp(est.loca[1])*100-100,1),"% lower finding probability than release at Kalvø")
paste("Release at sandvig or Sønderballe Hoved is not significantly different from release at Kalvø")



# month effect confidence interval
est.month <- as.list(sdr, "Est")$logitMonth.diff
sd.month <- as.list(sdr, "Std")$logitMonth.diff

paste("Release in April CI = [",round(exp(est.month[1]-1.96*sd.month[1]),3),",",round(exp(est.month[1]+1.96*sd.month[1]),3),"]")
paste("Release in May CI = [",round(exp(est.month[2]-1.96*sd.month[2]),3),",",round(exp(est.month[2]+1.96*sd.month[2]),3),"]")
paste("Release in June CI = [",round(exp(est.month[3]-1.96*sd.month[3]),3),",",round(exp(est.month[3]+1.96*sd.month[3]),3),"]")



paste("Release in April has",round(exp(est.month[1])*100-100,1),"% higher odds of recovery than release in March")
paste("Release in May has",round(exp(est.month[2])*100-100,1),"% higher odds of recovery than release in March")
paste("Release in June has",round(exp(est.month[3])*100-100,1),"% higher finding probability than release in March, overlapping CI")
paste("Release in June is not significantly different from release in March")
print("No significant difference between April, May, and June")
#OBS måske undersøg om fodringsforsøg kan relateres til udsætningsmåned.
# Fodring i juni
#####

# Final model. ADDITIVE effects with refinding efficiency TRUE
#####

# create size classes
tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350

tagging.exp$season[tagging.exp$month %in% c("marts","juni")] <- "less nesting"
tagging.exp$season[tagging.exp$month %in% c("april","maj")] <- "peak nesting"
tagging.exp$season <- factor(tagging.exp$season,levels = c("less nesting" ,"peak nesting" ))


dexp <- aggregate(found~PIT+species+length.class+season+year,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class+season+year,data=tagging.exp,FUN = length)$found

dFeed <- aggregate(found~PIT,data=feeding.exp,FUN = sum)
dFeed$n <- aggregate(found~PIT,data=feeding.exp,FUN = length)$found

dat <- as.list(dexp)
dat$feed.exp <- dFeed

str(dexp)
par <- list(logitP_cod = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logitTag.diff = 0,logitYear.diff = 0,logitSeason.diff = 0,
            logitRefind.eff = 0)


nll <- function(par){
  getAll(par, dat)

  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  jnll <- -sum(dbinom(feed.exp$found,feed.exp$n,plogis(c(logitRefind.eff,logitRefind.eff+logitTag.diff)),log=TRUE))
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- rep(0:1)[unique(PIT)==PIT[i]]
    season.idx <- rep(0:1)[unique(season)==season[i]]
    year.idx <- rep(0:1)[unique(year)==year[i]]
    
    diff.PIT <- 0
    diff.season <- 0
    diff.year <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    if(season.idx>0) diff.season <- logitSeason.diff[season.idx]
    if(year.idx>0) diff.year <- logitYear.diff
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.year+diff.season)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.year+diff.season)
    }
  }
  codLow22 <- log(plogis(logitP_cod)/plogis(logitRefind.eff))
  codLow24 <- log(plogis(logitP_cod+logitYear.diff)/plogis(logitRefind.eff))
  codHigh22 <- log(plogis(logitP_cod+logitSeason.diff)/plogis(logitRefind.eff))
  codHigh24 <- log(plogis(logitP_cod+logitSeason.diff+logitYear.diff)/plogis(logitRefind.eff))
  
  flounderLow22 <- log(plogis(logitP_flounder)/plogis(logitRefind.eff))
  flounderLow24 <- log(plogis(logitP_flounder+logitYear.diff)/plogis(logitRefind.eff))
  flounderHigh22 <- log(plogis(logitP_flounder+logitSeason.diff)/plogis(logitRefind.eff))
  flounderHigh24 <- log(plogis(logitP_flounder+logitSeason.diff+logitYear.diff)/plogis(logitRefind.eff))
  
  ADREPORT(codLow22)
  ADREPORT(codLow24)
  ADREPORT(codHigh22)
  ADREPORT(codHigh24)
  ADREPORT(flounderLow22)
  ADREPORT(flounderLow24)
  ADREPORT(flounderHigh22)
  ADREPORT(flounderHigh24)
  
  jnll <- jnll -sum(dbinom(found,n,p,log=TRUE))
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

#cod
est.codLow22 <- as.list(sdr, "Est",report=TRUE)$codLow22
sd.codLow22 <- as.list(sdr, "Std",report=TRUE)$codLow22
est.codLow24 <- as.list(sdr, "Est",report=TRUE)$codLow24
sd.codLow24 <- as.list(sdr, "Std",report=TRUE)$codLow24
est.codHigh22 <- as.list(sdr, "Est",report=TRUE)$codHigh22
sd.codHigh22 <- as.list(sdr, "Std",report=TRUE)$codHigh22
est.codHigh24 <- as.list(sdr, "Est",report=TRUE)$codHigh24
sd.codHigh24 <- as.list(sdr, "Std",report=TRUE)$codHigh24
#flounder
est.flounderLow22 <- as.list(sdr, "Est",report=TRUE)$flounderLow22
sd.flounderLow22 <- as.list(sdr, "Std",report=TRUE)$flounderLow22
est.flounderLow24 <- as.list(sdr, "Est",report=TRUE)$flounderLow24
sd.flounderLow24 <- as.list(sdr, "Std",report=TRUE)$flounderLow24
est.flounderHigh22 <- as.list(sdr, "Est",report=TRUE)$flounderHigh22
sd.flounderHigh22 <- as.list(sdr, "Std",report=TRUE)$flounderHigh22
est.flounderHigh24 <- as.list(sdr, "Est",report=TRUE)$flounderHigh24
sd.flounderHigh24 <- as.list(sdr, "Std",report=TRUE)$flounderHigh24

df <- data.frame(species = rep(c("cod","flounder"),each=5*4),
                 season=rep(rep(c("less nesting" ,"peak nesting"),each=5),4),
                 length_class=c(rep(factor(c("11-15","16-20","21-25","26-30",">30"),levels=c("11-15","16-20","21-25","26-30",">30")),4),
                                rep(factor(c("10-12","13-15","16-18","19-21",">21"),levels=c("10-12","13-15","16-18","19-21",">21")),4)),
                 year = rep(rep(c(2022,2024),each=10),2),est = c(est.codLow22,est.codHigh22,est.codLow24,est.codHigh24,
                                                                 est.flounderLow22,est.flounderHigh22,est.flounderLow24,est.flounderHigh24),
                 se = c(sd.codLow22,sd.codHigh22,sd.codLow24,sd.codHigh24,
                         sd.flounderLow22,sd.flounderHigh22,sd.flounderLow24,sd.flounderHigh24))


# OBS tænk over størrelsesklasser of evt. +gruppe
# PIT size confidence interval
est.PIT <- as.list(sdr, "Est")$logitTag.diff
sd.PIT <- as.list(sdr, "Std")$logitTag.diff

z <- est.PIT / sd.PIT
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Pit size CI = [",round(exp(est.PIT-2*sd.PIT),2),",",round(exp(est.PIT+2*sd.PIT),2),"]")

paste("23 mm tags have",round(exp(est.PIT)*100-100,1),"% higher odds of recovery than 14 mm tags")

# Year effect confidence interval
est.year <- as.list(sdr, "Est")$logitYear.diff
sd.year <- as.list(sdr, "Std")$logitYear.diff

z <- est.year / sd.year
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("year effect CI = [",round(exp(est.year-2*sd.year),2),",",round(exp(est.year+2*sd.year),2),"]")

paste("2024 has",-round(exp(est.year)*100-100,1),"% lower odds of recovery than 2022")


# season effect confidence interval
est.season <- as.list(sdr, "Est")$logitSeason.diff
sd.season <- as.list(sdr, "Std")$logitSeason.diff

z <- est.season / sd.season
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Release in low season CI = [",round(exp(est.season[1]-1.96*sd.season[1]),3),",",round(exp(est.season[1]+1.96*sd.season[1]),3),"]")



paste("Release in low season has",-round(exp(est.season[1])*100-100,1),"% lower odds of recovery than release in the high season")

#OBS måske undersøg om fodringsforsøg kan relateres til udsætningsmåned.
# Fodring i juni
#####

# Feeding experiment - estimate refinding efficiency
#####

# create size classes

dFeed <- aggregate(found~PIT,data=feeding.exp,FUN = sum)
dFeed$n <- aggregate(found~PIT,data=feeding.exp,FUN = length)$found

dat <- list(feed.exp = dFeed)

par <- list(logitRefind.eff = 0)


nll <- function(par){
  getAll(par, dat)
  
  jnll <- -sum(dbinom(feed.exp$found,feed.exp$n,plogis(c(logitRefind.eff,logitRefind.eff+0.2782587)),log=TRUE))
  
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr
#####


Refinding_eff <- plogis(-0.852338)


# Final model. ADDITIVE effects
#####

# create size classes
tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350

tagging.exp$season[tagging.exp$month %in% c("marts","juni")] <- "low season"
tagging.exp$season[tagging.exp$month %in% c("april","maj")] <- "high season"
tagging.exp$season <- factor(tagging.exp$season,levels = c("low season" ,"high season" ))


dexp <- aggregate(found~PIT+species+length.class+season+year,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class+season+year,data=tagging.exp,FUN = length)$found

dat <- as.list(dexp)

par <- list(logitP_cod = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logitTag.diff = 0,logitYear.diff = 0,logitSeason.diff = 0)


nll <- function(par){
  getAll(par, dat)
  
  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- rep(0:1)[unique(PIT)==PIT[i]]
    season.idx <- rep(0:1)[unique(season)==season[i]]
    year.idx <- rep(0:1)[unique(year)==year[i]]
    
    diff.PIT <- 0
    diff.season <- 0
    diff.year <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    if(season.idx>0) diff.season <- logitSeason.diff[season.idx]
    if(year.idx>0) diff.year <- logitYear.diff
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.year+diff.season)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.year+diff.season)
    }
  }
  #codLow22 <- logitP_cod+diff.PIT+diff.year+diff.season+logitRecap.eff
  
  
  
  
  jnll <- -sum(dbinom(found,n,p,log=TRUE))
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

paste("AIC = ",round(2*length(opt$par)+2*obj$fn(opt$par)))

# OBS tænk over størrelsesklasser of evt. +gruppe
# PIT size confidence interval
est.PIT <- as.list(sdr, "Est")$logitTag.diff
sd.PIT <- as.list(sdr, "Std")$logitTag.diff

z <- est.PIT / sd.PIT
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Pit size CI = [",round(exp(est.PIT-2*sd.PIT),2),",",round(exp(est.PIT+2*sd.PIT),2),"]")

paste("23 mm tags have",round(exp(est.PIT)*100-100,1),"% higher odds of recovery than 14 mm tags")

# Year effect confidence interval
est.year <- as.list(sdr, "Est")$logitYear.diff
sd.year <- as.list(sdr, "Std")$logitYear.diff

z <- est.year / sd.year
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("year effect CI = [",round(exp(est.year-2*sd.year),2),",",round(exp(est.year+2*sd.year),2),"]")

paste("2024 has",-round(exp(est.year)*100-100,1),"% lower odds of recovery than 2022")


# season effect confidence interval
est.season <- as.list(sdr, "Est")$logitSeason.diff
sd.season <- as.list(sdr, "Std")$logitSeason.diff

paste("Release in low season CI = [",round(exp(est.season[1]-1.96*sd.season[1]),3),",",round(exp(est.season[1]+1.96*sd.season[1]),3),"]")



paste("Release in low season has",-round(exp(est.season[1])*100-100,1),"% lower odds of recovery than release in the high season")

#OBS måske undersøg om fodringsforsøg kan relateres til udsætningsmåned.
# Fodring i juni
#####

# Mean recapture with refinding eff. ADDITIVE effects
#####

# create size classes
tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350

tagging.exp$season[tagging.exp$month %in% c("marts","juni")] <- "low season"
tagging.exp$season[tagging.exp$month %in% c("april","maj")] <- "high season"

dexp <- aggregate(found~PIT+species+length.class,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class,data=tagging.exp,FUN = length)$found

dFeed <- aggregate(found~PIT,data=feeding.exp,FUN = sum)
dFeed$n <- aggregate(found~PIT,data=feeding.exp,FUN = length)$found

#dexp$length.class[dexp$length.class<150] <- 150

dat <- as.list(dexp)
dat$feed.exp <- dFeed

par <- list(logitP_cod = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logitTag.diff = 0,logitRecap.eff = 0)


nll <- function(par){
  getAll(par, dat)
  
  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  jnll <- -sum(dbinom(feed.exp$found,feed.exp$n,plogis(c(logitRecap.eff,logitRecap.eff+logitTag.diff)),log=TRUE))
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- rep(0:1)[unique(PIT)==PIT[i]]
   
    diff.PIT <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT)
    }
  }
  #codLow22 <- logitP_cod+diff.PIT+diff.year+diff.season+logitRecap.eff
  
  
  
  
  jnll <- jnll -sum(dbinom(found,n,p*plogis(logitRecap.eff),log=TRUE))
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

paste("AIC = ",round(2*length(opt$par)+2*obj$fn(opt$par)))

# OBS tænk over størrelsesklasser of evt. +gruppe
# PIT size confidence interval
est.PIT <- as.list(sdr, "Est")$logitTag.diff
sd.PIT <- as.list(sdr, "Std")$logitTag.diff

z <- est.PIT / sd.PIT
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Pit size CI = [",round(exp(est.PIT-2*sd.PIT),2),",",round(exp(est.PIT+2*sd.PIT),2),"]")

paste("23 mm tags have",round(exp(est.PIT)*100-100,1),"% higher odds of recovery than 14 mm tags")

# Year effect confidence interval
est.year <- as.list(sdr, "Est")$logitYear.diff
sd.year <- as.list(sdr, "Std")$logitYear.diff

z <- est.year / sd.year
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("year effect CI = [",round(exp(est.year-2*sd.year),2),",",round(exp(est.year+2*sd.year),2),"]")

paste("2024 has",-round(exp(est.year)*100-100,1),"% lower odds of recovery than 2022")


# season effect confidence interval
est.season <- as.list(sdr, "Est")$logitSeason.diff
sd.season <- as.list(sdr, "Std")$logitSeason.diff

paste("Release in low season CI = [",round(exp(est.season[1]-1.96*sd.season[1]),3),",",round(exp(est.season[1]+1.96*sd.season[1]),3),"]")



paste("Release in low season has",-round(exp(est.season[1])*100-100,1),"% lower odds of recovery than release in the high season")

#OBS måske undersøg om fodringsforsøg kan relateres til udsætningsmåned.
# Fodring i juni
#####


#plots
#####
# plot function and CI for cod

p1 <- ggplot(df %>% filter(species=="cod" & year==2022),aes(x = length_class, y = exp(est),fill=season)) +
  geom_bar(stat = "identity", position="dodge", width = 0.5) +scale_fill_manual(values=c("less nesting"="brown2","peak nesting"="orange"))+
  geom_errorbar(aes(ymin = exp(est-se), ymax = exp(est+se)), 
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Cod 2022")+
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.5)+
  labs(y = "eating probability", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1.3),expand=c(0,0))+
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


p2 <- ggplot(df %>% filter(species=="cod" & year==2024),aes(x = length_class, y = exp(est),fill=season)) +
  geom_bar(stat = "identity", position="dodge", width = 0.5) +scale_fill_manual(values=c("less nesting"="brown2","peak nesting"="orange"))+
  geom_errorbar(aes(ymin = exp(est-se), ymax = exp(est+se)),
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Cod 2024")+
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.5)+
  labs(y = "eating probability", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1.3),expand=c(0,0))+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = "n",
                    legend.title = element_blank(),
                    plot.title = element_text(hjust = 0.5))

p3 <- ggplot(df %>% filter(species=="flounder" & year==2022),aes(x = length_class, y = exp(est),fill=season)) +
  geom_bar(stat = "identity", position="dodge", width = 0.5) +scale_fill_manual(values=c("less nesting"="brown2","peak nesting"="orange"))+
  geom_errorbar(aes(ymin = exp(est-se), ymax = exp(est+se)), 
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Flounder 2022")+
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.5)+
  labs(y = "eating probability", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1.3),expand=c(0,0))+
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

p4 <- ggplot(df %>% filter(species=="flounder" & year==2024),aes(x = length_class, y = exp(est),fill=season)) +
  geom_bar(stat = "identity", position="dodge", width = 0.5) +scale_fill_manual(values=c("less nesting"="brown2","peak nesting"="orange"))+
  geom_errorbar(aes(ymin = exp(est-se), ymax = exp(est+se)),
                width = 0.3, linewidth = 0.8,position = position_dodge(width = 0.5)) +ggtitle("Flounder 2024")+
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.5)+
  labs(y = "eating probability", x = "fish length [cm]") + scale_y_continuous(limits = c(0,1.3),expand=c(0,0))+
  theme_bw()+ theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    text=element_text(color="black", size=16,family="serif"),
                    legend.position = c(0.8, 0.88),
                    legend.title = element_blank(),
                    plot.title = element_text(hjust = 0.5))

p <- ggarrange(p1,p3,p2,p4,ncol=2)
p
ggsave("predation.png",p,height=8,width=10)

#####



tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350

dexp <- aggregate(found~PIT+species+length.class+month+release.loca+year,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class+month+release.loca+year,data=tagging.exp,FUN = length)$found



dexp$length.class <- factor(dexp$length.class)
dexp$lcs <- paste(dexp$length.class,dexp$species) # length class and species index
cod <- dexp %>% filter(species=="torsk")
flounder <- dexp %>% filter(species=="skrubbe")

model <- glm(cbind(found, n - found) ~ lcs*year+month+PIT, data = dexp, family = binomial)
summary(model)

flounder_model <- glm(cbind(found, n - found) ~ length.class+year+month+PIT, data = flounder, family = binomial)
summary(flounder_model)

logLik <- -opt$objective  # TMB returns -log-likelihood, so negate it
k <- length(opt$par)      # Number of fixed effects
AIC <- 2 * k - 2 * logLik

# Other models (e.g., double sigmoidal) and test
#####
dexp$found[dexp$PIT=="23 mm"] <- dexp$found[dexp$PIT=="23 mm"]/1.256275
dexp$found <- dexp$found/0.2881289

dFeed$found/dFeed$n

dexp$found[dexp$PIT=="14 mm"] <- dexp$found[dexp$PIT=="14 mm"]/0.2375000
dexp$found[dexp$PIT=="23 mm"] <- dexp$found[dexp$PIT=="23 mm"]/0.4516129


t


df <- data.frame(p = dexp$found/dexp$n, size=dexp$length.class,species=dexp$species,n=dexp$n,
                 p_cor = dexp$found/dexp$n*dexp$n)
ag <- aggregate(p_cor~size+species,data=df,FUN=sum)
ag$n <- aggregate(n~size+species,data=df,FUN=sum)$n
ag$w.m <- ag$p_cor/ag$n 
ag$p <- aggregate(p~size+species,data=df,FUN=mean)$p

(weighted.mean_C <- ag$w.m[ag$species=="torsk"])
(weighted.mean_F <- ag$w.m[ag$species=="skrubbe"])

plot(ag$size[ag$species=="torsk"],weighted.mean_C,ylim=c(0,1),pch=19)
plot(ag$size[ag$species=="skrubbe"],weighted.mean_F,ylim=c(0,1),pch=19)


est.cod
est.flounder

cod <- dexp %>% filter(species=="torsk")
flounder <- dexp %>% filter(species=="skrubbe")
table(cod$PIT)
table(flounder$PIT)

# data from feeding experiment
dFeed <- aggregate(found~PIT,data=feeding.exp,FUN = sum)
dFeed$n <- aggregate(found~PIT,data=feeding.exp,FUN = length)$found
# create size classes
tagging.exp$length <- ceiling(tagging.exp$length/10)*10

dexp <- aggregate(found~PIT+species+length,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length,data=tagging.exp,FUN = length)$found



dat <- list(feed.exp=dFeed,tag.exp=dexp %>% filter(species=="torsk"))

double_sigmoid <- function(x, A, M1, M2,B1, B2) {
  A / (1 + exp(-(x - M1) / B1)) - A / (1 + exp(-(x - M2) / B2))
}


par <- list(logA=log(0.3),logM1=log(130),logM2.gt.M1=log(120),logR1=log(10),logR2=log(5))

# logitP = logit(p), p is probability of being eaten if you are released close to Hopsø
# logitP14 = logit(p14), p14 is the probability of detecting an eaten 14mm tag in Hopsø
# logP23diff._1 = log(p23diff-1). p23diff > 1 is a factor that describes how much better 23mm tags are than 14mm.
# p23 = p14*p23diff

nll <- function(par){
  getAll(par, dat)
  
  A <- exp(logA)
  M1 <- exp(logM1)
  M2 <- exp(logM1)+exp(logM2.gt.M1)
  r1 <- exp(logR1)
  r2 <- exp(logR2)
  
  jnll <- 0
  for(i in 1:length(tag.exp$n)){

    found <- tag.exp$found[i]
    n <- tag.exp$n[i]
    l <- tag.exp$length[i]
    
    p <- double_sigmoid(l,A,M1,M2,r1,r2)+10^-9
    jnll <- jnll-dbinom(found,n,p,log=TRUE)
  }
  x <- 50:350
  logFit <- log(double_sigmoid(x,A,M1,M2,r1,r2)+10^-9)
  ADREPORT(logFit)
  jnll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

# plot function and CI for cod
est <- as.list(sdr, "Est", report=TRUE)$logFit
sd <- as.list(sdr, "Std", report=TRUE)$logFit

CI.up <- exp(est  +2*sd) #upper CI
CI.low <- exp(est  -2*sd) #lower CI
c(CI.low,CI.up)

plot(50:350,exp(est),ylim = c(0,1),xlab='length [mm]',
     ylab = 'probability of eaten',main='Cod',type='l',lwd=2)
lines(50:350,CI.low,col='red')
lines(50:350,CI.up,col='red')
points(d$length,d$found/d$n,pch=19,col='darkblue')
#####

t23 <- dexp %>% filter(PIT=="23 mm")
t14 <- dexp %>% filter(PIT=="14 mm")

sum(t23$found)/sum(t23$n)
sum(t14$found)/sum(t14$n)

# test if seasonality is similar between years
#####

# create size classes
#tagging.exp$found[which(tagging.exp$colony!="skarvkoloni-Hopsø")] <- 0

tagging.exp$length.class[tagging.exp$species=="torsk"] <- ceiling(tagging.exp$length[tagging.exp$species=="torsk"]/50)*50
tagging.exp$length.class[tagging.exp$species=="skrubbe"] <- ceiling(tagging.exp$length[tagging.exp$species=="skrubbe"]/30)*30

tagging.exp$length.class[which(tagging.exp$length.class>210 & tagging.exp$species=="skrubbe")] <- 240
tagging.exp$length.class[which(tagging.exp$length.class>300 & tagging.exp$species=="torsk")] <- 350

tagging.exp <- tagging.exp %>% filter(year==2024)

dexp <- aggregate(found~PIT+species+length.class+month+release.loca,data=tagging.exp,FUN = sum)
dexp$n <- aggregate(found~PIT+species+length.class+month+release.loca,data=tagging.exp,FUN = length)$found

dat <- as.list(dexp)
#dat <- as.list(dexp %>% filter(species=="torsk" & PIT=="14 mm"))
str(dexp)
par <- list(logitP_cod = rep(0,length(unique(dexp$length.class[dexp$species=="torsk"]))),
            logitP_flounder = rep(0,length(unique(dexp$length.class[dexp$species=="skrubbe"]))),
            logitTag.diff = 0,logitLoca.diff=0,logitMonth.diff = c(0,0,0))


nll <- function(par){
  getAll(par, dat)
  
  c_sizes <- unique(dexp$length.class[dexp$species=="torsk"])
  f_sizes <- unique(dexp$length.class[dexp$species=="skrubbe"])
  
  p <- rep(0,length(n))
  for(i in 1:length(n)){
    PIT.idx <- rep(0:1)[unique(PIT)==PIT[i]]
    month.idx <- rep(0:3)[unique(month)==month[i]]
#    year.idx <- rep(0:1)[unique(year)==year[i]]
    loca.idx <- rep(0:2)[unique(release.loca)==release.loca[i]]
    
    diff.PIT <- 0
    diff.month <- 0
   # diff.year <- 0
    diff.loca <- 0
    
    if(PIT.idx>0) diff.PIT <- logitTag.diff
    if(month.idx>0) diff.month <- logitMonth.diff[month.idx]
  #  if(year.idx>0) diff.year <- logitYear.diff
    if(loca.idx>0) diff.loca <- logitLoca.diff
    
    if(species[i]=='torsk'){
      p[i] <- plogis(logitP_cod[c_sizes==length.class[i]]+diff.PIT+diff.loca+diff.month)
    }
    if(species[i]=='skrubbe'){
      p[i] <- plogis(logitP_flounder[f_sizes==length.class[i]]+diff.PIT+diff.loca+diff.month)
    }
  }
  
  nll <- -sum(dbinom(found,n,p,log=TRUE))
  nll
}
obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr


# OBS tænk over størrelsesklasser of evt. +gruppe
# PIT size confidence interval
est.PIT <- as.list(sdr, "Est")$logitTag.diff
sd.PIT <- as.list(sdr, "Std")$logitTag.diff

z <- est.PIT / sd.PIT
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Pit size CI = [",round(exp(est.PIT-2*sd.PIT),2),",",round(exp(est.PIT+2*sd.PIT),2),"]")

paste("23 mm tags have",round(exp(est.PIT)*100-100,1),"% higher odds of recovery than 14 mm tags")

# Year effect confidence interval
est.year <- as.list(sdr, "Est")$logitYear.diff
sd.year <- as.list(sdr, "Std")$logitYear.diff

z <- est.year / sd.year
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("year effect CI = [",round(exp(est.year-2*sd.year),2),",",round(exp(est.year+2*sd.year),2),"]")

paste("2024 has",-round(exp(est.year)*100-100,1),"% lower odds of recovery than 2022")

# location effect confidence interval
est.loca <- as.list(sdr, "Est")$logitLoca.diff
sd.loca <- as.list(sdr, "Std")$logitLoca.diff

z <- est.loca / sd.loca
(p_val <- 2 * (1 - pnorm(abs(z))))

paste("Release at sandvig or Sønderballe Hoved CI = [",round(exp(est.loca[1]-2*sd.loca[1]),2),",",round(exp(est.loca[1]+2*sd.loca[1]),2),"]")

#paste("Release at Sandvig or Sønderballe Hoved has",-round(exp(est.loca[1])*100-100,1),"% lower finding probability than release at Kalvø")
paste("Release at sandvig or Sønderballe Hoved is not significantly different from release at Kalvø")



# month effect confidence interval
est.month <- as.list(sdr, "Est")$logitMonth.diff
sd.month <- as.list(sdr, "Std")$logitMonth.diff

paste("Release in April CI = [",round(exp(est.month[1]-1.96*sd.month[1]),3),",",round(exp(est.month[1]+1.96*sd.month[1]),3),"]")
paste("Release in May CI = [",round(exp(est.month[2]-1.96*sd.month[2]),3),",",round(exp(est.month[2]+1.96*sd.month[2]),3),"]")
paste("Release in June CI = [",round(exp(est.month[3]-1.96*sd.month[3]),3),",",round(exp(est.month[3]+1.96*sd.month[3]),3),"]")



paste("Release in April has",round(exp(est.month[1])*100-100,1),"% higher odds of recovery than release in March")
paste("Release in May has",round(exp(est.month[2])*100-100,1),"% higher odds of recovery than release in March")
paste("Release in June has",round(exp(est.month[3])*100-100,1),"% higher finding probability than release in March, overlapping CI")
paste("Release in June is not significantly different from release in March")
print("No significant difference between April, May, and June")
#OBS måske undersøg om fodringsforsøg kan relateres til udsætningsmåned.
# Fodring i juni
#####

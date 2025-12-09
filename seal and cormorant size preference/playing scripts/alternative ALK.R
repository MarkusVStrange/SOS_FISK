# Run "Prepare data" from the 'Prepare DATRAS.R' file
rm(list=setdiff(ls(),c('ca','hl')))

# Estimates age distribution
ca.ag <- aggregate(CANoAtLngt~LngtClass+species+Age+Year+Quarter,data = ca,FUN=sum)
hl.ag <- aggregate(HLNoAtLngt~LngtClass+species+Year+Quarter,data = hl,FUN=sum)
hl.ag <- hl.ag %>% filter(LngtClass>5)
ca.ag$Age[which(ca.ag$species=="cod" & ca.ag$Age>7)] <- 7
ca.ag$Age[which(ca.ag$species=="herring" & ca.ag$Age>8)] <- 8
ca.ag$Age[which(ca.ag$species %in% c("flounder","plaice","dab") & ca.ag$Age>10)] <- 10
ca.ag$yqID <- paste(ca.ag$Year,ca.ag$Quarter)
ca.ag <- ca.ag %>% filter(Quarter %in% c(1,4))

hl_N <- (ca.ag %>% select(Quarter,Year,LngtClass,species,Age,CANoAtLngt))[1,]
names(hl_N)[6] <- "N"
# All species
yq.s <- unique(ca.ag$yqID)
#i=12
for(i in 1:length(yq.s)){
  year.i <-  as.numeric(substr(yq.s[i],1,4))
  quarter.i <- as.numeric(substr(yq.s[i],6,7))
  ca.yq <- ca.ag %>% filter(Year==year.i & Quarter==quarter.i)
  hl.yq <- hl.ag %>% filter(Year==year.i & Quarter==quarter.i)
  
  ca.comb <- aggregate(CANoAtLngt~LngtClass+species,data=ca.yq,FUN=sum)
  ca.yq$prop <- ca.yq$CANoAtLngt/ca.comb$CANoAtLngt[match(paste(ca.yq$LngtClass,ca.yq$species),
                                                          paste(ca.comb$LngtClass,ca.comb$species))]
  ca.yq$totHL <- hl.yq$HLNoAtLngt[match(paste(ca.yq$LngtClass,ca.yq$species),
                                        paste(hl.yq$LngtClass,hl.yq$species))]
  ca.yq$N <- ca.yq$prop*ca.yq$totHL
  ca.yq <- ca.yq %>% filter(!is.na(N))
  #Cod
  #####
  lc <-0:120
  pdfC <- matrix(0,ncol=8,nrow=length(lc))
  m_minus1 <- 0
  sd_minus1 <- 0
  ca.cod <- ca.yq %>% filter(species=="cod")
  for(j in 0:7){
    N <- ca.cod$N[ca.cod$Age==j]
    if(length(N)==0) next
    l.obs <- ca.cod$LngtClass[ca.cod$Age==j]
    m <- sum(l.obs*N)/sum(N)
    s <- sqrt(sum(N*(l.obs-m)^2)/sum(N))
    if(m<m_minus1) m<- m_minus1
    if(s<sd_minus1) s<- sd_minus1
    if(is.na(sd(l.obs)) | sd(l.obs)==0) s<- sd_minus1
    #pdfC[,j+1] <- dnorm(lc,mean=m,sd=s)*sum(ca.cod$CANoAtLngt[ca.cod$Age==j])
    # Convert to log-normal parameters (meanlog, sdlog)
    meanlog <- log(m^2 / sqrt(s^2 + m^2))
    sdlog   <- sqrt(log(1 + (s^2 / m^2)))
    pdfC[, j+1] <- dlnorm(lc, meanlog = meanlog, sdlog = sdlog) * sum(N)# Evaluate log-normal density
    m_minus1 <- m
    sd_minus1 <- s
    #print(m)
    #print(s)
}
  if(sum(pdfC)!=0){ # If no age data use alk from previous year
    alkC <- pdfC/rowSums(pdfC)
  }
  hl.cod <- hl.yq %>% filter(species=="cod")
  hl_NC <- (ca.yq %>% select(Quarter,Year,LngtClass,species,Age,N))[1,]
  for(j in 1:length(unique(hl.cod$LngtClass))){
    l.idx <- lc==sort(unique(hl.cod$LngtClass))[j]
    hl.idx <- hl.cod$LngtClass==lc[l.idx]
    N.i <- hl.cod$HLNoAtLngt[hl.idx]
    N.age <- N.i*alkC[l.idx,]
    hl.i <- data.frame(Quarter=hl.cod$Quarter[hl.idx],
                       Year=hl.cod$Year[hl.idx],
                       LngtClass=lc[l.idx],species=hl.cod$species[hl.idx],
                       Age=0:7,N=N.age)
    hl_NC <- rbind(hl_NC,hl.i)
  }
  hl_NC <- hl_NC[-1,]
  #####
  
  #Herring
  #####
  lh <-0:50
  pdfH <- matrix(0,ncol=9,nrow=length(lh))
  m_minus1 <- 0
  sd_minus1 <- 0
  ca.herring <- ca.yq %>% filter(species=="herring")
  for(j in 0:8){
    N <- ca.herring$N[ca.herring$Age==j]
    if(length(N)==0) next
    l.obs <- ca.herring$LngtClass[ca.herring$Age==j]
    m <- sum(l.obs*N)/sum(N)
    s <- sqrt(sum(N*(l.obs-m)^2)/sum(N))
    if(m<m_minus1) m<- m_minus1
    if(s<sd_minus1) s<- sd_minus1
    if(is.na(sd(l.obs)) | sd(l.obs)==0) s<- sd_minus1
    #pdfH[,j+1] <- dnorm(lh,mean=m,sd=s)*sum(ca.herring$CANoAtLngt[ca.herring$Age==j])
    # Convert to log-normal parameters (meanlog, sdlog)
    meanlog <- log(m^2 / sqrt(s^2 + m^2))
    sdlog   <- sqrt(log(1 + (s^2 / m^2)))
    pdfH[, j+1] <- dlnorm(lh, meanlog = meanlog, sdlog = sdlog) * sum(N)# Evaluate log-normal density
    m_minus1 <- m
    sd_minus1 <- s
    #print(m)
    #print(s)
  }  
  if(sum(pdfH)!=0){ # If no age data use alk from previous year
    alkH <- pdfH/rowSums(pdfH)
  }
  hl.herring <- hl.yq %>% filter(species=="herring")
  hl_NH <- (ca.yq %>% select(Quarter,Year,LngtClass,species,Age,N))[1,]
  for(j in 1:length(unique(hl.herring$LngtClass))){
    l.idx <- lh==sort(unique(hl.herring$LngtClass))[j]
    hl.idx <- hl.herring$LngtClass==lh[l.idx]
    N.i <- hl.herring$HLNoAtLngt[hl.idx]
    N.age <- N.i*alkH[l.idx,]
    hl.i <- data.frame(Quarter=hl.herring$Quarter[hl.idx],
                       Year=hl.herring$Year[hl.idx],
                       LngtClass=lh[l.idx],species=hl.herring$species[hl.idx],
                       Age=0:8,N=N.age)
    hl_NH <- rbind(hl_NH,hl.i)
  }
  hl_NH <- hl_NH[-1,]
  #####
  
  #flounder
  #####
  if(year.i==1991){
    ca.yq <- ca.ag %>% filter(Year==1992 & Quarter==quarter.i)
    hl.92 <- hl.ag %>% filter(Year==1992 & Quarter==quarter.i)
    
    ca.comb <- aggregate(CANoAtLngt~LngtClass+species,data=ca.yq,FUN=sum)
    ca.yq$prop <- ca.yq$CANoAtLngt/ca.comb$CANoAtLngt[match(paste(ca.yq$LngtClass,ca.yq$species),
                                                            paste(ca.comb$LngtClass,ca.comb$species))]
    ca.yq$totHL <- hl.92$HLNoAtLngt[match(paste(ca.yq$LngtClass,ca.yq$species),
                                          paste(hl.92$LngtClass,hl.92$species))]
    ca.yq$N <- ca.yq$prop*ca.yq$totHL
    ca.yq <- ca.yq %>% filter(!is.na(N))
  }
  lf <-0:60
  pdfF <- matrix(0,ncol=11,nrow=length(lf))
  m_minus1 <- 0
  sd_minus1 <- 0
  ca.flounder <- ca.yq %>% filter(species=="flounder")
  for(j in 0:10){
    N <- ca.flounder$N[ca.flounder$Age==j]
    if(length(N)<=1) next
    l.obs <- ca.flounder$LngtClass[ca.flounder$Age==j]
    m <- sum(l.obs*N)/sum(N)
    s <- sqrt(sum(N*(l.obs-m)^2)/sum(N))
    if(m<m_minus1) m<- m_minus1
    if(s<sd_minus1) s<- sd_minus1
    if(is.na(sd(l.obs)) | sd(l.obs)==0) s<- sd_minus1
    #pdfF[,j+1] <- dnorm(lf,mean=m,sd=s)*sum(ca.flounder$CANoAtLngt[ca.flounder$Age==j])
    # Convert to log-normal parameters (meanlog, sdlog)
    meanlog <- log(m^2 / sqrt(s^2 + m^2))
    sdlog   <- sqrt(log(1 + (s^2 / m^2)))
    pdfF[, j+1] <- dlnorm(lf, meanlog = meanlog, sdlog = sdlog) * sum(N)# Evaluate log-normal density
    m_minus1 <- m
    sd_minus1 <- s
    #print(m)
    #print(s)
  }
  if(sum(pdfF)!=0){ # If no age data use alk from previous year
    alkF <- pdfF/rowSums(pdfF)
  }
  hl.flounder <- hl.yq %>% filter(species=="flounder")
  hl_NF <- (ca.yq %>% select(Quarter,Year,LngtClass,species,Age,N))[1,]
  for(j in 1:length(unique(hl.flounder$LngtClass))){
    l.idx <- lf==sort(unique(hl.flounder$LngtClass))[j]
    hl.idx <- hl.flounder$LngtClass==lf[l.idx]
    N.i <- hl.flounder$HLNoAtLngt[hl.idx]
    N.age <- N.i*alkF[l.idx,]
    hl.i <- data.frame(Quarter=quarter.i,
                       Year=year.i,
                       LngtClass=lf[l.idx],species=hl.flounder$species[hl.idx],
                       Age=0:10,N=N.age)
    hl_NF <- rbind(hl_NF,hl.i)
  }
  hl_NF <- hl_NF[-1,]
  #####
  
  #plaice
  #####
  if(year.i==1991){
    ca.yq <- ca.ag %>% filter(Year==1994 & Quarter==quarter.i)
    hl.94 <- hl.ag %>% filter(Year==1994 & Quarter==quarter.i)
    
    ca.comb <- aggregate(CANoAtLngt~LngtClass+species,data=ca.yq,FUN=sum)
    ca.yq$prop <- ca.yq$CANoAtLngt/ca.comb$CANoAtLngt[match(paste(ca.yq$LngtClass,ca.yq$species),
                                                            paste(ca.comb$LngtClass,ca.comb$species))]
    ca.yq$totHL <- hl.94$HLNoAtLngt[match(paste(ca.yq$LngtClass,ca.yq$species),
                                          paste(hl.94$LngtClass,hl.94$species))]
    ca.yq$N <- ca.yq$prop*ca.yq$totHL
    ca.yq <- ca.yq %>% filter(!is.na(N))
  }
  lp <-0:65
  pdfP <- matrix(0,ncol=11,nrow=length(lp))
  m_minus1 <- 0
  sd_minus1 <- 0
  ca.plaice <- ca.yq %>% filter(species=="plaice")
  for(j in 0:10){
    N <- ca.plaice$N[ca.plaice$Age==j]
    if(length(N)==0) next
    l.obs <- ca.plaice$LngtClass[ca.plaice$Age==j]
    m <- sum(l.obs*N)/sum(N)
    s <- sqrt(sum(N*(l.obs-m)^2)/sum(N))
    if(m<m_minus1) m<- m_minus1
    if(s<sd_minus1) s<- sd_minus1
    if(is.na(sd(l.obs)) | sd(l.obs)==0) s<- sd_minus1
    #pdfP[,j+1] <- dnorm(lp,mean=m,sd=s)*sum(ca.plaice$CANoAtLngt[ca.plaice$Age==j])
    # Convert to log-normal parameters (meanlog, sdlog)
    meanlog <- log(m^2 / sqrt(s^2 + m^2))
    sdlog   <- sqrt(log(1 + (s^2 / m^2)))
    pdfP[, j+1] <- dlnorm(lp, meanlog = meanlog, sdlog = sdlog) * sum(N)# Evaluate log-normal density
    m_minus1 <- m
    sd_minus1 <- s
    #print(m)
    #print(s)
  }
  if(sum(pdfP)!=0){ # If no age data use alk from previous year
    alkP <- pdfP/rowSums(pdfP)
  }
  
  hl.plaice <- hl.yq %>% filter(species=="plaice")
  hl_NP <- (ca.yq %>% select(Quarter,Year,LngtClass,species,Age,N))[1,]
  for(j in 1:length(unique(hl.plaice$LngtClass))){
    l.idx <- lp==sort(unique(hl.plaice$LngtClass))[j]
    hl.idx <- hl.plaice$LngtClass==lp[l.idx]
    N.i <- hl.plaice$HLNoAtLngt[hl.idx]
    N.age <- N.i*alkP[l.idx,]
    hl.i <- data.frame(Quarter=hl.plaice$Quarter[hl.idx],
                       Year=hl.plaice$Year[hl.idx],
                       LngtClass=lp[l.idx],species=hl.plaice$species[hl.idx],
                       Age=0:10,N=N.age)
    hl_NP <- rbind(hl_NP,hl.i)
  }
  hl_NP <- hl_NP[-1,]
  #####
  
  #dab
  #####
  if(year.i==1991){
    ca.yq <- ca.ag %>% filter(Year==1994 & Quarter==quarter.i)
    hl.94 <- hl.ag %>% filter(Year==1994 & Quarter==quarter.i)
    
    ca.comb <- aggregate(CANoAtLngt~LngtClass+species,data=ca.yq,FUN=sum)
    ca.yq$prop <- ca.yq$CANoAtLngt/ca.comb$CANoAtLngt[match(paste(ca.yq$LngtClass,ca.yq$species),
                                                            paste(ca.comb$LngtClass,ca.comb$species))]
    ca.yq$totHL <- hl.94$HLNoAtLngt[match(paste(ca.yq$LngtClass,ca.yq$species),
                                          paste(hl.94$LngtClass,hl.94$species))]
    ca.yq$N <- ca.yq$prop*ca.yq$totHL
    ca.yq <- ca.yq %>% filter(!is.na(N))
  }
  ld <-0:50
  pdfD <- matrix(0,ncol=11,nrow=length(ld))
  m_minus1 <- 0
  sd_minus1 <- 0
  ca.dab <- ca.yq %>% filter(species=="dab")
  for(j in 0:10){
    N <- ca.dab$N[ca.dab$Age==j]
    if(length(N)==0) next
    l.obs <- ca.dab$LngtClass[ca.dab$Age==j]
    m <- sum(l.obs*N)/sum(N)
    s <- sqrt(sum(N*(l.obs-m)^2)/sum(N))
    if(m<m_minus1) m<- m_minus1
    if(s<sd_minus1) s<- sd_minus1
    if(is.na(sd(l.obs)) | sd(l.obs)==0) s<- sd_minus1
    #pdfD[,j+1] <- dnorm(ld,mean=m,sd=s)*sum(ca.dab$CANoAtLngt[ca.dab$Age==j])
    # Convert to log-normal parameters (meanlog, sdlog)
    meanlog <- log(m^2 / sqrt(s^2 + m^2))
    sdlog   <- sqrt(log(1 + (s^2 / m^2)))
    pdfD[, j+1] <- dlnorm(ld, meanlog = meanlog, sdlog = sdlog) * sum(N)# Evaluate log-normal density
    m_minus1 <- m
    sd_minus1 <- s
    #print(m)
    #print(s)

  }
  if(sum(pdfD)!=0){ # If no age data use alk from previous year
    alkD <- pdfD/rowSums(pdfD)
  }
  hl.dab <- hl.yq %>% filter(species=="dab")
  hl_ND <- (ca.yq %>% select(Quarter,Year,LngtClass,species,Age,N))[1,]
  for(j in 1:length(unique(hl.dab$LngtClass))){
    l.idx <- ld==sort(unique(hl.dab$LngtClass))[j]
    hl.idx <- hl.dab$LngtClass==ld[l.idx]
    N.i <- hl.dab$HLNoAtLngt[hl.idx]
    N.age <- N.i*alkD[l.idx,]
    hl.i <- data.frame(Quarter=hl.dab$Quarter[hl.idx],
                       Year=hl.dab$Year[hl.idx],
                       LngtClass=ld[l.idx],species=hl.dab$species[hl.idx],
                       Age=0:10,N=N.age)
    hl_ND <- rbind(hl_ND,hl.i)
  }
  hl_ND <- hl_ND[-1,]
  #####
  
  hl_N <- rbind(hl_N,hl_NC,hl_NH,hl_NF,hl_NP,hl_ND)
  print(year.i)
}

qday <- aggregate(jday~Quarter+Year,data=ca,FUN=mean) # average sampling day pr. quarter

hl_N$jday <- qday$jday[match(paste(hl_N$Year,hl_N$Quarter),
                         paste(qday$Year,qday$Quarter))]
hl_N$jday <- round(hl_N$jday)/365

haulDura <- aggregate(HaulDur~haulID+Quarter+Year,data=hl,FUN = mean)
hd <- aggregate(HaulDur~Quarter+Year,data=haulDura,FUN = sum)
hl_N <- hl_N %>% left_join(hd)
hl_N$CPUE <- hl_N$N/hl_N$HaulDur*60 # CPUE in ind. pr. hour


rm(list=setdiff(ls(),c('ca.ag','hl.ag','hl_N')))

yr <- 2005
q <- 1
spec <- "plaice"

# Plots

ca.yq <- ca.ag %>% filter(Year==yr & Quarter==q & species==spec)
hl.yq <- hl.ag %>% filter(Year==yr & Quarter==q & species==spec)

ca.comb <- aggregate(CANoAtLngt~LngtClass,data=ca.yq,FUN=sum)
ca.yq$prop <- ca.yq$CANoAtLngt/ca.comb$CANoAtLngt[match(ca.yq$LngtClass,ca.comb$LngtClass)]
ca.yq$totHL <- hl.yq$HLNoAtLngt[match(ca.yq$LngtClass,hl.yq$LngtClass)]
ca.yq$N <- ca.yq$prop*ca.yq$totHL
ca.yq <- ca.yq %>% filter(!is.na(N))

ggplot(ca.yq, aes(x = LngtClass, y = log(N), fill = factor(Age))) +
  geom_col() +xlim(0,50)+
  labs(x = "Length [cm]", y = "log(ind.)", fill = "Age [years]") +
  theme_minimal()

#ggplot(ca.yq, aes(x = LngtClass, y = (N), fill = factor(Age))) +
#  geom_col() +xlim(0,50)+facet_wrap(~Age)+
#  labs(x = "Length [cm]", y = "log(ind.)", fill = "Age [years]") +
#  theme_minimal()


ggplot(hl_N  %>% filter(Year==yr & Quarter==q & species==spec), aes(x = LngtClass, y = log(CPUE), fill = factor(Age))) +
  geom_col() +xlim(0,50)+ylim(0,2)+
  labs(x = "Length [cm]", y = "log(ind.)", fill = "Age [years]") +
  theme_minimal()


N_sums <- aggregate(CPUE~LngtClass,data=hl_N   %>% filter(Year==yr & Quarter==q & species==spec),FUN=sum)
hl_N$tot <- N_sums$CPUE[match(hl_N$LngtClass,N_sums$LngtClass)]
ggplot(hl_N  %>% filter(Year==yr & Quarter==q & species==spec), aes(x = LngtClass, y = CPUE/tot, fill = factor(Age))) +
  geom_col() +xlim(0,50)+
  labs(x = "Length [cm]", y = "log(ind.)", fill = "Age [years]") +
  theme_minimal()

rm(list=setdiff(ls(),c('hl_N')))

plaice98 <- hl_N %>% filter(Year==1998 & species=="plaice")
hl98 <- hl.
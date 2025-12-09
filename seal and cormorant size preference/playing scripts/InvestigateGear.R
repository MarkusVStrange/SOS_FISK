library(LaplacesDemon)
library(surveyIndex)
library(RTMB)

# working directory to coefficient-file
wd_coef <- "C:/Users/mavast/Documents/GitHub/SOS data/"

dAll <- readRDS(paste(wd_coef,"DatrasExchange.R",sep=""))
WB <- st_read(paste(wd_coef,"shapefiles/ICES_areas.shp",sep=""))
plot(WB %>% filter(ICES_SUB %in% 22:24))
WB <- as(WB,"Spatial")

dAll <- addSpatialData(dAll,WB)

# Cod ALK QI
#####
d<-subset(dAll, Species=="Gadus morhua",ICES_SUB %in% 22:24)
dAll<-NULL; gc(); ## lose dAll because it takes up a lot of memory
d<-addSpectrum(d,by=1)

ca <- d[[1]]
hh <- d[[2]]
hl <- d[[3]]

hl <- aggregate(HLNoAtLngt~haul.id+LngtCm,data=hl,FUN=sum)

ca$lengthID <- paste(ca$haul.id,ca$LngtCm)
hl$lengthID <- paste(hl$haul.id,hl$LngtCm)

ca_cod <- left_join(ca,hl %>% select(lengthID,HLNoAtLngt))
ca_cod$cohort <- as.numeric(as.character(ca_cod$Year))-ca_cod$Age
ca_cod <- ca_cod %>% filter(!is.na(HLNoAtLngt))
hh$Jday <- yday(make_date(year=hh$Year,month=hh$Month,day=hh$Day))/yday(make_date(year=hh$Year,month=12,day=31))
ca_cod <- left_join(ca_cod,hh %>% select(haul.id,Jday))
ca_cod <- ca_cod %>% filter(!is.na(Jday))
ca_cod <- ca_cod %>% filter(!is.na(Age))
#ca_cod <- ca_cod %>% filter(!(Age<1 & LngtCm>30))

ca_cod <- ca_cod %>% filter(Age<10)
H20 <- ca_cod %>% filter(Gear=="H20")
FOT <- ca_cod %>% filter(Gear=="FOT")
GOV <- ca_cod %>% filter(Gear=="GOV")
GRT <- ca_cod %>% filter(Gear=="GRT")
SON <- ca_cod %>% filter(Gear=="SON")
TVL <- ca_cod %>% filter(Gear=="TVL")
TVS <- ca_cod %>% filter(Gear=="TVS")


ca_cod$aq <- paste(ca_cod$Age,ca_cod$Quarter)

H20_mean <- numeric(length(unique(ca_cod$aq)))
FOT_mean <- numeric(length(unique(ca_cod$aq)))
GOV_mean <- numeric(length(unique(ca_cod$aq)))
GRT_mean <- numeric(length(unique(ca_cod$aq)))
SON_mean <- numeric(length(unique(ca_cod$aq)))
TVL_mean <- numeric(length(unique(ca_cod$aq)))
TVS_mean <- numeric(length(unique(ca_cod$aq)))

aqs <- sort(unique(ca_cod$aq))
for(i in 1:length(aqs)){
  Q <- substr(aqs[i],3,3)
  Age <- as.numeric(as.character(substr(aqs[i],1,1)))
  H20_mean[i] <- sum(H20$LngtCm[which(H20$Age==Age & H20$Quarter==Q)]*H20$HLNoAtLngt[which(H20$Age==Age & H20$Quarter==Q)],na.rm = T)/
    sum(H20$HLNoAtLngt[which(H20$Age==Age & H20$Quarter==Q)])
  FOT_mean[i] <- sum(FOT$LngtCm[which(FOT$Age==Age & FOT$Quarter==Q)]*FOT$HLNoAtLngt[which(FOT$Age==Age & FOT$Quarter==Q)],na.rm = T)/
    sum(FOT$HLNoAtLngt[which(FOT$Age==Age & FOT$Quarter==Q)])
  GOV_mean[i] <- sum(GOV$LngtCm[which(GOV$Age==Age & GOV$Quarter==Q)]*GOV$HLNoAtLngt[which(GOV$Age==Age & GOV$Quarter==Q)],na.rm = T)/
    sum(GOV$HLNoAtLngt[which(GOV$Age==Age & GOV$Quarter==Q)])
  GRT_mean[i] <- sum(GRT$LngtCm[which(GRT$Age==Age & GRT$Quarter==Q)]*GRT$HLNoAtLngt[which(GRT$Age==Age & GRT$Quarter==Q)],na.rm = T)/
    sum(GRT$HLNoAtLngt[which(GRT$Age==Age & GRT$Quarter==Q)])
  SON_mean[i] <- sum(SON$LngtCm[which(SON$Age==Age & SON$Quarter==Q)]*SON$HLNoAtLngt[which(SON$Age==Age & SON$Quarter==Q)],na.rm = T)/
    sum(SON$HLNoAtLngt[which(SON$Age==Age & SON$Quarter==Q)])
  TVL_mean[i] <- sum(TVL$LngtCm[which(TVL$Age==Age & TVL$Quarter==Q)]*TVL$HLNoAtLngt[which(TVL$Age==Age & TVL$Quarter==Q)],na.rm = T)/
    sum(TVL$HLNoAtLngt[which(TVL$Age==Age & TVL$Quarter==Q)])
  TVS_mean[i] <- sum(TVS$LngtCm[which(TVS$Age==Age & TVS$Quarter==Q)]*TVS$HLNoAtLngt[which(TVS$Age==Age & TVS$Quarter==Q)],na.rm = T)/
    sum(TVS$HLNoAtLngt[which(TVS$Age==Age & TVS$Quarter==Q)])
}

A <- rep(sort(unique(ca_cod$Age)),each=2)+rep(c(0.125,0.875),10)

plot(A,H20_mean,pch=19,cex=1.5,col="black",ylim=c(0,120))
points(A,FOT_mean,pch=19,cex=1.5,col="red")
points(A,GOV_mean,pch=19,cex=1.5,col="orange")
points(A,GRT_mean,pch=19,cex=1.5,col="green")
points(A,SON_mean,pch=19,cex=1.5,col="yellow")
points(A,TVL_mean,pch=19,cex=1.5,col="purple")
points(A,TVS_mean,pch=19,cex=1.5,col="blue")
legend(1,100,c("H20","FOT","GOV","GRT","SON","TVL","TVS"),pch=19,col=c("black","red","orange","green","yellow","purple","blue"))

library(scales)
ca_cod$gear_counter <- 1
Gear_freq <- aggregate(gear_counter~Gear+Year,data=ca_cod,FUN=sum)
Year_fish <- aggregate(gear_counter~Year,data=ca_cod,FUN=sum)
names(Year_fish) <- c("Year","n")
Gear_freq <- left_join(Gear_freq,Year_fish)
Gear_freq$Year <- as.numeric(as.character(Gear_freq$Year))
Gear_freq$Fraq <- Gear_freq$gear_counter/Gear_freq$n

ggplot(Gear_freq, aes(x = Year, y = Fraq, fill = Gear)) +
  geom_bar(stat = "identity", position = "stack", color = "black",width=1) +
  ylab("Diet proportion") +
  ggtitle("Grey Seal food") +
  xlab("Year") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    text = element_text(color = "black", size = 14, family = "serif"),
    legend.title = element_blank()
  )

############
# Problems with some of the German data - 2015 and 2018 Age = 0
###########
DE <- ca_cod %>% filter(Country=="DE" & Year==2015 & Age==0 & Quarter==4)
hist(DE$LngtCm)


t <- ca_cod %>% filter(Age==0 & Year==2016)
ggplot(t, aes(x = LngtCm, y = HLNoAtLngt)) +
  geom_bar(stat = "identity", color = "black",width=1) +
  ylab("Frequency") +
  ggtitle("Age 0 of the 2016 cohort") +
  xlab("Length [cm]") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    text = element_text(color = "black", size = 14, family = "serif"),
    legend.title = element_blank()
  )

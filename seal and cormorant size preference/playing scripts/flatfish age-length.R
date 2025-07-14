source('prepare DATRAS.R')
rm(list=setdiff(ls(),c('hl_N')))


hl_N$cohort <- hl_N$Year-hl_N$Age
table(hl_N$cohort)

# make a 7+ age group for all species
hl_N$Age[hl_N$Age>7] <- 7

N.ag <- aggregate(N_age~LngtClass+species,data=hl_N,FUN = mean)

ggplot(data=N.ag,aes(x=LngtClass,y=N_age,fill = species))+
  geom_bar(stat = "identity", width = 1)+
  facet_wrap(~species)+
  labs(x='Age [years]',y=~paste('CPUE [#  ',h^-1,year^-1,']'))+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position =  'none')


N.year <- aggregate(N_age~species+Year,data=hl_N,FUN = sum)

ggplot(data=N.year %>% filter(Year %in% 1992:1994),aes(x=factor(Year),y=N_age,color = species))+
  geom_point()+
  facet_wrap(~species)+
  labs(x='year',y=~paste('CPUE [#  ',h^-1,']'))+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position =  'none')

N.age <- aggregate(N_age~Age+species,data=hl_N,FUN = mean)

ggplot(data=N.age,aes(x=Age,y=N_age,fill = species))+
  geom_bar(stat = "identity", width = 1)+
  facet_wrap(~species)+
  labs(x='Age [years]',y=~paste('CPUE [#  ',h^-1,']'))+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position =  'none')


cohort88 <- hl_N %>% filter(cohort==2009)
cohort88_ag <- aggregate(N_age~Age+Quarter+species,data=cohort88,FUN = sum)

cohort88.Q1 <- cohort88_ag %>% filter(Quarter==1)


ggplot(data=cohort88.Q1,aes(x=Age,y=log(N_age),color = species))+
  geom_point(pch=19,cex=2)+
  facet_wrap(~species)+
  labs(x='Age [years]',y=~paste('CPUE [#  ',h^-1,']'))+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position =  'none')


table(hl_N$cohort) 

# estimate mortality, Z (M+F), 1985-2014
hl_N.ag <-  aggregate(N_age~Age+Quarter+species+cohort,data=hl_N,FUN = sum)
cohort <- 1985:2014
d <- hl_N.ag %>% filter(cohort %in% cohort & Quarter %in% c(1,4) & Age>2)
table(paste(d$Quarter,d$cohort))
M_dab.Q1 <- rep(0,length(cohort))
M_flounder.Q1 <- rep(0,length(cohort))
M_plaice.Q1 <- rep(0,length(cohort))
M_dab.Q4 <- rep(0,length(cohort))
M_flounder.Q4 <- rep(0,length(cohort))
M_plaice.Q4 <- rep(0,length(cohort))

est_N <- hl_N.ag[1,]
est_N[1,] <- NA


for(i in 1:length(cohort)){
  cohort.i <- cohort[i]
  dada <- d %>% filter(cohort==cohort.i)
  
  fit.Q1 <- coef(lm(log(N_age)~Age,data=dada %>% filter(Quarter==1 & species=="dab")))[2]
  summary(fit.Q1)
  fit.Q1 <- lm(log(N_age)~Age,data=dada %>% filter(Quarter==1 & species=="plaice"))
  summary(fit.Q1)
  fit.Q1 <- lm(log(N_age)~Age,data=dada %>% filter(Quarter==1 & species=="flounder"))
  summary(fit.Q1)
  
  
  M_dab.Q1[i] <- coef(lm(log(N_age)~Age,data=dada %>% filter(Quarter==1 & species=="dab")))[2]
  M_dab.Q4[i] <- coef(lm(log(N_age)~Age,data=dada %>% filter(Quarter==4 & species=="dab")))[2]
  
  M_plaice.Q1[i] <- coef(lm(log(N_age)~Age,data=dada %>% filter(Quarter==1 & species=="plaice")))[2]
  M_plaice.Q4[i] <- coef(lm(log(N_age)~Age,data=dada %>% filter(Quarter==4 & species=="plaice")))[2]
  
  M_flounder.Q1[i] <- coef(lm(log(N_age)~Age,data=dada %>% filter(Quarter==1 & species=="flounder")))[2]
  M_flounder.Q4[i] <- coef(lm(log(N_age)~Age,data=dada %>% filter(Quarter==4 & species=="flounder")))[2]
  
  d.i <- hl_N.ag[1:36,]
  d.i$Age <- rep(0:5,6)
  d.i$Quarter <- rep(rep(c(1,4),each=6),3)
  d.i$species <- rep(c("dab","plaice","flouder"),each=12)
  d.i$cohort <- cohort.i
  d.i$N_age <- 0
}


par(mfrow=c(1,3))
plot(cohort,M_dab.Q1,ylim = c(-2,0),main="Dab",ylab="mortality coefficient")
points(cohort,M_dab.Q4,pch=19)
plot(cohort,M_plaice.Q1,ylim = c(-2,0),main="Plaice",ylab="")
points(cohort,M_plaice.Q4,pch=19)
plot(cohort,M_flounder.Q1,ylim = c(-2,0),main="Flounder",ylab="")
points(cohort,M_flounder.Q4,pch=19)

par(mfrow=c(1,3))
plot(cohort,M_dab.Q1-M_dab.Q4,ylim = c(-1,1),main="Dab",ylab="Q1 and Q4 difference",pch=19)
plot(cohort,M_plaice.Q1-M_plaice.Q4,ylim = c(-1,1),main="Plaice",ylab="",pch=19)
plot(cohort,M_flounder.Q1-M_flounder.Q4,ylim = c(-1,1),main="Flounder",ylab="",pch=19)

mean(M_dab.Q1-M_dab.Q4)



N.Y_Q <- aggregate(N_age~species+Year+Quarter+Age,data=hl_N,FUN = sum)

dab <- N.Y_Q %>% filter(species=="dab")
plaice <- N.Y_Q %>% filter(species=="plaice")
flounder <- N.Y_Q %>% filter(species=="flounder")
table(dab$Age)
table(plaice$Age)
table(flounder$Age)





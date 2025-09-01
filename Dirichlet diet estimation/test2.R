data_wd <- paste(dirname(dirname(getwd())),"/SOS data/",sep="") # data working directory

sealFood <- read.table(paste(data_wd,"sealFood_samplings.csv",sep=""),header=TRUE,sep=';')
sealFood1 <- read.table(paste(data_wd,"sealFood_month.csv",sep=""),header=TRUE,sep=';')
sealFood2 <- read.table("Sealfood_sim.csv",header=TRUE,sep=';')


species <- c("cod","herring","flounder","plaice","dab")
par(mfrow=c(2,3))
for (i in 1:5){
  plot(sealFood$B_index[sealFood$species==species[i]],
       sealFood1$biomass[sealFood1$species==species[i]],xlab = "dnorm",ylab="rnorm")
  lines(c(0,10^8),c(0,10^8),lwd=2,col="red")
}


cormFood <- read.table(paste(data_wd,"cormorantFood_samplings.csv",sep=""),header=TRUE,sep=';')
cormFood1 <- read.table(paste(data_wd,"cormorantFood_month.csv",sep=""),header=TRUE,sep=';')
cormFood2 <- read.table("cormorantFood_month.csv",header=TRUE,sep=';')


species <- c("cod","herring","flounder","plaice","dab")
par(mfrow=c(2,3))
for (i in 1:5){
  plot(cormFood$B_index[cormFood$species==species[i]],
       cormFood1$biomass[cormFood1$species==species[i]],
       xlab = "dnorm",ylab="rnorm",
       xlim=c(0,max(c(cormFood$B_index[cormFood$species==species[i]],
                      cormFood1$biomass[cormFood1$species==species[i]]))),
       ylim=c(0,max(c(cormFood$B_index[cormFood$species==species[i]],
                      cormFood1$biomass[cormFood1$species==species[i]]))),
       main=species[i])
  lines(c(0,10^9),c(0,10^9),lwd=2,col="red")
}

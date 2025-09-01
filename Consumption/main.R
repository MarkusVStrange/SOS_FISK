# packages required
library(ggplot2)
library(dplyr)
library(egg)
library(smsR)
# read in functions from other scripts
source('getPredators.R')
source('getDiets.R')
source('getEnergyBudget.R')
source('getConsumption.R')
source('getConsPlot.R')
source('getDietPlot.R')

# predator abundance
Predators <- getPredators(years = 1985:2023,species = c('cormorant','grey seal'))
# predator diets
Diets <- getDiets(Predators,method='data') # takes a little while
# how much fish (g) the predators eat of selected prey fish
g_eaten <- getEnergyBudget(Diets,method='model')
# the total consumption (tonnes) of selected prey fish by the predators
consumption <- getConsumption(g_eaten,Predators)
# plot the consumption compared to the fishery
getConsPlot(prey='cod',consumption,type="total")
# plot the diets
getDietPlot(predators=c('cormorant','grey seal'),Diets,year=c(1993)) # this is weight-based diets and n denote number of excrements


p <- getConsPlot(prey='cod',consumption,type="total")
p
#ggsave("plots/consumption.png",p,width=7,height=6)

source('getConsInd.R')
df <- getConsInd(consumption)
write.table(df,"consumed_ind.csv",row.names = FALSE,sep=';')


gseal <- Predators %>% filter(species=="grey seal" & quarter=="Q1")
plot(gseal$year,gseal$N_preds,type='l',lwd=3,xlab="Year",ylab="Grey seals",ylim=c(0,2000))

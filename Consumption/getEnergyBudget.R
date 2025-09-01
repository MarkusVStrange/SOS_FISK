########################
# a function to get biomass fish eaten by a predator, based on bioenergetics
# 
#' @param Diets = a dataframe containing information of predators and their diets
#' 
#' 
#' @return 
#' a dataframe containing diets in g/predator/day in years specified in "Diets"
#########################

getEnergyBudget <- function(Diets,method){
  E.budget <- data.frame(year = NA,quarter=NA,prey=NA,g_eaten=NA,predator=NA)
  data_wd <- paste(dirname(dirname(getwd())),"/SOS data/",sep="") # data working directory
  
  # Cormorants
  #####
  if(method=='data' & 'cormorant' %in% unique(Diets$predator)){
    # cormorant diet
    corm.diet <- Diets %>% filter(predator=='cormorant')
    
    # energy requirement for a cormorant
    corm_energy_demand <- 2094 # kJ/day/bird (Keller & Visser, 1999) - Winter value for P. carbo sinensis. 
    # there are different factors affecting energy expenditure: temperature, breeding (early breeding low; late high), migratory behaviour
    
    # relevant prey species
    energy_density <- data.frame(species = c("cod","eelpout","flatfish","goby","herring","other","sandeel","sculpin"),
                                 kJ_pr._g_Ww = c(4.9,3.4,5.7,5.56,mean(c(6.21,9.73)),mean(c(4.9,5.7,mean(c(6.21,9.73)),4.2,5.93,3.22,3.4)),5.93,3.22)) # kJ/g Ww
    # Atlantic cod (mean from Mårtenssen et al., 1996) 
    # eelpout (Eder & Lewis 2005)
    # flatfish (Pleuronectiformes from Spitz et al., 2010)
    # goby (Jane Behrens data)
    # Atlantic herring (mean from Mårtenssen et al., 1996). Assuming a 1:1 distribution of immature (6.21 kJ/g) and mature (9.73 kJ/g)
    # sandeel (Karlsen & Andersen, 2011)
    # sculpin (great sculpin from Anthony et al., 2000)
    # other (mean of the other species)
    
    # calculate the biomass eaten as the energy demand for a bird divided by the energy density for prey, distributed out according to the weight-based diet
    corm.diet$g_eaten <- corm.diet$diet*corm_energy_demand/energy_density$kJ_pr._g_Ww #g/day/bird
    
    corm.diet <- corm.diet[,names(E.budget)]

    E.budget <- rbind(E.budget,corm.diet)
    rm(list=setdiff(ls(),c('E.budget','Diets','method')))
  } 
  
  
  if(method=='model' & 'cormorant' %in% unique(Diets$predator)){
    diet_pred <- read.table(paste(data_wd,"pred_diet.csv",sep=""),header=TRUE,sep=';')

    # cormorant diet
    corm.diet <- Diets %>% filter(predator=='cormorant')
    energy_weight <- aggregate(diet~prey,data=corm.diet %>% filter(!(prey %in% c("cod","flatfish"))),FUN=mean)
    energy_weight$weight <- energy_weight$diet/sum(energy_weight$diet) 
    corm.diet <- diet_pred %>% filter(predator=='cormorant')
    yrs_constant <- 1985:1990
    constant_diet <- data.frame(prey=rep(unique(corm.diet$prey),length(yrs_constant)*4),
                                quarter=rep(rep(c("Q1","Q2","Q3","Q4"),each=3),length(yrs_constant)),
                                year=rep(yrs_constant,each=3*4),predator="cormorant",
                                diet=rep(corm.diet$diet[1:12],length(yrs_constant)))
    corm.diet <- rbind(constant_diet,corm.diet)
    
    # energy requirement for a cormorant
    corm_energy_demand <- 2094 # kJ/day/bird (Keller & Visser, 1999) - Winter value for P. carbo sinensis. 
    # there are different factors affecting energy expenditure: temperature, breeding (early breeding low; late high), migratory behaviour
    
    # relevant prey species
    energy_density <- data.frame(species = c("cod","eelpout","flatfish","goby","herring","other","sandeel","sculpin"),
                                 kJ_pr._g_Ww = c(4.9,3.4,5.7,5.56,mean(c(6.21,9.73)),mean(c(4.9,5.7,mean(c(6.21,9.73)),4.2,5.93,3.22,3.4)),5.93,3.22)) # kJ/g Ww
    # Atlantic cod (mean from Mårtenssen et al., 1996) 
    # eelpout (Eder & Lewis 2005)
    # flatfish (Pleuronectiformes from Spitz et al., 2010)
    # goby (Jane Behrens data)
    # Atlantic herring (mean from Mårtenssen et al., 1996). Assuming a 1:1 distribution of immature (6.21 kJ/g) and mature (9.73 kJ/g)
    # sandeel (Karlsen & Andersen, 2011)
    # sculpin (great sculpin from Anthony et al., 2000)
    # other (mean of the other species)
    energy_density <- c(energy_density$kJ_pr._g_Ww[1],energy_density$kJ_pr._g_Ww[3],
                        sum(energy_density$kJ_pr._g_Ww[-c(1,3)]*energy_weight$weight))
    # calculate the biomass eaten as the energy demand for a bird divided by the energy density for prey, distributed out according to the weight-based diet
    corm.diet$g_eaten <- corm.diet$diet*corm_energy_demand/energy_density #g/day/bird
    
    corm.diet <- corm.diet[,names(E.budget)]
    
    E.budget <- rbind(E.budget,corm.diet)
    rm(list=setdiff(ls(),c('E.budget','Diets','method','data_wd')))
  }
  #####
  
  
  
  # Grey seals
  #####
  if(method=='data' & 'grey seal' %in% unique(Diets$predator)){
    # grey seal diet
    g.seal.diet <- Diets %>% filter(predator=='grey seal')
    
    # Calculate bady mass for a grey seal - parameters from Hauksson 2007
    # 1 = male,2=female
    x0 <- -0.59
    A_inf <- c(279.2,164.1)
    a <- c(0.116,0.332) 
    b <- c(0.716,0.849)
    x <- 6 # assuming the average age of grey seal in WB is 6 years
    WB_gseal_mass <- A_inf*(1-exp(-a*(x-x0)))*b
    metabolic_exponent <- 0.75 # no reference
    energy_correction <- (mean(WB_gseal_mass)^metabolic_exponent)/(mean(A_inf)^metabolic_exponent)    
    # energy requirement for a grey seal
    sparling <- energy_correction*(c(25.5,25.5,20.1,42.4)*1000+ c(25.5,25.5,20.1,42.4)*1000*c(1.23,1.23,1.23,1.07))/2
    #energy_req_pr_seal_pr_day <- sparling # Quarterly. kJ/day/seal (sparling and smout, 2003). Assuming a 1:1 sex ratio and that all observed seals are adults
    energy_req_pr_seal_pr_day <- rep(mean(sparling),4)      # Averaged over quarters  kJ/day/seal (sparling and smout, 2003). Assuming a 1:1 sex ratio and that all observed seals are adults
    
    #energy_req_pr_seal_pr_day.AdultF <- c(25.5,25.5,20.1,42.4)*1000 # Quarterly. kJ/day/seal (sparling and smout, 2003). Adult females
    #energy_req_pr_seal_pr_day.AdultM <- c(25.5,25.5,20.1,42.4)*1000*c(1.23,1.23,1.23,1.07) # Quarterly. kJ/day/seal (sparling and smout, 2003). Adult males
    
    #energy_req_pr_seal_pr_day.JuvF <- c(25.5,25.5,20.1,42.4)*1000 # Quarterly. kJ/day/seal (sparling and smout, 2003). Adult females
    #energy_req_pr_seal_pr_day.AdultM <- c(25.5,25.5,20.1,42.4)*1000*c(1.23,1.23,1.23,1.07) # Quarterly. kJ/day/seal (sparling and smout, 2003). Adult males
    
    # relevant prey species
    energy_density <- data.frame(species = c("cod","flatfish","goby","herring","other","sandeel","sprat","whiting"),
                                 kJ_pr._g_Ww = c(4.9,5.7,5.56,mean(c(6.21,9.73)),mean(c(4.9,5.7,5.56,mean(c(6.21,9.73)),5.93,6.5,4.2)),5.93,6.5,4.2)) # kJ/g Ww
    # Atlantic cod (mean from Mårtenssen et al., 1996) 
    # flatfish (Pleuronectiformes from Spitz et al., 2010)
    # goby (Jane Behrens data)
    # Atlantic herring (mean from Mårtenssen et al., 1996). Assuming a 1:1 distribution of immature (6.21 kJ/g) and mature (9.73 kJ/g)
    # other (mean of the other species) 
    # sandeel (Karlsen & Andersen, 2011)
    # sprat (Spitz et al., 2010)
    # whiting (Pedersen 6 Hislop, 2001)
    
    # put the quarterly energy requirements into the dataframe
    g.seal.diet$E_req <- rep(energy_req_pr_seal_pr_day,each=length(unique(g.seal.diet$prey)))
    # calculate the biomass eaten as the energy demand for a seal divided by the energy density for prey, distributed out according to the weight-based diet
    g.seal.diet$g_eaten <- g.seal.diet$diet*g.seal.diet$E_req/energy_density$kJ_pr._g_Ww #g/day/bird
    
    g.seal.diet <- g.seal.diet[,names(E.budget)]
    
    E.budget <- rbind(E.budget,g.seal.diet)
    rm(list=setdiff(ls(),c('E.budget','Diets','method','data_wd')))
  }
  
  
  if(method=='model' & 'grey seal' %in% unique(Diets$predator)){
    diet_pred <- read.table(paste(data_wd,"pred_diet.csv",sep=""),header=TRUE,sep=';')

    # grey seal diet
    g.seal.diet <- Diets %>% filter(predator=='grey seal')
    energy_weight <- aggregate(diet~prey,data=g.seal.diet %>% filter(!(prey %in% c("cod","flatfish"))),FUN=mean)
    energy_weight$weight <- energy_weight$diet/sum(energy_weight$diet) 
    g.seal.diet <- diet_pred %>% filter(predator=='grey seal')
    yrs_constant <- 1985:1990
    constant_diet <- data.frame(prey=rep(unique(g.seal.diet$prey),length(yrs_constant)*4),
                                quarter=rep(rep(c("Q1","Q2","Q3","Q4"),each=3),length(yrs_constant)),
                                year=rep(yrs_constant,each=3*4),predator="grey seal",
                                diet=rep(g.seal.diet$diet[1:12],length(yrs_constant)))
    g.seal.diet <- rbind(constant_diet,g.seal.diet)
    
    # Calculate bady mass for a grey seal - parameters from Hauksson 2007
    # 1 = male,2=female
    x0 <- -0.59
    A_inf <- c(279.2,164.1)
    a <- c(0.116,0.332) 
    b <- c(0.716,0.849)
    x <- 6 # assuming the average age of grey seal in WB is 6 years
    WB_gseal_mass <- A_inf*(1-exp(-a*(x-x0)))*b
    metabolic_exponent <- 0.75 # no reference
    energy_correction <- (mean(WB_gseal_mass)^metabolic_exponent)/(mean(A_inf)^metabolic_exponent)    
    # energy requirement for a grey seal
    sparling <- energy_correction*(c(25.5,25.5,20.1,42.4)*1000+ c(25.5,25.5,20.1,42.4)*1000*c(1.23,1.23,1.23,1.07))/2
    #energy_req_pr_seal_pr_day <- sparling # Quarterly. kJ/day/seal (sparling and smout, 2003). Assuming a 1:1 sex ratio and that all observed seals are adults
    energy_req_pr_seal_pr_day <- mean(sparling)      # Averaged over quarters  kJ/day/seal (sparling and smout, 2003). Assuming a 1:1 sex ratio and that all observed seals are adults
    
    #energy_req_pr_seal_pr_day.AdultF <- c(25.5,25.5,20.1,42.4)*1000 # Quarterly. kJ/day/seal (sparling and smout, 2003). Adult females
    #energy_req_pr_seal_pr_day.AdultM <- c(25.5,25.5,20.1,42.4)*1000*c(1.23,1.23,1.23,1.07) # Quarterly. kJ/day/seal (sparling and smout, 2003). Adult males
    
    #energy_req_pr_seal_pr_day.JuvF <- c(25.5,25.5,20.1,42.4)*1000 # Quarterly. kJ/day/seal (sparling and smout, 2003). Adult females
    #energy_req_pr_seal_pr_day.AdultM <- c(25.5,25.5,20.1,42.4)*1000*c(1.23,1.23,1.23,1.07) # Quarterly. kJ/day/seal (sparling and smout, 2003). Adult males
    
    # relevant prey species
    energy_density <- data.frame(species = c("cod","flatfish","goby","herring","other","sandeel","sprat","whiting"),
                                 kJ_pr._g_Ww = c(4.9,5.7,5.56,mean(c(6.21,9.73)),mean(c(4.9,5.7,5.56,mean(c(6.21,9.73)),5.93,6.5,4.2)),5.93,6.5,4.2)) # kJ/g Ww
    # Atlantic cod (mean from Mårtenssen et al., 1996) 
    # flatfish (Pleuronectiformes from Spitz et al., 2010)
    # goby (Jane Behrens data)
    # Atlantic herring (mean from Mårtenssen et al., 1996). Assuming a 1:1 distribution of immature (6.21 kJ/g) and mature (9.73 kJ/g)
    # other (mean of the other species) 
    # sandeel (Karlsen & Andersen, 2011)
    # sprat (Spitz et al., 2010)
    # whiting (Pedersen 6 Hislop, 2001)
    energy_density <- c(energy_density$kJ_pr._g_Ww[1],energy_density$kJ_pr._g_Ww[2],
                        sum(energy_density$kJ_pr._g_Ww[-c(1,2)]*energy_weight$weight))
    # put the quarterly energy requirements into the dataframe
    #g.seal.diet$E_req <- rep(energy_req_pr_seal_pr_day,each=length(unique(g.seal.diet$prey)))
    # calculate the biomass eaten as the energy demand for a seal divided by the energy density for prey, distributed out according to the weight-based diet
    g.seal.diet$g_eaten <- g.seal.diet$diet*energy_req_pr_seal_pr_day/energy_density #g/day/bird
    
    g.seal.diet <- g.seal.diet[,names(E.budget)]
    
    E.budget <- rbind(E.budget,g.seal.diet)
    rm(list=setdiff(ls(),c('E.budget','Diets','method')))
  }
  #####
  
  E.budget <- E.budget[-1,]
  E.budget
}

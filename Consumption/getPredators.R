########################
# a function to get the number of predators hunting in a region of interest
# 
#' @param years = a vector with the years to include, e.g. c(1985:2023)
#' @param species = a vector with species names, e.g. c('cormorant','grey seal','harbour seal','harbour porpois')
#' @param area = a vector with the ICES sub.div.s to include, e.g. c(22,23,24)
#' 
#' @return 
#' a dataframe containing number of predators, for each species, year, quarter, and area
#########################

getPredators <- function(years,species,area){
  data_wd <- paste(dirname(dirname(getwd())),"/SOS data/",sep="") # data working directory
  
  # Cormorants
  if ('cormorant' %in% species){
    # define the annual succession of cormorants based on 2014 data
    dk_corms2014 <- data.frame(month = factor(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
                                                "Sep","Oct","Nov","Dec"),levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
                                                                                  "Sep","Oct","Nov","Dec")),
                               cormorants = c(24087,30547,70579,101331,107012,123352,
                                              143633,155047,132985,89269,45036,30866
                               ))
    
    # Find the quarter means
    dk_corms2014.Q <- data.frame(Quarter = c('Q1','Q2','Q3','Q4'),
                                 cormorants = c(mean(dk_corms2014$cormorants[1:3]),
                                                mean(dk_corms2014$cormorants[4:6]),
                                                mean(dk_corms2014$cormorants[7:9]),
                                                mean(dk_corms2014$cormorants[10:12])))
    # Read in the number of cormorants nests
    nests <- read.table(paste(data_wd,"cormorant nests in Denmark.csv",sep=""),sep=';',header=TRUE)
    nests <- nests %>% filter(year %in% years)
    
    # Find the proportion of nests to birds in 2014
    corm.ann.succ. <- dk_corms2014.Q$cormorants/nests$Nests[nests$year==2014] #cormorants pr. nest
    # Extrapolate into the whole time series of nests. Dimensions: [quarter, year]
    dk_corms <- corm.ann.succ. %o% nests$Nests
    WB_cor <- 1-matrix(rep(nests$non.WB.feeders,each=4),nrow=4)  # Disregard cormorants in Northern Kattegat and the Jutland west coast
    dk_WB.corms <- dk_corms*WB_cor
    
    
    Preds <- data.frame(year = rep(nests$year,each=4),
                        quarter = rep(c('Q1','Q2','Q3','Q4'),length(nests$year)),
                        species = rep('cormorant',4*length(nests$year)),
                        N_preds = as.vector(dk_WB.corms))
  }
  
  # Grey seals
  if ('grey seal' %in% species){
    # read in grey seal counts from ICES sub.div. 22:24. From Ander Galatius pers. comm.
    WB_greyseals <- data.frame(year = 2011:2024,
                               count = c(612,405,653,554,802,838,1069,1246,615,1245,1311,
                                         1917,1800,2512))
    WB_greyseals <- WB_greyseals %>% filter(WB_greyseals$year %in% years)
    
    # Include in the predator df. Only one counting pr. year. Assuming representative for all quarters
    Preds <- rbind(Preds,data.frame(year = rep(WB_greyseals$year,each=4),
                                    quarter = rep(c('Q1','Q2','Q3','Q4'),length(WB_greyseals$year)),
                                    species = rep('grey seal',4*length(WB_greyseals$year)),
                                    N_preds = rep(WB_greyseals$count,each=4)))
    
  }
  Preds
}

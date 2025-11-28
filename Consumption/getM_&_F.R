########################
# a function to get mortality from predators (M) and fishery (F)
# 
#' @param consumption = a dataframe with year, prey, predator, and consumption (biomass)
#' 
#' @return 
#' a dataframe containing yearly mortality parameters for each predator and the fishery
#########################
library(reshape)
getMortality <- function(consumption){
  data_wd <- "C:/Users/mavast/Documents/GitHub/SOS data/"
  
  rel_age <- readRDS(paste(data_wd,"pred_ages.R",sep=""))
  
  M <- data.frame(year=NA,age=NA,predator="O",M=NA) 
  
  if("cormorant" %in% consumption$predator){
    # Quarterly cormorant consumption of cod post 1990
    cod.cons <- consumption %>% filter(predator=="cormorant" & prey=="cod" & year>1990)
    # relative biomass-ind. key for all ages. unit: 1000s ind. / tonne
    cod.pred <- rel_age$cod_corm
    
    ind.cons <- cod.pred
    ind.cons[,-(1:2)] <- ind.cons[,-(1:2)]*cod.cons$consumption_t # consumed ind. unit: 1000s ind.
    #stock assessment numbers
    sms_nis <- readRDS(paste(data_wd,"sms_mortality.RDS",sep=""))
    N <- sms_nis$N %>% filter(years %in% ind.cons$year)
    
    # reshape quarterly consumed ind. matrix to sum over a whole year
    x<- melt(ind.cons[,-(1:2)])
    names(x) <- c("Age","ind.")
    x$quarter <- 1:4
    x$year <- rep(cod.pred$year,8)
    # compute to annual number of individuals consumed
    d <- aggregate(ind.~Age+year,data=x,FUN=sum)
    
    #calculate natural mortality based on ind. consumption
    M_c <- data.frame(year=d$year,age=N$ages,
                      predator="danish cormorants",
                      M = -log(1-pmin(0.9,d$ind./N$N)))
    
    M <- rbind(M,M_c)
    
  }
  if("grey seal" %in% consumption$predator){
    # Quarterly cormorant consumption of cod post 1990
    cod.cons <- consumption %>% filter(predator=="grey seal" & prey=="cod" & year>1990)
    # relative biomass-ind. key for all ages. unit: 1000s ind. / tonne
    cod.pred <- rel_age$cod_gseal %>% filter(year %in% cod.cons$year)
    
    ind.cons <- cod.pred
    ind.cons[,-(1:2)] <- ind.cons[,-(1:2)]*cod.cons$consumption_t # consumed ind. unit: 1000s ind.
    #stock assessment numbers
    sms_nis <- readRDS(paste(data_wd,"sms_mortality.RDS",sep=""))
    N <- sms_nis$N %>% filter(years %in% ind.cons$year)
    
    # reshape quarterly consumed ind. matrix to sum over a whole year
    x<- melt(ind.cons[,-(1:2)])
    names(x) <- c("Age","ind.")
    x$quarter <- 1:4
    x$year <- rep(cod.pred$year,8)
    # compute to annual number of individuals consumed
    d <- aggregate(ind.~Age+year,data=x,FUN=sum)
    
    #calculate natural mortality based on ind. consumption
    M_gs <- data.frame(year=d$year,age=N$ages,
                       predator="grey seals",
                       M = -log(1-pmin(0.9,d$ind./N$N)))
    M <- rbind(M,M_gs)
  }
  M <- M[-1,]
  
  ggplot(data=M) +
    geom_line(aes(x = year, y = M,color=predator),linewidth=1.5)+
    ggtitle("Cod in the western Baltic")+facet_wrap(~age,scales ="free_y")+
    scale_color_manual(values = c("cod stock" = "black","danish cormorants" = "darkblue",
                                  "grey seals"="darkorange","fisheries"="darkred"),name="")+
    xlab("year")+
    theme_bw()+ theme(axis.line = element_line(color='black'),
                      plot.background = element_blank(),
                      text=element_text(color="black", size=18,family="serif"),
                      legend.position = c(0.8,0.15),
                      legend.title = element_blank())
}








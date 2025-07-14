########################
# a function to plot the diets of predators
# 
#' @param predator = the predator species to plot, e.g., 'cormorant' 
#' @param Diet = a dataframe with diet data
#' @param year = the year from which to plot the diet, e.g., 1997
#' 
#' 
#' @return 
#' a plot of diets of a specific species in each quarter
#########################

getDietPlot <- function(predators,diets,year){
  library(egg)
  diets$prey <- factor(diets$prey,levels = c('cod','flatfish','herring','goby','sandeel',
                                             'sprat','whiting','eelpout','sculpin','other'))
  if('cormorant' %in% predators){
    p.corm <- list()
    for(i in 1:length(year)){
      corm.diet <- diets[which(diets$year==year[i] & diets$predator=='cormorant'),]
      lg.pos <- "top"
      
      if(i>1)lg.pos <- "n"
      
      p.corm[[i]] <- ggplot(data=corm.diet) +
        geom_bar( aes(x = quarter, y = 100*diet, fill = prey),
                 stat = "identity",position="stack",alpha=0.6,width=0.7)+xlab("quarter")+
        geom_text(aes(x=quarter,y=105,label = paste("n =",n)))+
        ggtitle(paste('Cormorants',year[i]))+
        scale_y_continuous(name = 'diet [%]',limits = c(0,105)) +
        theme_bw()+ theme(axis.line = element_line(color='black'),
                          plot.background = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.grid.major = element_blank(),
                          text=element_text(color="black", size=16,family="serif"),
                          legend.position = lg.pos,
                          legend.title = element_blank())
    }
  }
  
  if('grey seal' %in% predators){
    p.gseal <- list()
    lg.pos <- "top"
    
    for(i in 1:length(year)){
      if(i>1)lg.pos <- "n"
      g.seal.diet <- diets[which(diets$year==year[i] & diets$predator=='grey seal'),]
      p.gseal[[i]] <- ggplot(data=g.seal.diet) +
        geom_bar( aes(x = quarter, y = 100*diet, fill = prey),
                  stat = "identity",position="stack",alpha=0.6,width=0.7)+xlab("quarter")+
        geom_text(aes(x=quarter,y=105,label = paste("n =",n)))+
        ggtitle(paste('Grey seals',year[i]))+
        scale_y_continuous(name = 'diet [%]',limits = c(0,105)) +
        theme_bw()+ theme(axis.line = element_line(color='black'),
                          plot.background = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.grid.major = element_blank(),
                          text=element_text(color="black", size=16,family="serif"),
                          legend.position = lg.pos,
                          legend.title = element_blank())
    }
    
  }
  if(length(year)==1) ggarrange(p.corm[[1]], p.gseal[[1]],ncol=2)
  if(length(year)==2) ggarrange(p.corm[[1]], p.gseal[[1]],p.corm[[2]],p.gseal[[2]],ncol=2)
  if(length(year)==3) ggarrange(p.corm[[1]], p.gseal[[1]],p.corm[[2]],p.gseal[[2]],p.corm[[3]],p.gseal[[3]],ncol=2)
  
            
}




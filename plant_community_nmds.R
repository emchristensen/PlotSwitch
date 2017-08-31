# are plant communities different on different treatment types just before the 2015 switch?

library(vegan)
library(dplyr)


# =============================================================================================
# Functions

plot_nmds = function(treatments,nmds_out) {
  # treatment info
  treat=treatments$treatment
  
  treatmentcolors = data.frame(treatment=c('C','E','X'),colour=c('green','blue','red'))
  treatmentcolors = merge(treatments,treatmentcolors)
  treatmentcolors = treatmentcolors[order(treatmentcolors$plot),]
  
  ordiplot(nmds_out,type="n",xlim=c(-.5,.5),ylim=c(-.5,.5))
  orditorp(nmds_out,display="sites",col=as.character(treatmentcolors$colour),air=0.01,cex=1.25)
  legend(-.55,.5, c("C","E","X"), cex=0.8, 
         col=c("green","blue","red"), pch=15:15)
  ordihull(nmds_out, treat, display="si",lty=1, col="green", show.groups="C")
  ordihull(nmds_out, treat, display="si", lty=1, col="blue", show.groups="E")
  ordihull(nmds_out, treat, display="si", lty=1, col="red", show.groups="X")
  
}
 
#'
#'
#' @param plantdat data frame of plant quadrat data: Portal_plant_quadrats.csv
#' @param plantsp data frame of plant species data: Portal_plant_species.csv
#' @param plant_duration Annual, Perennial, or blank
#' @param year year of interest
#' @param season winter or summer
#'
#'
get_plant_data = function(plantdat, plantsp, yr, seas, plant_duration=c('Annual','Perennial')) {

  spcodes = filter(plantsp,commonname!='Unknown',duration %in% plant_duration)
  
  # plant census right before flip: remove stake 17 beacuse plot 24 doesn't have one, also remove unknown species
  censusdat = filter(plantdat,season==seas,year==yr,quadrat!=17,species %in% spcodes$speciescode)
  
  census_total = aggregate(censusdat$abundance,by=list(plot = censusdat$plot,species=censusdat$species),FUN=sum)
  census_table = reshape(census_total, idvar='plot', timevar='species',direction='wide')
  census_table[is.na(census_table)] = 0
  census_table = census_table[order(census_table$plot),]
  rownames(census_table) = census_table$plot
  
  return(census_table)
}

# ==========================================================================================

plantdat = read.csv('C:/Users/EC/Desktop/git/PortalData/Plants/Portal_plant_quadrats.csv')
plantsp = read.csv('C:/Users/EC/Desktop/git/PortalData/Plants/Portal_plant_species.csv')

# just for curiosity, look at total abundance by census
totalplant = aggregate(plantdat$abundance,by=list(year = plantdat$year, season=plantdat$season),FUN=sum,na.rm=T)

treatments = data.frame(treatment=c('C','C','E','C','X','E',
                                    'X','C','C','X','C','C',
                                    'E','C','E','X','C','E',
                                    'E','E','E','C','X','X'),
                        plot=seq(24))

#treatments after switch
treatments2 = data.frame(treatment=c('X','E','E','C','C','C',
                                     'C','E','X','X','C','X',
                                     'C','C','E','X','C','C',
                                     'E','E','E','E','X','C'),
                         plot=seq(24))
treatments_combined = data.frame(treatment = c('CX','CE','EE','CC','XC','EC',
                                               'XC','CE','CX','XX','CC','CX',
                                               'EC','CC','EE','XX','CC','EC',
                                               'EE','EE','EE','CE','XX','XC'),
                                 plot=seq(24))

plant_table = get_plant_data(plantdat,plantsp,2008,'winter','Annual')

# distance -- I took this stuff from code Baiser gave me, I think...
plantdist = vegdist(plant_table[-1],'bray')

# the actual nmds calculation.  k=2 to restrict to 2 dimensions
nmds_out<-metaMDS(plantdist,k=2, trace=T)

#stressplot(nmdsrat)

plot_nmds(treatments,nmds_out)

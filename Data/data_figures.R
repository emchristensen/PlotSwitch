# code for visualizing data
library(ggplot2)




# Dipodomys numbers
dipoN = read.csv('Data/Dipo_counts.csv')

dipoN$censusdate = as.Date(dipoN$censusdate)

ggplot(dipoN, aes(x = censusdate, y =n, colour = treatment)) +
  geom_point() +
  geom_smooth(method = 'loess', se = TRUE) +
  scale_colour_brewer(type = 'qual', palette = 'Dark2') +
  theme(legend.position = 'top') +
  geom_vline(xintercept=as.Date('2015-04-10'))

dipo_before = dplyr::filter(dipoN,censusdate < as.Date('2015-04-01'))
dipo_after = dplyr::filter(dipoN,censusdate > as.Date('2015-04-01'))
aggregate(dipo_before$n,by=list(treatment=dipo_before$treatment),FUN=mean)
aggregate(dipo_after$n,by=list(treatment=dipo_after$treatment),FUN=mean)

# # Species richness
# sprich = read.csv('Data/SpeciesRichness.csv')
# sprich$censusdate = as.Date(sprich$censusdate)
# 
# ggplot(sprich, aes(x = censusdate, y =n, colour = treatment)) +
#   geom_point() +
#   geom_smooth(method = 'loess', se = TRUE) +
#   scale_colour_brewer(type = 'qual', palette = 'Dark2') +
#   theme(legend.position = 'top') +
#   geom_vline(xintercept=as.Date('2015-04-10'))
# 
# sprich_before = dplyr::filter(sprich,censusdate < as.Date('2015-04-01'))
# sprich_after = dplyr::filter(sprich,censusdate > as.Date('2015-04-01'))
# aggregate(sprich_before$n,by=list(treatment=sprich_before$treatment),FUN=mean)
# aggregate(sprich_after$n,by=list(treatment=sprich_after$treatment),FUN=mean)


# Total rodent energy
energy = read.csv('Data/TotalCommunityEnergy.csv')
energy$censusdate = as.Date(energy$censusdate)

ggplot(energy, aes(x = censusdate, y =n, colour = treatment)) +
  geom_point() +
  geom_smooth(method = 'loess', se = TRUE) +
  scale_colour_brewer(type = 'qual', palette = 'Dark2') +
  theme(legend.position = 'top') +
  geom_vline(xintercept=as.Date('2015-04-10'))

energy_before = dplyr::filter(energy,censusdate < as.Date('2015-04-01'))
energy_after = dplyr::filter(energy,censusdate > as.Date('2015-04-01'))
aggregate(energy_before$n,by=list(treatment=energy_before$treatment),FUN=mean)
aggregate(energy_after$n,by=list(treatment=energy_after$treatment),FUN=mean)


# Small granivores
smgran = read.csv('Data/SmallGranivores.csv')
smgran$censusdate = as.Date(smgran$censusdate)
ggplot(smgran, aes(x = censusdate, y =n, colour = treatment)) +
  geom_point() +
  geom_smooth(method = 'loess', se = TRUE) +
  scale_colour_brewer(type = 'qual', palette = 'Dark2') +
  theme(legend.position = 'top') +
  geom_vline(xintercept=as.Date('2015-04-10'))

smgran_before = dplyr::filter(smgran,censusdate < as.Date('2015-04-01'))
smgran_after = dplyr::filter(smgran,censusdate > as.Date('2015-04-01'))
aggregate(smgran_before$n,by=list(treatment=smgran_before$treatment),FUN=mean)
aggregate(smgran_after$n,by=list(treatment=smgran_after$treatment),FUN=mean)

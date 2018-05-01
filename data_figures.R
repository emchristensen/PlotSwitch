# code for visualizing data
library(ggplot2)




# Dipodomys numbers
dipoN = read.csv('Dipo_counts.csv')

dipoN$censusdate = as.Date(dipoN$censusdate)

ggplot(dipoN, aes(x = censusdate, y =dipos, colour = treatment)) +
  geom_point() +
  geom_smooth(method = 'loess', se = TRUE) +
  scale_colour_brewer(type = 'qual', palette = 'Dark2') +
  theme(legend.position = 'top') +
  geom_vline(xintercept=as.Date('2015-04-10'))


# Species richness
sprich = read.csv('SpeciesRichness.csv')
sprich$censusdate = as.Date(sprich$censusdate)

ggplot(sprich, aes(x = censusdate, y =n, colour = treatment)) +
  geom_point() +
  geom_smooth(method = 'loess', se = TRUE) +
  scale_colour_brewer(type = 'qual', palette = 'Dark2') +
  theme(legend.position = 'top') +
  geom_vline(xintercept=as.Date('2015-04-10'))


# Total rodent energy
energy = read.csv('TotalCommunityEnergy.csv')
energy$censusdate = as.Date(energy$censusdate)

ggplot(energy, aes(x = censusdate, y =energy, colour = treatment)) +
  geom_point() +
  geom_smooth(method = 'loess', se = TRUE) +
  scale_colour_brewer(type = 'qual', palette = 'Dark2') +
  theme(legend.position = 'top') +
  geom_vline(xintercept=as.Date('2015-04-10'))

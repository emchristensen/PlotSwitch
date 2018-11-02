# example of permanova from:
# https://rpubs.com/collnell/manova scroll down to permanova

library(vegan)
library(dplyr)

# bird example from rpubs ----
birds<-read.csv('https://raw.githubusercontent.com/collnell/lab-demo/master/bird_by_fg.csv')
bird.matrix<-as.matrix(birds[,3:9])##response variables in a sample x species matrix
# 1. transform or standaradize data to minimize influence of most abundant groups
bird.mat <- sqrt(bird.matrix)
bird.prop <- decostand(bird.matrix, method = 'total')
# 2. calculate ecological distance
bird.dist <- vegdist(bird.mat, method='bray')
# 3. perMANOVA
bird.div <- adonis2(bird.dist~DIVERSITY, data = birds, permutations = 999, method = 'bray', strata = 'PLOT')
bird.div
# 4. multivariate dispersion: average distance to group centroid. used as a measure of multivariate beta diversity
dispersion <- betadisper(bird.dist, group=birds$DIVERSITY)
permutest(dispersion)
plot(dispersion, hull=F, ellipse=T)
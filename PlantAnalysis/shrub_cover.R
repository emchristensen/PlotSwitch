library(portalr)
library(dplyr)
library(ggplot2)

# check out shrub cover data
shrubs = shrub_cover(path = 'repo', type = "Shrubs", plots = "all",
                        unknowns = FALSE, correct_sp = TRUE)

# filter out year 2015 (data was collected differently), total by plot
total_byplot = shrubs %>%
  filter(year!=2015) %>%
  filter(species == 'acac cons',plot==3) %>%
  group_by(year,plot,treatment) %>%
  summarise(totalcover = sum(cover))


ggplot(total_byplot,aes(x=year,y=totalcover,color=treatment)) +
  geom_point()


# =============================================================================
# comparison by treatment in 2016 ----
# just the plots used in this chapter
oldtreatments = data.frame(plot=c(4,11,14,17,6,13,18,5,7,24),oldtreat=c(rep('control',4),rep('exclosure',3),rep('removal',3)))
# all the plots
oldtreatments = data.frame(plot=c(4,11,14,17,2,8,22,1,9,12,6,13,18,3,15,19,20,21,5,7,24,10,16,23),
                           oldtreat=c(rep('control',10),rep('exclosure',8),rep('removal',6)))

shrubcover = merge(total_byplot,oldtreatments) %>% filter(year==2016)

cbPalette <- c( "#e19c02","#999999", "#56B4E9", "#0072B2", "#D55E00", "#F0E442", "#009E73", "#CC79A7")

shrubbox = ggplot(shrubcover,aes(x=oldtreat,y=totalcover, colour=oldtreat)) +
  geom_boxplot() +
  geom_point() +
  ggtitle('Kangaroo rat abundances 2013-2015') +
  scale_colour_manual(values=cbPalette) +
  labs(y = 'Count', x = NULL) +
  scale_x_discrete(breaks=c("control","exclosure","removal"),
                   labels=c("long-term\n control", "kangaroo rat\n removal", "rodent\n removal")) +
  theme(legend.position = 'none')
shrubbox


# check if data within factor levels is normal
shrub_control = dplyr::filter(shrubcover,oldtreat=='control')
shrub_excl = dplyr::filter(shrubcover,treatment=='exclosure')
shrub_rem = dplyr::filter(shrubcover,treatment=='removal')

shapiro.test(shrub_rem$totalcover)
# the removals are normal but the other two are not
# need a kruskal-wallis test
kruskal.test(totalcover~oldtreat, data=shrubcover)
pairwise.wilcox.test(shrubcover$totalcover, shrubcover$oldtreat,
                     p.adjust.method = "bonferroni")
# no differences between treatments

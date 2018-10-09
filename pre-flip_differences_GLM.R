# comparing values (number of dipos, number of non-dipo species, and total energy) between treatments, before treatment change in 2015

library(dplyr)
library(multcomp)
library(ggplot2)

# dipos -----
df1 = read.csv('Dipo_counts.csv') %>% filter(numericdate<16.543)
# plot data
ggplot(df1,aes(x=censusdate,y=n,colour=treatment)) + geom_jitter(width=.1,height=.2)

# glm of krat abundance
dipo.glm= glm(n ~ treatment+numericdate, data = df1, family = poisson())
summary(dipo.glm)
confint(dipo.glm)

# glht (general linear hypothesis)
summary(glht(dipo.glm, linfct=mcp(treatment='Tukey')))

# small granivores ----
df2 = read.csv('SmallGranivores.csv') %>% filter(numericdate<16.543)
ggplot(df2,aes(x=censusdate,y=n,colour=treatment)) + geom_jitter(width=.1,height=.2)

sm.glm = glm(n ~ treatment+numericdate, data = df2, family = poisson())
summary(sm.glm)
confint(sm.glm)

summary(glht(sm.glm, linfct=mcp(treatment='Tukey')))

# energy ---- 
df3 = read.csv('TotalCommunityEnergy.csv') %>% filter(numericdate<16.543)
ggplot(df3,aes(x=censusdate,y=n,colour=treatment)) +geom_jitter(width=.1,height=.2)

en.glm = glm(n ~ treatment+numericdate, data = df3, family = gaussian)
summary(en.glm)     
confint(en.glm)

summary(glht(en.glm, linfct=mcp(treatment='Tukey')))

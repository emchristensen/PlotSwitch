# exploring rare/transient species
library(dplyr)

dat = get_data(startdate = "2013-03-11",
               min_num_plots = 21,
               treatments=c('CC','CE','EC','EE','XC','CX','XX'))

# small granivores: c('BA','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO')
# look at each individually to identify core vs transient
BA = filter(dat,species=='BA')
PB = filter(dat,species=='PB')
PE = filter(dat,species=='PE')
PF = filter(dat,species=='PF')
PH = filter(dat,species=='PH')
PI = filter(dat,species=='PI')
PL = filter(dat,species=='PL')
PM = filter(dat,species=='PM')
PP = filter(dat,species=='PP')
RF = filter(dat,species=='RF')
RM = filter(dat,species=='RM')
RO = filter(dat,species=='RO')
DM = filter(dat,species=='DM')
DO = filter(dat,species=='DO')
DS = filter(dat,species=='DS')

ggplot(data=PB,aes(x=censusdate,y=abundance,colour=before_after)) +
  geom_point() +
  geom_smooth()


# core == species present in >=33% of time steps at site-level
testsp = DO

test = aggregate(testsp$abundance,by=list(testsp$censusdate),FUN=sum)
length(test$x[test$x>0])/length(test$x)
plot(test)

#' as of 7/2018:
#'  BA = .41
#'  PB = 1
#'  PE = .98
#'  PF = .54
#'  PH = .11
#'  PI = 0
#'  PL = .20
#'  PM = .33
#'  PP = .96
#'  RF = .09
#'  RM = .61
#'  RO = .33
#'  
#'  DM = 1
#'  DO = .96
#'  DS = .30

library(dplyr)


dat = read.csv('C:/Users/EC/Desktop/git/PortalData/Rodents/Portal_rodent.csv',stringsAsFactors = F)


# what is average apparent life span of dipo species for the year before the switch
dipo_before = filter(dat,period>=424,period<=435,species %in% c('DM')) 

# find tag numbers for animals newly tagged in first 3 months
new_before = filter(dipo_before,note2=='*',period %in% c(424:426),prevrt=='',prevlet=='',plot %in% c(4,11,14,17)) %>% select(tag)

RT_before = data.frame()
for (indiv in new_before$tag) {
  encounters = filter(dipo_before,tag==indiv)
  residence = max(encounters$period) - min(encounters$period)
  moved = 0
  if (length(unique(encounters$plot))>1) {moved=1}
  RT_before = rbind(RT_before,c(residence,moved))
}
names(RT_before) = c('res_time','moved')
hist(RT_before[,1])  

# ==================================================================
# long term controls after
dipo_after = filter(dat,period>=436,period<=447,species=='DM')

new_cc_after = filter(dipo_after,note2=='*',period %in% c(436:438),prevrt=='',prevlet=='',plot %in% c(4,11,14,17)) %>% select(tag)

RT_after = data.frame()
for (indiv in new_cc_after$tag) {
  encounters = filter(dipo_after,tag==indiv)
  residence = max(encounters$period) - min(encounters$period)
  moved=0
  if (length(unique(encounters$plot))>1) {moved=1}
  RT_after = rbind(RT_after,c(residence,moved))
}
names(RT_after) = c('res_time','moved')
hist(RT_after$res_time)

# ==================================================================
# exclosures to controls
new_xc_after = filter(dipo_after,note2=='*',period %in% c(436:438),prevrt=='',prevlet=='',plot %in% c(5,7,24)) %>% select(tag)

RT_after_xc = data.frame()
for (indiv in new_xc_after$tag) {
  encounters = filter(dipo_after,tag==indiv)
  residence = max(encounters$period) - min(encounters$period)
  moved=0
  if (length(unique(encounters$plot))>1) {moved=1}
  RT_after_xc = rbind(RT_after_xc,c(residence,moved))
}
names(RT_after_xc) = c('res_time','moved')
hist(RT_after_xc$res_time)
}

new_ec_after = filter(dipo_after,note2=='*',period %in% c(436:438),prevrt=='',prevlet=='',plot %in% c(6,13,18)) %>% select(tag)

RT_after_ec = data.frame()
for (indiv in new_ec_after$tag) {
  encounters = filter(dipo_after,tag==indiv)
  residence = max(encounters$period) - min(encounters$period)
  moved=0
  if (length(unique(encounters$plot))>1) {moved=1}
  RT_after_ec = rbind(RT_after_ec,c(residence,moved))
}
names(RT_after_ec) = c('res_time','moved')
hist(RT_after_ec$res_time)

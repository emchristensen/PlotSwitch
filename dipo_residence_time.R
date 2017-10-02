
library(dplyr)




#' @param tags list of unique tag numbers
#' @param dat data frame in which to look for captures
#' @return returns data frame with tag, first_cap, last_cap [dates], and moved [0/1 indicates whether the animal was captured on more than one plot]
first_last_cap = function(tags,dat) {
  df = c()
  for (indiv in tags) {
    encounters = filter(dat, tag==indiv)
    first = filter(encounters,period==min(encounters$period))
    first_cap = as.Date(paste(first$year,first$month,first$day,sep='-'))
    last = filter(encounters,period==max(encounters$period))
    last_cap = as.Date(paste(last$year,last$month,last$day,sep='-'))
    moved = 0
    if (length(unique(encounters$plot))>1) {moved = 1}
    df = rbind(df,c(indiv,first_cap,last_cap,moved))
  }
  df = data.frame(df,stringsAsFactors = F)
  names(df) = c('tag','first_cap','last_cap','moved')
  df$span = as.numeric(df$last_cap) - as.numeric(df$first_cap)
  return(df)
}



# what is average apparent life span of dipo species for the year before the switch
dat = read.csv('C:/Users/EC/Desktop/git/PortalData/Rodents/Portal_rodent.csv',stringsAsFactors = F)
newdipodat = filter(dat, period>0, species=='DM', note2=='*', prevrt=='', prevlet=='')

# find tag numbers for animals newly tagged in March-May
new_before_cc = filter(newdipodat, period %in% c(425:427), plot %in% c(4,11,14,17))

RT_before = first_last_cap(new_before_cc$tag,dat)
hist(RT_before$span)  

# ==================================================================
# long term controls after

new_cc_after = filter(newdipodat, period %in% c(437:439), plot %in% c(4,11,14,17))
RT_after_cc = first_last_cap(new_cc_after$tag,dat)

hist(RT_after_cc$span)


# ==================================================================
# exclosures to controls
new_xc_after = filter(newdipodat, period %in% c(437:439), plot %in% c(5,7,24))
RT_after_xc = first_last_cap(new_xc_after$tag,dat)

hist(RT_after_xc$span)


new_ec_after = filter(newdipodat, period %in% c(437:439), plot %in% c(6,13,18))
RT_after_ec = first_last_cap(new_ec_after$tag,dat)

hist(RT_after_ec$span)


# ==================================================================
# compare metrics

# how many newly tagged animals in span
length(RT_before$span)
length(RT_after_cc$span)
length(RT_after_xc$span)
length(RT_after_ec$span)

# average life span of animals captured more than once
filter(RT_before,span>0) %>% select(span) %>% colMeans()
filter(RT_after_cc,span>0) %>% select(span) %>% colMeans()
filter(RT_after_xc,span>0) %>% select(span) %>% colMeans()
filter(RT_after_ec,span>0) %>% select(span) %>% colMeans()

# percent of animals captured more than once
length(RT_before$span[RT_before$span>0])/length(RT_before$span)
length(RT_after_cc$span[RT_after_cc$span>0])/length(RT_after_cc$span)
length(RT_after_xc$span[RT_after_xc$span>0])/length(RT_after_xc$span)
length(RT_after_ec$span[RT_after_ec$span>0])/length(RT_after_ec$span)

# percent of animals that moved plot
length(RT_before$moved[RT_before$moved==1])/length(RT_before$moved)
length(RT_after_cc$moved[RT_after_cc$moved==1])/length(RT_after_cc$moved)
length(RT_after_xc$moved[RT_after_xc$moved==1])/length(RT_after_xc$moved)
length(RT_after_ec$moved[RT_after_ec$moved==1])/length(RT_after_ec$moved)


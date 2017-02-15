

# ==================================================================================
#  This function will take any data frame that has columns for "wgt" and "species" and 
# adds a column for "energy"
# If wgt is NA for any row, it will calculate average energy for that species

metabolic_energy_from_wgt = function(rdat) {
  # calculate metabolic energy from weight
  rdat$energy = 5.69*rdat$wgt**.75
  # average metabolic energy by species
  spmeans = aggregate(rdat$energy[!is.na(rdat$energy)],by=list(species=rdat$species[!is.na(rdat$energy)]),FUN=mean)
  # fill gaps in weight measurements with means
  for (n in seq(length(rdat$energy))) {
    if (is.na(rdat$energy[n])) {
      sp = rdat$species[n]
      fill = spmeans$x[spmeans$species == sp]
      rdat$energy[n] = fill
    }
  }
  return(rdat)
}
# ======================================================================================

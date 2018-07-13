# PlotSwitch
Code for my project on the effects of changing rodent treatments on plots at Portal in 2015. Specifically, I'm looking for evidence of priority effects in slowing the re-colonization of kangaroo rats on patches from which they were excluded, but were opened to them through changing the experimental treatments. 

## Data
  * __create_data_series.R__ code for creating time series data csvs
  * __data_functions.R__ functions used in create_data_series.R
  * __data_figures.R__ code to plot data
  * __Dipo_counts.csv__ counts by plot and sampling time
  * __SpeciesRichness.csv__ species richness by plot and sampling time
  * __TotalCommunityEnergy.csv__ total metabolic energy of all species on a plot combined
  * __SmallGranivores.csv__ counts of small granivore species (Dipo competitors); counts by plot and sampling time

## Model exploration and selection
  * __rodent-analysis.R__ GAM analysis of dipodomys count data
  * __rodent-analysis-energy.R__ GAM analysis of energy data
  * __rodent-analysis-richness.R__ GAM analysis of species richness
  * __rodent-analysis-smallgranivores.R__ GAM analysis of small granivore count data
  
## Final analysis and figures
  * __rodent-GAM-analyses.R__ final models for each analysis, with code for creating figures
  * __analysis_functions.R__ functions used by rodent-GAM-analyses.R

## Title: R script to run functions for *SEIR* stochastic compartmental model with gillespie algorithm  
## Author: Tamika Lunn, Griffith University
## Version: V22, created 5th May 2020

rm(list=ls())
### Set values for visualisation:
threshold = 10 #set threshold value for successful/failed outbreak 
modelrun <- "SEIR_stoch outputs_PO_p50" #change to whichever model set you are creating figures for (name as specified in file directory)

setwd(paste(rprojroot::find_rstudio_root_file(),"/Code/Full code to repliate results",sep=""))
source ("stoch_helperfunctions_V21.R")
# SIR.R0: Calculates average R0 for the roost, from the SIR model, based on the B matrix.
#    - Input = n, Ntree, Bmatrix, gamma
#    - Translation = number of tree groups (single value); string of values showing numbers of bats in each tree group; Beta matrix of within-group and between-group transmission coefficients; rate of recovery
# SIR.Bcrit: Calculates critical B where R0=1, from the SIR model
#    - Input = n, Ntree, scaled.data, gamma
#    - Translation = number of tree groups (single value); string of values showing numbers of bats in each tree group; matrix of distance between groups scaled by set radius value; rate of recovery
#    - Important: scaled.data needs to contain 1 values (not zero values) for within-tree distances (i.e along the diagonal, or for trees closer than the given radius)
# SIR.B: Overall within-roost B value given an average R0roost value, from the SIR model
#    - Input = n, Ntree, scaled.data, gamma, R0roost
#    - Translation = number of tree groups (single value); string of values showing numbers of bats in each tree group; matrix of distance between groups scaled by set radius value; rate of recovery; average R0 of roost
#    - Important: scaled.data needs to contain 1 values (not zero values) for within-tree distances (i.e along the diagonal, or for trees closer than the given radius)
# successful: return a list of successful epidemics, from simend output. Threshold= 10 infected individuals by end of simulation
#    - Input = simend.output, threshold
#    - Translation = values extracted from the end of simualtion (1 dataframe per set of simulations); single value of how many total infected individuals per simulation should classify an outbreak as successful
# mergelistA and mergelistB: merges dataframes in list, based on a pre-defined subset, and returns the merged dataframes as a list. A includes line to remove extra I categories, B does not remove I categories
#    - A Input = data, series.params, param, subset.ids, parameter.values
#    - A Translation = summary dataset of simulation output (format = list of dataframes, each dataframe is a set of simulations for one parameter set); dataframe containing parameter values of each simulation set; name of the parameter that is used to subset the data; vector of position values of data to be subsetted; dataframe of values of parameters in subset from series.params
#    - B Input = data, series.params, param, parameter.values
#    - B Translation = summary dataset of simulation output (format = list of dataframes, each dataframe is a set of simulations for one parameter set); dataframe containing parameter values of each simulation set; name of the parameter that is used to subset the data; dataframe of values of parameters in subset from series.params
# subsetdata.identifiers and subsetdata.data: subset dataframe from dataframe, based on a pre-defined subset. Return either the identifying values of subset, or all data in subset
#    - Input = values, param, dataframe, identifier
#    - Translation = values (numeric) to be matched in the dataframe, name of parameter for values of interest, dataframe to extract values from, identifier of interest (e.g. id or position)
# Numextract and NumextractB: extract number values from strings. Note that Numextract B is specifically for extracting subplot values from single digit strings, and re-formatting with the 000 format
#    - Input = string
#    - Translation = string to extract values from
# read.all: reads in data from a list of csv files and attachs corresponding parameter values from series.params dataframe. read.subset reads in all data and then deletes those not in subset. Two id values are given to make sure that parameter values are assigned to the correct dataframe: idseries=id taken from the file name, id=id taken from series.params dataframe. Position (all) or newposition (subset) are assigned for selecting dataframes in list
#    - Input = series, series.params
#    - Translation = string of file names, dataframe with parameter values of stochastic simulations
# read.subset: reads in data from a list of csv files and attachs corresponding parameter values from series.params dataframe. read.subset reads in all data and then deletes those not in subset. Two id values are given to make sure that parameter values are assigned to the correct dataframe: idseries=id taken from the file name, id=id taken from series.params dataframe. Position (all) or newposition (subset) are assigned for selecting dataframes in list
#    - Input = series, subset, series.params, series.params.subset
#    - Translation = string of file names, identifying values of subset (output from subsetdata.identifiers), dataframe with parameter values of stochastic simulations, subset of dataframe with parameter values of stochastic simulations (output from subsetdata.data)
# addmaxdist: merges simulation dataframes with maxdistance dataframe (maximum distance of spread, per index case). Loops through entire list of dataframes, and merges dataframes one at a time by index and plotid
#    - Input = data, Maxdist
#    - Translation = list of dataframes being used, and a dataframe containing fartherst neighbour distances per plot, per index tree (i.e. max distance of spread per index tree). Data and maxdistance dataframes must have matching plotid values (to distinguist site and subplot containing trees) and index values (id of individual trees) 
# findloc: finds the position of value in dataframe that most closely matches a given, single value
#    - Input = dataframe, value
#    - Translation = dataframe to select value from, value of interest
# extractvalues: 
#    - Input = modstrs, colnames, structure
#    - Translation = names of model structures in the format "SIRns|SIR", colnames for table in the format c("modstr","site")- note that these should be in the same order as in the file name, dataframe with additional tree structure information (e.g Ctree nodes)

### Libraries & helping functions:
library(plyr) #for ddply 
library(anchors) #for replace.value
library(HH) #for horizontal stacked bar charts
library(dplyr)
library(mosaic) #for more efficient ifelse with mutate
library(knitr)
library(kableExtra)
library(gridExtra) #for grid.arrange
library(data.table) #for rbindlist
library(directlabels) #labels for ggplot2
library(spatstat) #for generating random spp for heat maps
library(MALDIquant) #for matching values
library(hrbrthemes) #for viridis colour pallette 
library(viridis)
library(gam) #calculate gam smoothed splines for plotting 

##------------------------------------------------------------------------------------------------------
##---------------------------------------Read in parameter values---------------------------------------
##------------------------------------------------------------------------------------------------------
### Create table of parameter values
modstrs <- "SEIRns|SEIR"  
colnames <- c("modstr","Bfun","site","id","subplot","n","Beta","gamma", "delta","N","Ntot","radius","theta","threshold","nsims","R0roost")  

setwd(paste(rprojroot::find_rstudio_root_file(),"/Data/Processed",sep=""))
structure <- read.csv("plot groups.csv", row.names=1, check.names=FALSE) #categories from Ctree analysis
structure$structure <- ifelse(structure$node == 1, "sparse",
                              ifelse(structure$node == 2, "dense",
                                     "intermediate"))
setwd(paste(rprojroot::find_rstudio_root_file(),"/Output/", modelrun, "/Compartment summary discrete",sep="")) #select any folder with simulation outputs  
series.params <- extractvalues(modstrs, colnames, structure)
series.params <- rename(series.params, replace = c("plotID" = "plotid")) #re-name to avoid issues with grep finding capitalized I

### Read-in information on max distance per index tree, to be merged with dataframes later
setwd(paste(rprojroot::find_rstudio_root_file(),"/Data/Processed",sep=""))
Maxdist <- read.csv("max-distance-from-index.csv", row.names=1, check.names=FALSE)
Maxdist <- rename(Maxdist, replace = c("plotID" = "plotid")) #re-name for merge

##------------------------------------------------------------------------------------------------------
##-------------------------------------Read in simulation summaries-------------------------------------
##------------------------------------------------------------------------------------------------------

unique(series.params$Ntot)
values <- c(288, 2880, 4320) #Select values of interest for plotting
subset <- subsetdata.identifiers(values, "Ntot", series.params, "position")
series.params.subset <- subsetdata.data(values, "Ntot", series.params, "position")

### Read in simulation summaries
## Proportion in each class - takes a bit longer
setwd(paste(rprojroot::find_rstudio_root_file(),"/Output/", modelrun, "/Compartment summary discrete",sep=""))
series <- list.files() #need for every read-in, because file names are different
summary.output.discrete <- read.subset(series, subset, series.params, series.params.subset) 

## Epidemic peak of each simulation
setwd(paste(rprojroot::find_rstudio_root_file(),"/Output/", modelrun, "/Epidemic peak",sep=""))
series <- list.files()
epipeak <- read.subset(series, subset, series.params, series.params.subset)

## Proportion and timing of extinction
## Winter
setwd(paste(rprojroot::find_rstudio_root_file(),"/Output/", modelrun, "/Observed extinction - winter",sep=""))
series <- list.files()
ext.wnt.cats <- read.subset(series, subset, series.params, series.params.subset)
## Year
setwd(paste(rprojroot::find_rstudio_root_file(),"/Output/", modelrun, "/Observed extinction - year",sep=""))
series <- list.files()
ext.yr.cats <- read.subset(series, subset, series.params, series.params.subset)

## Expected proportion of extinction & R0 values, per tree 
setwd(paste(rprojroot::find_rstudio_root_file(),"/Output/", modelrun, "/Expected extinction",sep=""))
series <- list.files()
R0tree.df <- read.subset(series, subset, series.params, series.params.subset)

## Proportion in each class
setwd(paste(rprojroot::find_rstudio_root_file(),"/Output/", modelrun, "/End simulation",sep=""))
series <- list.files()
simend.output <- read.subset(series, subset, series.params, series.params.subset)

## Distance of spread from index case
setwd(paste(rprojroot::find_rstudio_root_file(),"/Output/", modelrun, "/Distance discrete",sep=""))
series <- list.files()
Igroups.discrete <- read.subset(series, subset, series.params, series.params.subset) 
Igroups.discrete <- addmaxdist(Igroups.discrete, Maxdist, "time") #merge dataframes with information on maximum spread distance possible per index tree

##------------------------------------------------------------------------------------------------------
##---------------------------------visualise simulation summaries---------------------------------------
##------------------------------------------------------------------------------------------------------
setwd(paste(rprojroot::find_rstudio_root_file(),"/Output/Figures_", modelrun,sep=""))

### Define model structure for input into all plot functions:
spatmodstr <- "SEIR" #to subset runvalues by spatial model structure  
nonspatmodstr <- "SEIRns" #to subset runvalues by spatial model structure  

### Merge by Beta (combine all theta):
unique(series.params$Beta)
param.values <- unique(series.params$Beta) #select values. Or unique(series.params$Beta)
subset.param <- series.params.subset[series.params.subset$Beta %in% param.values,"newposition"]
param <- "Beta"

## Check numbers of successful cases:
simend.output.data <- success.filter(simend.output, simend.output, threshold) #filter out successful outbreaks
simend.output.merge <- mergelistB(simend.output.data, series.params.subset, param, param.values) #create merged datafiles from chosen subset
for (i in 1:length(simend.output.merge)) { simend.output.merge[[i]]$Group <- ifelse(simend.output.merge[[i]]$modstr==spatmodstr, paste(simend.output.merge[[i]]$structure, sep=""), paste("x homogenous", sep="")) } #create new collumn "group" so that spatial/non-spatial boxes can be presented side-by-side. Code says if model structure is spatial, paste tree structure label. if model structure is non-spatial, paste "homogenous"
successoverview <- ddply(simend.output.merge[[1]], c("modstr","structure", "R0roost", "Beta", "theta", "Ntot", "nsims"), summarise,
                         count.n =sum(!is.na(.n)),
                         mininf = min(max.infections), #check that only simulations with >threshold are included in average
                         average.ext.time = mean(end.time)) %>%
  mutate(prop.success = count.n/nsims) #note that this code doesn't show zero for combinations with zero rows
write.csv(successoverview, "successoverview.csv")

#------------------------Plot proportion and timing of extinction------------------------#
## (Check bars are consistent with successoverview)
## Winter
ext.wnt.filename <- "extinction-winter_structure with fixed Ntot and Beta values_SEIR"  
ext.wnt.cats.merge <- mergelistB(ext.wnt.cats, series.params.subset, param, param.values) #create merged datafiles from chosen subset
for (i in 1:length(ext.wnt.cats.merge)) { ext.wnt.cats.merge[[i]]$Group <- ifelse(ext.wnt.cats.merge[[i]]$modstr==spatmodstr, paste(ext.wnt.cats.merge[[i]]$structure, sep=""), paste("x homogenous", sep="")) } #create new collumn "group" so that spatial/non-spatial boxes can be presented side-by-side. Code says if model structure is spatial, paste tree structure label. if model structure is non-spatial, paste "homogenous"
collabels <- c("Not extinct by day 84","Extinct by day 84","Extinct by day 70","Extinct by day 56","Extinct by day 42","Extinct by day 28","Extinct by day 14", "Outbreak failed")
plot.ext.group(ext.wnt.cats.merge, "theta", "Ntot", "Group", ext.wnt.filename, collabels) #spatial models. Change category labels inside function if plotting over year period (will need to read in year version of data too)
## Year
ext.yr.filename <- "extinction-year_structure with fixed Ntot and Beta values_SEIR"  
ext.yr.cats.merge <- mergelistB(ext.yr.cats, series.params.subset, param, param.values) #create merged datafiles from chosen subset
for (i in 1:length(ext.yr.cats.merge)) { ext.yr.cats.merge[[i]]$Group <- ifelse(ext.yr.cats.merge[[i]]$modstr==spatmodstr, paste(ext.yr.cats.merge[[i]]$structure, sep=""), paste("x homogenous", sep="")) } #create new collumn "group" so that spatial/non-spatial boxes can be presented side-by-side. Code says if model structure is spatial, paste tree structure label. if model structure is non-spatial, paste "homogenous"
collabels <- c("Not extinct by day 360","Extinct by day 360","Extinct by day 300","Extinct by day 240","Extinct by day 180","Extinct by day 120","Extinct by day 60", "Outbreak failed")
plot.ext.group(ext.yr.cats.merge, "theta", "Ntot", "Group", ext.yr.filename, collabels) #spatial models. Change category labels inside function if plotting over year period (will need to read in year version of data too)
#----------------------------------------------------------------------------------------#

#------------------------------Plot spread from index case-------------------------------#
Igroups.data.discrete <- success.filter(Igroups.discrete, simend.output, threshold) #filter out successful outbreaks
Ispread.merge.discrete <- mergelistB(Igroups.data.discrete, series.params.subset, param, param.values) #create merged datafiles from chosen subset

## Extract median and 95% IQR
Ispread.merge.discrete.med <- list()
for (i in 1:length(Ispread.merge.discrete)) {
  Ispread.merge.discrete.med[[i]] <- ddply(Ispread.merge.discrete[[i]], c("modstr","structure", "plotid", "R0roost", "Beta", "theta", "Ntot", "time"), summarise, 
                                         median.spread = median(distance/maxdist),
                                         LIQR = Liqr(distance/maxdist),
                                         UIQR = Uiqr(distance/maxdist)
  )
} #this is the median and IQR across the simulations

## Plot
spread.filename <- "spread-discrete_structure with fixed Ntot and Beta values_SEIR"  
ylab <- "(Normalised) distance from index case"
xlab <- "Time"
mainlab <- "Spread of Infection"
line_ribbon_plot(Ispread.merge.discrete.med, spatmodstr, nonspatmodstr, "median.spread", "time", "LIQR", "UIQR", "structure", "theta", "Ntot", spread.filename, ylab, xlab, mainlab, 90, 1)
#----------------------------------------------------------------------------------------#

#-----------------------------Plot total infecteds over time-----------------------------#
summary.output.data.discrete <- success.filter(summary.output.discrete, simend.output, threshold) #filter out successful outbreaks
summary.output.merge.discrete <- mergelistB(summary.output.data.discrete, series.params.subset, param, param.values) #create merged datafiles from chosen subset

## Extract median and 95% IQR
summary.output.merge.discrete.med <- list()
for (i in 1:length(summary.output.merge.discrete)) {
  summary.output.merge.discrete.med[[i]] <- ddply(summary.output.merge.discrete[[i]], c("modstr","structure", "plotid", "R0roost", "Beta", "theta", "Ntot", "time"), summarise, 
                                         median.I.prop = median(I.prop),
                                         LIQR = Liqr(I.prop), 
                                         UIQR = Uiqr(I.prop)
  )
} #this is the median and IQR across the simulations

## Plot
inf.filename.wnt <- "infection-discrete-winter_structure with fixed Ntot and Beta values_SEIR"  
inf.filename.yr <- "infection-discrete-year_structure with fixed Ntot and Beta values_SEIR"  
ylab <- "Proportion of infected bats"
xlab <- "Time"
mainlab <- "Infections over time"
max(summary.output.merge.discrete.med[[1]]$median.I.prop) #check max y value
line_ribbon_plot(summary.output.merge.discrete.med, spatmodstr, nonspatmodstr, "median.I.prop", "time", "LIQR", "UIQR", "structure", "theta", "Ntot", inf.filename.wnt, ylab, xlab, mainlab, 90, 0.1)
line_ribbon_plot(summary.output.merge.discrete.med, spatmodstr, nonspatmodstr, "median.I.prop", "time", "LIQR", "UIQR", "structure", "theta", "Ntot", inf.filename.yr, ylab, xlab, mainlab, 365, 0.1)
#----------------------------------------------------------------------------------------#

#--------------------Box plots to visualise character of epidemic peak-------------------#
epipeak.data <- success.filter(epipeak, simend.output, threshold) #filter out successful outbreaks
epipeak.merge <- mergelistB(epipeak.data, series.params.subset, param, param.values) #create merged datafiles from chosen subset
for (i in 1:length(epipeak.merge)) { epipeak.merge[[i]]$Group <- ifelse(epipeak.merge[[i]]$modstr==spatmodstr, paste(epipeak.merge[[i]]$structure, sep=""), paste("x homogenous", sep="")) } #create new collumn "group" so that spatial/non-spatial boxes can be presented side-by-side. Code says if model structure is spatial, paste tree structure label. if model structure is non-spatial, paste "homogenous"

## Time of epidemic peak:
time.box.filename <- "time to epidemic peak_structure with fixed Ntot and Beta values_SEIR"  
ylab <- "Time of epidemic peak"  
xlab <- "Tree structure"
mainlab <- "Time of epidemic peak"
max(epipeak.merge[[1]]$time) #check max y value
box_plot_combined(epipeak.merge, spatmodstr, nonspatmodstr, "time", "Group", "theta", "Ntot", time.box.filename, ylab, xlab, mainlab, 300) #spatial and non-spatial models

## Magnitude of epidemic peak:
mag.box.filename <- "magnitude of epidemic peak_structure with fixed Ntot and Beta values_SEIR"  
ylab <- "Proportion of infections at epidemic peak"  
xlab <- "Tree structure"
mainlab <- "Magnitude of epidemic peak"
max(epipeak.merge[[1]]$I.prop) #check max y value
box_plot_combined(epipeak.merge, spatmodstr, nonspatmodstr, "I.prop", "Group", "theta", "Ntot", mag.box.filename, ylab, xlab, mainlab, 0.1) #spatial models only #the epidemic peak is the time when the %I (first) reaches it's highest. The magnitude boxplot should also be %I. Can't be new infections, because if a group reach their peak later, the total cumulative new infections will be higher by default

## Duration of epidemic peak:
summary.output.data <- success.filter(summary.output.discrete, simend.output, threshold) #filter out successful outbreaks
summary.output.merge <- mergelistB(summary.output.data, series.params.subset, param, param.values) #create merged datafiles from chosen subset
for (i in 1:length(summary.output.merge)) { summary.output.merge[[i]]$Group <- ifelse(summary.output.merge[[i]]$modstr==spatmodstr, paste(summary.output.merge[[i]]$structure, sep=""), paste("x homogenous", sep="")) } #create new collumn "group" so that spatial/non-spatial boxes can be presented side-by-side. Code says if model structure is spatial, paste tree structure label. if model structure is non-spatial, paste "homogenous"

epiduration <- list()
for (i in 1:length(summary.output.merge)) {
  datasubset <- summary.output.merge[[i]] #save separately to be compatable with ddply group + non-flexible data call in function. Data name in function is "datasubset"
  epiduration[[i]] <- ddply(datasubset, c("modstr","Group", "plotid", "R0roost", "Beta", "theta", "Ntot", ".n"), summarise,
                            LIQR = Liqr(I.prop),
                            UIQR = Uiqr(I.prop), 
                            LIQR.time = Liqr(time),
                            UIQR.time = Uiqr(time),
                            duration = UIQR.time-LIQR.time
  )
}

duration.box.filename <- "duration of epidemic peak_structure with fixed Ntot and Beta values_SEIR"  
ylab <- "Duration of epidemic peak (time in days)"  
xlab <- "Tree structure"
mainlab <- "Duration of epidemic peak"
max(epiduration[[1]]$duration) #check max y value
box_plot_combined(epiduration, spatmodstr, nonspatmodstr, "duration", "Group", "theta", "Ntot", duration.box.filename, ylab, xlab, mainlab, 200)


#------------------------Successful outbreaks by %trees affected-------------------------#
#View(summary.output.merge.discrete) #Relies on merge output from infection vs time plot. Need to add Group with: 
for (i in 1:length(summary.output.merge.discrete)) { summary.output.merge.discrete[[i]]$Group <- ifelse(summary.output.merge.discrete[[i]]$modstr==spatmodstr, paste(summary.output.merge.discrete[[i]]$structure, sep=""), paste("x homogenous", sep="")) } 

## Extract median and 95% IQR
summary.output.merge.discrete.med.proptrees <- list()
for (i in 1:length(summary.output.merge.discrete)) {
  summary.output.merge.discrete.med.proptrees [[i]] <- ddply(summary.output.merge.discrete[[i]], c("modstr","Group", "Beta", "theta", "Ntot", "time"), summarise, 
                                                  median.I.group.prop = median(I.group.prop),
                                                  LIQR = Liqr(I.group.prop),
                                                  UIQR = Uiqr(I.group.prop)
  )
} #this is the median and IQR across the simulations

## Plot
inftree.disc.filename.wnt <- "By trees infected-discrete time-winter_structure with fixed Ntot and Beta values_SEIR"  
inftree.disc.filename.yr <- "By trees infected-discrete time-year_structure with fixed Ntot and Beta values_SEIR"  
ylab <- "Proportion of infected trees"
xlab <- "Time"
mainlab <- "Tree infections over time"
line_ribbon_plot(summary.output.merge.discrete.med.proptrees, spatmodstr, nonspatmodstr, "median.I.group.prop", "time", "LIQR", "UIQR", "Group", "theta", "Ntot", inftree.disc.filename.wnt, ylab, xlab, mainlab, 90, 1.05)
line_ribbon_plot(summary.output.merge.discrete.med.proptrees, spatmodstr, nonspatmodstr, "median.I.group.prop", "time", "LIQR", "UIQR", "Group", "theta", "Ntot", inftree.disc.filename.yr, ylab, xlab, mainlab, 360, 1.05)
# -------------------------------------------- #

##------------------------------------------------------------------------------------------------------
##-----------------------Save summaries per dataset and variable combination----------------------------
##------------------------------------------------------------------------------------------------------
## The following code is written so that true averages can be calculated across groupings of combinations (e.g. spatial vs non-spatial, including all different Ntots in the average)
## The code is written so that values calculated per row (e.g. time of simulation end) will be averaged with the proper demominator. i.e. value/number of rows, not average value/number of groups)
## For values that are calculated once per group (e.g. proportion of all simulaions that were successful) the true denominator is the the number of groups

# -------------Time to extinction------------- #
## Average per bar
filename <- "SEIR_extinction time_by bar"  
grouplist <- c("modstr","structure", "plotid", "R0roost", "Beta", "theta", "Ntot", "nsims")
Nsuccessul <- ext.time.summary(simend.output.merge, grouplist, filename) #this code gives number of successful outbreaks, to calculate proportion of successful outbreaks below. If the origional dataframe has more than 1 list, the way this is saved will need to change to accomadate multiple lists

## Average per tree structure
filename <- "SEIR_extinction time_by tree"  
grouplist <- c("modstr","structure", "plotid", "Beta", "theta", "nsims")
ext.time.summary(simend.output.merge, grouplist, filename) 

## Average per abundance
filename <- "SEIR_extinction time_by Ntot"  
grouplist <- c("Ntot", "Beta", "theta", "nsims")
ext.time.summary(simend.output.merge, grouplist, filename) 

## Average per model structure
filename <- "SEIR_extinction time_by modstr"  
grouplist <- c("modstr", "Beta", "theta", "nsims")
ext.time.summary(simend.output.merge, grouplist, filename) 

## Average per model structure & population
filename <- "SEIR_extinction time_by modstr & Ntot"  
grouplist <- c("modstr", "Ntot", "Beta", "theta", "nsims")
ext.time.summary(simend.output.merge, grouplist, filename) 
# -------------------------------------------- #

# -----Proportion of successful outbreaks----- #
## Value per bar
Nsuccessul$prop.successful <- Nsuccessul$count.n/Nsuccessul$nsims #smallest level for the proportion of successful outbreaks (because calculated from all simulations per simulation run)

## Average per tree structure
filename <- "SEIR_prop extinction_by tree.csv"  
grouplist <- c("modstr","structure", "plotid", "Beta", "theta", "nsims")
success.summary(Nsuccessul, grouplist, filename) 

## Average per abundance
filename <- "SEIR_prop extinction_by Ntot.csv"  
grouplist <- c("Ntot", "Beta", "theta", "nsims")
success.summary(Nsuccessul, grouplist, filename) 

## Average per model structure 
filename <- "SEIR_prop extinction_by modstr.csv"  
grouplist <- c("modstr", "Beta", "theta", "nsims")
success.summary(Nsuccessul, grouplist, filename) 

## Average per model structure & population
filename <- "SEIR_prop extinction_by modstr & Ntot.csv"  
grouplist <- c("modstr", "Ntot", "Beta", "theta", "nsims")
success.summary(Nsuccessul, grouplist, filename) 
# -------------------------------------------- #

# -----Magnitude and time of epidemic peak---- #
## (note - repace tree structure with Group to get both spatial and non-spatial comparisons)
## Average per box
filename <- "SEIR_epidemic mag and time_by bar"  
grouplist <- c("modstr","Group", "R0roost", "Beta", "theta", "Ntot", "nsims")
epi.MagTime.summary(epipeak.merge, grouplist, filename) #this code gives number of successful outbreaks, to calculate proportion of successful outbreaks below. If the origional dataframe has more than 1 list, the way this is saved will need to change to accomadate multiple lists

## Average per Group structure
filename <- "SEIR_epidemic mag and time_by Group"  
grouplist <- c("Group", "Beta", "theta", "nsims")
epi.MagTime.summary(epipeak.merge, grouplist, filename) 

## Average per abundance
filename <- "SEIR_epidemic mag and time_by Ntot"  
grouplist <- c("Ntot", "Beta", "theta", "nsims")
epi.MagTime.summary(epipeak.merge, grouplist, filename) 

## Average per model structure
filename <- "SEIR_epidemic mag and time_by modstr"  
grouplist <- c("modstr", "Beta", "theta", "nsims")
epi.MagTime.summary(epipeak.merge, grouplist, filename) 

## Average per model structure & population
filename <- "SEIR_epidemic mag and time_by modstr & Ntot"  
grouplist <- c("modstr", "Ntot", "Beta", "theta", "nsims")
epi.MagTime.summary(epipeak.merge, grouplist, filename) 
# -------------------------------------------- #

# ----------Duration of epidemic peak--------- #
## Calculate duration of epidemic peak per .n
epiduration.n <- list()
for (i in 1:length(summary.output.merge)) {
  datasubset <- summary.output.merge[[i]] #save separately to be compatable with ddply group + non-flexible data call in function. Data name in function is "datasubset"
  epiduration.n[[i]] <- ddply(datasubset, c(".n","modstr","Group", "plotid", "R0roost", "Beta", "theta", "Ntot", ".n"), summarise,
                            LIQR = Liqr(I.prop),
                            UIQR = Uiqr(I.prop), 
                            LIQR.time = Liqr(time),
                            UIQR.time = Uiqr(time),
                            duration = UIQR.time-LIQR.time
  )
}

## (note - repace tree structure with Group to get both spatial and non-spatial comparisons)
## Average per box
filename <- "SEIR_epidemic duration_by box"  
grouplist <- c("modstr","Group", "R0roost", "Beta", "theta", "Ntot")
epi.duration.summary(epiduration.n, grouplist, filename) #this code gives number of successful outbreaks, to calculate proportion of successful outbreaks below. If the origional dataframe has more than 1 list, the way this is saved will need to change to accomadate multiple lists

## Average per Group structure
filename <- "SEIR_epidemic duration_by Group"  
grouplist <- c("Group", "plotid", "Beta", "theta")
epi.duration.summary(epiduration.n, grouplist, filename) 

## Average per abundance
filename <- "SEIR_epidemic duration_by Ntot"  
grouplist <- c("Ntot", "Beta", "theta")
epi.duration.summary(epiduration.n, grouplist, filename) 

## Average per model structure
filename <- "SEIR_epidemic duration_by modstr"  
grouplist <- c("modstr", "Beta", "theta")
epi.duration.summary(epiduration.n, grouplist, filename) 

## Average per model structure & population
filename <- "SEIR_epidemic duration_by modstr & Ntot"  
grouplist <- c("modstr", "Ntot", "Beta", "theta")
epi.duration.summary(epiduration.n, grouplist, filename)
# -------------------------------------------- #

# -------- Days at 100% trees infected-------- #
Itreesummary <- list()
for (i in 1:length(summary.output.merge.discrete.med.proptrees)) {
  datasubset <- summary.output.merge.discrete.med.proptrees[[i]] #save separately to be compatable with ddply group + non-flexible data call in function. Data name in function is "datasubset"
  Itreesummary[[i]] <- ddply(datasubset, c("modstr","Group", "Beta", "theta", "Ntot"), summarise,
                             max.Igroup.prop = max(median.I.group.prop), 
                             days.at.max = sum(!is.na(median.I.group.prop)[median.I.group.prop==max(median.I.group.prop)])
                             
  )
  write.csv(Itreesummary, paste("SEIR_days with max infected trees_list ", i, ".csv",sep="")) #save summary for later reference
}
# -------------------------------------------- #

# ----------- Table of R0 per tree ----------- #
R0tree.merge <- mergelistB(R0tree.df, series.params.subset, param, param.values) #create merged datafiles from chosen subset
R0tree.table <- list()
for (i in 1:length(R0tree.merge)) {
  R0tree.table[[i]] <- ddply(R0tree.merge[[i]], c("modstr","structure", "R0roost", "Beta", "theta", "Ntot", "nsims", "index"), summarise,
                             R0tree.full = mean(R0tree),
                             R0tree.fround = round(R0tree, digits=2))
  write.csv(R0tree.table, paste("R0 per tree_list ", i, ".csv",sep="")) #save summary for later reference
}
# -------------------------------------------- #

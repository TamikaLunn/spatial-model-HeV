## Title: R script defining helping functions for stochastic compartmental models with gillespie algorithm
## Author: Tamika Lunn, Griffith University
## Version: V20, created 5th May 2020

## Install libraries
library(tidyverse) 
library(ggplot2) 
library(cowplot) 
library(reshape2) 
library(igraph) 
library(data.table) 
library(R.utils) 

## Functions:
##----------------------------------------------------------------------------------------------------
##------------------------------------------------General---------------------------------------------
##----------------------------------------------------------------------------------------------------
## The function(s) (SIR.R0, SIR.Bcrit, SIR.B) contain equations for calculating, respectively: Average R0 for the roost, based on the B matrix; Critical B where R0=1; Overall within-roost B value given an average R0roost value. Note that scaled.data needs to contain 1 values (not zero values) for within-tree distances (i.e along the diagonal, or for trees closer than the given radius)
## The function (SIRsummary) calculates the sum and proportion of individuals in the roost in each category. Length of data is .n*time
## The function (epidemicpeak) return data at the time of the epidemic peak (first time %I reaches it's max value). Length of data is .n
## The function (extinction) calculates the proportion and timing of extinction from simulation output. Threshold= 10 infected individuals by end of simulation
## The function (expected) calculates and returns the expected proportion of extinction, R0 values and centrality (closeness) values, per tree, in a dataframe
## The function (simend) extracts values from end of simualtion, and calculates the proportion of trees that were infected over the simulation. Use maxinfections to create list of successful simulations
## The functin (successful) return a list of successful epidemics. Threshold= 10 infected individuals by end of simulation
## The function (groups) return infections per group over simulation @ each simulation (note - not a tally of infections) and the max distance of infection from the index case at each time step
## The function(s) (mergelistA and mergelistB) merges dataframes in list, based on a pre-defined subset, and returns the merged dataframes as a list. A includes line to remove extra I categories, B does not remove I categories
## The function(s) (subsetdata.identifiers and subsetdata.data) subset dataframe from dataframe, based on a pre-defined subset. Returns either the identifying values of subset, or all data in subset
## The function (addmaxdist) merges simulation dataframes with maxdistance dataframe (maximum distance of spread, per index case). Loops through entire list of dataframes, and merges dataframes one at a time by index and plotid
## The function(s) (Numextract and NumextractB) extract number values from strings. Note that Numextract B is specifically for extracting subplot values from single digit strings, and re-formatting with the 000 format
## The function(s) (read.all and read.subset) reads in data from a list of csv files and attachs corresponding parameter values from series.params dataframe. read.subset reads in all data and then deletes those not in subset. Two id values are given to make sure that parameter values are assigned to the correct dataframe: idseries=id taken from the file name, id=id taken from series.params dataframe. Position (all) or newposition (subset) are assigned for selecting dataframes in list
## The function (findloc) finds the position of value in dataframe that most closely matches a given, single value
## The function(s) (read.DATA, read.site, read.subplot) reads in distance matrix data as a list for stochastic simulations. List should be the length of the number of different data matrices being read in. The function is repeated three times to return the different elements (DATA, site and subplot) from the function, each as a list
## The function(s) (setRunvalue.fixedBeta, setRunvalue.fixedR0roost.initial, setRunvalue.fixedR0roost) create a parameter table (data frame) for running values using either fixed B or fixed R0roost. R0roost is split into two functions. The first (setRunvalue.fixedR0roost.initial) is a stepping stone, the second (setRunvalue.fixedR0roost) is the one to use in the main program
## The function(s) (setmatrix_) returns a list of beta matrices to be used in the main program. This should be the length of nrow of the parameter table, and with each dataframe in the list corresponding to the row in the parameter dataframe
## The function(s) (rescaleData.SCALED.DATA, rescaleData.SCALED.DATA.ALT, rescaleData.SCALED.DATA.DISTANCE) re-scales the data by radius and returns a list of dataframes to be used in the main program. This should be the length of nrow of the parameter table, and with each dataframe in the list corresponding to the row in the parameter dataframe. The function is repeated three times to return the different elements (SCALED.DATA, SCALED.DATA.ALT, SCALED.DATA.DISTANCE) from the function
##     *Note* naming convention throughout, DATA = list of all dataframes, Data or runvalues$Data = the position of the dataframe in the list, data = single df containing the actual data. Ditto caps convention for BETA/Beta/B
## The function (extractvalues) Creates a table of parameter values. Specifically, it: Extracts parameter values from file names, adds information on tree structure (merge from Ctree outputs), and Adjusts numeric/factor/character categories if needed, and sets the order by position
## The function (discrete_dataframe)
## The function (thin_sim) thins output to a single simulation


##---------------------------------Function(s): thin output ---------------------------------
# sim is an output array
# save_t is a vector of time points when the values must be saved
thin_sim <- function(sim,save_t)
{
  sim_t <- sim[,"time"]
  n_t <- length(sim_t)
  # Trim save_t after the last simulation event.
  save__t <- save_t[1:min(which(save_t > sim_t[n_t]))]
  # Select the rows that need to be saved
  select <- sapply(save__t, function(t){max(which(sim_t<=t))})
  out <- sim[select,]
  # Update time stamps
  out[,"time"] <- save__t
  return(out)  
}

##---------------------------------Function(s): equations for SIR models---------------------------------
# R0:
SIR.R0tree <- function(n, Ntree, Bmatrix, gamma) {
R0tree <- sapply(1:n,function(i) {sum(B[i,]*Ntree[i])/gamma}) 
R0roost <- sum(Ntree*R0tree)/sum(Ntree) 
return(R0tree)
}

SIR.R0roost <- function(n, Ntree, Bmatrix, gamma) {
  R0tree <- sapply(1:n,function(i) {sum(B[i,]*Ntree[i])/gamma}) 
  R0roost <- sum(Ntree*R0tree)/sum(Ntree) 
  return(R0roost)
}

# Critical Beta from gravity model:
SIR.Bcrit.gravity <- function(n, Ntree, scaled.data.alt, gamma, theta) {
tree <- sapply(1:n,function(i)  { sum(Ntree[i]/(theta*scaled.data.alt[i,])^2) }) #inf values when divided by zero
Broost <- gamma*(sum(Ntree)/sum(Ntree*tree))
return(Broost)
}

# Beta from R0 from gravity model:
SIR.B.gravity <- function(n, Ntree, scaled.data.alt, gamma, R0roost, theta) {
tree <- sapply(1:n,function(i)  { sum(Ntree[i]/(theta*scaled.data.alt[i,])^2) }) #inf values when divided by zero
Broost <- (gamma*R0roost)*(sum(Ntree)/sum(Ntree*tree))  
return(Broost)
}

# Critical Beta from plus-one model:
SIR.Bcrit.plusone <- function(n, Ntree, data.long, gamma, theta) {
  tree <- sapply(1:n,function(i)  { sum(Ntree[i]/(data.long[i,]/theta)+1) })
  Broost <- gamma*(sum(Ntree)/sum(Ntree*tree))
  return(Broost)
}

# Beta from R0 from plus-one model:
SIR.B.plusone <- function(n, Ntree, data.long, gamma, R0roost, theta) {
  tree <- sapply(1:n,function(i)  { sum(Ntree[i]/(data.long[i,]/theta)+1) }) #inf values when divided by zero
  Broost <- (gamma*R0roost)*(sum(Ntree)/sum(Ntree*tree))  
  return(Broost)
}

# Beta from R0 for a well mixed population:
SIR.B.homo <- function(Ntot, gamma, R0roost) {
  Broost <- (R0roost*gamma)/Ntot
  return(Broost)
}

##---------------------------------Function: Calculate sum and proportion of individuals in S,I,R groups
SIRsummary <- function(simdat, Ntree) {
#compile tree groups from each compartment
  S.cols <- grep("S", colnames(simdat))
  I.cols <- grep("I", colnames(simdat))
  R.cols <- grep("R", colnames(simdat))
#calculate sum and proportion in each compartment per time step
  summary.output <- simdat %>%
    dplyr::mutate(
    S.tot = rowSums(simdat[,S.cols]),
    I.tot = rowSums(simdat[,I.cols]),
    R.tot = rowSums(simdat[,R.cols]),
    S.prop = S.tot/sum(Ntree),
    I.prop = I.tot/sum(Ntree),
    R.prop = R.tot/sum(Ntree),
    S.group.prop = rowSums(!isZero(simdat)[,S.cols])/length(S.cols), #proportion of trees with at least 1 susceptible bat
    I.group.prop = rowSums(!isZero(simdat)[,I.cols])/length(I.cols),
    R.group.prop = rowSums(!isZero(simdat)[,R.cols])/length(R.cols)) 
#remove groups to speed up writing and reading data
return(summary.output[,c(".n", "time", "index", "infections", "S.tot", "I.tot", "R.tot", "S.prop", "I.prop", "R.prop", "S.group.prop", "I.group.prop", "R.group.prop")])
}

##---------------------------------Function: Calculate sum and proportion of individuals in S,E,I,R groups
SEIRsummary <- function(simdat, Ntree) {
  #compile tree groups from each compartment
  S.cols <- grep("S", colnames(simdat))
  E.cols <- grep("E", colnames(simdat))
  I.cols <- grep("I", colnames(simdat))
  R.cols <- grep("R", colnames(simdat))
  #calculate sum and proportion in each compartment per time step
  summary.output <- simdat %>%
    dplyr::mutate(
      S.tot = rowSums(simdat[,S.cols]),
      E.tot = rowSums(simdat[,E.cols]),
      I.tot = rowSums(simdat[,I.cols]),
      R.tot = rowSums(simdat[,R.cols]),
      S.prop = S.tot/sum(Ntree),
      E.prop = E.tot/sum(Ntree),
      I.prop = I.tot/sum(Ntree),
      R.prop = R.tot/sum(Ntree),
      S.group.prop = rowSums(!isZero(simdat)[,S.cols])/length(S.cols), #proportion of trees with at least 1 susceptible bat
      E.group.prop = rowSums(!isZero(simdat)[,E.cols])/length(E.cols),
      I.group.prop = rowSums(!isZero(simdat)[,I.cols])/length(I.cols),
      R.group.prop = rowSums(!isZero(simdat)[,R.cols])/length(R.cols)) 
  #remove groups to speed up writing and reading data
  return(summary.output[,c(".n", "time", "index", "infections", "S.tot", "I.tot", "R.tot", "S.prop", "I.prop", "R.prop", "S.group.prop", "E.group.prop", "I.group.prop", "R.group.prop")])
  }

##---------------------------------Function: return data at time of epidemic peak per simulation---------------------------------
epidemicpeak <- function(summary.output) {
  epidemicpeak.output <- summary.output %>%
    group_by(.n,index) %>%
    arrange(time) %>% #ensures it takes the first max only
    slice(which.max(I.prop)) %>%
    as.data.frame(x)
  return(epidemicpeak.output[,c(".n", "time", "index", "infections", "I.prop", "I.group.prop")])
}

##---------------------------------Function: calculate proportion and timing of extinction---------------------------------
extinction <- function(simdat, threshold, Max, By) { #for winter period set: Max=84, By=14; for year period set: Max=360, By=60
  ext.cats <- simdat %>%
  ddply(c(".n", "index"), summarise,
        max = max(infections), #number of infections at end of each simulation
        end.time = max(time)) %>% #time at end of simulation
  mutate(ext.time = findInterval(end.time, seq(0, Max, by=By), rightmost.closed = TRUE)) %>% #categorised end times 
  ddply(c("index"), summarise, 
        sim.fail = sum(!is.na(max)[max<=threshold])/sum(!is.na(max)), #proportion of simulations that failed
        ext.1 = sum(!is.na(ext.time)[ext.time==1 & max>threshold])/sum(!is.na(ext.time)), #proportion of all simulations that went extinct inside category 1. "Extinct by" category 2
        ext.2 = sum(!is.na(ext.time)[ext.time==2 & max>threshold])/sum(!is.na(ext.time)),
        ext.3 = sum(!is.na(ext.time)[ext.time==3 & max>threshold])/sum(!is.na(ext.time)),
        ext.4 = sum(!is.na(ext.time)[ext.time==4 & max>threshold])/sum(!is.na(ext.time)),
        ext.5 = sum(!is.na(ext.time)[ext.time==5 & max>threshold])/sum(!is.na(ext.time)),
        ext.6 = sum(!is.na(ext.time)[ext.time==6 & max>threshold])/sum(!is.na(ext.time)),
        ext.7 = sum(!is.na(ext.time)[ext.time==7 & max>threshold])/sum(!is.na(ext.time))
  ) %>% #keep labelling as ext.1-7 (instead of ext.14 ext.28) to make ordering easier later on
  melt(id.vars = c("index"), measure.vars = c("ext.7", "ext.6", "ext.5", "ext.4", "ext.3", "ext.2","ext.1", "sim.fail"),
       variable.name = c("sim.category"), value.name="prop.sim")
  return(as.data.frame(ext.cats))
}

##---------------------------------Function: compile expected proportion of extinction, R0 values and centrality (closeness) values---------------------------------
expected <- function(data, n, R0tree) {
# create network & calculate centrality scores
  wtmatrix <- as.matrix(data)
  nodes <- as.numeric(colnames(wtmatrix))
  wtlist <- subset(melt(wtmatrix), value!=0)
  colnames(wtlist)[colnames(wtlist)=="Var1"] <- "from"
  colnames(wtlist)[colnames(wtlist)=="Var2"] <- "to"
  colnames(wtlist)[colnames(wtlist)=="value"] <- "weight"
  wtlinks <- as.matrix(wtlist)
  network <- graph_from_data_frame(d = wtlinks, vertices = nodes, directed = FALSE)
  centrality <- closeness(network, weights=wtlist$weight, normalized=T) #closeness = centrality based on distance to others in the network. Normalization is performed by multiplying the raw closeness by n-1, where n is the number of vertices in the graph
# calculate expected extinction values from R0 tree values  
  expect.tree <- array(dim=c(1,n))
  colnames(expect.tree) <- 1:n
  expect.tree[1,] <- (1/R0tree)^1
# arrange into data frame 
  R0tree.df <- array(dim=c(n,4))
  colnames(R0tree.df) <- c("index","R0tree", "centrality","expect.tree")
  R0tree.df[,1] <- 1:n
  R0tree.df[,2] <- R0tree
  R0tree.df[,3] <- centrality
  R0tree.df[,4] <- expect.tree
  R0tree.df <- as.data.frame(R0tree.df)
  R0tree.df <- R0tree.df %>%
    mutate(expect.tree = replace(expect.tree, expect.tree > 1, 1)) #assume expected proportions are 1 when R0 is greater than 1.00 
  return(R0tree.df)
}

##---------------------------------Function: extract values from end of simulation. Use to create list of successful simulations---------------------------------
simend <- function(simdat) {
  simend.output.start <- simdat %>%
    ddply(c(".n", "index"), summarise,
          max.infections = max(infections), #number of infections at end of each simulation
          end.time = max(time)) #time at end of simulation
  ## Calculate total proportion of trees infected during simulations:
  ##select tree groups from each I compartment
  I.cols <- grep("I", colnames(simdat))
  ##create dataframe to store values in
  I.trees <- data.frame(matrix(ncol = 2, nrow = length(unique(simdat$.n))))
  colnames(I.trees) <- c(".n", "total.prop.infected.trees")
  for (i in c(unique(simdat$.n))) {
    rows <- which(simdat$.n==i)
    ##identify which trees became infected in the simulation (1 = infected, 0 = never infected)
    infected.trees <- ifelse(colSums(simdat[c(rows),I.cols])==0,0,1)
    total.prop.infected.trees <- sum(infected.trees)/length(infected.trees)
    I.trees[i,".n"] <- i
    I.trees[i,"total.prop.infected.trees"] <- total.prop.infected.trees
  }
  ## Merge with total proportion of infected trees 
  simend.output <- merge(simend.output.start, I.trees, by=".n")
  return(simend.output)
  }

##---------------------------------Function: return list of successful epidemics, given a defined threshold---------------------------------
## When using this function in a loop (looped by the length of the dataset), need to make sure that the length of simend.output matches the length of the dataset. If not, the returned dataframe will only be a subset of the data you really want
## But because .n isn't unique to each combination, it needs to be non-merged datasets that are filtered (a successful .n in one simulation set might not be a successful .n in another, but if you run this on merged datasets all matching .n will be pulled)
successful <- function(simend.output, threshold) {
  successful.output <- simend.output %>%
    filter(max.infections>threshold)
    return(as.data.frame(successful.output))
} #Then: filter(.n %in% successful.output$.n)

nonsuccessful <- function(simend.output, threshold) {
  nonsuccessful.output <- simend.output %>%
    filter(max.infections<=threshold)
  return(as.data.frame(nonsuccessful.output))
} 


##---------------------------------Function: return infections per group over simulation @ each simulation 
groups <- function(simdat, scaled.data, data) {
  #select just I groups to return
  Is <- simdat[,c(grep("I", colnames(simdat)))]
  other <- simdat[,c(".n", "time", "index", "infections")]
  Igroups <- cbind(other, Is) #potential bug, if cbind doesn't preserve row order
  #if Igroup value is >0 fill with data or scaled.data value of that tree from the index case
  scdistance <- ifelse(Igroups[,c(grep("I", colnames(Igroups)))]>0,scaled.data[Igroups$index,],0)
  distance <- ifelse(Igroups[,c(grep("I", colnames(Igroups)))]>0,data[Igroups$index,],0)
  #take the max distance value for each time step (row) and attach it to the df
  Igroups[, "scdistance"] <- apply(scdistance, 1, max)
  Igroups[, "distance"] <- apply(distance, 1, max)
  #remove groups to speed up writing and reading data
  return(Igroups[,c(".n", "time", "index", "infections","scdistance", "distance")])
  #return(Igroups)
}

##---------------------------------Function: merge dataframes in list, based on a pre-defined subset---------------------------------
## A includes line to remove extra I categories, B does not remove I categories
mergelistA <- function(data, series.params, param, subset.ids, parameter.values) { #data = df of summarised output data, series.params = df of parameter values, param = col name of parameter of interest (in series.params), subset.ids = vector with the ids of dfs with the desired parameter value, parameter.values = vector of the desired parameter values
  newdata <- data #to avoid over-riding raw data
  # Remove all I groups for merge (if present. Be cautious - will remove all colnames with a capitalized I)
  for (i in c(subset.ids)) { #repeat through idenfidied dfs
    newdata[[i]] <- newdata[[i]][,-c(grep("I", colnames(newdata[[i]])))] 
  } 
  # Merge dataframes that have same subset values
  merge <- list() #create blank list to save output into
  k <- 0 #reset k before loop
  for (i in c(parameter.values)) { #repeat through identified parameter values
    k <- k+1
    subset <- series.params[which(series.params[param]==i),]$newposition #identify which dfs have the specified parameter values, and repeat the merge through these parameter values one by one
    merge[[k]] <- rbindlist(newdata[newposition=subset[c(1:length(subset))]]) #working from inside out: subset contains the id/newposition values of the dfs to be merged. subset[] idenfities which df should be merged (1:nth position in subset). Newdata[] extracts the identified df from the list of dfs. rbindlist binds all of the extracted dfs together by row
  } 
  return(merge)  
} 


mergelistB <- function(data, series.params, param, parameter.values) { #data = df of summarised output data, series.params = df of parameter values, param = col name of parameter of interest (in series.params), subset.ids = vector with the ids of dfs with the desired parameter value, parameter.values = vector of the desired parameter values
  newdata <- data #to avoid over-riding raw data
# Merge dataframes that have same subset values
  merge <- list() #create blank list to save output into
  k <- 0 #re-set k before loop
  for (i in c(parameter.values)) { #repeat through identified parameter values
    k <- k+1 #for whatever reason, code will no longer save values like list[[1.5]]. Have introduced k into loop so that merges will save separately
    subset <- series.params[which(series.params[param]==i),]$newposition #identify which dfs have the specified parameter values, and repeat the merge through these parameter values one by one
    merge[[k]] <- rbindlist(newdata[newposition=subset[c(1:length(subset))]]) #working from inside out: subset contains the id/newposition values of the dfs to be merged. subset[] idenfities which df should be merged (1:nth position in subset). Newdata[] extracts the identified df from the list of dfs based on id. rbindlist binds all of the extracted dfs together by row
  } 
  return(merge)  
}

##---------------------------------Function: subset dataframe from dataframe, based on a pre-defined subset. Return either the identifying values of subset, or all data in subset---------------------------------
subsetdata.identifiers <- function(values, param, dataframe, identifier) {
  subset <- dataframe[dataframe[[param]] %in% values,identifier]
  return(subset)
}

subsetdata.data <- function(values, param, dataframe, identifier) {
  subset <- dataframe[dataframe[[param]] %in% values,identifier]
  dataframe.subset <- dataframe[dataframe[[param]] %in% values,]
  dataframe.subset$newposition <- 1:nrow(dataframe.subset) 
  return(dataframe.subset)
}

##---------------------------------Function: extract number values from strings-------------------------------
Numextract <- function(string){
  unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
}

## specifically for extracting subplot values from single digit strings, and re-formatting with the 000 format
NumextractB <- function(string){
  val <- unlist(regmatches(string,gregexpr("[[:digit:]]+",string)))
  newval <- paste("00", val, sep="") 
  newval <- revalue(newval, c("0010"="010"))
  return(newval)
}
  
##------------------------------Function: read in data (all or subset) and attach parameter values-------------------
read.all <- function(series, series.params) {
  data <- list()
  for (i in 1:length(series)) { #read in all data and attach id values
    data[[i]] <- read_csv(series[[i]])
    data[[i]]$idseries <- Numextract(str_extract(series[[i]],"id=...")) #take from real series values, so have a point of reference to check. Others are taken from series.params because this is much easier 
    data[[i]]$id <- series.params[which(series.params$position==i),"id"] #should be the same as idseries
    data[[i]]$modstr <- series.params[which(series.params$position==i),"modstr"]
    data[[i]]$site <- series.params[which(series.params$position==i),"site"]
    data[[i]]$subplot <- series.params[which(series.params$position==i),"subplot"]
    data[[i]]$plotposition <- series.params[which(series.params$position==i),"plot"]
    data[[i]]$plotid <- series.params[which(series.params$position==i),"plotid"]
    data[[i]]$n <- series.params[which(series.params$position==i),"n"]
    data[[i]]$N <- series.params[which(series.params$position==i),"N"]
    data[[i]]$Ntot <- series.params[which(series.params$position==i),"Ntot"]
    data[[i]]$Beta <- series.params[which(series.params$position==i),"Beta"]
    data[[i]]$gamma <- series.params[which(series.params$position==i),"gamma"]
    data[[i]]$radius <- series.params[which(series.params$position==i),"radius"]
    data[[i]]$theta <- series.params[which(series.params$position==i),"theta"]
    data[[i]]$threshold <- series.params[which(series.params$position==i),"threshold"]
    data[[i]]$nsims <- series.params[which(series.params$position==i),"nsims"]
    data[[i]]$R0roost <- series.params[which(series.params$position==i),"R0roost"]
    data[[i]]$structure <- series.params[which(series.params$position==i),"structure"]
    data[[i]]$meandist <- series.params[which(series.params$position==i),"meandist"]
    #data[[i]]$maxdist <- series.params[which(series.params$position==i),"maxdist"] #needs to be maxdist unique to index...
    data[[i]]$position <- series.params[which(series.params$position==i),"position"]
  }
  return(data)
 } 

read.subset <- function(series, subset, series.params, series.params.subset) {
  data <- list()
  for (i in c(subset)) { #read in all data and attach id values. Attach ids from file name and parameter value df to make sure loop is working correctly
    data[[i]] <- read_csv(series[[i]])
    data[[i]]$idseries <- Numextract(str_extract(series[[i]],"id=...")) #take from real series values, so have a point of reference to check. Others are taken from series.params because this is much easier 
    data[[i]]$id <- series.params[which(series.params$position==i),"id"] #should be the same as idseries. And yes, should be series.params, not series.params subset, because at this stage is going through entire list and just skipping those not selected
    data[[i]]$modstr <- series.params[which(series.params$position==i),"modstr"]
    data[[i]]$site <- series.params[which(series.params$position==i),"site"]
    data[[i]]$subplot <- series.params[which(series.params$position==i),"subplot"]
    data[[i]]$plotid <- series.params[which(series.params$position==i),"plotid"]
    data[[i]]$plotposition <- series.params[which(series.params$position==i),"plot"]
    data[[i]]$n <- series.params[which(series.params$position==i),"n"]
    data[[i]]$N <- series.params[which(series.params$position==i),"N"]
    data[[i]]$Ntot <- series.params[which(series.params$position==i),"Ntot"]
    data[[i]]$Beta <- series.params[which(series.params$position==i),"Beta"]
    data[[i]]$gamma <- series.params[which(series.params$position==i),"gamma"]
    data[[i]]$radius <- series.params[which(series.params$position==i),"radius"]
    data[[i]]$theta <- series.params[which(series.params$position==i),"theta"]
    data[[i]]$threshold <- series.params[which(series.params$position==i),"threshold"]
    data[[i]]$nsims <- series.params[which(series.params$position==i),"nsims"]
    data[[i]]$R0roost <- series.params[which(series.params$position==i),"R0roost"]
    data[[i]]$structure <- series.params[which(series.params$position==i),"structure"]
    data[[i]]$meandist <- series.params[which(series.params$position==i),"meandist"]
    #data[[i]]$maxdist <- series.params[which(series.params$position==i),"maxdist"]
    data[[i]]$position <- series.params[which(series.params$position==i),"position"]
  }
  data <- data[sapply(data, function(x) length(x)[1]) > 0] #remove dfs with no information (those not in subset)
  for (i in 1:length(data)) { #attach new position value of subset
    data[[i]]$newposition <- series.params.subset[which(series.params.subset$newposition==i),"newposition"]
  }
  return(data)
}

read.subset.rand <- function(series, subset, series.params, series.params.subset) {
  data <- list()
  for (i in c(subset)) { #read in all data and attach id values. Attach ids from file name and parameter value df to make sure loop is working correctly
    data[[i]] <- read_csv(series[[i]])
    data[[i]]$idseries <- Numextract(str_extract(series[[i]],"id=...")) #take from real series values, so have a point of reference to check. Others are taken from series.params because this is much easier 
    data[[i]]$id <- series.params[which(series.params$position==i),"id"] #should be the same as idseries. And yes, should be series.params, not series.params subset, because at this stage is going through entire list and just skipping those not selected
    data[[i]]$modstr <- series.params[which(series.params$position==i),"modstr"]
    #data[[i]]$site <- series.params[which(series.params$position==i),"site"]
    data[[i]]$subplot <- series.params[which(series.params$position==i),"subplot"]
    data[[i]]$plotid <- series.params[which(series.params$position==i),"plotid"]
    #data[[i]]$plotposition <- series.params[which(series.params$position==i),"plot"]
    data[[i]]$n <- series.params[which(series.params$position==i),"n"]
    data[[i]]$N <- series.params[which(series.params$position==i),"N"]
    data[[i]]$Ntot <- series.params[which(series.params$position==i),"Ntot"]
    data[[i]]$Beta <- series.params[which(series.params$position==i),"Beta"]
    data[[i]]$gamma <- series.params[which(series.params$position==i),"gamma"]
    data[[i]]$radius <- series.params[which(series.params$position==i),"radius"]
    data[[i]]$theta <- series.params[which(series.params$position==i),"theta"]
    data[[i]]$threshold <- series.params[which(series.params$position==i),"threshold"]
    data[[i]]$nsims <- series.params[which(series.params$position==i),"nsims"]
    data[[i]]$R0roost <- series.params[which(series.params$position==i),"R0roost"]
    #data[[i]]$structure <- series.params[which(series.params$position==i),"structure"]
    #data[[i]]$meandist <- series.params[which(series.params$position==i),"meandist"]
    #data[[i]]$maxdist <- series.params[which(series.params$position==i),"maxdist"]
    data[[i]]$position <- series.params[which(series.params$position==i),"position"]
  }
  data <- data[sapply(data, function(x) length(x)[1]) > 0] #remove dfs with no information (those not in subset)
  for (i in 1:length(data)) { #attach new position value of subset
    data[[i]]$newposition <- series.params.subset[which(series.params.subset$newposition==i),"newposition"]
  }
  return(data)
}

##---------------------------Function: Add information on maxdistance (of spread) per df-------------------
addmaxdist <- function(data, Maxdist, ordername) {
  for (i in 1:length(data)) {
    data[[i]] <- merge(Maxdist[,c("index", "plotid", "maxdist")],data[[i]], by=c("plotid","index"))
    data[[i]] <- data[[i]][order(data[[i]]$.n,data[[i]][[ordername]]),] #to origional order, pre-merge 
  }
  return(data)
  }
  
##---------------------------Function: find position of value in dataframe that most closely matched a given value-------------------
findloc <- function(dataframe, value) {
  find1 <- abs(dataframe-value)
  find2 <- min(find1)
  findloc <- which(find1==find2)
  findval<- nntreeplots[findloc]
  return(findloc)
}

findvalue.singledf <- function(dataframe, collumn, matchvalue, returnvalue) {
  find1 <- abs(dataframe[collumn]-matchvalue)
  find2 <- min(find1)
  findloc <- which(find1==find2)
  findval<- dataframe[findloc,returnvalue]
  return(findval)
}

## This function is specific for finding the times inside other functions like ddply, when data needs to be flexible and can't be defined rigidly by the function
findvalue.singledf.duration <- function(collumn, matchvalue, returnvalue) {
  find1 <- abs(datasubset[collumn]-matchvalue)
  find2 <- min(find1)
  findloc <- which(find1==find2)
  findval<- datasubset[findloc,returnvalue]
  return(findval)
}


##---------------------------Function: calculate lower and upper interquartile ranges---------------------------
Liqr<-function(x) { 
  return(round(quantile(x,0.25,na.rm=TRUE),10)) 
}
Uiqr<-function(x) { 
  return(round(quantile(x,0.75,na.rm=TRUE),10)) 
}

##---------------------------Function: read in distance matrix data for stochastic simulations-------------------
## note: function is repeated three times to return 1) DATA; 2) site and 3) subplot 
## to edit this code, make edits on one and copy-paste, only changing return() at end of function

read.DATA <- function(path) {
  setwd(path)
  series <- list.files() #save filenames in working directory as a vector of characters
  DATA <- list() #create blank list to store values from loop
  site <- list() #create blank list to store values from loop
  subplot <- list() #create blank list to store values from loop
  for (i in 1:length(series)) { #use filenames in working directory to create list of dfs:
    DATA[[i]] <- read.csv(series[[i]], row.names=1, check.names=FALSE)
    site[[i]] <- str_extract(series[[i]],"DAVO|DBUR|DCAN|DCLU|DLIS|DRED|DTOW|DSUN")
    subplot[[i]] <- NumextractB(series[[i]]) #NumextractB returns whole numbers
  }
  return(DATA)
}

read.site <- function(path) {
  setwd(path)
  series <- list.files() #save filenames in working directory as a vector of characters
  DATA <- list() #create blank list to store values from loop
  site <- list() #create blank list to store values from loop
  subplot <- list() #create blank list to store values from loop
  for (i in 1:length(series)) { #use filenames in working directory to create list of dfs:
    DATA[[i]] <- read.csv(series[[i]], row.names=1, check.names=FALSE)
    site[[i]] <- str_extract(series[[i]],"DAVO|DBUR|DCAN|DCLU|DLIS|DRED|DTOW|DSUN")
    subplot[[i]] <- NumextractB(series[[i]]) #NumextractB returns whole numbers
  }
  return(site)
}

read.subplot <- function(path) {
  setwd(path)
  series <- list.files() #save filenames in working directory as a vector of characters
  DATA <- list() #create blank list to store values from loop
  site <- list() #create blank list to store values from loop
  subplot <- list() #create blank list to store values from loop
  for (i in 1:length(series)) { #use filenames in working directory to create list of dfs:
    DATA[[i]] <- read.csv(series[[i]], row.names=1, check.names=FALSE)
    site[[i]] <- str_extract(series[[i]],"DAVO|DBUR|DCAN|DCLU|DLIS|DRED|DTOW|DSUN")
    subplot[[i]] <- NumextractB(series[[i]]) #NumextractB returns whole numbers
  }
  return(subplot)
}

read.n <- function(path) {
  setwd(path)
  series <- list.files() #save filenames in working directory as a vector of characters
  n <- list() #create blank list to store values from loop
  for (i in 1:length(series)) { #use filenames in working directory to create list of dfs:
    n[[i]] <- NumextractB(series[[i]]) #NumextractB returns whole numbers
  }
  return(n)
}

##---------------------------Function: create parameter table for running values with either fixed B or fixed R0roost-------------------
## note: function is repeated twice to return parameter table for 1) fixed Beta values (setRunvalue.fixedBeta); 2) fixed R0roost values (setRunvalue.fixedR0roost.initial)
## to edit this code, make edits on one and copy-paste, only changing "return() at end of function"Beta=Beta" or "R0roost=R0roost" at expand.grid
## there is a second step for R0roost values: R0roost has an additional chunk of code for calculating B values, AFTER data is scaled. This is encompased within (setRunvalue.fixedR0roost)
## Consider making functions for fixed N or Ntot only, in additon to combined as is here

setRunvalue.fixedBeta <- function(paramvalues, site, subplot, DATA) {
  ### Create expanded df of parameter values. One for fixed N (A), one for fixed Ntot (B)
  runvaluesA <- expand.grid(paramvalues[-2]) #expand grid to create table with all combinations of parameters, without Ntot
  runvaluesB <- expand.grid(paramvalues[-1]) #expand grid to create table with all combinations of parameters, without N
  ## Fill data frame with complementary info
  ## Fixed N:
  for (i in 1:length(c(unique(runvaluesA$Data)))) { #Add information about site & subplot - repeat function for each value of data
    subset <- as.numeric(row.names(runvaluesA[which(runvaluesA$Data==i),])) #subset row names of data with specified value
    runvaluesA[c(subset),"site"] <- site[[i]] #fill in site and subplot information of corresponding data subset
    runvaluesA[c(subset),"subplot"] <- subplot[[i]]
    runvaluesA[c(subset),"n"] <- nrow(DATA[[i]])
  }
  runvaluesA$Ntot <- runvaluesA$N*runvaluesA$n
  ## Fixed Ntot:
  for (i in 1:length(c(unique(runvaluesB$Data)))) { #Add information about site & subplot - repeat function for each value of data
    subset <- as.numeric(row.names(runvaluesB[which(runvaluesB$Data==i),])) #subset row names of data with specified value
    runvaluesB[c(subset),"site"] <- site[[i]] #fill in site and subplot information of corresponding data subset
    runvaluesB[c(subset),"subplot"] <- subplot[[i]]
    runvaluesB[c(subset),"n"] <- nrow(DATA[[i]])
  }
  runvaluesB$N <- runvaluesB$Ntot/runvaluesB$n
  ## Merge fixed N and Ntot
  runvalues <- rbind(runvaluesA, runvaluesB)
  ## Add info cont.
  runvalues$id <- 1:nrow(runvalues)
  runvalues$plotid <- paste(runvalues$site,runvalues$subplot, sep="")
  return(runvalues)
} 

## The first function (setRunvalue.fixedR0roost.initial) is to create parameter table without Beta values to use to re-scale data. The second function (setRunvalue.fixedR0roost) is then used to give the final parameter table with calculated Beta values
setRunvalue.fixedR0roost.initial <- function(paramvalues, site, subplot, DATA) {
  ### Create expanded df of parameter values. One for fixed N (A), one for fixed Ntot (B)
  runvaluesA <- expand.grid(paramvalues[-2]) #expand grid to create table with all combinations of parameters, without Ntot
  runvaluesB <- expand.grid(paramvalues[-1]) #expand grid to create table with all combinations of parameters, without N
  ## Fill data frame with complementary info
  ## Fixed N:
  for (i in 1:length(c(unique(runvaluesA$Data)))) { #Add information about site & subplot - repeat function for each value of data
    subset <- as.numeric(row.names(runvaluesA[which(runvaluesA$Data==i),])) #subset row names of data with specified value
    runvaluesA[c(subset),"site"] <- site[[i]] #fill in site and subplot information of corresponding data subset
    runvaluesA[c(subset),"subplot"] <- subplot[[i]]
    runvaluesA[c(subset),"n"] <- nrow(DATA[[i]])
  }
  runvaluesA$Ntot <- runvaluesA$N*runvaluesA$n
  ## Fixed Ntot:
  for (i in 1:length(c(unique(runvaluesB$Data)))) { #Add information about site & subplot - repeat function for each value of data
    subset <- as.numeric(row.names(runvaluesB[which(runvaluesB$Data==i),])) #subset row names of data with specified value
    runvaluesB[c(subset),"site"] <- site[[i]] #fill in site and subplot information of corresponding data subset
    runvaluesB[c(subset),"subplot"] <- subplot[[i]]
    runvaluesB[c(subset),"n"] <- nrow(DATA[[i]])
  }
  runvaluesB$N <- runvaluesB$Ntot/runvaluesB$n
  ## Merge fixed N and Ntot
  runvalues <- rbind(runvaluesA, runvaluesB)
  ## Add info cont.
  runvalues$id <- 1:nrow(runvalues)
  runvalues$plotid <- paste(runvalues$site,runvalues$subplot, sep="")
  return(runvalues)
}

## Second function finalised the parameter table with Beta values
## Note, the beta function within may need to be changed depending on model structure being used... 
setRunvalue.fixedR0roost <- function(runvalues.initial, SCALED.DATA.ALT, fun) {
  ### calculate beta from R0:
  runvalues.initial$Beta <- NA #need to create collumn in advance of loop, otherwise won't fill Broost into correct positions (will repeat along multiple of i)
  for (i in 1:nrow(runvalues.initial)) { #note that this will need to be run separately for different SEIR+ model structure types. For loop may need to select for the relevant model structure
    ## read values from relevant row for input into equation
    scaled.data.alt <- SCALED.DATA.ALT[[runvalues.initial$id[i]]] #inside [[]] is the value of the parameter/variable given at the row position in runvalues. This is then used to subset the correct df from list SCALED.DATA.ALT. Note that id should be in the correct order here - even though the data was split to calculate between spatial and non-spatial models, the gaps were maintained and filled.I've checked that the scaled.data.alt order is maintained throughout, and that calculated B values match the param values 
    n <- nrow(scaled.data.alt)
    Ntree <- rep(runvalues.initial$N[i],runvalues.initial$n[i]) 
    R0roost <- runvalues.initial$R0roost[i]
    gamma <- runvalues.initial$gamma[i]
    ## calculate beta values from R0
    Broost <- fun(n, Ntree, scaled.data.alt, gamma, R0roost) #fun is function for B values from R0tree, variable depending on model structure. This is for meta-population model without birth and death
    ## fill calculated B values into parameter data frame
    runvalues.initial$Beta[i] <- Broost
  }
  return(runvalues.initial)
}

##---------------------------Function(s): Set beta matrices for data---------------------------
## Gravity model
setBmatrix_gravity <- function(runvalues, SCALED.DATA, theta) {
BETA <- list()
for (i in 1:nrow(runvalues)) { #same for spatial and non-spatial models if distance matrices are set
  #fix values inside loop, per run of loop
  scaled.data <- SCALED.DATA[[runvalues$id[i]]] #need to use id, not just Data, because they're different between spat and non-spat
  Beta <- runvalues$Beta[i]
  theta <- runvalues$theta[i]
  #set Beta matrix
  Btemp <- Beta/((theta*scaled.data)^2) #rate of transmission between groups, scaled by tree radius (gravity model)
  BETA[[i]] <- ifelse(Btemp=="Inf",Beta,Btemp) #rate of transmission within groups
}
return(BETA)
}

## "Plus-one" model
setBmatrix_plusone <- function(runvalues, DATA.long, theta) {
  BETA <- list()
  for (i in 1:nrow(runvalues)) { #same for spatial and non-spatial models if distance matrices are set
    #fix values inside loop, per run of loop
    data.numbered <- DATA.long[[runvalues$id[i]]] #need to use id, not just Data, because they're different between spat and non-spat
    Beta <- runvalues$Beta[i]
    theta <- runvalues$theta[i]
    #set Beta matrix
    BETA[[i]] <- Beta/((data.numbered/theta)+1) #rate of transmission - don't need extra step for within-groups, because Beta/1 is Beta
  }
  return(BETA)
}

##---------------------------Function: re-scale data by radius-------------------
## note - function is repeated four times to return 1) DATA; 2) SCALED.DATA; 3) SCALED.DATA.ALT; 4)SCALED.DATA.DISTANCE 
## DATA.long: Raw distance matrices with only tree numbers re-named to start from 1
## DATA.CORRECTED: #data with between-tree distance threshold applied (spatia models) or with all values changed to zero (non-spatial). This should be main input into setBmatrix of plus-one model
## SCALED.DATA: matrices to be used for setBmatrix of gravity model. Zero values will be replaced with Beta value (within-in group Beta). Have used zero values simply because they are easier to multiply by Beta and then consistently replace. All diagonal values, and trees less than specified radius, will have within-group value for Beta
## SCALED.DATA.ALT: matrices to be used for B equations of gravity models. Replace zero values with 1, so that within-group transmission is included in B calculation (i.e multiplied by 1 instead of 0). 
## SCALED.DATA.DISTANCE: matrices to be used for evaluating spread of virus - do not apply between-tree distance threshold
## to edit this code, make edits on one and copy-paste, only changing return() at end of function
## there is lots that could be trimmed for each (e.g. for DATA, all changes after re-naming), but have kept all the same to make it easier to keep edits consistent

## --- ##
rescaleData.DATA <- function(DATA, runvalues, modstrname, modstrnsname) {
  smods<- runvalues$id[which(runvalues$modstr==modstrname)] #id values of spatial models. 
  nsmods<- runvalues$id[which(runvalues$modstr==modstrnsname)] #id values of non-spatial models. 
  DATA.long <- list() #data with tree numbers re-labelled from 1
  DATA.CORRECTED <- list() #data matrix edited so that trees with distance < radius have zero distance (spatial) or all values are zero (non-spatial). This should be main input into setBmatrix of plus-one model
  SCALED.DATA <- list() #scaled distance matrix with threshold correction applied. For input into setBmatrix of gravity model
  SCALED.DATA.ALT <- list() #scaled distance matrix, with threshold correction applied, but with 1 values for beta and R0 equations
  SCALED.DATA.DISTANCE <- list() #for input into distance spread plot
  
  ## Spatial models:
  for (i in c(smods)) { #run loop on spatial models
    ## fix values inside loop, per run of loop
    data <- DATA[[runvalues$Data[i]]]
    radius <- runvalues$radius[i]
    ## Re-label tree numbers from 1:
    n <- nrow(data)
    colnames(data) <- 1:n
    rownames(data) <- 1:n
    data <- as.matrix(data,ncol=sqrt(length(data)))
    DATA.long[[i]] <- data #over-ride data with re-numbered rows and collumns
    ## Manually change distance values of close trees that should be considered as 1 tree:
    data.corrected <- ifelse(data<radius,0,data)
    DATA.CORRECTED[[i]] <- data.corrected
    ## Scale data matrix by radius (on uncorrected data, for absoulute distance moved):
    scaled.data.distance <- data/radius 
    SCALED.DATA.DISTANCE[[i]] <- scaled.data.distance 
    ## Scale data matrix by radius (on corrected data):
    scaled.data <- data.corrected/radius
    SCALED.DATA[[i]] <- scaled.data #scaled distance with threshold
    ## For gravity model, replace zero values with 1 values:
    scaled.data.alt <- ifelse(scaled.data==0,1,scaled.data) 
    SCALED.DATA.ALT[[i]] <- scaled.data.alt #scaled distance with threshold, but with 1 values for beta and R0 equations
  }
  
  ## Non-spatial models:
  for (i in c(nsmods)) { #run loop on non-spatial models
    ## fix values inside loop, per run of loop
    data <- DATA[[runvalues$Data[i]]]
    radius <- runvalues$radius[i]
    ## re-label tree numbers from 1
    n <- nrow(data)
    colnames(data) <- 1:n
    rownames(data) <- 1:n
    data <- as.matrix(data,ncol=sqrt(length(data)))
    DATA.long[[i]] <- data #over-ride data with re-numbered rows and collumns
    ## Manually change all distance values to zero, because distance doesn't matter in non-spatial models
    data.corrected <- data*0
    DATA.CORRECTED[[i]] <- data.corrected
    ## Scale data matrix by radius (on uncorrected data, for absoulute distance moved):
    scaled.data.distance <- data/radius 
    SCALED.DATA.DISTANCE[[i]] <- scaled.data.distance 
    ## Scale data matrix by radius (on corrected data):
    scaled.data <- data*0 #data*0 is just to preserve data matrix structure
    SCALED.DATA[[i]] <- scaled.data #all values=0 for non-spatial model so that distance not a factor in B matrix (zero values are replaced with the within-group Beta value)
    ## For gravity model, replace zero values with 1 values:
    scaled.data.alt <- (data*0) + 1
    SCALED.DATA.ALT[[i]] <- scaled.data.alt #all values=1 for non-spatial model so that distance not a factor in B and R0 equations
  }
  return(DATA.long)
}

## --- ##
rescaleData.DATA.COR <- function(DATA, runvalues, modstrname, modstrnsname) {
  smods<- runvalues$id[which(runvalues$modstr==modstrname)] #id values of spatial models. 
  nsmods<- runvalues$id[which(runvalues$modstr==modstrnsname)] #id values of non-spatial models. 
  DATA.long <- list() #data with tree numbers re-labelled from 1
  DATA.CORRECTED <- list() #data matrix edited so that trees with distance < radius have zero distance (spatial) or all values are zero (non-spatial). This should be main input into setBmatrix of plus-one model
  SCALED.DATA <- list() #scaled distance matrix with threshold correction applied. For input into setBmatrix of gravity model
  SCALED.DATA.ALT <- list() #scaled distance matrix, with threshold correction applied, but with 1 values for beta and R0 equations
  SCALED.DATA.DISTANCE <- list() #for input into distance spread plot
  
  ## Spatial models:
  for (i in c(smods)) { #run loop on spatial models
    ## fix values inside loop, per run of loop
    data <- DATA[[runvalues$Data[i]]]
    radius <- runvalues$radius[i]
    ## Re-label tree numbers from 1:
    n <- nrow(data)
    colnames(data) <- 1:n
    rownames(data) <- 1:n
    data <- as.matrix(data,ncol=sqrt(length(data)))
    DATA.long[[i]] <- data #over-ride data with re-numbered rows and collumns
    ## Manually change distance values of close trees that should be considered as 1 tree:
    data.corrected <- ifelse(data<radius,0,data)
    DATA.CORRECTED[[i]] <- data.corrected
    ## Scale data matrix by radius (on uncorrected data, for absoulute distance moved):
    scaled.data.distance <- data/radius 
    SCALED.DATA.DISTANCE[[i]] <- scaled.data.distance 
    ## Scale data matrix by radius (on corrected data):
    scaled.data <- data.corrected/radius
    SCALED.DATA[[i]] <- scaled.data #scaled distance with threshold
    ## For gravity model, replace zero values with 1 values:
    scaled.data.alt <- ifelse(scaled.data==0,1,scaled.data) 
    SCALED.DATA.ALT[[i]] <- scaled.data.alt #scaled distance with threshold, but with 1 values for beta and R0 equations
  }
  
  ## Non-spatial models:
  for (i in c(nsmods)) { #run loop on non-spatial models
    ## fix values inside loop, per run of loop
    data <- DATA[[runvalues$Data[i]]]
    radius <- runvalues$radius[i]
    ## re-label tree numbers from 1
    n <- nrow(data)
    colnames(data) <- 1:n
    rownames(data) <- 1:n
    data <- as.matrix(data,ncol=sqrt(length(data)))
    DATA.long[[i]] <- data #over-ride data with re-numbered rows and collumns
    ## Manually change all distance values to zero, because distance doesn't matter in non-spatial models
    data.corrected <- data*0
    DATA.CORRECTED[[i]] <- data.corrected
    ## Scale data matrix by radius (on uncorrected data, for absoulute distance moved):
    scaled.data.distance <- data/radius 
    SCALED.DATA.DISTANCE[[i]] <- scaled.data.distance 
    ## Scale data matrix by radius (on corrected data):
    scaled.data <- data*0 #data*0 is just to preserve data matrix structure
    SCALED.DATA[[i]] <- scaled.data #all values=0 for non-spatial model so that distance not a factor in B matrix (zero values are replaced with the within-group Beta value)
    ## For gravity model, replace zero values with 1 values:
    scaled.data.alt <- (data*0) + 1
    SCALED.DATA.ALT[[i]] <- scaled.data.alt #all values=1 for non-spatial model so that distance not a factor in B and R0 equations
  }
  return(DATA.CORRECTED)
}

## --- ##
rescaleData.SCALED.DATA <- function(DATA, runvalues, modstrname, modstrnsname) {
  smods<- runvalues$id[which(runvalues$modstr==modstrname)] #id values of spatial models. 
  nsmods<- runvalues$id[which(runvalues$modstr==modstrnsname)] #id values of non-spatial models. 
  DATA.long <- list() #data with tree numbers re-labelled from 1
  DATA.CORRECTED <- list() #data matrix edited so that trees with distance < radius have zero distance (spatial) or all values are zero (non-spatial). This should be main input into setBmatrix of plus-one model
  SCALED.DATA <- list() #scaled distance matrix with threshold correction applied. For input into setBmatrix of gravity model
  SCALED.DATA.ALT <- list() #scaled distance matrix, with threshold correction applied, but with 1 values for beta and R0 equations
  SCALED.DATA.DISTANCE <- list() #for input into distance spread plot
  
  ## Spatial models:
  for (i in c(smods)) { #run loop on spatial models
    ## fix values inside loop, per run of loop
    data <- DATA[[runvalues$Data[i]]]
    radius <- runvalues$radius[i]
    ## Re-label tree numbers from 1:
    n <- nrow(data)
    colnames(data) <- 1:n
    rownames(data) <- 1:n
    data <- as.matrix(data,ncol=sqrt(length(data)))
    DATA.long[[i]] <- data #over-ride data with re-numbered rows and collumns
    ## Manually change distance values of close trees that should be considered as 1 tree:
    data.corrected <- ifelse(data<radius,0,data)
    DATA.CORRECTED[[i]] <- data.corrected
    ## Scale data matrix by radius (on uncorrected data, for absoulute distance moved):
    scaled.data.distance <- data/radius 
    SCALED.DATA.DISTANCE[[i]] <- scaled.data.distance 
    ## Scale data matrix by radius (on corrected data):
    scaled.data <- data.corrected/radius
    SCALED.DATA[[i]] <- scaled.data #scaled distance with threshold
    ## For gravity model, replace zero values with 1 values:
    scaled.data.alt <- ifelse(scaled.data==0,1,scaled.data) 
    SCALED.DATA.ALT[[i]] <- scaled.data.alt #scaled distance with threshold, but with 1 values for beta and R0 equations
  }
  
  ## Non-spatial models:
  for (i in c(nsmods)) { #run loop on non-spatial models
    ## fix values inside loop, per run of loop
    data <- DATA[[runvalues$Data[i]]]
    radius <- runvalues$radius[i]
    ## re-label tree numbers from 1
    n <- nrow(data)
    colnames(data) <- 1:n
    rownames(data) <- 1:n
    data <- as.matrix(data,ncol=sqrt(length(data)))
    DATA.long[[i]] <- data #over-ride data with re-numbered rows and collumns
    ## Manually change all distance values to zero, because distance doesn't matter in non-spatial models
    data.corrected <- data*0
    DATA.CORRECTED[[i]] <- data.corrected
    ## Scale data matrix by radius (on uncorrected data, for absoulute distance moved):
    scaled.data.distance <- data/radius 
    SCALED.DATA.DISTANCE[[i]] <- scaled.data.distance 
    ## Scale data matrix by radius (on corrected data):
    scaled.data <- data*0 #data*0 is just to preserve data matrix structure
    SCALED.DATA[[i]] <- scaled.data #all values=0 for non-spatial model so that distance not a factor in B matrix (zero values are replaced with the within-group Beta value)
    ## For gravity model, replace zero values with 1 values:
    scaled.data.alt <- (data*0) + 1
    SCALED.DATA.ALT[[i]] <- scaled.data.alt #all values=1 for non-spatial model so that distance not a factor in B and R0 equations
  }
  return(SCALED.DATA)
}

## --- ##
rescaleData.SCALED.DATA.ALT <- function(DATA, runvalues, modstrname, modstrnsname) {
  smods<- runvalues$id[which(runvalues$modstr==modstrname)] #id values of spatial models. 
  nsmods<- runvalues$id[which(runvalues$modstr==modstrnsname)] #id values of non-spatial models. 
  DATA.long <- list() #data with tree numbers re-labelled from 1
  DATA.CORRECTED <- list() #data matrix edited so that trees with distance < radius have zero distance (spatial) or all values are zero (non-spatial). This should be main input into setBmatrix of plus-one model
  SCALED.DATA <- list() #scaled distance matrix with threshold correction applied. For input into setBmatrix of gravity model
  SCALED.DATA.ALT <- list() #scaled distance matrix, with threshold correction applied, but with 1 values for beta and R0 equations
  SCALED.DATA.DISTANCE <- list() #for input into distance spread plot
  
  ## Spatial models:
  for (i in c(smods)) { #run loop on spatial models
    ## fix values inside loop, per run of loop
    data <- DATA[[runvalues$Data[i]]]
    radius <- runvalues$radius[i]
    ## Re-label tree numbers from 1:
    n <- nrow(data)
    colnames(data) <- 1:n
    rownames(data) <- 1:n
    data <- as.matrix(data,ncol=sqrt(length(data)))
    DATA.long[[i]] <- data #over-ride data with re-numbered rows and collumns
    ## Manually change distance values of close trees that should be considered as 1 tree:
    data.corrected <- ifelse(data<radius,0,data)
    DATA.CORRECTED[[i]] <- data.corrected
    ## Scale data matrix by radius (on uncorrected data, for absoulute distance moved):
    scaled.data.distance <- data/radius 
    SCALED.DATA.DISTANCE[[i]] <- scaled.data.distance 
    ## Scale data matrix by radius (on corrected data):
    scaled.data <- data.corrected/radius
    SCALED.DATA[[i]] <- scaled.data #scaled distance with threshold
    ## For gravity model, replace zero values with 1 values:
    scaled.data.alt <- ifelse(scaled.data==0,1,scaled.data) 
    SCALED.DATA.ALT[[i]] <- scaled.data.alt #scaled distance with threshold, but with 1 values for beta and R0 equations
  }
  
  ## Non-spatial models:
  for (i in c(nsmods)) { #run loop on non-spatial models
    ## fix values inside loop, per run of loop
    data <- DATA[[runvalues$Data[i]]]
    radius <- runvalues$radius[i]
    ## re-label tree numbers from 1
    n <- nrow(data)
    colnames(data) <- 1:n
    rownames(data) <- 1:n
    data <- as.matrix(data,ncol=sqrt(length(data)))
    DATA.long[[i]] <- data #over-ride data with re-numbered rows and collumns
    ## Manually change all distance values to zero, because distance doesn't matter in non-spatial models
    data.corrected <- data*0
    DATA.CORRECTED[[i]] <- data.corrected
    ## Scale data matrix by radius (on uncorrected data, for absoulute distance moved):
    scaled.data.distance <- data/radius 
    SCALED.DATA.DISTANCE[[i]] <- scaled.data.distance 
    ## Scale data matrix by radius (on corrected data):
    scaled.data <- data*0 #data*0 is just to preserve data matrix structure
    SCALED.DATA[[i]] <- scaled.data #all values=0 for non-spatial model so that distance not a factor in B matrix (zero values are replaced with the within-group Beta value)
    ## For gravity model, replace zero values with 1 values:
    scaled.data.alt <- (data*0) + 1
    SCALED.DATA.ALT[[i]] <- scaled.data.alt #all values=1 for non-spatial model so that distance not a factor in B and R0 equations
  }
  return(SCALED.DATA.ALT)
}

## --- ##
rescaleData.SCALED.DATA.DISTANCE <- function(DATA, runvalues, modstrname, modstrnsname) {
  smods<- runvalues$id[which(runvalues$modstr==modstrname)] #id values of spatial models. 
  nsmods<- runvalues$id[which(runvalues$modstr==modstrnsname)] #id values of non-spatial models. 
  DATA.long <- list() #data with tree numbers re-labelled from 1
  DATA.CORRECTED <- list() #data matrix edited so that trees with distance < radius have zero distance (spatial) or all values are zero (non-spatial). This should be main input into setBmatrix of plus-one model
  SCALED.DATA <- list() #scaled distance matrix with threshold correction applied. For input into setBmatrix of gravity model
  SCALED.DATA.ALT <- list() #scaled distance matrix, with threshold correction applied, but with 1 values for beta and R0 equations
  SCALED.DATA.DISTANCE <- list() #for input into distance spread plot
  
  ## Spatial models:
  for (i in c(smods)) { #run loop on spatial models
    ## fix values inside loop, per run of loop
    data <- DATA[[runvalues$Data[i]]]
    radius <- runvalues$radius[i]
    ## Re-label tree numbers from 1:
    n <- nrow(data)
    colnames(data) <- 1:n
    rownames(data) <- 1:n
    data <- as.matrix(data,ncol=sqrt(length(data)))
    DATA.long[[i]] <- data #over-ride data with re-numbered rows and collumns
    ## Manually change distance values of close trees that should be considered as 1 tree:
    data.corrected <- ifelse(data<radius,0,data)
    DATA.CORRECTED[[i]] <- data.corrected
    ## Scale data matrix by radius (on uncorrected data, for absoulute distance moved):
    scaled.data.distance <- data/radius 
    SCALED.DATA.DISTANCE[[i]] <- scaled.data.distance 
    ## Scale data matrix by radius (on corrected data):
    scaled.data <- data.corrected/radius
    SCALED.DATA[[i]] <- scaled.data #scaled distance with threshold
    ## For gravity model, replace zero values with 1 values:
    scaled.data.alt <- ifelse(scaled.data==0,1,scaled.data) 
    SCALED.DATA.ALT[[i]] <- scaled.data.alt #scaled distance with threshold, but with 1 values for beta and R0 equations
  }
  
  ## Non-spatial models:
  for (i in c(nsmods)) { #run loop on non-spatial models
    ## fix values inside loop, per run of loop
    data <- DATA[[runvalues$Data[i]]]
    radius <- runvalues$radius[i]
    ## re-label tree numbers from 1
    n <- nrow(data)
    colnames(data) <- 1:n
    rownames(data) <- 1:n
    data <- as.matrix(data,ncol=sqrt(length(data)))
    DATA.long[[i]] <- data #over-ride data with re-numbered rows and collumns
    ## Manually change all distance values to zero, because distance doesn't matter in non-spatial models
    data.corrected <- data*0
    DATA.CORRECTED[[i]] <- data.corrected
    ## Scale data matrix by radius (on uncorrected data, for absoulute distance moved):
    scaled.data.distance <- data/radius 
    SCALED.DATA.DISTANCE[[i]] <- scaled.data.distance 
    ## Scale data matrix by radius (on corrected data):
    scaled.data <- data*0 #data*0 is just to preserve data matrix structure
    SCALED.DATA[[i]] <- scaled.data #all values=0 for non-spatial model so that distance not a factor in B matrix (zero values are replaced with the within-group Beta value)
    ## For gravity model, replace zero values with 1 values:
    scaled.data.alt <- (data*0) + 1
    SCALED.DATA.ALT[[i]] <- scaled.data.alt #all values=1 for non-spatial model so that distance not a factor in B and R0 equations
  }
  return(SCALED.DATA.DISTANCE)
}

##---------------------------Function: Create table of parameter values: Extract parameter values from file names, add information on tree structure (merge from Ctree outputs), and Adjust numeric/factor/character categories if needed, and set order
extractvalues <- function(modstrs, colnames, structure) {
  ## Extract parameter values from file names
  series <- list.files()
  series.params1 <- as.data.frame(t(sapply(1:length(series),function(i) {Numextract(series[[i]])})))
  series.params2 <- as.data.frame(t(t(sapply(1:length(series),function(i) {str_extract(series[[i]],"DAVO|DBUR|DCAN|DCLU|DLIS|DRED|DTOW|DSUN")}))))
  series.params3 <- as.data.frame(t(t(sapply(1:length(series),function(i) {str_extract(series[[i]],modstrs)})))) #be cautious here - will match to first partial match
  series.params4 <- as.data.frame(t(t(sapply(1:length(series),function(i) {str_extract(series[[i]],"PO|G")}))))
  series.params <- cbind(series.params3, series.params4, series.params2, series.params1)
  colnames(series.params) <- c(colnames)
  series.params$position <- 1:nrow(series.params) #add position number now - risk that additional manipulation will change order (e.g. merge with additional datasets)
  series.params$plotID <- paste(series.params$site, series.params$subplot, sep="")
  ## Add information on tree structure (merge from Ctree outputs)
  series.params <- merge(structure[,c("plotID","structure", "meandist")],series.params, by="plotID") #re-orders series.params here
  ## Adjust numeric/factor/character categories if needed, and set order
  #summary(series.params)
  series.params$id <- as.numeric(as.character(series.params$id)) #as.character wrap to stop R from changing values
  series.params$plotID <- as.factor(series.params$plotID)
  series.params$Bfun <- as.factor(series.params$Bfun)
  series.params$n <- as.numeric(as.character(series.params$n)) #as.character wrap to stop R from changing values
  series.params$Beta <- as.numeric(as.character(series.params$Beta)) #as.character wrap to stop R from changing values
  series.params$gamma <- as.numeric(as.character(series.params$gamma)) #as.character wrap to stop R from changing values
  series.params$theta <- as.numeric(as.character(series.params$theta))
  series.params$N <- as.numeric(as.character(series.params$N)) #as.character wrap to stop R from changing values
  series.params$Ntot <- as.numeric(as.character(series.params$Ntot)) #as.character wrap to stop R from changing values
  series.params$radius <- as.numeric(as.character(series.params$radius)) #as.character wrap to stop R from changing values
  series.params$threshold <- as.numeric(as.character(series.params$threshold)) #as.character wrap to stop R from changing values
  series.params$nsims <- as.numeric(as.character(series.params$nsims)) #as.character wrap to stop R from changing values
  series.params$R0roost <- as.numeric(as.character(series.params$R0roost)) #as.character wrap to stop R from changing values
  series.params$structure <- as.factor(series.params$structure)
  
  series.params <- series.params[order(series.params$position),] #re-order again to match with series, for read-ing in data
  return(series.params)
}

extractvalues.rand <- function(modstrs, colnames) {
  ## Extract parameter values from file names
  series <- list.files()
  series.params1 <- as.data.frame(t(sapply(1:length(series),function(i) {Numextract(series[[i]])})))
  series.params2 <- as.data.frame(t(t(sapply(1:length(series),function(i) {str_extract(series[[i]],"DAVO|DBUR|DCAN|DCLU|DLIS|DRED|DTOW|DSUN")}))))
  series.params3 <- as.data.frame(t(t(sapply(1:length(series),function(i) {str_extract(series[[i]],modstrs)})))) #be cautious here - will match to first partial match
  series.params4 <- as.data.frame(t(t(sapply(1:length(series),function(i) {str_extract(series[[i]],"PO|G")}))))
  series.params <- cbind(series.params3, series.params4, series.params2, series.params1)
  colnames(series.params) <- c(colnames)
  series.params$position <- 1:nrow(series.params) #add position number now - risk that additional manipulation will change order (e.g. merge with additional datasets)
  series.params$plotID <- paste(series.params$site, series.params$subplot, sep="")
  ## Adjust numeric/factor/character categories if needed, and set order
  #summary(series.params)
  series.params$id <- as.numeric(as.character(series.params$id)) #as.character wrap to stop R from changing values
  series.params$plotID <- as.factor(series.params$plotID)
  series.params$Bfun <- as.factor(series.params$Bfun)
  series.params$n <- as.numeric(as.character(series.params$n)) #as.character wrap to stop R from changing values
  series.params$Beta <- as.numeric(as.character(series.params$Beta)) #as.character wrap to stop R from changing values
  series.params$gamma <- as.numeric(as.character(series.params$gamma)) #as.character wrap to stop R from changing values
  series.params$theta <- as.numeric(as.character(series.params$theta))
  series.params$N <- as.numeric(as.character(series.params$N)) #as.character wrap to stop R from changing values
  series.params$Ntot <- as.numeric(as.character(series.params$Ntot)) #as.character wrap to stop R from changing values
  series.params$radius <- as.numeric(as.character(series.params$radius)) #as.character wrap to stop R from changing values
  series.params$threshold <- as.numeric(as.character(series.params$threshold)) #as.character wrap to stop R from changing values
  series.params$nsims <- as.numeric(as.character(series.params$nsims)) #as.character wrap to stop R from changing values
  series.params$R0roost <- as.numeric(as.character(series.params$R0roost)) #as.character wrap to stop R from changing values
  
  series.params <- series.params[order(series.params$position),] #re-order again to match with series, for read-ing in data
  return(series.params)
}

##--------------------Function: calculate and save summaries per dataset and variable combination-------------------
## Average extinction times for grouping
ext.time.summary <- function(dataframe, grouplist, filename) {
  for (i in 1:length(dataframe)) {
    dataframeI <- dataframe[[i]]
    save <- ddply(dataframeI, grouplist, summarise,
                  count.n =sum(!is.na(.n)), #gives the number of successful outbreaks per grouping. Also gives the number of rows, because single row per .n
                  mininf = min(max.infections), #check that only simulations with >threshold are included in average
                  mean.end.time = mean(end.time), #average extinction time in specified grouping
                  median.end.time = median(end.time),
                  LIQR.end.time = Liqr(end.time), 
                  UIQR.end.time = Uiqr(end.time)) 
    write.csv(save, paste(filename, "_list ", i, ".csv", sep = "")) #save summary for later reference
    return(save)
  }
}

## Average proportion of successful outbreaks
success.summary <- function(dataframe, grouplist, filename) {
  save <- ddply(dataframe, grouplist, summarise,
                sum.successful = sum(count.n), #gives the number of successful outbreaks per grouping. Also gives the number of rows, because single row per .n
                denominator = sum(is.na(count.n)),
                mean.prop.successful = mean(prop.successful), 
                LIQR.prop.successful = Liqr(prop.successful),
                UIQR.prop.successful = Uiqr(prop.successful)) #average extinction time in specified grouping
  write.csv(save, paste(filename, "_list ", i, ".csv", sep = "")) #save summary for later reference
  return(save)
}
  
## Average magnitude and time of epidemic peak for grouping
epi.MagTime.summary <- function(dataframe, grouplist, filename) {
  for (i in 1:length(dataframe)) {
    dataframeI <- dataframe[[i]]
    save <- ddply(dataframeI, grouplist, summarise,
                  mean.epi.Iprop = mean(I.prop),
                  LIQR.epi.Iprop = Liqr(I.prop),
                  UIQR.epi.Iprop = Uiqr(I.prop),
                  mean.epi.infections = mean(infections),
                  mean.epi.time = mean(time),
                  LIQR.epi.time = Liqr(time),
                  UIQR.epi.time = Uiqr(time),
                  denominator =sum(!is.na(.n))) 
    write.csv(save, paste(filename, "_list ", i, ".csv", sep = "")) #save summary for later reference
    return(save)
  }
}

## Average duration of epidemic peak for grouping
epi.duration.summary <- function(dataframe, grouplist, filename) {
  for (i in 1:length(dataframe)) {
    dataframeI <- dataframe[[i]]
    save <- ddply(dataframeI, grouplist, summarise,
                  mean.epi.duration = mean(duration),
                  LIQR.epi.duration = Liqr(duration),
                  UIQR.epi.duration = Uiqr(duration),
                  denominator =sum(!is.na(.n))) 
    write.csv(save, paste(filename, "_list ", i, ".csv", sep = "")) #save summary for later reference
    return(save)
  }
}
##----------------------------------------------------------------------------------------------------
##------------------------------------------------Plotting--------------------------------------------
##----------------------------------------------------------------------------------------------------

##--------------------Function: filter out non-successful outbreaks before plotting-------------------
success.filter <- function(dataset, simend.output, threshold) {
success.data <- list()
for (i in 1:length(dataset)) {
  success.list <- successful(simend.output[[i]], threshold)
  success.data[[i]] <- dataset[[i]][dataset[[i]]$.n %in% success.list$.n,] #run for successful outbreaks only
}
return(success.data)
}

##---------------------------Function: Plot proportion and timing of extinction by group (e.g. structure or site) and facet by R0 or Beta ---------
## Note - possible trap here because this code takes the END time of simulations, which could be when the virus becomes extinct OR when the simulation was forced to end. Make sure all simulations are allowed to run to the max time specified in this code
## Spatial plots:
plot.ext.spat <- function(data, spatmodstr, fixedparam, subsetparam, group, filename, collabels) {
plots.spat <- list()
plots.spat.group <- list()
plots.spat.single <- list()
data.sub <- list()
dataspat <- list()
for (i in 1:length(data)) {
  ## prep data:
  data.sub[[i]] <- data[[i]] #run for all data (successful and non-successful)
  #data.sub <- data[[i]][data[[i]]$.n %in% success.list$.n,] #run for successful outbreaks only
  dataspat[[i]] <- data.sub[[i]][which(data.sub[[i]]$modstr==spatmodstr),] %>% #select spatial data and arrange in plotting order
    mutate(sim.category = factor(sim.category,levels=c(paste("ext",7:1,sep="."),"sim.fail"))) %>%
    arrange(desc(sim.category)) 
  ## plot:
  plots.spat[[i]] <- dataspat[[i]] %>%
    ggplot(aes(fill=sim.category, y=prop.sim, x=dataspat[[i]][[group]])) + #bars by group (e.g. structure or site)
    geom_bar(position="fill", stat="identity")+
    labs(y="% simulations", x="Group", fill="Outcome of simulation") +
    ggtitle("Outbreak success - heterogeneous mixing") +
    theme_bw() +
    background_grid("none") +
    scale_fill_brewer(palette="RdYlBu",labels=c(collabels)) +
    facet_grid(as.formula(paste(fixedparam,"~",subsetparam))) 
    ## save each plot with element sizes for either single plot, or grouped plot:
  plots.spat.group[[i]] <- plots.spat[[i]] +
    theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.title = element_text(size=8)) #adjust legend and text size as needed
  plots.spat.single[[i]] <- plots.spat[[i]] +
    theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(1, "cm"), legend.text = element_text(size=10), legend.title = element_text(size=12), axis.title = element_text(size=10)) #adjust legend and text size as needed
  #png(filename = paste(filename, "_spatial", "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=1140) #square with 1 row
  png(filename = paste(filename, "_spatial", "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=3250) #square with 4 rows
  plot(plots.spat.single[[i]])
  dev.off()
  ## save rds - workaround for return(plots)
  #saveRDS(plots.spat.group[[i]], file = paste(filename, "_spatial", "_PLOT ", i,".rds", sep = ""))
}
#return(plots.spat) #This won't work with facet_grid(dataspat[[i]][[fixedparam]]~dataspat[[i]][[subsetparam]]). Get: "Error in `$<-.data.frame`(x, name, value) : replacement has 168 rows, data has 176"
}

## Non-spatial plots
plot.ext.nonspat <- function(data, nonspatmodstr, fixedparam, subsetparam, group, filename, collabels) {
  plots.nonspat <- list()
  plots.nonspat.group <- list()
  plots.nonspat.single <- list()
  data.sub <- list()
  datanonspat <- list()
  for (i in 1:length(data)) {
    ## prep data:
    data.sub[[i]] <- data[[i]] #run for all data (successful and non-successful)
    datanonspat[[i]] <- data.sub[[i]][which(data.sub[[i]]$modstr==nonspatmodstr),] %>% #select spatial data and arrange in plotting order
      mutate(sim.category = factor(sim.category,levels=c(paste("ext",7:1,sep="."),"sim.fail"))) %>%
      arrange(desc(sim.category)) 
    ## plot:
    plots.nonspat[[i]] <- datanonspat[[i]] %>%
      ggplot(aes(fill=sim.category, y=prop.sim, x=datanonspat[[i]][[group]])) + #bars by group (e.g. structure or site)
      geom_bar(position="fill", stat="identity")+
      labs(y="% simulations", x="Group", fill="Outcome of simulation") +
      ggtitle("Outbreak success - homogeneous mixing") +
      theme_bw() +
      background_grid("none") +
      scale_fill_brewer(palette="RdYlBu",labels=c(collabels)) +
      facet_grid(as.formula(paste(fixedparam,"~",subsetparam)))
    ## save each plot with element sizes for either single plot, or grouped plot:
    plots.nonspat.group[[i]] <- plots.nonspat[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.title = element_text(size=8)) #adjust legend and text size as needed
    plots.nonspat.single[[i]] <- plots.nonspat[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(1, "cm"), legend.text = element_text(size=10), legend.title = element_text(size=12), axis.title = element_text(size=10)) #adjust legend and text size as needed
    #png(filename = paste(filename, "_non spatial", "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=1140) #square with 1 row
    png(filename = paste(filename, "_non spatial", "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=3250) #square with 4 rows
    plot(plots.nonspat.single[[i]])
    dev.off()
    ## save rds - workaround for return(plots)
    #saveRDS(plots.nonspat.group[[i]], file = paste(filename, "_non spatial", "_PLOT ", i,".rds", sep = ""))
  }
  #return(plots.nonspat) #This won't work with facet_grid(datanonspat[[i]][[fixedparam]]~datanonspat[[i]][[subsetparam]]). Get: "Error in `$<-.data.frame`(x, name, value) : replacement has 168 rows, data has 176"
}

## Combined, split by group
plot.ext.group <- function(data, fixedparam, subsetparam, group, filename, collabels) {
  plots.bar <- list()
  plots.bar.group <- list()
  plots.bar.single <- list()
  data.sub <- list()
  for (i in 1:length(data)) {
    ## prep data:
    data.sub[[i]] <- data[[i]] %>% #arrange in plotting order
      mutate(sim.category = factor(sim.category,levels=c(paste("ext",7:1,sep="."),"sim.fail"))) %>%
      arrange(desc(sim.category)) 
    ## plot:
    plots.bar[[i]] <- data.sub[[i]] %>%
      ggplot(aes(fill=sim.category, y=prop.sim, x=data.sub[[i]][[group]])) + #bars by group (e.g. structure or site)
      geom_bar(position="fill", stat="identity")+
      labs(y="Proportion of simulations", x="Structure", fill="Outcome of simulation") +
      ggtitle("Outbreak success") +
      theme_bw() +
      background_grid("none") +
      scale_fill_brewer(palette="RdYlBu",labels=c(collabels)) +
      facet_grid(as.formula(paste(fixedparam,"~",subsetparam))) 
    ## save each plot with element sizes for either single plot, or grouped plot:
    plots.bar.group[[i]] <- plots.bar[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.title = element_text(size=8)) #adjust legend and text size as needed
    plots.bar.single[[i]] <- plots.bar[[i]] +
      #theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(1, "cm"), legend.text = element_text(size=10), legend.title = element_text(size=12), axis.title = element_text(size=10)) #adjust legend and text size as needed
      theme(axis.text.x = element_text(angle = 45, hjust=0.9, size=15),legend.key.size = unit(1, "cm"), legend.text = element_text(size=15), legend.title = element_text(size=17),axis.text.y = element_text(size=15), axis.title = element_text(size=17), strip.text.x = element_text(size = 15),strip.text.y = element_text(size = 15)) #adjust legend and text size as needed
    png(filename = paste(filename, "_spatial", "_PLOT ", i, ".png", sep = ""), res = 300, width=3050, height=3700) #full page, square with 5 rows
    #png(filename = paste(filename, "_spatial", "_PLOT ", i, ".png", sep = ""), res = 300, width=2850, height=2250)
    plot(plots.bar.single[[i]])
    dev.off()
    }
}

##---------------------------Function: Plot spread from index case by group (e.g. structure or site) and facet by R0 or Beta ---------
line_ribbon_plot <- function(dataset, spatmodstr, nonspatmodstr, yname, xname, Lvarname, Uvarname, group, fixedparam, subsetparam, filename, ylab, xlab, mainlab, xlim, ylim) {
  spread.plots <- list() #note - stat_smooth will average all values not specified as separate under the colour aes (e.g. if plotting by structure, will combine multiple sites of the same structure)
  spread.plots.single <- list()
  spread.plots.group <- list()
  data <- list()
  dataspat <- list()
  datanspat <- list()
  for (i in 1:length(dataset)) {
    data[[i]] <- dataset[[i]] #run for all data (successful and non-successful)
    dataspat[[i]] <- data[[i]][which(data[[i]]$modstr==spatmodstr),] #select spatial data
    datanspat[[i]] <- data[[i]][which(data[[i]]$modstr==nonspatmodstr),] #select non-spatial data
    spread.plots[[i]] <- ggplot() +
      
      #spatial (solid line, solid labels)
      geom_line(data=dataspat[[i]], aes(y=dataspat[[i]][[yname]], x=dataspat[[i]][[xname]], colour=dataspat[[i]][[group]]), size=0.5, linetype="solid") + #geom_line is not a single, good looking line - needs to have only one combination set plotted to be this. Is likely not ok to loess smooth my already smoothed lines... 
      geom_ribbon(data=dataspat[[i]], aes(ymin=dataspat[[i]][[Lvarname]], ymax=dataspat[[i]][[Uvarname]], x=dataspat[[i]][[xname]], group=dataspat[[i]][[group]]), alpha=0.5, fill = "grey70") +
      
      #non spatial (dashed line, transparant labels)
      geom_line(data=datanspat[[i]], aes(y=datanspat[[i]][[yname]], x=datanspat[[i]][[xname]], colour=datanspat[[i]][[group]]), size=0.5, linetype="dashed") + #geom_line is not a single, good looking line - needs to have only one combination set plotted to be this. Is likely not ok to loess smooth my already smoothed lines... 
      geom_ribbon(data=datanspat[[i]], aes(ymin=datanspat[[i]][[Lvarname]], ymax=datanspat[[i]][[Uvarname]], x=datanspat[[i]][[xname]], group=datanspat[[i]][[group]]), alpha=0.5, fill = "grey70") +
      
      #rest of plot
      coord_cartesian(xlim = c(0, xlim), ylim = c(0, ylim)) + #max possible non-normalised distance in a plot would be 28.3. Max normalised is 1
      labs(y=ylab, x=xlab, main=mainlab, colour=group) +
      ggtitle(mainlab) +
      theme_bw() +
      facet_grid(as.formula(paste(fixedparam,"~",subsetparam))) + #drop=FALSE meant to make ggplot draw the empty facets
      background_grid("none")
    
    ## save each plot with element sizes for either single plot, or grouped plot:
    spread.plots.group[[i]] <- spread.plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.title = element_text(size=8)) #adjust legend and text size as needed
    spread.plots.single[[i]] <- spread.plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9, size=15),legend.key.size = unit(1, "cm"), legend.text = element_text(size=15), legend.title = element_text(size=17),axis.text.y = element_text(size=15), axis.title = element_text(size=17), strip.text.x = element_text(size = 15),strip.text.y = element_text(size = 15)) #adjust legend and text size as needed
      #theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(1, "cm"), legend.text = element_text(size=10), legend.title = element_text(size=12), axis.title = element_text(size=10)) #adjust legend and text size as needed
    #png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=1080) #square with 1 row
    #png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=3250) #square with 4 rows
    png(filename = paste(filename, "_spatial", "_PLOT ", i, ".png", sep = ""), res = 300, width=2850, height=3700)
    plot(spread.plots.single[[i]])
    dev.off()
    ## save rds - workaround for return(plots) - takes a while with large simulations
    #saveRDS(spread.plots.group[[i]], file = paste(filename, "_PLOT ", i, ".rds", sep = ""))
  }
  #return(spread.plots) #This won't work with facet_grid
}

## Plot with inbuild gam
line_ribbon_plot_inbuilt <- function(dataset, spatmodstr, nonspatmodstr, yname, xname, group, fixedparam, subsetparam, filename, ylab, xlab, mainlab, xlim, ylim) {
  spread.plots <- list() #note - stat_smooth will average all values not specified as separate under the colour aes (e.g. if plotting by structure, will combine multiple sites of the same structure)
  spread.plots.single <- list()
  spread.plots.group <- list()
  data <- list()
  dataspat <- list()
  datanspat <- list()
  for (i in 1:length(dataset)) {
    data[[i]] <- dataset[[i]] #run for all data (successful and non-successful)
    dataspat[[i]] <- data[[i]][which(data[[i]]$modstr==spatmodstr),] #select spatial data
    datanspat[[i]] <- data[[i]][which(data[[i]]$modstr==nonspatmodstr),] #select non-spatial data
    spread.plots[[i]] <- ggplot() +
      stat_smooth(data=dataspat[[i]], aes(y=dataspat[[i]][[yname]], x=dataspat[[i]][[xname]], colour=dataspat[[i]][[group]]), method = "gam", formula = y ~s(x), size=0.5, linetype="solid") +
      stat_smooth(data=datanspat[[i]], aes(y=datanspat[[i]][[yname]], x=datanspat[[i]][[xname]], colour=datanspat[[i]][[group]]), method = "gam", formula = y ~s(x), size=0.5, linetype="dashed") +
      #rest of plot
      coord_cartesian(xlim = c(0, xlim), ylim = c(0, ylim)) + #max possible non-normalised distance in a plot would be 28.3. Max normalised is 1
      labs(y=ylab, x=xlab, main=mainlab, colour=group) +
      ggtitle(mainlab) +
      theme_bw() +
      facet_grid(as.formula(paste(fixedparam,"~",subsetparam))) + #drop=FALSE meant to make ggplot draw the empty facets
      background_grid("none")
    ## save each plot with element sizes for either single plot, or grouped plot:
    spread.plots.group[[i]] <- spread.plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.title = element_text(size=8)) #adjust legend and text size as needed
    spread.plots.single[[i]] <- spread.plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(1, "cm"), legend.text = element_text(size=10), legend.title = element_text(size=12), axis.title = element_text(size=10)) #adjust legend and text size as needed
    #png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=2000)
    png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=1080) #square
    plot(spread.plots.single[[i]])
    dev.off()
    ## save rds - workaround for return(plots) - takes a while with large simulations
    #saveRDS(spread.plots.group[[i]], file = paste(filename, "_PLOT ", i, ".rds", sep = ""))
  }
  #return(spread.plots) #This won't work with facet_grid
}

## Line plot for visualisation of individual simulations (run per modstr) - without continuous gams (just sim lines and discrete median + IQR)
line_ribbon_plot_.n <- function(SimData, DiscData, modstrSub, Simyname, Simxname, Discyname, Discxname, DiscLvarname, DiscUvarname, group, fixedparam, subsetparam, filename, ylab, xlab, mainlab, xlim, ylim) {
  SimDataSub <- list()
  DiscDataSub <- list()
  plots <- list() 
  plots.single <- list()
  plots.group <- list()
  
  for (i in 1:length(SimData)) {
    ##subset datasets by spatial structure
    SimDataSub[[i]] <- SimData[[i]][which(SimData[[i]]$modstr==modstrSub),] 
    DiscDataSub[[i]] <- DiscData[[i]][which(DiscData[[i]]$modstr==modstrSub),] 
    
    plots[[i]] <- ggplot() +
      ##individual simulations
      geom_line(data=SimDataSub[[i]], aes(y=SimDataSub[[i]][[Simyname]], x=SimDataSub[[i]][[Simxname]], group=SimDataSub[[i]][[group]]), colour="gray", size=0.1, linetype="solid") + 
      ##median discrete line
      geom_line(data=DiscDataSub[[i]], aes(y=DiscDataSub[[i]][[Discyname]], x=DiscDataSub[[i]][[Discxname]]), colour="mediumpurple4", size=0.5, linetype="solid") + 
      geom_ribbon(data=DiscDataSub[[i]], aes(ymin=DiscDataSub[[i]][[DiscLvarname]], ymax=DiscDataSub[[i]][[DiscUvarname]], x=DiscDataSub[[i]][[Discxname]]), alpha=0.5, fill = "mediumpurple4") +
      ##rest of plot
      coord_cartesian(xlim = c(0, xlim), ylim = c(0, ylim)) + 
      labs(y=ylab, x=xlab, colour=group) +
      ggtitle(paste(mainlab, " - ",modstrSub, sep="")) +
      theme_bw() +
      facet_grid(as.formula(paste(fixedparam,"~",subsetparam))) + 
      background_grid("none")
    ## save each plot with element sizes for either single plot, or grouped plot:
    plots.group[[i]] <- plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.title = element_text(size=8)) #adjust legend and text size as needed
    plots.single[[i]] <- plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(1, "cm"), legend.text = element_text(size=10), legend.title = element_text(size=12), axis.title = element_text(size=10)) #adjust legend and text size as needed
    png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=2000)
    #png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=3180) 
    plot(plots.single[[i]])
    dev.off()
    ## save rds - workaround for return(plots) - takes a while with large simulations
    #saveRDS(plots.group[[i]], file = paste(filename, "_PLOT ", i, ".rds", sep = ""))
  }
  #return(plots) #This won't work with facet_grid
}

## Line plot for visualisation of individual simulations (run per modstr) - with continuous gams
line_ribbon_plot_.n_withcont <- function(SimData, ContData, DiscData, modstrSub, Simyname, Simxname, Contyname, Contxname, ContLvarname, ContUvarname, Discyname, Discxname, DiscLvarname, DiscUvarname, group, fixedparam, subsetparam, filename, ylab, xlab, mainlab, xlim, ylim) {
  SimDataSub <- list()
  ContDataSub <- list()
  DiscDataSub <- list()
  plots <- list() 
  plots.single <- list()
  plots.group <- list()
  
  for (i in 1:length(SimData)) {
    ##subset datasets by spatial structure
    SimDataSub[[i]] <- SimData[[i]][which(SimData[[i]]$modstr==modstrSub),] 
    ContDataSub[[i]] <- ContData[[i]][which(ContData[[i]]$modstr==modstrSub),] 
    DiscDataSub[[i]] <- DiscData[[i]][which(DiscData[[i]]$modstr==modstrSub),] 
    
    plots[[i]] <- ggplot() +
      ##individual simulations
      geom_line(data=SimDataSub[[i]], aes(y=SimDataSub[[i]][[Simyname]], x=SimDataSub[[i]][[Simxname]], group=SimDataSub[[i]][[group]]), colour="gray", size=0.1, linetype="solid") + 
      ##fitted continuous line
      geom_line(data=ContDataSub[[i]], aes(y=ContDataSub[[i]][[Contyname]], x=ContDataSub[[i]][[Contxname]]), colour="black", size=0.5, linetype="solid") + 
      geom_ribbon(data=ContDataSub[[i]], aes(ymin=(ContDataSub[[i]][[Contyname]]-ContDataSub[[i]][[ContLvarname]]), ymax=(ContDataSub[[i]][[Contyname]]+ContDataSub[[i]][[ContUvarname]]), x=ContDataSub[[i]][[Contxname]]), alpha=0.5, fill = "grey70") +
      ##median discrete line
      geom_line(data=DiscDataSub[[i]], aes(y=DiscDataSub[[i]][[Discyname]], x=DiscDataSub[[i]][[Discxname]]), colour="mediumpurple4", size=0.5, linetype="solid") + 
      geom_ribbon(data=DiscDataSub[[i]], aes(ymin=DiscDataSub[[i]][[DiscLvarname]], ymax=DiscDataSub[[i]][[DiscUvarname]], x=DiscDataSub[[i]][[Discxname]]), alpha=0.5, fill = "mediumpurple4") +
      ##rest of plot
      coord_cartesian(xlim = c(0, xlim), ylim = c(0, ylim)) + 
      labs(y=ylab, x=xlab, colour=group) +
      ggtitle(paste(mainlab, " - ",modstrSub, sep="")) +
      theme_bw() +
      facet_grid(as.formula(paste(fixedparam,"~",subsetparam))) + 
      background_grid("none")
    ## save each plot with element sizes for either single plot, or grouped plot:
    plots.group[[i]] <- plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.title = element_text(size=8)) #adjust legend and text size as needed
    plots.single[[i]] <- plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(1, "cm"), legend.text = element_text(size=10), legend.title = element_text(size=12), axis.title = element_text(size=10)) #adjust legend and text size as needed
    #png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=2000)
    png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=3180) #square
    plot(plots.single[[i]])
    dev.off()
    ## save rds - workaround for return(plots) - takes a while with large simulations
    #saveRDS(plots.group[[i]], file = paste(filename, "_PLOT ", i, ".rds", sep = ""))
  }
  #return(plots) #This won't work with facet_grid
}

##---------------------------Function: Plot characteristics of epidemic peak (timing, magnitude, duration) by group (e.g. structure or site) and facet by R0 or Beta ---------
box_plot <- function(dataset, spatmodstr, nonspatmodstr, yname, group, fixedparam, subsetparam, filename, ylab, xlab, mainlab, ylim) {
  epipeak.plots <- list() 
  epipeak.plots.single <- list() 
  epipeak.plots.group <- list() 
  data <- list()
  dataspat <- list()
  datanspat <- list()
  for (i in 1:length(dataset)) {
    data[[i]] <- dataset[[i]]
    dataspat[[i]] <- data[[i]][which(data[[i]]$modstr==spatmodstr),] #select spatial data
    datanspat[[i]] <- data[[i]][which(data[[i]]$modstr==nonspatmodstr),] #select non-spatial data
    epipeak.plots[[i]] <- ggplot() +
      
      #spatial (solid line, solid labels)
      geom_boxplot(data=dataspat[[i]], aes(y=dataspat[[i]][[yname]], x=dataspat[[i]][[group]], colour=dataspat[[i]][[group]])) + 
      #geom_boxplot(data=dataspat[[i]], aes(y=dataspat[[i]][[yname]], x=dataspat[[i]][[group]], colour=dataspat[[i]][[group]],fill=dataspat[[i]][[group]], alpha=0.1)) + #ideally, would like both the line and fill to be transperant, but can't figure this?
      
      #non spatial (dashed line, transparant labels)
      #geom_boxplot(data=datanspat[[i]], aes(y=datanspat[[i]][[yname]], x=datanspat[[i]][[group]], colour=datanspat[[i]][[group]],fill=datanspat[[i]][[group]], alpha=0.01)) + 
      
      #rest of plot
      coord_cartesian(ylim = c(0, ylim)) + 
      labs(y=ylab, x=xlab, main=mainlab, colour=group, fill=NULL, alpha=NULL) +
      ggtitle(mainlab) +
      theme_bw() +
      facet_grid(as.formula(paste(fixedparam,"~",subsetparam))) + #drop=FALSE meant to make ggplot draw the empty facets
      background_grid("none")
    
    ## save each plot with element sizes for either single plot, or grouped plot:
    epipeak.plots.group[[i]] <- epipeak.plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.title = element_text(size=8)) #adjust legend and text size as needed
    epipeak.plots.single[[i]] <- epipeak.plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(1, "cm"), legend.text = element_text(size=10), legend.title = element_text(size=12), axis.title = element_text(size=10)) #adjust legend and text size as needed
    #png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=2500, height=1300) #square with 1 row
    png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=2500, height=2725) #square with 4 rows
    plot(epipeak.plots.single[[i]])
    dev.off()
    ## save rds - workaround for return(plots)
    #saveRDS(epipeak.plots.group[[i]], file = paste(filename, "_PLOT ", i, ".rds", sep = ""))
  }
  #return(epipeak.plots) #This won't work with facet_grid
}

box_plot_combined <- function(dataset, spatmodstr, nonspatmodstr, yname, group, fixedparam, subsetparam, filename, ylab, xlab, mainlab, ylim) {
  epipeak.plots <- list() 
  epipeak.plots.single <- list() 
  epipeak.plots.group <- list() 
  data <- list()
  dataspat <- list()
  datanspat <- list()
  for (i in 1:length(dataset)) {
    data[[i]] <- dataset[[i]]
    epipeak.plots[[i]] <- ggplot() +
      geom_boxplot(data=data[[i]], aes(y=data[[i]][[yname]], x=data[[i]][[group]], colour=data[[i]][[group]])) + #can't have a piping operator for this because I need to choose groups from a named dataframe
      coord_cartesian(ylim = c(0, ylim)) + 
      labs(y=ylab, x=xlab, main=mainlab, colour=group, fill=NULL, alpha=NULL) +
      ggtitle(mainlab) +
      theme_bw() +
      facet_grid(as.formula(paste(fixedparam,"~",subsetparam))) + #drop=FALSE meant to make ggplot draw the empty facets
      background_grid("none")
    
    ## save each plot with element sizes for either single plot, or grouped plot:
    epipeak.plots.group[[i]] <- epipeak.plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.title = element_text(size=8)) #adjust legend and text size as needed
    epipeak.plots.single[[i]] <- epipeak.plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9, size=15),legend.key.size = unit(1, "cm"), legend.text = element_text(size=15), legend.title = element_text(size=19),axis.text.y = element_text(size=15), axis.title = element_text(size=19), strip.text.x = element_text(size = 15),strip.text.y = element_text(size = 15)) #adjust legend and text size as needed
    #png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=2000)
    #png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=2500, height=1053) #square #1300 for box plots with 3 boxes, 1053 for 4 boxes
    #png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=2500, height=2725) #square with 4 rows
    png(filename = paste(filename, "_spatial", "_PLOT ", i, ".png", sep = ""), res = 300, width=2900, height=3900) #square with 5 rows
    plot(epipeak.plots.single[[i]])
    dev.off()
    ## save rds - workaround for return(plots)
    #saveRDS(epipeak.plots.group[[i]], file = paste(filename, "_PLOT ", i, ".rds", sep = ""))
  }
  #return(epipeak.plots) #This won't work with facet_grid
}

##---------------------------Function: Plot proportion of outbreaks by proportion of trees infected (bar plot)
## Across entire simulation (simend)
total_infectedtrees_bar_plot <- function(data, fixedparam, subsetparam, group, filename) {
  plots <- list()
  plots.group <- list()
  plots.single <- list()
  data.sub <- list()
  for (i in 1:length(data)) {
    ## prep data:
    data.sub <- data[[i]]
    data.binned <- ddply(data.sub, c("modstr", "Ntot", "Beta","theta", "nsims", "Group", "newposition"), summarise,
                         count.n = sum(!is.na(.n)), #gives the number of successful outbreaks per grouping
                         bin.1 = sum(!is.na(.n)[total.prop.infected.trees<=0.25])/sum(!is.na(.n)), 
                         bin.2 = sum(!is.na(.n)[total.prop.infected.trees>0.25 & total.prop.infected.trees<=0.50])/sum(!is.na(.n)),
                         bin.3 = sum(!is.na(.n)[total.prop.infected.trees>0.50 & total.prop.infected.trees<=0.75])/sum(!is.na(.n)),
                         bin.4 = sum(!is.na(.n)[total.prop.infected.trees>0.75 & total.prop.infected.trees<=1.0])/sum(!is.na(.n))) %>%
      melt(id.vars = c("modstr", "Ntot", "Beta","theta", "nsims", "Group", "newposition"), measure.vars = c("bin.1", "bin.2", "bin.3", "bin.4"),
           variable.name = c("bin.category"), value.name="prop.sim")
    ## plot:
    plots[[i]] <- data.binned %>%
      ggplot() + #bars by group (e.g. structure or site)
      geom_bar(aes(fill=data.binned[[group]], y=prop.sim, x=bin.category), position="dodge", stat="identity")+
      labs(y="% successful simulations", x="Proportion of trees infected in outbreak", fill="Tree structure", main="Tree breakdown of successful outbreaks - entire outbreak") +
      ggtitle("Tree breakdown of successful outbreaks - entire outbreak") +
      theme_bw() +
      background_grid("none") +
      scale_x_discrete(labels=c("bin.1" = "0-25%", "bin.2" = "25-50%","bin.3" = "50-75%", "bin.4"="75-100%"))+
      #scale_fill_brewer(palette="RdYlBu",labels=c("Dense", "Intermediate", "Sparse", "Homogenous")) + 
      facet_grid(as.formula(paste(fixedparam,"~",subsetparam))) 
    ## save each plot with element sizes for either single plot, or grouped plot:
    plots.group[[i]] <- plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.title = element_text(size=8)) #adjust legend and text size as needed
    plots.single[[i]] <- plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(1, "cm"), legend.text = element_text(size=10), legend.title = element_text(size=12), axis.title = element_text(size=10)) #adjust legend and text size as needed
    #png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=2000)
    #png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=1140) #square with 1 row
    png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=3250) #square with 4 rows
    plot(plots.single[[i]])
    dev.off()
    ## save rds - workaround for return(plots)
    #saveRDS(plots.group[[i]], file = paste(filename, "_spatial", "_PLOT ", i,".rds", sep = ""))
  }
  #return(plots.spat) #This won't work with facet_grid(dataspat[[i]][[fixedparam]]~dataspat[[i]][[subsetparam]]). Get: "Error in `$<-.data.frame`(x, name, value) : replacement has 168 rows, data has 176"
}

## At the epidemic peak only (epipeak)
## Same as above, but for data from the epidemic peak. To change, copy and paste, and change variable name from total.prop.infected.trees to I.group.prop
epipeak_infectedtrees_bar_plot <- function(data, fixedparam, subsetparam, group, filename) {
  plots <- list()
  plots.group <- list()
  plots.single <- list()
  data.sub <- list()
  for (i in 1:length(data)) {
    ## prep data:
    data.sub <- data[[i]]
    data.binned <- ddply(data.sub, c("modstr", "Ntot", "Beta","theta", "nsims", "Group", "newposition"), summarise,
                         count.n = sum(!is.na(.n)), #gives the number of successful outbreaks per grouping
                         bin.1 = sum(!is.na(.n)[I.group.prop<=0.25])/sum(!is.na(.n)), 
                         bin.2 = sum(!is.na(.n)[I.group.prop>0.25 & I.group.prop<=0.50])/sum(!is.na(.n)),
                         bin.3 = sum(!is.na(.n)[I.group.prop>0.50 & I.group.prop<=0.75])/sum(!is.na(.n)),
                         bin.4 = sum(!is.na(.n)[I.group.prop>0.75 & I.group.prop<=1.0])/sum(!is.na(.n))) %>%
      melt(id.vars = c("modstr", "Ntot", "Beta","theta", "nsims", "Group", "newposition"), measure.vars = c("bin.1", "bin.2", "bin.3", "bin.4"),
           variable.name = c("bin.category"), value.name="prop.sim")
    ## plot:
    plots[[i]] <- data.binned %>%
      ggplot() + #bars by group (e.g. structure or site)
      geom_bar(aes(fill=data.binned[[group]], y=prop.sim, x=bin.category), position="dodge", stat="identity")+
      labs(y="% successful simulations", x="Proportion of trees infected in outbreak", fill="Tree structure", main="Tree breakdown of successful outbreaks - entire outbreak") +
      ggtitle("Tree breakdown of successful outbreaks - at epidemic peak") +
      theme_bw() +
      background_grid("none") +
      scale_x_discrete(labels=c("bin.1" = "0-25%", "bin.2" = "25-50%","bin.3" = "50-75%", "bin.4"="75-100%"))+
      #scale_fill_brewer(palette="RdYlBu",labels=c("Dense", "Intermediate", "Sparse", "Homogenous")) + 
      facet_grid(as.formula(paste(fixedparam,"~",subsetparam))) 
    ## save each plot with element sizes for either single plot, or grouped plot:
    plots.group[[i]] <- plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.title = element_text(size=8)) #adjust legend and text size as needed
    plots.single[[i]] <- plots[[i]] +
      theme(axis.text.x = element_text(angle = 45, hjust=0.9),legend.key.size = unit(1, "cm"), legend.text = element_text(size=10), legend.title = element_text(size=12), axis.title = element_text(size=10)) #adjust legend and text size as needed
    #png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=2000)
    #png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=1140) #square with 1 row
    png(filename = paste(filename, "_PLOT ", i, ".png", sep = ""), res = 300, width=3000, height=3250) #square with 4 rows
    plot(plots.single[[i]])
    dev.off()
    ## save rds - workaround for return(plots)
    #saveRDS(plots.group[[i]], file = paste(filename, "_spatial", "_PLOT ", i,".rds", sep = ""))
  }
  #return(plots.spat) #This won't work with facet_grid(dataspat[[i]][[fixedparam]]~dataspat[[i]][[subsetparam]]). Get: "Error in `$<-.data.frame`(x, name, value) : replacement has 168 rows, data has 176"
}

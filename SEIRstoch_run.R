## Title: R script to run functions for *SEIR* stochastic compartmental model with gillespie algorithm
## Author: Tamika Lunn, Griffith University
## Version: V22, created 5th May 2020

rm(list=ls())

####----------------------------------------------------------------------------------------------------
##---------------------------------Call to functions----------------------------------------------------
####----------------------------------------------------------------------------------------------------
source ("stoch_gillespie_V14.R")
source ("Transition matrices fun_V1.R")
source ("stoch_helperfunctions_V20.R")

####-------------------------------------------------------------------------------------------------------
##---------------------------------Create df of parameter values to be run---------------------------------
####----------------------------------------------------------------------------------------------------

### Read in all data to be used
## path <- set path to a folder with matrix data only. Code reads in all files by name
DATA <- read.DATA(path)
site <- read.site(path)
subplot <- read.subplot(path)

### Set values & create skeleton of data frame with parameter values to be run
Data <- 1:length(DATA) #number of dfs in list (i.e. number of roost structures being compared)
radius <- 3
theta <- c(0.1, 0.5, 1, 2, 10) 
gamma <- 1/7 #infectious period. 7 days used in Plowright 2011, 7 days used in Wang 2013, 16 days reported in Halpin et al. 2011 (RNA in urine on days 3-19). 
#sigma <- 1/1500
delta <- 1/6 #6 days to excretion post innoculation - reported in Halpin et al. 2011 *** 
nsims <- 500
threshold <- 10 
N <- c(5) #variable total abundance (set number per tree)
Ntot <- c(288, 2880, 4320) #set total abundance (variable number per tree)
modstr <- c("SEIR", "SEIRns") #***
prev <- 0.50 #start population with X seroprevalence

### Select functions for input into code #***
spatmodstr <- "SEIR" #to subset runvalues by spatial model structure ***
nonspatmodstr <- "SEIRns" #to subset runvalues by spatial model structure ***
transitions.fun <- transitions.SEIR.fun #function for defining transition matrix ***
fun.R0tree <- SIR.R0tree #function for calcuating R0 values per index tree. Don't need to change for SIRS (no impact of S on infections)
fun.R0roost <- SIR.R0roost #function for calcuating R0 values averaged across roost. Don't need to change for SIRS (no impact of S on infections)
fun.index.simulation <- SEIR.index.simulation #function for choosing index case per simulation set. Don't need to change for SIRS (no impact of S on SIR groups)
fun.summary <- SEIRsummary #function for creating summary data from simulation output. Don't need to change for SIRS (no impact of S on SIR groups)
uni.simulation <- SEIR.uni.simulation.time #function for defining when to stop simulation ***
stoch.fun <- SEIRstoch #function for main program ***
## Select between gravity model or alternate model:
fun.B <- SIR.B.gravity #function for calculating B value from R0 value. Don't need to change for SIRS (no impact of S on infections)
#fun.B <- SIR.B.plusone
name <- "G" #gravity
#name <- "PO" #plus-one

#------calculate fixed Beta-------#
R0roost <- 1/prev #calculate R0 when equalibrium
paramvalues <- list(N=N, Ntot=Ntot, Data=Data, radius=radius, theta=theta, R0roost=R0roost, gamma=gamma, delta=delta, nsims=nsims, threshold=threshold, modstr=modstr) ##note - need to specify N as first and Ntot as second elements in the list, for indexing within function ***
runvalues <- setRunvalue.fixedBeta(paramvalues, site, subplot, DATA)
DATA.long <- rescaleData.DATA(DATA, runvalues, spatmodstr, nonspatmodstr) 
DATA.COR <- rescaleData.DATA.COR(DATA, runvalues, spatmodstr, nonspatmodstr) 
SCALED.DATA <- rescaleData.SCALED.DATA(DATA, runvalues, spatmodstr, nonspatmodstr)
SCALED.DATA.ALT <- rescaleData.SCALED.DATA.ALT(DATA, runvalues, spatmodstr, nonspatmodstr) 
SCALED.DATA.DISTANCE <- rescaleData.SCALED.DATA.DISTANCE(DATA, runvalues, spatmodstr, nonspatmodstr)
runvalues$Beta <- NA
for (i in 1:length(runvalues$id)) { 
  runvalues$Beta[i] <- SIR.B.homo(Ntot[length(Ntot)/2], runvalues$gamma[i], runvalues$R0roost[i]) #Set Ntot as Ntot[i] if the goal is to keep R0 the same between population size comparisons. Set Ntot as Ntot[length(Ntot)/2] if you want a single (middle) Beta value to apply to all population sizes
} #calculate Beta after and outside of the setRunvalue, because want to keep the calculated Beta specific to the R0*Ntot*site combination. If calculate before/outside, will end up with 3 Betas per combination 
runvalues$R0roost <- NA
#----------------------------------#

### Set beta matrices for data
BETA <- setBmatrix_gravity(runvalues, SCALED.DATA, theta)
#BETA <- setBmatrix_plusone(runvalues, DATA.COR, theta)

####----------------------------------------------------------------------------------------------------
##---------------------------------Repeat simulations with param values---------------------------------
####----------------------------------------------------------------------------------------------------
for (i in 1:length(runvalues$id)) { # start loop! Repeat through each row ("id") in runvalues:

  ##---------------------------------Define conditions for model---------------------------------
  ## Set values based on values in runvalues
  data <- DATA.long[[i]] 
  scaled.data <- SCALED.DATA.DISTANCE[[i]] #used for spread function
  scaled.data.alt <- SCALED.DATA.ALT[[i]] #used for spread function
  radius <- runvalues$radius[i]
  theta <- runvalues$theta[i]
  Beta <- runvalues$Beta[i]
  B <- BETA[[runvalues$id[i]]] #sub in pre-defined B matrix (either spatial or not)
  R0roost <- runvalues$R0roost[i]
  gamma <- runvalues$gamma[i]
  #sigma <- runvalues$sigma[i] #***
  delta <- runvalues$delta[i] #***
  N <- runvalues$N[i]
  n <- nrow(data)
  #nsims <- runvalues$nsims[i]
  threshold <- runvalues$threshold[i]
  id <- runvalues$id[i]
  plotid <- runvalues$plotid[i]
  Ntot <- runvalues$Ntot[i]
  modstr <- runvalues$modstr[i]
  
  ## Set initial conditions
  params <- list(n=n, beta=B,gamma=gamma, delta=delta) #***
  Ntree <- rep(N,n) #rep specified number of bats per tree group
  S=round(Ntree*(1-prev), digits = 0) #note - this will introduce a bug if you want EXACTLY prev. Have included because code works in whole individuals
  E=rep(0,n) #***
  I=rep(0,n)
  R=Ntree-S
  start <- c(S=S, E=E, I=I, R=R,time=0, index=0, infections=0) #initial conditions ***
  
  ## Calculate R0 of individual trees
  R0tree <- fun.R0tree(n, Ntree, B, gamma) 
  R0roost <- fun.R0roost (n, Ntree, B, gamma)
  runvalues$R0roost[i] <- fun.R0roost (n, Ntree, B, gamma)
  
  ## Transitions
  transitions <- transitions.fun(n)
  
  ## Define conditions for model
  maxstep = 500000 #doesn't influence ending of simulation when .time function used, but it does need to be > the number of rows expected
  endtime = 365
  
  ##---------------------------------Run program---------------------------------
  simdat <- rdply(
    nsims,
    uni.simulation(fun.index.simulation(start),params,transitions,stoch.fun, maxstep, endtime))
  head(simdat)
  
  ##---------------------------------Create summary output---------------------------------
  ## Sum and proportion in each class
  summary.output.discrete <- fun.summary(simdat, Ntree)
  
  ## Output at epidemic peak of each simulation
  epipeak <- epidemicpeak(summary.output.discrete)
  
  ## Proportion and timing of extinction
  ext.wnt <- extinction(simdat, threshold, 84, 14) #winter period
  ext.yr <- extinction(simdat, threshold, 360, 60) #year period
  
  ## Expected proportion of extinction & R0 values, per tree
  R0tree.df <- expected(data, n, R0tree)
  
  ## Values at end of simualtion
  simend.output <- simend(simdat)
  
  ## Spread from index case
  Igroups.discrete <- groups(simdat, scaled.data, data)
  
  ### Frequency of epidemic peaks = max(simdat$infections)
  maxvals <- simdat %>%
    ddply(c(".n", "index"), summarise,
          max = max(infections), 
          end.time = max(time)
    )
  breaks <- 200
  
  ## Plots
  histplot <- ggplot(maxvals, aes(x=max)) + 
    geom_histogram(bins=breaks) +
    #xlim(c(0,max(maxvals$max+1))) +
    labs(y="Frequency of epidemic size", x="Epidemic size", title="All simulations combined") +
    theme_bw() +
    background_grid("none")
  
  ## Save counts per break
  length <- length(hist(maxvals$max, breaks=breaks, plot=FALSE)$count)
  histcount <- array(dim=c(length,2))
  colnames(histcount) <- c("breaks", "count")
  #data.frame(histcount)
  histcount[,1] <- hist(maxvals$max, breaks=breaks, plot=FALSE)$breaks[1:length] #these read as 1 to 2, 2 to 3, 3 to 4 ... n-1 (to n)
  histcount[,2] <- hist(maxvals$max, breaks=breaks, plot=FALSE)$count
  
  ##---------------------------------Save outputs---------------------------------
  param.labels <- paste('id=',id,' ',modstr,' ', name,' ',plotid,' n=',n,' B=',format(round(Beta,digits=6), scientific=FALSE),' g=',format(round(gamma,digits=4), scientific=FALSE), 'd=',format(round(delta,digits=4), scientific=FALSE),' N=',N,' Nt=',Ntot,' rad=',radius,  ' tht=',theta,' thrsh=',threshold,' nsims=',nsims,' RO=',round(R0roost, digits=3),sep='') #***
  
  ## Set a working directory
  #filename<-paste(param.labels,'.csv',sep='')
  #write.csv(simdat, file=filename, row.names=FALSE)
  
  ## Set a working directory
  #filename<-paste(param.labels,'.csv',sep='')
  #write.csv(summary.output, file=filename, row.names=FALSE)
  
  ## Set a working directory
  filename<-paste(param.labels,'.csv',sep='')
  write.csv(summary.output.discrete, file=filename, row.names=FALSE)
  
  ## Set a working directory
  filename<-paste(param.labels,'.csv',sep='')
  write.csv(epipeak, file=filename, row.names=FALSE)
  
  ## Set a working directory
  filename<-paste(param.labels,'.csv',sep='')
  write.csv(ext.wnt, file=filename, row.names=FALSE)
  
  ## Set a working directory
  filename<-paste(param.labels,'.csv',sep='')
  write.csv(ext.yr, file=filename, row.names=FALSE)
  
  ## Set a working directory
  filename<-paste(param.labels,'.csv',sep='')
  write.csv(R0tree.df, file=filename, row.names=FALSE)
  
  ## Set a working directory
  filename<-paste(param.labels,'.csv',sep='')
  write.csv(simend.output, file=filename, row.names=FALSE)
  
  ## Set a working directory
  #filename<-paste(param.labels,'.csv',sep='')
  #write.csv(Igroups, file=filename, row.names=FALSE)
  
  ## Set a working directory
  filename<-paste(param.labels,'.csv',sep='')
  write.csv(Igroups.discrete, file=filename, row.names=FALSE)
  
  ## Set a working directory
  filename <- paste(param.labels,'.csv',sep='')
  write.csv(histcount, file=filename, row.names=FALSE)
  png(filename = paste(param.labels,'.png',sep=''), res = 300, width=3000, height=3000) #square
  plot(histplot)
  dev.off()
  
} # end loop!

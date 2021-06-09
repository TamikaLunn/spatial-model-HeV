## Title: R script defining functions for stochastic compartmental models with gillespie algorithm
## Author: Tamika Lunn, Griffith University
## Version: V14, created 15th April 2020

##---------------------------------Model specification---------------------------------
## Code overview:
## Code simulates disease dynamics in example roost structure types.   
## Models are applied to real tree location data (subplots).   
## This R file specifies a stochastic within-roost compartmental meta-population model
## Separate compartmental functions ("SIRstoch", "SIRSstoch" etc.) are written. These can each be subbed into the functions to run the model

## Functions:
## The function (SIRstoch) defines the model structure: specifies parameters, variables, rates, and event selection for the SIR model
## The function (uni.simulation) specifies the bounds of the simulation, and stores simulated data in a matrix. The simulation runs until either the maxstep is reaches, or there are no more infectives
## The function (SIR.index.simulation) randomly chooses the location of index case per simulation, and resets the number of bats in initial S, I and R groups per tree for the SIR model

## Model specifications:
## Transmission is density-dependent
## Transmission between tree-groups (denoted i and j) is a function of distance
## Bij is the rate at which susceptible individuals in the ith tree group and infectious individuals in the jth tree group come into effective contact per unit time
## Si are the number of susceptible individuals in the ith tree group
## Ij are the number of infectious individuals in the jth tree group
## Yi(t) denotes the force of infection in a single tree per unit time (Bij*Ij, summed over all values of j)
## model structures have no birth/death or immigration/emmigration - is a completely closed roost

## Install libraries
library(deSolve)
library(tidyverse) 
library(plyr) #for running simulations
library(reshape2) #for melt

##---------------------------------Function: define the SIR model structure---------------------------------
SIRstoch <- function(state, params, transitions) { 
  ## Define parameters
  n <- params$n
  B <- params$beta
  gamma <- params$gamma
  ## Define states       
  S <- state[1:n]
  I <- state[(n+1):(2*n)]
  R <- state[(2*n+1):(3*n)]
  N <- S+I+R
  ## Define rates       
  rates <- c(
    infection = colSums(B*I)*S, #note: infection needs to be first rate specified for infections tally
    recovery=gamma*I #returns a vector of infection and recovery rate per group
  )
  ## Define event selection and rate
  total.rate <- sum(rates) #total rate of events
  tau <- rexp(n=1,rate=total.rate) #waiting time (note exponentially distributed random events)
  event <- sample.int(length(rates),size=1,prob=rates/total.rate) ##now which event occurs. Code says: from the list of rates, pick 1 event, given probabilities of each rate. The chosen event is index location (1:2*n) in list of rates
  state+c(transitions[event,],tau,0,ifelse(event<length(rates[grep("infection", names(rates))]),1,0)) # Choose the selected event from the tranition matrix and add to state. Note, don't need to specify which tree in state, as the transmission matrix specifies events per tree. +tau for cumulative time. + 0 for index location as doesn't change for each step. ifelse to tally new infection events (adds 1 if event choice includes "infection")
}

##---------------------------------Function: define the SEIR model structure---------------------------------
SEIRstoch <- function(state, params, transitions) { 
  ## Define parameters
  n <- params$n
  B <- params$beta
  gamma <- params$gamma
  delta <- params$delta
  ## Define states       
  S <- state[1:n]
  E <- state[(n+1):(2*n)]
  I <- state[(2*n+1):(3*n)]
  R <- state[(3*n+1):(4*n)]
  N <- S+E+I+R
  ## Define rates       
  rates <- c(
    exposure = colSums(B*I)*S,
    infection = delta*E,
    recovery=gamma*I #returns a vector of infection and recovery rate per group
  )
  ## Define event selection and rate
  total.rate <- sum(rates) #total rate of events
  tau <- rexp(n=1,rate=total.rate) #waiting time (note exponentially distributed random events)
  event <- sample.int(length(rates),size=1,prob=rates/total.rate) ##now which event occurs. Code says: from the list of rates, pick 1 event, given probabilities of each rate. The chosen event is index location (1:2*n) in list of rates
  state+c(transitions[event,],tau,0,ifelse(event>n & event<=n+length(rates[grep("infection", names(rates))]),1,0)) # Choose the selected event from the tranition matrix and add to state. Note, don't need to specify which tree in state, as the transmission matrix specifies events per tree. +tau for cumulative time. + 0 for index location as doesn't change for each step. ifelse to tally new infection events (adds 1 if event choice includes "infection")
}

##---------------------------------Function: define bounds of simulation---------------------------------
## Stop when maxstep reached or infecteds fall to zero
SIR.uni.simulation <- function (state, params, transitions, stoch.fun, maxstep) {
  output <- array(dim=c(maxstep+1,length(state)))
  colnames(output) <- names(state)
  output[1,] <- state #add current values of state to first row
  k <- 1 #k will determine which row subsequent state values get added to. In loop, k+1 steps the program down by rows
  
  ## loop until either k > maxstep or there are no more infectives
  while ((k <= maxstep) && sum(state[grep("I", names(state))]) > 0) { #grep to pick any collumn of infectives from state (named vector). Run loop until maxstep is reached, or no more infectious individuals 
    k <- k+1 
    output[k,] <- state <- stoch.fun(state,params,transitions) #add outcome of stochastic program (transition changes) to appropriate row
  }
  thin_output <- thin_sim(output[1:k,], 0:endtime)
  print(max(thin_output[,"time"]))
  return(as.data.frame(thin_output))
}

SEIR.uni.simulation <- function (state, params, transitions, stoch.fun, maxstep) {
  output <- array(dim=c(maxstep+1,length(state)))
  colnames(output) <- names(state)
  output[1,] <- state #add current values of state to first row
  k <- 1 #k will determine which row subsequent state values get added to. In loop, k+1 steps the program down by rows
  
  ## loop until either k > maxstep or there are no more infectives
  while ((k <= maxstep) && sum(state[grep("I", names(state))] + state[grep("E", names(state))]) > 0) { #grep to pick any collumn of exposed or infectives from state (named vector). Run loop until maxstep is reached, or no more exposed or infectious individuals 
    k <- k+1 
    output[k,] <- state <- stoch.fun(state,params,transitions) #add outcome of stochastic program (transition changes) to appropriate row
  }
  thin_output <- thin_sim(output[1:k,], 0:endtime)
  print(max(thin_output[,"time"]))
  return(as.data.frame(thin_output))
}

## Alternatively stop when specified time is reached, or infecteds fall to zero
SIR.uni.simulation.time <- function (state, params, transitions, stoch.fun, maxstep, endtime) {
  output <- array(dim=c(maxstep+1,length(state)))
  colnames(output) <- names(state)
  output[1,] <- state #add current values of state to first row
  k <- 1 #k will determine which row subsequent state values get added to. In loop, k+1 steps the program down by rows
  
  ## loop until either time > endtime or there are no more infectives
  while ((state["time"] <= endtime) && sum(state[grep("I", names(state))]) > 0) { #grep to pick any collumn of infectives from state (named vector). Run loop until maxstep is reached, or no more infectious individuals 
    k <- k+1 
    output[k,] <- state <- stoch.fun(state,params,transitions) #add outcome of stochastic program (transition changes) to appropriate row
  }
  thin_output <- thin_sim(output[1:k,], 0:endtime)
  print(max(thin_output[,"time"]))
  return(as.data.frame(thin_output))
}

SEIR.uni.simulation.time <- function (state, params, transitions, stoch.fun, maxstep, endtime) {
  output <- array(dim=c(maxstep+1,length(state)))
  colnames(output) <- names(state)
  output[1,] <- state #add current values of state to first row
  k <- 1 #k will determine which row subsequent state values get added to. In loop, k+1 steps the program down by rows
  
  ## loop until either time > endtime or there are no more infectives
  while ((state["time"] <= endtime) && sum(state[grep("I", names(state))] + state[grep("E", names(state))]) > 0) { #grep to pick any collumn of exposed or infectives from state (named vector). Run loop until maxstep is reached, or no more exposed or infectious individuals 
    k <- k+1 
    output[k,] <- state <- stoch.fun(state,params,transitions) #add outcome of stochastic program (transition changes) to appropriate row
  }
  thin_output <- thin_sim(output[1:k,], 0:endtime)
  print(max(thin_output[,"time"]))
  return(as.data.frame(thin_output))
}

##---------------------------------Function: choose location of index case per simulation---------------------------------
SIR.index.simulation <- function(start){{
  index <- sample(1:n, 1)
  S[index]=S[index]-1
  I[index]=I[index]+1 #1 infected bat in randomly chosen tree
  xstart <- c(S=S, I=I, R=R,time=0, index=index, infections=1)
}
  xstart
}

SEIR.index.simulation <- function(start){{
  index <- sample(1:n, 1)
  S[index]=S[index]-1
  I[index]=I[index]+1 #1 infected bat in randomly chosen tree
  xstart <- c(S=S,E=E, I=I, R=R,time=0, index=index, infections=1)
}
  xstart
}


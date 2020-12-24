## Title: R script defining transition matrices for stochastic compartmental models with gillespie algorithm
## Author: Tamika Lunn, Griffith University
## Version: V1, created 24th Feb 2020

### Functions:
## Transmission matrix specifies what happens to each state for each event. 
## Events and state are split by group number. So e.g. an infection in group 2 ("Inf2") only changes the states in group 2 ("S2", "I2", "R2"). 
## Code fills in a blank matrix by collumns. Each line in code corresponds to event changes for individual states (S(E)IR). Numbers in rep(c()) can be read in blocks of 3 (each block of 3 representing an event), with the middle number being the change in state with the event and the zeros either side controlling the spacing 
## Transitions coded below include:
##    - transitions.SIR.fun
##    - transitions.SIRS.fun
##    - transitions.SEIR.fun
##    - transitions.SEIRS.fun

transitions.SIR.fun <- function(n) {
  transitions.SIR <- matrix(ncol = 3*n, nrow = 2*n) #Set blank matrix (ncol= #states*n, nrow= #events*n)
  position = 0
  for(counter in 1:n){
    transitions.SIR[,counter]<-rep(c(0,-1,0,0,0,0), times = c(position,1,n-counter,position,1,n-counter)) #fill S states by collumn (through events)
    transitions.SIR[,counter+n]<-rep(c(0,1,0,0,-1,0), times = c(position,1,n-counter,position,1,n-counter)) #fill I states by collumn (through events)
    transitions.SIR[,counter+n*2]<-rep(c(0,0,0,0,1,0), times = c(position,1,n-counter,position,1,n-counter)) #fill R states by collumn (through events)
    position <- position+1
  }
  return(transitions.SIR)
}

transitions.SIRS.fun <- function(n) {
  transitions.SIRS <- matrix(ncol = 3*n, nrow = 3*n) #Set blank matrix (ncol= #states*n, nrow= #events*n)
  position = 0
  for(counter in 1:n){
    transitions.SIRS[,counter]<-rep(c(0,-1,0,0,0,0,0,1,0), times = c(position,1,n-counter,position,1,n-counter,position,1,n-counter)) #fill S states by collumn (through events)
    transitions.SIRS[,counter+n]<-rep(c(0,1,0,0,-1,0,0,0,0), times = c(position,1,n-counter,position,1,n-counter,position,1,n-counter)) #fill I states by collumn (through events)
    transitions.SIRS[,counter+n*2]<-rep(c(0,0,0,0,1,0,0,-1,0), times = c(position,1,n-counter,position,1,n-counter,position,1,n-counter)) #fill R states by collumn (through events)
  position <- position+1
  } #fill values
  return(transitions.SIRS)
}

transitions.SEIR.fun <- function(n) {
  transitions.SEIR <- matrix(ncol = 4*n, nrow = 3*n) #Set blank matrix (ncol= #states*n, nrow= #events*n)
  position = 0
  for(counter in 1:n){
    transitions.SEIR[,counter]<-rep(c(0,-1,0,0,0,0,0,0,0), times = c(position,1,n-counter,position,1,n-counter,position,1,n-counter)) #fill S states by collumn (through events)
    transitions.SEIR[,counter+n]<-rep(c(0,1,0,0,-1,0,0,0,0), times = c(position,1,n-counter,position,1,n-counter,position,1,n-counter)) #fill E states by collumn (through events)
    transitions.SEIR[,counter+n*2]<-rep(c(0,0,0,0,1,0,0,-1,0), times = c(position,1,n-counter,position,1,n-counter,position,1,n-counter)) #fill I states by collumn (through events)
    transitions.SEIR[,counter+n*3]<-rep(c(0,0,0,0,0,0,0,1,0), times = c(position,1,n-counter,position,1,n-counter,position,1,n-counter)) #fill R states by collumn (through events)
  position <- position+1
  } #fill values
  return(transitions.SEIR)
}

transitions.SEIRS.fun <- function(n) { 
  transitions.SEIRS <- matrix(ncol = 4*n, nrow = 4*n) #Set blank matrix (ncol= #states*n, nrow= #events*n)
  position = 0
  for(counter in 1:n){
    transitions.SEIRS[,counter]<-rep(c(0,-1,0,0,0,0,0,0,0,0,1,0), times = c(position,1,n-counter,position,1,n-counter,position,1,n-counter,position,1,n-counter)) #fill S states by collumn (through events) (rep(Exp*3, Inf*3, Rec*3, Waning*3)) #infection before exposure to match "new infections" algorithm in main program
    transitions.SEIRS[,counter+n]<-rep(c(0,1,0,0,-1,0,0,0,0,0,0,0), times = c(position,1,n-counter,position,1,n-counter,position,1,n-counter,position,1,n-counter)) #fill E states by collumn (through events)
    transitions.SEIRS[,counter+n*2]<-rep(c(0,0,0,0,1,0,0,-1,0,0,0,0), times = c(position,1,n-counter,position,1,n-counter,position,1,n-counter,position,1,n-counter)) #fill I states by collumn (through events)
    transitions.SEIRS[,counter+n*3]<-rep(c(0,0,0,0,0,0,0,1,0,0,-1,0), times = c(position,1,n-counter,position,1,n-counter,position,1,n-counter,position,1,n-counter)) #fill R states by collumn (through events)
  position <- position+1
  } #fill values
  return(transitions.SEIRS)
}
# spatial-model-HeV
Data and code repository for Lunn et al (2021) 'Roost tree structure is an important driver of Hendra virus infection within communally roosting Pteropus species'

Three datafiles are provided, these are distance matrices (meters) of roosting trees. The set give an example of a i) sparse tree distribution, ii) intermediate tree distribution and iii) dense tree distribution. Visuals of the tree structures are also provided

Code is organised into 5 files:
- stoch_gillespie.R : Defines functions for stochastic compartmental models with gillespie algorithm (model structure, bounds of simulation and location of index case per simulation)
- stoch_helperfunctions.R : Miscellaneous functions for the running of code (e.g. compiling simuations, plotting output, reading lists)
- Transition matrices fun.R : Defines transition matrices for model
- SIRstoch_run.R : Performs SIR stochastic simulation
- SEIRstoch_run.R : Performs SEIR stochastic simulation


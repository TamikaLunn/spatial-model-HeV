# spatial-model-HeV
Data and code repository for Lunn et al (2021) 'Spatial dynamics of pathogen transmission in communally roosting species: impacts of changing habitats on bat-virus dynamics'

Code is organised into an RProject containing 8 files:
- stoch_gillespie_VX.R : Defines functions for stochastic compartmental models with gillespie algorithm (model structure, bounds of simulation and location of index case per simulation)
- stoch_helperfunctions_VX.R : Miscellaneous functions for the running of code (e.g. compiling simuations, plotting output, reading lists)
- Transition matrices fun_VX.R : Defines transition matrices for model. An explanation of transition matrices is given in the word document 'Transmission matrix structures.docx'
- SIRstoch_run_VX.R : Performs SIR stochastic simulation
- SEIRstoch_run_VX.R : Performs SEIR stochastic simulation
- SIRstoch_stats_VX.R : Code to plot SIR model output
- SEIRstoch_stats_VX.R : Code to plot SEIR model output
- Distribution of PW distances.R : Calculates meta-data of pairwise distances between trees in the stand
- ** Important to load the code as the RProject so that root folders in setwd are correct **

Data is organised into an RProject containing 6 files:
- /Raw : Gives data matricies (meters) used in models. Naming links with other projects using this data. DTOW = sparse tree distrubution; DCLU = intermediate tree distrubution; DTOW = dense tree distrubution. Visuals of the stand structures are also provided in the root folder
- /Processed : Contains extra meta-data used to investigate patterns in the data (including pairwise distances between all trees, maximum pairwise distance per tree stand)

Output is organised into an RProject containing many files:
- Figures_SIR_stoch outputs_PO_p50 gives the figures used in the paper
- SIR_stoch outputs_PO_p50 gives simulation outputs that were analysed in the paper, for the SIR model structure (SIR), the main function (PO) (i.e. not the gravity function, G), at a starting seroprevalence of 50 (p50)
- SEIR_stoch outputs_PO_p50 gives simulation outputs for the equivalent SEIR model, analysed in the paper
- New run provides a blank template of the output folder directory to match saving system in the code

# spatial-model-HeV
Data and code repository for Lunn et al (2021) 'Spatial dynamics of pathogen transmission in communally roosting species: impacts of changing habitats on bat-virus dynamics'

To explore how infection dynamics are influenced by heterogeneity in stand structure, we applied spatially explicit and stochastic compartmental models to three empirical examples of flying-fox roost stand structures, representing sparse, intermediate, and dense stand structures respectively. Generation of the spatial model structure needs input of a distance matrix between tree-groups, and specification of how transmission is expected to relate to distance. 

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
- /Raw : Gives data matricies (meters) used in models. Naming links with other projects using this data. DTOW = sparse tree distrubution; DCLU = intermediate tree distrubution; DTOW = dense tree distrubution, which include the distances between 4 trees (tree #69-72), 32 trees (tree #90-#121), and 72 trees (tree #61-#132) in subplots, respectively. The data is structured as a standard distance matrix, where trees are named on the first row and first collumn. Visuals of the stand structures are also provided in the root folder
- /Processed : Contains extra meta-data used to investigate patterns in the data (including pairwise distances between all trees, maximum pairwise distance per tree stand)

Output is organised into an RProject containing many files:
- Figures_SIR_stoch outputs_PO_p50 gives the figures used in the paper
- SIR_stoch outputs_PO_p50 gives simulation outputs that were analysed in the paper, for the SIR model structure (SIR), the main function (PO) (i.e. not the gravity function, G), at a starting seroprevalence of 50 (p50)
- SEIR_stoch outputs_PO_p50 gives simulation outputs for the equivalent SEIR model, analysed in the paper
- New run provides a blank template of the output folder directory to match saving system in the code

Brief overview of methods used to generate distance matrices:

The distance matrices used in the manuscript are provided. These represent a tree stand structure, where the spatial arrangement of all overstory, canopy and midstory trees were mapped in a subplot (20x20 meters each) using an ultrasound distance instrument (Vertex Hypsometer, Haglöf Sweden). We did not map trees or shrubs in the understory as these are no suitable roosting habitat. Trees were mapped and tagged using tree survey methods described in the “Ausplots Forest Monitoring Network, Large Tree Survey Protocol” (https://portal.tern.org.au/tern-ausplots-forest-2012-2015/21755). Briefly, subplots were georeferenced at one corner. Distances were measured from the N/S or E/W subplot boundaries using an ultrasound distance instrument (Vertex Hypsometer, Haglöf Sweden, accurate to 10-30 cm) along the defined orientation bearing. Trees within the subplot were then mapped with the X-Y coordinate in relation to the georeferenced corner (0,0). To achieve maximum accuracy with the Vertex Hypsometer, only distances of up to 10 meters were recorded. If a tree was greater than 10 meters from the west/south origin (0 meter) subplot boundary, the tree was measured from the opposite (20 meter) subplot boundary, and the measured distance subtracted from 20 to give the distance from the origin boundary. Each tree was individually tagged and assigned a crown class following definitions in the Ausplots survey protocol. This approach allowed for precise spatial mapping of trees, with locations of trees within subplots accurate to 10-30 cm.

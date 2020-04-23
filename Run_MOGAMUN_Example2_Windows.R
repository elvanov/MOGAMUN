# load functions and main code
source("MOGAMUN_FUN.R")
source("MOGAMUN.R")

# mandatory function to setup needed parameters
setUpMetaParameters()

# setup the 8 parameters that are needed for the evolution process:
#    Generations - number of generations to run (default = 500)
#    PopSize - number of subnetworks in the population (default = 100)
#    MinNumberOfNodesPerIndividual - minimum size of the subnetworks (default = 15)
#    MaxNumberOfNodesPerIndividual - maximum size of the subnetworks (default = 50)
#    CrossoverRate - rate for the crossover (default = 0.8, i.e., 80%)
#    MutationRate - rate for the mutation (default = 0.1, i.e., 10%)
#    JaccardSimilarityThreshold - subnetworks over this Jaccard similarity threshold are considered as duplicated (default = 30)
#    TournamentSize - size of the tournament (default = 2)
#    MyObjectiveNames - list containing the names of the objectives (default =  c("AverageNodesScore", "Density"))
setUpEvolutionProcessParameters(Generations=3)

# define the data to be loaded for the run:
#    DEPath - full path to the differential expression results file (in CSV format). This file must contain It must contain at least the columns "gene" with the gene names, and ("PValue" or "FDR"). It can also contain "logFC"
#    NodesScoresPath - full path to the file containing the nodes scores (in CSV format). It must contain two columns: "gene" and "nodescore". NOTE. If no file is available, it wll be automatically generated
#    NetworkLayersDir - path of the folder that contains the networks that will be the layers of the multiplex
#    Layers - string of numbers, where the numbers correspond to the first character of the name of the network files
loadData(DEPath = "C:/Users/User/Documents/My_Experiment/Differential_Expression_Results.csv",
         NodesScoresPath = "C:/Users/User/Documents/My_Experiment/Genes_With_Nodes_Scores.csv",
         NetworkLayersDir = "C:/Users/User/Documents/Layers_Multiplex/",
         Layers = "123")

# runs mogamun and outputs the results in the specified folder in resultsDir
mogamun(resultsDir = "C:/Users/User/Documents/My_Experiment/MOGAMUN_Results/")
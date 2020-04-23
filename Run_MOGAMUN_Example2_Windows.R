source("MOGAMUN_FUN.R")
source("MOGAMUN.R")
setUpMetaParameters()
setUpEvolutionProcessParameters(Generations=3)
loadData(DEPath = "C:/Users/User/Documents/My_Experiment/Differential_Expression_Results.csv",
         NodesScoresPath = "C:/Users/User/Documents/My_Experiment/Genes_With_Nodes_Scores.csv",
         NetworkLayersDir = "C:/Users/User/Documents/Layers_Multiplex/",
         LayersToUse <- "123")
mogamun(resultsDir = "C:/Users/User/Documents/My_Experiment/MOGAMUN_Results/")
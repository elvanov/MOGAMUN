source("MOGAMUN_FUN.R")
source("MOGAMUN.R")
setUpMetaParameters()
setUpEvolutionProcessParameters(Generations=2)
loadData(DEPath = "~/My_Experiment/Differential_Expression_Results.csv",
         NodesScoresPath = "~/My_Experiment/Genes_With_Nodes_Scores.csv",
         NetworkLayersDir = "~/Layers_Multiplex/",
         Layers <- "123")
mogamun(resultsDir = "~/My_Experiment/MOGAMUN_Results/")

# change accordingly to the folder containing the unzipped ExampleFiles/ folder
setwd('.')

# init parameters. Check ?mogamun.init for usage
mogamun.init(Generations=5)

# load data files. Check ?mogamun.load.data for usage
mogamun.load.data(DifferentialExpressionPath = "DifferentialExpressionData/Banerji2017.csv",
         NodesScoresPath = "DifferentialExpressionData/Banerji2017_NodesScore",
         NetworkLayersDir = "LayersMultiplex/",
         Layers = "123")

# run. Check ?mogamun.run for usage 
mogamun.run(resultsDir = "MOGAMUN_Results/")

# ?mogamun.init, ?mogamun.load.data, ?mogamun.run to open the man pages.
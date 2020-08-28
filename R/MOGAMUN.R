# Multi objective Genetic algorithm to look for differentially expressed 
# subnetworks in a multiplex network
# Integer representation (nodes' IDs), size of the individuals:
# between 15 and 50 CONNECTED nodes, tournament selection,
# customized crossover and mutation operators (to never lose connectivity),
# elitist generational replacement


#' @title mogamun_init
#'
#' @description initialize evolution parameters
#'
#' @param Generations number of generations to run (default = 500)
#' @param PopSize number of subnetworks in the population (default = 100)
#' @param MinSize minimum size (no. of nodes) of the subnetworks (default = 15)
#' @param MaxSize maximum size (no. of nodes) of the subnetworks (default = 50)
#' @param CrossoverRate rate for the crossover (default = 0.8)
#' @param MutationRate rate for the mutation (default = 0.1)
#' @param JaccardSimilarityThreshold subnetworks over this Jaccard similarity 
#' threshold are considered as duplicated (default = 30)
#' @param TournamentSize size of the tournament (default = 2)
#' @param Measure measure to calculate the nodes scores and to determine which 
#' genes are differentially expressed 
#' (possible values PValue and FDR, default = FDR)
#' @param ThresholdDEG threshold to consider a gene as significantly 
#' differerentially expressed. Note: if there is a logFC available, it is also 
#' considered |logFC|>1  (default = 0.05)
#' @param MaxNumberOfAttempts maximum number of attempts to find compatible 
#' parents (default = 3) 
#'
#' @return EvolutionParameters
#'
#' @examples
#' EvolutionParameters <- 
#'     mogamun_init(
#'         Generations = 1,
#'         PopSize = 10,
#'         MinSize = 15,
#'         MaxSize = 50,
#'         CrossoverRate = 0.8,
#'         MutationRate = 0.1,
#'         JaccardSimilarityThreshold = 30,
#'         TournamentSize = 2,
#'         Measure = "FDR",
#'         ThresholdDEG = 0.05,
#'         MaxNumberOfAttempts = 3
#'     )
#'
#' @export
#' @import doParallel igraph stringr RUnit
#' @importFrom devtools install_github
#' @importFrom foreach `%dopar%` foreach
#' @importFrom utils write.table read.table combn read.csv write.csv 
#' @importFrom stats runif qnorm
#' @importFrom RCy3 cytoscapePing createNetworkFromDataFrames loadTableData 
#' @importFrom RCy3 setEdgeColorMapping setNodeColorMapping 
#' @importFrom RCy3 setNodeBorderColorMapping setNodeBorderWidthMapping 
#' @importFrom RCy3 layoutNetwork createColumnFilter createSubnetwork 
#' @importFrom RCy3 closeSession saveSession
#' @importFrom graphics boxplot plot legend
#' @importFrom grDevices svg dev.off rainbow
mogamun_init <- function(Generations = 500, PopSize = 100,
    MinSize = 15, MaxSize = 50,
    CrossoverRate = 0.8, MutationRate = 0.1, JaccardSimilarityThreshold = 30,
    TournamentSize = 2, Measure = "FDR", ThresholdDEG = 0.05,
    MaxNumberOfAttempts = 3) {

    pkgname <- c("doParallel", "igraph", "stringr")
    
    for (p in pkgname) {
        require(p, quietly = TRUE, character.only = TRUE) || 
            stop("Package '", p, "' not found")
    }
    
    # determines the parameters that will be used for the evolution
    EvolutionParameters <- list(
        Generations = Generations, 
        PopSize = PopSize, 
        MinSize = MinSize,
        MaxSize = MaxSize,
        CrossoverRate = CrossoverRate,
        MutationRate = MutationRate, 
        JaccardSimilarityThreshold = JaccardSimilarityThreshold,
        TournamentSize = TournamentSize,
        ObjectiveNames = c("AverageNodesScore", "Density"),
        Measure = Measure,
        ThresholdDEG = ThresholdDEG,
        MaxNumberOfAttempts = MaxNumberOfAttempts 
    )
    
    return(EvolutionParameters)    
}

#' @title mogamun_load_data
#'
#' @description Load the data to process
#'
#' @param EvolutionParameters evolution paramenters returned by mogamun_init()
#' @param DifferentialExpressionPath full path to the differential expression 
#' results file (in CSV format). This file must contain at least the columns 
#' "gene" with the gene names, and ("PValue" or "FDR"). 
#' It can also contain "logFC"
#' @param NodesScoresPath full path to an existing CSV file containing the 
#' nodes scores (columns "gene" and "nodescore"). NOTE. If the file does not 
#' exist, MOGAMUN will generate it in the provided path with the specified name
#' @param NetworkLayersDir path of the folder that contains the networks that 
#' will be the layers of the multiplex. NOTE. Each file must start with a 
#' different digit
#' @param Layers string of numbers, where the numbers correspond to the first 
#' character of the name of the network files (e.g. "123" builds a multiplex 
#' with layers 1, 2, and 3)
#'
#' @return List with the data to process
#'
#' @examples
#' DEGPath <- system.file("extdata/DE/Sample_DE.csv", package = "MOGAMUN")
#' NodesScoresPath <- 
#'     system.file("extdata/DE/Sample_NodesScore.csv", package = "MOGAMUN")
#' LayersPath <- 
#'     paste0(system.file("extdata/LayersMultiplex", package = "MOGAMUN"), "/")
#' EvolutionParameters <- mogamun_init(Generations = 1, PopSize = 10)
#' LoadedData <- 
#'     mogamun_load_data(
#'         EvolutionParameters = EvolutionParameters,
#'         DifferentialExpressionPath = DEGPath,
#'         NodesScoresPath = NodesScoresPath,
#'         NetworkLayersDir = LayersPath,
#'         Layers = "23"
#'     )
#' @export

mogamun_load_data <- function(EvolutionParameters, DifferentialExpressionPath, 
    NodesScoresPath, NetworkLayersDir, Layers) {
    Measure <- EvolutionParameters$Measure # "FDR" or "PValue"
    ThresholdDEG <- EvolutionParameters$ThresholdDEG # threshold for DEG
    DifferentialExpressionPath <- DifferentialExpressionPath # path for DE res
    NodesScoresPath <- NodesScoresPath # full path to the nodes score
    NetworkLayersDir <- NetworkLayersDir # folder containing the networks 
    LayersToUse <- Layers # layers to use to build the multiplex
    DE_results <- data.frame(read.csv(DifferentialExpressionPath)) # load DE
    DE_results <- RemoveDuplicates_DE_Results(DE_results) # remove dup entries
    DEG <- DE_results[DE_results$FDR < ThresholdDEG, ] # get list of DEG
    
    # verify existence of log(fold change) and consider it for the DEG
    if ("logFC" %in% colnames(DEG)) { DEG <- DEG[abs(DEG$logFC) > 1, ] }
    
    # read the file names of the networks for the current experiment
    Files <- list.files(NetworkLayersDir, pattern = 
                paste0("^[", LayersToUse, "]_"), full.names = TRUE)
    
    if ( ! file.exists(NodesScoresPath) ) { # if no nodes scores file exists
        # calculate the nodes scores for all the genes in DE analysis results
        NodesScores <- as.numeric(GetNodesScoresOfListOfGenes(DE_results, 
            as.character(DE_results$gene), Measure))
        
        # data frame of genes and scores. NOTE. Genes not in the list have 0
        GenesWithNS <- data.frame("gene" = as.character(DE_results$gene), 
            "nodescore" = NodesScores)
        write.csv(GenesWithNS, file = NodesScoresPath, row.names = FALSE)
    } else {
        GenesWithNS <- 
            data.frame(read.csv(NodesScoresPath, stringsAsFactors = FALSE))
    }
    
    Multiplex <- GenerateMultiplexNetwork(Files) # make the multiplex network
    Merged <- GenerateMergedNetwork(Files, Multiplex) # make the merged network 
    DensityPerLayerMultiplex <- unlist(lapply(Multiplex, graph.density)) # dens
    
    LoadedData <- c(EvolutionParameters, list(
        NetworkLayersDir = NetworkLayersDir, Layers = Layers, 
        DE_results = DE_results, DEG = DEG, GenesWithNodesScores = GenesWithNS,
        Multiplex = Multiplex, 
        DensityPerLayerMultiplex = DensityPerLayerMultiplex,
        Merged = Merged))
    
    return (LoadedData)
}

#' @title mogamun_run
#'
#' @description Run the algorithm with the specified values for the evolution 
#' parameters
#'
#' @param LoadedData list returned by mogamun_load_data()
#' @param Cores to run MOGAMUN in parallel on the given number of cores (in 
#' line with the number of physical processor cores) (default = 1)
#' @param NumberOfRunsToExecute number of runs (default = 1)
#' @param ResultsDir outputs the results in the specified folder
#'
#' @return None
#'
#' @examples
#' 
#' DEGPath <- system.file("extdata/DE/Sample_DE.csv", package = "MOGAMUN")
#' NodesScoresPath <- 
#'     system.file("extdata/DE/Sample_NodesScore.csv", package = "MOGAMUN")
#' LayersPath <- 
#'     paste0(system.file("extdata/LayersMultiplex", package = "MOGAMUN"), "/")
#' EvolutionParameters <- mogamun_init(Generations = 1, PopSize = 10)
#' LoadedData <- 
#'     mogamun_load_data(
#'         EvolutionParameters = EvolutionParameters,
#'         DifferentialExpressionPath = DEGPath,
#'         NodesScoresPath = NodesScoresPath,
#'         NetworkLayersDir = LayersPath,
#'         Layers = "23"
#'     )
#' ResultsDir <- paste0(system.file("SampleResults", package="MOGAMUN"), "/")
#' mogamun_run(
#'     LoadedData = LoadedData,
#'     Cores = 1,
#'     NumberOfRunsToExecute = 1,
#'     ResultsDir = ResultsDir
#' )
#' @export
mogamun_run <- function(LoadedData, Cores = 1, NumberOfRunsToExecute = 1,
    ResultsDir = '.') {
    registerDoParallel(cores = Cores) # in line with the no. of physical cores
    PopSize <- LoadedData$PopSize
    Generations <- LoadedData$Generations

    if (exists("LoadedData")) {
        Multiplex <- LoadedData$Multiplex
        ResultsPath <- paste0(ResultsDir, "/Experiment_", Sys.Date(), "/")
        dir.create(ResultsPath, recursive = TRUE)  # create result folder 
        BestIndsPath <- paste0(ResultsPath, "MOGAMUN_Results_") # path for res
        
        RunNumber <- 1
        
        # loop to execute the algorithm many times
        foreach(RunNumber = seq_len(NumberOfRunsToExecute)) %dopar% {
            BestIndsFile <- paste0(BestIndsPath, "_Run_", RunNumber, ".txt")
            MyInitPop <- GenerateInitialPop(PopSize, Multiplex, LoadedData) 
            FitnessData <- EvaluatePopulation(MyInitPop, Multiplex, LoadedData)
            Population <- data.frame("Individual" = I(MyInitPop), FitnessData) 
            
            # obtain ranking and crowding distances
            Population <- NonDomSort(
                PopulationToSort = Population, LoadedData = LoadedData )
            
            g <- 1  # initilizes the number of generation
            StatsGen <- data.frame(matrix(ncol = 3, nrow = 0))
            colnames(StatsGen) <- 
                c( "Generation", "BestAverageNodesScore", "BestDensity" )
            
            # evolution's loop for g generations or until all inds have rank = 1
            while (g <= Generations && !all(Population$Rank == 1)) {
                Population <- MakeNewPopulation(LoadedData, Population) 
                
                # add the best values for the two objective functions
                StatsGen[nrow(StatsGen) + 1, ] <- c(g, 
                    max(Population$AverageNodesScore), max(Population$Density))
                print(paste0("Run ", RunNumber, ". Gen. ", g, " completed"))
                g <- g + 1 # increments the generation
            }   
            # saves data in files
            write.csv(StatsGen, file = paste0(BestIndsPath,
                "StatisticsPerGeneration_Run", RunNumber, ".csv"), 
                row.names = FALSE)
            SaveFinalPop(BestIndsFile, Population, PopSize, Multiplex[[1]]) 
            print(paste0("FINISH TIME, RUN ", RunNumber, ": ", Sys.time()))
            gc()
        }        
    } else {if (!exists(LoadedData)) {print ("Missing parameter: LoadedData")}}
}



#' @title mogamun_postprocess
#'
#' @description Postprocess the results. This function: 
#' i) calculates the accumulated Pareto front, i.e. the individuals on the 
#' first Pareto front after re-ranking the results from multiple runs
#' (NOTE. If there is a single run, the result is the set of individuals in 
#' the first Pareto front),
#' ii) filters the networks to leave only the interactions between the genes 
#' that are included in the results, 
#' iii) generates some plots of interest, such as scatter plots and boxplots, 
#' and
#' iv) (optional) creates a Cytoscape file to visualize the results, merging 
#' the subnetworks with a Jaccard similarity coefficient superior to 
#' JaccardSimilarityThreshold 
#' (NOTE. Make sure to open Cytoscape if VisualizeInCytoscape is TRUE)
#'
#' @param ExperimentDir folder containing the results to be processed. It is 
#' the same folder specified as ResultsDir in mogamun_run
#' @param LoadedData list returned by mogamun_load_data()
#' @param JaccardSimilarityThreshold subnetworks over this Jaccard similarity 
#' threshold are merged in a single subnetwork
#' @param VisualizeInCytoscape TRUE if you wish to visualize the accumulated 
#' Pareto front in Cytoscape, FALSE otherwise
#'
#' @return None
#'
#' @examples
#' 
#' DEGPath <- system.file("extdata/DE/Sample_DE.csv", package = "MOGAMUN")
#' NodesScoresPath <- 
#'     system.file("extdata/DE/Sample_NodesScore.csv", package = "MOGAMUN")
#' LayersPath <- 
#'     paste0(system.file("extdata/LayersMultiplex", package = "MOGAMUN"), "/")
#' EvolutionParameters <- mogamun_init(Generations = 1, PopSize = 10)
#' LoadedData <- 
#'     mogamun_load_data(
#'         EvolutionParameters = EvolutionParameters,
#'         DifferentialExpressionPath = DEGPath,
#'         NodesScoresPath = NodesScoresPath,
#'         NetworkLayersDir = LayersPath,
#'         Layers = "23"
#'     )
#' ResultsDir <- paste0(system.file("SampleResults", package="MOGAMUN"), "/")
#' mogamun_run(
#'     LoadedData = LoadedData,
#'     Cores = 1,
#'     NumberOfRunsToExecute = 1,
#'     ResultsDir = ResultsDir
#' )
#' mogamun_postprocess(
#'     ExperimentDir = ResultsDir,
#'     LoadedData = LoadedData,
#'     JaccardSimilarityThreshold = 70,
#'     VisualizeInCytoscape = FALSE
#' ) 
#'
#' @export
mogamun_postprocess <- function(ExperimentDir = '.', LoadedData = LoadedData,
    JaccardSimilarityThreshold = 70, VisualizeInCytoscape = TRUE) {
    
    PostprocessResults(ExperimentDir = ExperimentDir, LoadedData = LoadedData, 
        JaccardSimilarityThreshold = JaccardSimilarityThreshold, 
        VisualizeInCytoscape = VisualizeInCytoscape)
}

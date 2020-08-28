# These are the unit test to check the functionality of the four main functions
# of MOGAMUN

test.mogamun_init <- function() {
    EvolutionParameters <- 
        mogamun_init(
            Generations = 1, PopSize = 10, MinSize = 10, MaxSize = 20, 
            CrossoverRate = 1, MutationRate = 1, 
            JaccardSimilarityThreshold = 30, TournamentSize = 2, 
            Measure = "FDR", ThresholdDEG = 0.05, MaxNumberOfAttempts = 1)
    
    # check length 
    checkEquals(length(EvolutionParameters), 12)
    
    # check values
    checkEquals(EvolutionParameters$Generations, 1)
    checkEquals(EvolutionParameters$PopSize, 10)
    checkEquals(EvolutionParameters$MinSize, 10)
    checkEquals(EvolutionParameters$MaxSize, 20)
    checkEquals(EvolutionParameters$CrossoverRate, 1)
    checkEquals(EvolutionParameters$MutationRate, 1)
    checkEquals(EvolutionParameters$JaccardSimilarityThreshold, 30)
    checkEquals(EvolutionParameters$TournamentSize, 2)
    checkEquals(EvolutionParameters$Measure, "FDR")
    checkEquals(EvolutionParameters$ThresholdDEG, 0.05)
    checkEquals(EvolutionParameters$MaxNumberOfAttempts, 1)
    checkEquals(
        EvolutionParameters$ObjectiveNames, c("AverageNodesScore", "Density"))
}


test.mogamun_load_data <- function() {
    DEGPath <- system.file("extdata/DE/Sample_DE.csv", package = "MOGAMUN")
    
    NodesScoresPath <- 
        system.file("extdata/DE/Sample_NodesScore.csv", package = "MOGAMUN")
    
    LayersPath <- 
        paste0(system.file("extdata/LayersMultiplex", package = "MOGAMUN"), "/")
    
    EvolutionParameters <- 
        mogamun_init(
            Generations = 1, PopSize = 10, MinSize = 10, MaxSize = 20, 
            CrossoverRate = 1, MutationRate = 1, MaxNumberOfAttempts = 1)
    
    LoadedData <- 
        mogamun_load_data(
            EvolutionParameters = EvolutionParameters,
            DifferentialExpressionPath = DEGPath, 
            NodesScoresPath = NodesScoresPath, NetworkLayersDir = LayersPath, 
            Layers = "23")
    
    # check length
    checkEquals(length(LoadedData), 20)
    
    # check that all the data in EvolutionParameters is in LoadedData 
    checkTrue(all(names(EvolutionParameters) %in% names(LoadedData)))
    
    # check that the values in EvolutionParameters are the same in LoadedData
    checkIdentical(
        LoadedData[names(LoadedData) %in% names(EvolutionParameters)], 
        EvolutionParameters)
 
    # check other values
    checkEquals(LoadedData$NetworkLayersDir, LayersPath)
    checkEquals(LoadedData$Layers, "23")
    checkTrue(all(c("gene", EvolutionParameters$Measure) %in% 
        names(LoadedData$DE_results)))
    
    # check types
    checkTrue(all(is.numeric(
        LoadedData$DE_results[[EvolutionParameters$Measure]])))
    checkTrue(all(LoadedData$DE_results[[EvolutionParameters$Measure]] < 1))
    checkTrue(all(is.character(LoadedData$DE_results$gene)))
    checkTrue(all(LoadedData$DEG[[EvolutionParameters$Measure]] < 
        EvolutionParameters$ThresholdDEG))
    checkTrue(all(abs(LoadedData$DEG$logFC) > 1))
    checkTrue(all(LoadedData$GenesWithNodesScores$gene %in% 
        LoadedData$DE_results$gene))
    checkTrue(all(LoadedData$GenesWithNodesScores$nodescore >= 0 && 
        LoadedData$GenesWithNodesScores$nodescore <=1))
    checkEquals(length(LoadedData$Multiplex), nchar(LoadedData$Layers))
    
    for (i in seq_len(length(LoadedData$Multiplex))) {
        checkTrue(is.igraph(LoadedData$Multiplex[[i]]))
    }
    
    checkEquals(length(LoadedData$DensityPerLayerMultiplex), 
        length(LoadedData$Multiplex))
    checkTrue(all(is.numeric(LoadedData$DensityPerLayerMultiplex)))
    checkTrue(is.igraph(LoadedData$Merged))
}

test.mogamun_run <- function() {
    DEGPath <- system.file("extdata/DE/Sample_DE.csv", package = "MOGAMUN")
    
    NodesScoresPath <- 
        system.file("extdata/DE/Sample_NodesScore.csv", package = "MOGAMUN")
    
    LayersPath <- 
        paste0(system.file("extdata/LayersMultiplex", package = "MOGAMUN"), "/")
    
    EvolutionParameters <- 
        mogamun_init(
            Generations = 1, PopSize = 10, MinSize = 10, MaxSize = 20, 
            CrossoverRate = 1, MutationRate = 1, MaxNumberOfAttempts = 1)
    
    LoadedData <- 
        mogamun_load_data(
            EvolutionParameters = EvolutionParameters,
            DifferentialExpressionPath = DEGPath, 
            NodesScoresPath = NodesScoresPath, NetworkLayersDir = LayersPath, 
            Layers = "23")
    
    set.seed(123)
    ResDir <- paste0(system.file("SampleResults", package="MOGAMUN"), "/")
    mogamun_run(LoadedData = LoadedData, ResultsDir = ResDir)    
    
    # check existence of results' folder and files 
    ResDirAll <- paste0(ResDir, "Experiment_", Sys.Date())
    checkTrue(dir.exists(ResDirAll))
    checkEquals(length(list.files(ResDirAll)), 2)
    Files <- list.files(ResDirAll, pattern = "MOGAMUN_Results_")
    checkIdentical(c("MOGAMUN_Results__Run_1.txt", 
        "MOGAMUN_Results_StatisticsPerGeneration_Run1.csv"), Files)
    SampleFiles <- list.files(ResDir, recursive = FALSE)

    # check that the file that was created corresponds to the expected output
    for (f in Files) {
        RealOutput <- readLines(paste0(ResDirAll, "/", f))
        ExpectedOutput <- readLines(paste0(ResDir, f))
        checkIdentical(RealOutput, ExpectedOutput)
    }
    
    # remove directory and files
#    unlink(ResDirAll, recursive = TRUE)
}

test.mogamun_postprocess <- function() {
    DEGPath <- system.file("extdata/DE/Sample_DE.csv", package = "MOGAMUN")
    
    NodesScoresPath <- 
        system.file("extdata/DE/Sample_NodesScore.csv", package = "MOGAMUN")
    
    LayersPath <- 
        paste0(system.file("extdata/LayersMultiplex", package = "MOGAMUN"), "/")
    
    EvolutionParameters <- 
        mogamun_init(
            Generations = 1, PopSize = 10, MinSize = 10, MaxSize = 20, 
            CrossoverRate = 1, MutationRate = 1, MaxNumberOfAttempts = 1)
    
    LoadedData <- 
        mogamun_load_data(
            EvolutionParameters = EvolutionParameters,
            DifferentialExpressionPath = DEGPath, 
            NodesScoresPath = NodesScoresPath, NetworkLayersDir = LayersPath, 
            Layers = "23")
    
    set.seed(123)
    ResDir <- paste0(system.file("SampleResults", package="MOGAMUN"), "/")
    mogamun_run(LoadedData = LoadedData, ResultsDir = ResDir)     
    
    mogamun_postprocess(
        ExperimentDir = ResDir,
        LoadedData = LoadedData,
        JaccardSimilarityThreshold = 70,
        VisualizeInCytoscape = FALSE
    ) 
    
    ResDirAll <- paste0(ResDir, "Experiment_", Sys.Date(), "/")
    
    # check existence of files
    ResFiles <- list.files(ResDirAll)
    ExpectedFiles <- 
        c(
            "2_Pathways_Sample_FILTERED.csv", 
            "3_Coexpression_Sample_FILTERED.csv", 
            "A_AccPF_CYTOSCAPE_JacSimT_70.csv", "A_AccPF_JacSimT_70.csv", 
            "A_ALL_FILTERED_INTERACTIONS_CYTOSCAPE.csv", 
            "A_ScatterPlot_AccPF.svg" 
        )
    checkTrue(all(ExpectedFiles %in% ResFiles))
    
    SampleFiles <- list.files(ResDir, recursive = FALSE)
    
    # check content of files
    for (f in SampleFiles) {
        RealOutput <- readLines(paste0(ResDirAll, "/", f))
        ExpectedOutput <- readLines(paste0(ResDir, f))
        checkIdentical(RealOutput, ExpectedOutput)
    }
    
    # remove directory and files
#    unlink(ResDirAll, recursive = TRUE)
}



# Multi objective Genetic algorithm to look for differentially expressed subnetworks in a multiplex network
# Integer representation (nodes' IDs), size of the individuals: 
# between 15 and 50 CONNECTED nodes, tournament selection,
# customized crossover and mutation operators (to never lose connectivity), 
# elitism of 1 and generational replacement, each node has a score associated 


#' @title mogamun.init
#'
#' @description initialize evolution parameters
#'
#' @param Generations number of generations to run (default = 500)
#' @param PopSize number of subnetworks in the population (default = 100)
#' @param MinNumberOfNodesPerIndividual minimum size of the subnetworks (default = 15)
#' @param MaxNumberOfNodesPerIndividual maximum size of the subnetworks (default = 50)
#' @param CrossoverRate rate for the crossover (default = 0.8)
#' @param MutationRate rate for the mutation (default = 0.1)
#' @param JaccardSimilarityThreshold subnetworks over this Jaccard similarity threshold are considered as duplicated (default = 30)
#' @param TournamentSize size of the tournament (default = 2)
#' @param ObjectiveNames list containing the names of the objectives (default =  c("AverageNodesScore", "Density"))
#' @param ThresholdSignificantlyDEGenes threshold to consider a gene as significantly differerentially expressed. Note: if there is a logFC available, it is also considered |logFC|>1  (default = 0.05)
#' @param MaxNumberOfAttempts maximum number of attempts to find compatible parents (default = 3) 
#' @param Measure measure to calculate the nodes scores and to determine which genes are differentially expressed (possible values PValue and FDR, default = FDR)
#' @param NumberOfRunsToExecute number of runs (default = 1)
#'
#' @return None
#'
#' @examples
#' mogamun.init(Generations = 5,
#'                PopSize = 100,
#'                MinNumberOfNodesPerIndividual = 15,
#'                MaxNumberOfNodesPerIndividual = 50,
#'                CrossoverRate = 0.8,
#'                MutationRate = 0.1,
#'                JaccardSimilarityThreshold = 30,
#'                TournamentSize = 2,
#'                ObjectiveNames = c("AverageNodesScore", "Density"),
#'                ThresholdSignificantlyDEGenes = 0.05,
#'                MaxNumberOfAttempts = 3,
#'                Measure = "FDR",
#'                NumberOfRunsToExecute = 1)
#'
#' @export
#' 
mogamun.init <- function(Generations = 500,
                         PopSize = 100,
                         MinNumberOfNodesPerIndividual = 15,
                         MaxNumberOfNodesPerIndividual = 50,
                         CrossoverRate = 0.8,
                         MutationRate = 0.1,
                         JaccardSimilarityThreshold = 30,
                         TournamentSize = 2,
                         ObjectiveNames = c("AverageNodesScore", "Density"),
                         ThresholdSignificantlyDEGenes = 0.05,
                         MaxNumberOfAttempts = 3,
                         Measure = "FDR",
                         NumberOfRunsToExecute = 1) {
     # determines the parameters that will be used for the evolution
     Generations <<- Generations # 500 # test with different values and check convergence
     PopSize <<- PopSize # test with different values and check if the fitness improves significantly with different sizes
     MinNumberOfNodesPerIndividual <<- MinNumberOfNodesPerIndividual
     MaxNumberOfNodesPerIndividual <<- MaxNumberOfNodesPerIndividual
     CrossoverRate <<- CrossoverRate
     MutationRate <<- MutationRate # test with values from 0.1 - 0.5
     JaccardSimilarityThreshold <<- JaccardSimilarityThreshold
     TournamentSize <<- TournamentSize
     ObjectiveNames <<- ObjectiveNames
     ThresholdSignificantlyDEGenes <<- ThresholdSignificantlyDEGenes
     MaxNumberOfAttempts <<- MaxNumberOfAttempts # this is used to find compatible parents and to create a valid size individual
     Measure <<- Measure # determines how to get the DE genes, either by 'PValue' or 'FDR'
     NumberOfRunsToExecute <<- NumberOfRunsToExecute
}



#' @title mogamun.load.data
#'
#' @description Load the data to process
#'
#' @param DifferentialExpressionPath - full path to the differential expression results file (in CSV format). This file must contain It must contain at least the columns "gene" with the gene names, and ("PValue" or "FDR"). It can also contain "logFC"
#' @param NodesScoresPath - full path to the file containing the nodes scores (in CSV format). It must contain two columns: "gene" and "nodescore". NOTE. If no file is available, it wll be automatically generated
#' @param NetworkLayersDir - path of the folder that contains the networks that will be the layers of the multiplex
#' @param Layers - string of numbers, where the numbers correspond to the first character of the name of the network files
#'
#' @return None
#'
#' @examples
#' mogamun.init()
#' mogamun.load.data(DifferentialExpressionPath =
#'                    "~/My_Experiment/Differential_Expression_Results.csv",
#'         NodesScoresPath = "~/My_Experiment/Genes_With_Nodes_Scores.csv",
#'         NetworkLayersDir = "~/My_Experiment/Layers_Multiplex/",
#'         Layers = "123")
#'
#' @export
mogamun.load.data <- function(DifferentialExpressionPath, NodesScoresPath, NetworkLayersDir, Layers) {
     # full path for the DE results
     # NOTE. It must contain at least the columns "gene" and ("PValue" or "FDR"). It can also contain "logFC"
     DifferentialExpressionPath <- DifferentialExpressionPath
     
     # full path to the file containing the nodes scores.
     # If the file pointed by the path does not exist, it will be created by mogamun,
     # Otherwise it will be used by the algorithm.
     # The file must contain 2 columns: "gene" and "nodescore"
     NodesScoresPath <- NodesScoresPath
     
     # Define the directory where the networks to be used as layer are located
     # Layers' filenames must start by a single digit identifier'
     NetworkLayersDir <- NetworkLayersDir
     
     # Layers to use in the experiment (first character of filename). "123" builds a multiplx with layers 1* 2* and 3*
     LayersToUse <- Layers
     
     # load DE analysis results
     DE_results <<- data.frame(read.csv(DifferentialExpressionPath))
     DE_results <<- RemoveDuplicates_DE_Results(DE_results)
     
     # get list of differentially expressed genes
     DifferentiallyExpressedGenes <<-
          DE_results[ abs(DE_results$logFC) > 1 &
                           DE_results$FDR < ThresholdSignificantlyDEGenes, ]
     
     # read the file names of the networks for the current experiment
     Files <- list.files(NetworkLayersDir, pattern = paste0("^[", LayersToUse, "]_"), full.names=T)
     
     # if no nodes scores file was given
     if ( ! file.exists(NodesScoresPath) ) {
          # calculates the nodes scores for all the genes in DE analysis results
          NodesScores <- as.numeric(GetNodesScoresOfListOfGenes(as.character(DE_results$gene), Measure))
          
          # makes a data frame with the gene names and their corresponding nodes scores
          # IMPORTANT!!! ANY gene not included in this list, has node score = 0
          GenesWithNodesScores <<- data.frame("gene" = as.character(DE_results$gene),  "nodescore" = NodesScores)
          
          write.csv(GenesWithNodesScores, file = NodesScoresPath, row.names = FALSE)
     } else {
          GenesWithNodesScores <<- data.frame(read.csv(NodesScoresPath, stringsAsFactors = FALSE))
     }     

     # generates the multiplex and merged network
     Multiplex <<- GenerateMultiplexNetwork(Files)
     Merged <<- GenerateMergedNetwork(Files)
     
     # calculate density per layer of the multiplex
     DensityPerLayerMultiplex <<- unlist(lapply(Multiplex, graph.density))
}


#' @title mogamun.run
#'
#' @description Run the algo
#'
#' @param resultsDir outputs the results in the specified folder
#'
#' @return None
#'
#' @examples
#' mogamun.init(Generations=0)
#' mogamun.load.data(DifferentialExpressionPath =
#'                    "~/My_Experiment/Differential_Expression_Results.csv",
#'         NodesScoresPath = "~/My_Experiment/Genes_With_Nodes_Scores.csv",
#'         NetworkLayersDir = "~/My_Experiment/Layers_Multiplex/",
#'         Layers = "123")
#' mogamun.run(resultsDir = '~/My_Experiment/MOGAMUN_Results/')
#'
#' @export
mogamun.run <- function(resultsDir = '.') {
     
     # uncomment the following two lines if executing multiple runs in parallel. Also uncomment the foreach of the main loop and comment the for
     # library(doParallel)
     # registerDoParallel(cores = 1) # should be in line with the number of physical processor cores
     
     # Directory where results are stored
     ResultsDir <- resultsDir
     
     # create result folder for new experiment
     ResultsPath <- paste0(ResultsDir, "Experiment_", Sys.Date(), "/")
     dir.create(ResultsPath, recursive=T)
     
     # defines the path and filename to store the results
     BestIndividualsPath <- paste0(ResultsPath, "MOGAMUN_Results_")
     
     ########################################################################################################################
     ############# --------------------------------- M A I N   L O O P -------------------------#############################
     ########################################################################################################################
     # loop to execute the algorithm many times
     # to measure the execution time: comment the following line if you are not benchmarking
     #ptime <- system.time(
     
     #foreach(RunNumber = 1:NumberOfRunsToExecute) %dopar% {
     for (RunNumber in 1:NumberOfRunsToExecute) {
          BestIndividualsFile <- paste0(BestIndividualsPath, "_Run_", RunNumber, ".csv")
          
          # ---------------------------------
          # ----- Create initial population of specified size
          # ---------------------------------
          MyInitialPopulation <- GenerateInitialPopulation(PopSize = PopSize, Multiplex = Multiplex)
          
          # --------------------------------
          # ----- Evaluate individuals
          # ---------------------------------
          FitnessData <- EvaluatePopulation(MyInitialPopulation, Multiplex, GenesWithNodesScores)
          
          # generate data frame with the individuals and their fitness
          Population <- data.frame("Individual" = I(MyInitialPopulation), FitnessData)
          
          # obtain ranking and crowding distances
          Population <- SortByNonDominationAndObtainCrowdingDistance(PopulationToSort = Population)
          
          
          # --------------------------------
          # ----- Evolution process starts here
          # ---------------------------------
          g <- 1  # initilizes the number of generation
          
          StatisticsPerGeneration <- data.frame(matrix(ncol = 3, nrow = 0))
          colnames(StatisticsPerGeneration) <- c(
               "GenerationNumber",
               "AverageNodesScoreOfBestIndividualInGeneration",
               "DensityOfBestIndividualInGeneration"
          )
          
          # evolution's loop for a given number of generations or untill all the individuals are in the first Pareto front
          while (g <= Generations && !all(Population$Rank == 1)) {
               MyNewPopulation <- vector("list", PopSize) # initialize an empty vector of the population size
               
               # loop to generate the new population. In each loop, 2 children are created
               for(i in seq(from=1, to=PopSize, by=2)) {
                    # initialize control variables
                    AttemptsToFindParents <- 0 # counter to control maximum number of attemps to find parents
                    KeepLooking <- TRUE # flag to keep looking for parents
                    
                    
                    while ( AttemptsToFindParents < MaxNumberOfAttempts & KeepLooking == TRUE ) {
                         # --------------------------------
                         # ----- Selection of parents (tournament)
                         # ---------------------------------
                         Parent1 <- TournamentSelection(TournamentSize, Population) # get parent 1
                         
                         ####### ---------------------------------------------------------------------
                         # with the new crossover, only "near" parents can mate, so the second parent
                         # has to be chosen with respect to the first one, therefore we need to filter
                         # the population first
                         ####### ---------------------------------------------------------------------
                         
                         # get the nodes' ids of parent 1, with respect to the big network
                         NodesIDsOfParent1 <- sort( unlist(Parent1$Individual) )
                         
                         # get the list of all the nodes IDs in parent 1 and their neighbors
                         AllNeighborsNodesOfParent1 <- GetNeighborsOfNodeList(NodesIDsOfParent1, Multiplex)
                         
                         # get the list of individuals in the population that contain at least one
                         # node from the previous list (nodes in parent 1 and their neighbors)
                         IndividualsInTheNeighborhoodOfParent1 <-
                              unlist( sapply( 1:PopSize, function(X) {
                                   if (length(intersect(unlist(Population[X,"Individual"]), AllNeighborsNodesOfParent1)) > 0){
                                        X
                                   } } ) )

                         # filter the original population to leave only those individuals that are near parent 1
                         PotentialIndividualsForParent2 <- Population[IndividualsInTheNeighborhoodOfParent1, ]
                         
                         # verify if parent 1 is in the list of potential individuals for parent 2
                         if ( rownames(Parent1) %in% rownames(PotentialIndividualsForParent2) ) {
                              # if this is the case, remove it from the list
                              IDofParent1InParent2List <- which(rownames(PotentialIndividualsForParent2) == rownames(Parent1))
                              PotentialIndividualsForParent2 <- PotentialIndividualsForParent2[ -IDofParent1InParent2List, ]
                         }
                         
                         # verify the length of the list of potential parent 2
                         if ( nrow(PotentialIndividualsForParent2) >= 2 ) {
                              
                              # perform the tournament to choose parent 2, on the filtered population
                              Parent2 <- TournamentSelection(TournamentSize, PotentialIndividualsForParent2) # get parent 2
                              KeepLooking <- FALSE
                         } else if (nrow(PotentialIndividualsForParent2) == 1) { # verify if there is 1 potential parent2
                              # if only one individual is compatible with parent 1, use it as parent 2
                              Parent2 <- PotentialIndividualsForParent2
                              KeepLooking <- FALSE
                         } else {
                              # increment the number of attempts to find the parents
                              AttemptsToFindParents <- AttemptsToFindParents + 1
                         }
                    }
                    
                    if (AttemptsToFindParents == MaxNumberOfAttempts) {
                         ### ***** ADD NEW GENERATED RANDOM INDIVIDUAL *****
                         print("Maximum number of attemps to find compatible parents reached. Adding two random individuals to the new population.")
                         
                         # generate two random individuals
                         Children <- GenerateInitialPopulation(PopSize = 2, Multiplex = Multiplex)
                    } else {
                         # --------------------------------
                         # ----- Crossover
                         # ---------------------------------
                         Children <- Crossover(Parent1, Parent2)
                    }
                    
                    # --------------------------------
                    # ----- Mutation
                    # ---------------------------------
                    Children <- Mutation(Children, Multiplex)
                    
                    # NOTE. Doing "MyNewPopulation[i] <- list(Children[[1]])" and "MyNewPopulation[i] <- Children[1]" is equivalent
                    MyNewPopulation[i] <- Children[1] # add individual to the population
                    MyNewPopulation[i+1] <- Children[2] # add individual to the population
               }
               
               # ----- Here we already have a whole new population
               # evaluate offspring
               FitnessData <- EvaluatePopulation(MyNewPopulation, Multiplex, GenesWithNodesScores)                    
                         
               # generate data frame with the individuals and their fitness
               NewPopulation <-
                    data.frame(
                         "Individual" = I(MyNewPopulation),
                         FitnessData,
                         Rank = rep(0, nrow(FitnessData)),
                         CrowdingDistance = rep(0, nrow(FitnessData))
                    )
               
               # --------------------------------
               # ----- Replacement
               # ---------------------------------
               NewPopulationForReplacement <- Replacement(Parents = Population, Children = NewPopulation)

               # replace old population
               Population <- NewPopulationForReplacement
               
               StatisticsPerGeneration[nrow(StatisticsPerGeneration) + 1, ] <-
                    c(g, # "GenerationNumber"
                      max(Population$AverageNodesScore), # "AverageNodesScoreOfBestIndividualInGeneration"
                      max(Population$Density) # "DensityOfBestIndividualInGeneration"
                    )
               
               print(paste0("Generation ", g, " completed"))
               
               g <- g + 1 # increments the generation
          }     ### end of      while (g <= Generations)
          
          # saves the statistics of the population
          write.csv(
               StatisticsPerGeneration,
               file = paste0(BestIndividualsPath,
                             "StatisticsPerGeneration_Run",
                             RunNumber,
                             ".csv"), row.names = FALSE
          )
          
          # saves the best individuals of the final population
          SaveTheBestIndividualsFromFinalPopulation(
               BestIndividualsFile,
               Population,
               PopSize,
               Multiplex[[1]] # this layer will be used to get the names of the corresponding genes, but
               # actually all the layers have the nodes in the same order, so it is indistinct
          )          

          print(paste0("FINISH TIME, RUN ", RunNumber, ": ", Sys.time()))
          
          gc()
     }
}

# Multi objective Genetic algorithm to look for differentially expressed subnetworks in a multiplex network
# Integer representation (nodes' IDs), size of the individuals: 
# between 15 and 50 CONNECTED nodes, tournament selection,
# customized crossover and mutation operators (to never lose connectivity), 
# elitism of 1 and generational replacement, each node has a z-score associated 

source('~/Doctorado/Code_files/MyGeneticAlgorithm/MyAlgorithm/MO_GA_MUN_PARALLEL_ALL_FUNCTIONS.R')

library(doParallel)
registerDoParallel(cores = 2) # should be in line with the number of physical processor cores

# in order to make the random functions NOT random at all
# for the purpose of benchmarking, comparing the reproducibility of results, etc.
# set.seed(1)

ThresholdSignificantlyDEGenes <- 0.05
#PathToSaveZScores <- "~/Doctorado/Data_files/"
PathToSaveZScores <- "~/Doctorado/Data_files/New_TCGA_Data/Downloaded_May_2019/"
# PathToSaveZScores <- "~/Doctorado/Data_files/Axis3_ResultsGA_MyProposal/Experiment_2019-01-23/"

MinNumberOfNodesPerIndividual <- 15
MaxNumberOfNodesPerIndividual <- 50
MaxNumberOfAttempts <- 3 # this is used to find compatible parents and to create a valid size individual
Measure <- "FDR" # determines how to get the DE genes, either by 'PValue' or 'FDR'

ExperimentToRun <- 5
NumberOfRunsToExecute <- 2
# NumberOfRunsToExecute <- 1


########################################################################################################
############# DEFINITION OF PARAMETERS FOR THE EVOLUTION PROCESS #######################################
########################################################################################################
# determines the parameters that will be used for the evolution
Generations <- 500 # test with different values and check convergence
PopSize <- 100 # test with different values and check if the fitness improves significantly with different sizes
#CrossoverRates <- c(0.8, 0.9, 1)
CrossoverRates <- 0.8
#MutationRates <- c(0.1, 0.2, 0.3) # test with values from 0.1 - 0.5
MutationRates <- 0.1 # test with values from 0.1 - 0.5
JaccardSimilarityThreshold <- 30
TournamentSize <- 2
MyObjectiveNames <- c("AverageZScore", "Density")

# definition of the experiments' names
MyExperimentsNames <- 
  c(
    "Exp1_Monoplex_PPI", 
    "Exp2_Monoplex_Pathways", 
    "Exp3_Monoplex_Coexp", 
    "Exp4_Monoplex_Merged",
    "Exp5_Multiplex_AvZSc1x"
  )

# DE_AnalysisResults_FileNameAndPath <-
#   "~/Doctorado/Data_files/DE_analysis_results_FIBROBLAST_LogFC_SingCorrected_WithDE_Pvalue_FDR.csv"

DE_AnalysisResults_FileNameAndPath <-
  "~/Doctorado/Data_files/New_TCGA_Data/Downloaded_May_2019/DE_edgeR_without_control_88_with_gene_names.csv"

# load DE analysis results
DE_results <- data.frame(read.csv(DE_AnalysisResults_FileNameAndPath))

DE_results <- RemoveDuplicates_DE_Results(DE_results)

# get list of differentially expressed genes
DifferentiallyExpressedGenes <- 
  DE_results[ abs(DE_results$logFC) > 1 & 
                DE_results[[Measure]] < ThresholdSignificantlyDEGenes, ]


# SET TO NULL WHEN NO Z-SCORE HAS BEEN CALCULATED
ZScores_FileAndPath <- "~/Doctorado/Data_files/New_TCGA_Data/Downloaded_May_2019/GenesWithZScores.csv"
# ZScores_FileAndPath <- NULL

# set working path
WorkingPath <- "~/Doctorado/Data_files/New_TCGA_Data/Downloaded_May_2019/"

# define the path where the networks to be used as layer are 
MyLayersPath <- paste0(WorkingPath, "LayersForTheMultiplex/")
# MyLayersPath <- "~/Doctorado/Data_files/LayersForTheMultiplex_SMALL_FOR_TESTING/"

LayersToUseForExperiment <- c(1:4, "123")

# read the file names of the networks for the current experiment 
Files <- list.files(MyLayersPath, pattern = paste0("^[", LayersToUseForExperiment[ExperimentToRun], "]_"))

# create folder for new experiment
WorkingPath <- paste0(WorkingPath, "Experiment_", Sys.time(), "/")
dir.create(WorkingPath)

# if no z-scores file was given
if ( is.null(ZScores_FileAndPath) ) {
  # calculates the z-scores for all the genes in DE analysis results
  ZScores <- as.numeric(GetZScoresOfListOfGenes(as.character(DE_results$gene), Measure))

  # makes a data frame with the gene names and their corresponding z-scores
  # IMPORTANT!!! ANY gene not included in this list, has Z-score = 0
  GenesWithZScores <- data.frame("gene" = as.character(DE_results$gene),  "zscore" = ZScores)

  write.csv(GenesWithZScores, file = paste0(PathToSaveZScores, "GenesWithZScores.csv"), row.names = FALSE)
} else {
  GenesWithZScores <- data.frame(read.csv(ZScores_FileAndPath, stringsAsFactors = FALSE))
}

# generates the multiplex and merged network 
Multiplex <- GenerateMultiplexNetwork(Files)
Merged <- GenerateMergedNetwork(Files)

# calculate density per layer of the multiplex
DensityPerLayerMultiplex <- unlist(lapply(Multiplex, graph.density))

for (CrossoverRate in CrossoverRates) {
  for (MutationRate in MutationRates) {

    # create folder with the experiment's name
    ExperimentPath <- 
      paste0(WorkingPath, MyExperimentsNames[ExperimentToRun], "_FDR_Crossover_", CrossoverRate, "_Mutation_", MutationRate, "/")
    dir.create(ExperimentPath)
    
    # defines the path and filename to store the results
    BestIndividualsPath <- 
      paste0(ExperimentPath, "Results_GeneticAlgorithm_Crossover_", CrossoverRate, "_Mutation_", MutationRate)
    
    ########################################################################################################################
    ############# --------------------------------- M A I N   L O O P -------------------------#############################
    ########################################################################################################################
    # loop to execute the algorithm many times
  # to measure the execution time: comment the following line if you are not benchmarking
  ptime <- system.time(
    foreach(RunNumber = 1:NumberOfRunsToExecute) %dopar% {
    # for (RunNumber in 1:NumberOfRunsToExecute) {
      BestIndividualsFile <- paste0(BestIndividualsPath, "_Run_", RunNumber, "_", Sys.time(), ".csv")
      
      # ---------------------------------
      # ----- Create initial population of specified size
      # ---------------------------------
      MyInitialPopulation <- GenerateInitialPopulation(PopSize = PopSize, Multiplex = Multiplex)
      
      # --------------------------------
      # ----- Evaluate individuals
      # ---------------------------------
      FitnessData <- EvaluatePopulation(MyInitialPopulation, Multiplex, GenesWithZScores)
    
      # generate data frame with the individuals and their fitness
      Population <- 
        data.frame("Individual" = I(MyInitialPopulation), FitnessData)
        
      # obtain ranking and crowding distances
      Population <- SortByNonDominationAndObtainCrowdingDistance(PopulationToSort = Population)
      
      # --------------------------------
      # ----- Evolution process starts here
      # ---------------------------------
      g <- 1  # initilizes the number of generation
      
      BestIndividualFitnessPerGeneration <- NULL
      StatisticsPerGeneration <- data.frame(matrix(ncol = 6, nrow = 0))
      colnames(StatisticsPerGeneration) <- c(
        "GenerationNumber", 
        "AverageIndividualsSizeInPopulation", 
        "AverageIndividualsFitnessInPopulation",
        "FitnessOfBestIndividualInGeneration",
        "AverageZScoreOfBestIndividualInGeneration",
        "DensityOfBestIndividualInGeneration"
        )
        
      DiversityInPopulation <- data.frame(
        Children = rep(0, Generations),
        CombinedPopulation = rep(0, Generations),
        AfterReplacement = rep(0, Generations)
      ) 
      
      # evolution's loop for a given number of generations or untill all the individuals are in the first Pareto front
      while (g <= Generations && !all(Population$Rank == 1)) {
        BestDensityBefore <- max(Population$Density)
        
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
        FitnessData <- EvaluatePopulation(MyNewPopulation, Multiplex, GenesWithZScores)
    
        # generate data frame with the individuals and their fitness
        NewPopulation <- 
          data.frame(
            "Individual" = I(MyNewPopulation), 
            FitnessData, 
            Rank = rep(0, nrow(FitnessData)), 
            CrowdingDistance = rep(0, nrow(FitnessData))
          )
    
        DiversityInPopulation$Children[g] <- nrow(unique(NewPopulation[, MyObjectiveNames]))
    
        WholePopulation <- rbind(Population[, MyObjectiveNames], NewPopulation[, MyObjectiveNames])
        
        DiversityInPopulation$CombinedPopulation[g] <- nrow(unique(WholePopulation[, MyObjectiveNames]))
        
        # --------------------------------
        # ----- Replacement
        # ---------------------------------
        NewPopulationForReplacement <- Replacement(Parents = Population, Children = NewPopulation)
        
        BestDensityAfter <- max(NewPopulationForReplacement$Density)
        
        stopifnot(BestDensityAfter >= BestDensityBefore)

        DiversityInPopulation$AfterReplacement[g] <- nrow(unique(NewPopulationForReplacement[, MyObjectiveNames]))
        
        # replace old population
        Population <- NewPopulationForReplacement
    
        # get the sizes of each individual in the population
        IndividualsSizesNP <- sapply(Population$Individual, function(X) { length(X) })
    
        # save best individual's fitness to check for convergence
        BestIndividualFitnessPerGeneration <- 
          rbind(BestIndividualFitnessPerGeneration, max(Population$Fitness))
        
        StatisticsPerGeneration[nrow(StatisticsPerGeneration) + 1, ] <- 
          c(g, # "GenerationNumber"
            mean(IndividualsSizesNP), # "AverageIndividualsSizeInPopulation"
            mean(Population$Fitness), #"AverageIndividualsFitnessInPopulation"
            max(Population$Fitness), # "FitnessOfBestIndividualInGeneration"
            max(Population$AverageZScore), # "AverageZScoreOfBestIndividualInGeneration"
            max(Population$Density) # "DensityOfBestIndividualInGeneration"
          )
          
        print(paste0("Generation ", g, ": ", sum(Population$Rank == 1), " individuals in first Pareto front, max. density = ", max(Population$Density), ", max. AverageZScore = ", max(Population$AverageZScore)))

        g <- g + 1 # increments the generation
      }     ### end of      while (g <= Generations)
          
    
      # saves the best individual fitness per generation in a file for the convergence plots 
      write.csv(
        BestIndividualFitnessPerGeneration, 
        file = paste0(BestIndividualsPath, 
                      "EvolutionOfBestIndividual_Run", 
                      RunNumber, 
                      "_", 
                      Sys.time(), ".csv")
      )
      
      # saves the statistics of the population 
      write.csv(
        StatisticsPerGeneration, 
        file = paste0(BestIndividualsPath, 
                      "StatisticsPerGeneration_Run", 
                      RunNumber, 
                      "_", 
                      Sys.time(), 
                      ".csv")
      )
      
      # saves the evolution of diversity of the population in a file for the convergence plots 
      write.csv(
        DiversityInPopulation, 
        file = paste0(BestIndividualsPath, 
                      "EvolutionOf_DIVERSITY_Run", 
                      RunNumber, 
                      "_", 
                      Sys.time(), ".csv")
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
      Sys.time()  
    }  
    
    # to print out exectution time: comment the 2 following lines if you are not benchmarking
  ) # end of ptime <- system.time( 
    print(ptime)
    
    gc()
  } # Mutation rates loop
} # Crossover rates loop


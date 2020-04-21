# Multi objective Genetic algorithm to look for differentially expressed subnetworks in a multiplex network
# Integer representation (nodes' IDs), size of the individuals: 
# between 15 and 50 CONNECTED nodes, tournament selection,
# customized crossover and mutation operators (to never lose connectivity), 
# elitism of 1 and generational replacement, each node has a node score associated 
library(igraph)
library(nsga2R)

#################### --------------------------------------------------------------- ####################
#################### --------------------------------------------------------------- ####################
#################### --------------------- ALL MY FUNCTIONS ------------------------ ####################
#################### --------------------------------------------------------------- ####################
#################### --------------------------------------------------------------- ####################


# definition of the function that filters the DE results, to remove all the duplicates 
# NOTE. This function had to be added becuase many reapeated genes with different expression values were found in the TCGA dataset
# INPUTS: DE_results: data frame with 5 columns (output of edgeR), {gene, logFC, logCPM, PValue, FDR}
# OUTPUT: DE analysis results, with only unique values
RemoveDuplicates_DE_Results <- function(DE_results) {
  # remove NAs
  DE_results <- DE_results[!is.na(DE_results$gene), ]
     
  dup <- unique(as.character(DE_results$gene[which(duplicated( as.character(DE_results$gene) ) )]))
  
  DE_results_filtered <- DE_results[!as.character(DE_results$gene) %in% dup, ]
  
  for (d in dup) {
    DupRows <- DE_results[as.character(DE_results$gene) == d, ]
    
    DE_results_filtered <- rbind(DE_results_filtered, DupRows[which.min(DupRows$FDR), ])
  }
  
  return (DE_results_filtered)
}


# definition of the function that generates the multiplex network 
# INPUTS: Files - list of file names that contain the biological networks data
# OUTPUT: Multiplex network (list of igraph objects)
GenerateMultiplexNetwork <- function(Files) {
  # declare empty list to store the multiplex
  Multiplex <- list()
  
  AllNodesToTake <- as.character(GetListOfAllNodesPresentInLayers(Files))
  
  # loop through all the layers to get the corresponding subnetwork from each of them
  for (LayerFile in Files) {
    # load layer
    Layer <- data.frame(read.table(paste0(MyLayersPath, LayerFile), header = FALSE, stringsAsFactors = FALSE))
    
    # create network with the current layer, making sure all the nodes have the same ID in every layer
    CurrentNetwork <- graph_from_data_frame(d = Layer, vertices = AllNodesToTake, directed = FALSE)
    
    # add subnetwork as a layer into the multiplex network
    Multiplex[[ length(Multiplex) + 1 ]] <- CurrentNetwork
  }  
  
  return (Multiplex)
}


# definition of the function that generates the merged network 
# INPUTS: Files - list of file names that contain the biological networks data
# OUTPUT: Merged network (igraph object)
GenerateMergedNetwork <- function(Files) {
  # declare empty data frame for the edges of the merged 
  Merged <- data.frame(V1 = character(), V2 = character(), Layer = character())

  # loop through all the layers to get the corresponding subnetwork from each of them
  for (LayerFile in Files) {
    # load layer
    Layer <- data.frame(read.table(paste0(MyLayersPath, LayerFile), header = FALSE, stringsAsFactors = FALSE), Layer = LayerFile)
    Merged <- rbind(Merged, Layer)
  }  

  # create network with the current layer, making sure all the nodes have the same ID in every layer
  MergedNetwork <- graph_from_data_frame(d = Merged, vertices = names(V(Multiplex[[1]])), directed = FALSE)
  
  return (MergedNetwork)
}

# definition of the function that gets the sorted list of nodes present in at least one of the layers
# INPUTS: Files - list of file names that contain the biological networks data
# OUTPUT: Sorted list of nodes
GetListOfAllNodesPresentInLayers <- function(Files) {
  # declare empty variable to store the nodes
  AllNodes <- NULL
  
  # loop through all the layers to get the corresponding subnetwork from each of them
  for (LayerFile in Files) {
    # load layer
    Layer <- data.frame(read.table(paste0(MyLayersPath, LayerFile), header = FALSE, stringsAsFactors = FALSE))
    
    AllNodes <- c(AllNodes, unique(union(Layer[, 1], Layer[, 2])))
  }
  
  # remove duplicates
  AllNodes <- unique(AllNodes)
  
  # order alphabetically
  AllNodes <- sort(AllNodes)
  
  return (AllNodes)
}



# definition of the function that calculates the individual node score for a list of genes
# INPUTS: ListOfGenes - list of gene names to calculate the node score
#         Measure - either consider 'PValue' or 'FDR' for the node score calculus
# OUTPUT: nodes scores
GetNodesScoresOfListOfGenes <- function(ListOfGenes, Measure) {
  
  ListOfGenes <- unique(ListOfGenes)
  
  # calculate nodes scores for the genes with formula: z=inverse CDF(1-p) if gene has an expression value and z=-Inf otherwise
  NodesScores <- as.numeric(sapply(
    ListOfGenes, 
    function(X) {
      ifelse (X %in% as.character(DE_results$gene), qnorm( 1 - DE_results[[Measure]][DE_results$gene == X] ), (-Inf))
    }
  ))
  
  # LEAVE THE +INF AND -INF UNTIL AFTER THE NORMALIZATION!!! THESE SHOULD ALWAYS BE 1 AND 0, RESPECTIVELY

  # normalizes nodes scores to be in the range [0-1], ignoring +Inf and -Inf values
  NodesScores <- 
    (NodesScores - min(NodesScores[which(is.finite(NodesScores))])) / 
    (max(NodesScores[which(is.finite(NodesScores))]) - min(NodesScores[which(is.finite(NodesScores))]))
  
  # replace +Inf with 1
  NodesScores[which(is.infinite(NodesScores) & NodesScores > 1)] <- 1
  
  # replace -Inf with 0
  NodesScores[which(is.infinite(NodesScores) & NodesScores < 1)] <- 0
  
  return(NodesScores)
}



# definition of the function to generate the initial population
# INPUTS: PopSize - population size, i.e. number of individuals to create
#         Multiplex - multiplex network (it is the search space)
# OUTPUT: Population of size PopSize with integer codification
GenerateInitialPopulation <- function (PopSize, Multiplex) {
  MyPopulation <- vector("list", PopSize) # initialize an empty vector of the population size
  for (i in 1:PopSize) { # loop from 1 to the total population size
    # randomly determine the size of the individual, i.e., the number of nodes it will contain 
    SizeOfIndividual <- sample(x = MinNumberOfNodesPerIndividual:MaxNumberOfNodesPerIndividual, size = 1)
    
    Root <- PickRandomRoot(MyMultiplex = Multiplex) 
    
    # use DFS (Depth First Search) from the chosen root, until the size of the individual is reached
    # NOTE. Only DFS will be used in the initial population because it allows to go further from the root
    
    # makes a random depth first search, i.e., the branches are shuffled before visiting them, in order
    # to generate different individuals each time
    DFS <-
      DFS_iterative_Multiplex(
        MyMultiplex = Multiplex,
        Root = Root,
        SizeOfIndividual = SizeOfIndividual
      )
    
    # add individual to the population
    MyPopulation[i] <- list(DFS)
  }
  return(MyPopulation) 
}


# definition of the function to perform a RANDOM depth first search of a specified length on a MULTIPLEX network.
# This means that the branch to be visted is randomnly picked from all the available ones
# INPUTS:
#         Multiplex - the multiplex network
#         Root - name of the node that will be used as root for the search 
#         SizeOfIndividual - desired length of the search
# OUTPUT: Individual with the desired length
DFS_iterative_Multiplex <- function (MyMultiplex, Root, SizeOfIndividual) {
  KeepLooking <- TRUE # flag to control the execution of the current function
  Attempts <- 0
  
  while(KeepLooking == TRUE) {
    Discovered <- NULL  # initialize empty list of already visited nodes 
    Result <- NULL  # initialize empty list of the resulting individual
    Stack <- Root  # initialize the stack of the pending nodes to visit with the root of the tree
    CurrentLayer <- NULL
    
    # loop while there are pending elements to visit and the desired size of the individual hasn't been met
    while(length(Stack) > 0 && length(Result) < SizeOfIndividual) {
      # if this is the first node to be picked (CurrentLayer == NULL) or if with 50% of probability we change of layer
      if (is.null(CurrentLayer) || runif(1, 0, 1) <= 0.5) {
        AvailableLayers <- 1:length(MyMultiplex) # build a list with the layers' numbers
        
        if (length(AvailableLayers) > 1) {
          AvailableLayers <- AvailableLayers[!AvailableLayers %in% CurrentLayer] # remove current layer
        }
        
        CurrentLayer <- ifelse(length(AvailableLayers > 0), sample(AvailableLayers, 1), CurrentLayer) # choose a new layer
      }
      
      Node <- Stack[length(Stack)]  # take a node from the stack to visit it
      Stack <- Stack[-length(Stack)] # remove node from stack
      
      # verify if the node hasn't been visited yet
      if (!(Node %in% Discovered)) {
        Discovered <- c(Discovered, Node) # add node to the list of visited 
        
        # --- NOTE. Any node can be disconnected in one layer or more, therefore, if we find a node with no neighbors
        # --------  we will change of layer, until finding neighbors or visiting all layers
        KeepGoing <- TRUE
        VisitedLayers <- NULL
        AvailableLayers <- 1:length(MyMultiplex)
        
        while (KeepGoing == TRUE) {
          # get the neighbors of the current node in current layer
          MyNeighbors <- names(neighbors(MyMultiplex[[CurrentLayer]], Node)) 
          
          # remove already discovered nodes
          MyNeighbors <- MyNeighbors[!MyNeighbors %in% Discovered]
          
          if (length(MyNeighbors) == 0) { 
            # add the current layer to the list of visited layers
            VisitedLayers <- c(VisitedLayers, CurrentLayer)
            
            # update the list of available layers to visit
            AvailableLayers <- AvailableLayers[!AvailableLayers %in% VisitedLayers]
            
            # verify if there are still available layers to visit
            if (length(AvailableLayers) > 0) {
              CurrentLayer <- AvailableLayers[sample(length(AvailableLayers), 1)]
            } else {
              KeepGoing <- FALSE # no neighbors could be found for this node
            }
          } else {
            KeepGoing <- FALSE
          }
        }
        
        MyNeighbors <- MyNeighbors[sample(length(MyNeighbors))] # shuffle neighbors to RANDOMLY pick a branch
        Stack <- c(Stack, MyNeighbors) # add shuffled neighbors to the stack
        Result <- c(Result, Node) # add node to the individual
      } 
    }
    
    # security check to make sure that the generated individual has the desired size, otherwise a new attempt will be done
    # with another root
    if (length(Result) != SizeOfIndividual) {
      # if no "good" root was found in the maximum number of attempts, generate a random individual with the original multiplex
      if (Attempts > MaxNumberOfAttempts) {
        MyMultiplex <- Multiplex # use the big original multiplex
        Root <- PickRandomRoot(MyMultiplex = MyMultiplex) # pick a random root with the big original multiplex
      } else { # if we still have attempts left
        Attempts <- Attempts + 1 # increment number of attempts
        Root <- PickRandomRoot(MyMultiplex = MyMultiplex) # pick a root with a potentially reduced version of the multiplex
      }
    } else { # if a good root was found and the individual has the desired size
      KeepLooking <- FALSE # deactivate flag to stop the search
    }
  }
  
  ##### the following conversion has to be done due to the fact that when a subnetwork is created,
  ##### the nodes get new IDs, according to the size of the new subnetwork. But the individuals
  ##### should always have IDs with respect to the global network, therefore the IDs need to be
  ##### "re-calculated"
  
  # get the names of the nodes in the local network with the local IDs
  Nodes <- names(V(MyMultiplex[[1]])[Result])
  
  # get the global IDs of the corresponding nodes
  NodesIDs <- GetNetworkIDofListOfNodes(Nodes, Multiplex[[1]])
  
  return (NodesIDs)
}


# Definition of the function to pick a random node to be the root of a search. It prefers DE genes over non-DE
# INPUT:  MyMultiplex - network containing all the available nodes 
# OUTPUT: Root node
PickRandomRoot <- function (MyMultiplex) {
  SearchSpaceGenes <- names( V(MyMultiplex[[1]]) )
  
  # get the list of all the DE genes (from the whole DE analysis test)
  DEGenesNames <- as.character(DifferentiallyExpressedGenes$gene)
  
  # get the list of DE genes in the search space 
  DEGenesInSearchSpace <- as.character(DifferentiallyExpressedGenes[DEGenesNames %in% SearchSpaceGenes, "gene"])
  
  # verify if there is at least one DE gene
  if (length(DEGenesInSearchSpace) > 0) {
    # randomly pick a DE gene from the list, to be the root of the tree (individual)
    Root <- DEGenesInSearchSpace[sample(1:length(DEGenesInSearchSpace), 1)]
  } else {
    # randomly pick any gene from the list of available genes
    Root <- SearchSpaceGenes[sample(1:length(SearchSpaceGenes), 1)]
  }
  
  return(Root)
}


# definition of the function to evaluate a whole population 
# INPUTS: MyPopulation - population to evaluate (only the codes of the individuals)
#         Multiplex - network to which the population belongs to
#         GenesWithNodesScores - nodes scores of the genes
# OUTPUT: Evaluation of the individual (fitness)
EvaluatePopulation <- function (MyPopulation, Multiplex, GenesWithNodesScores) {
  
  FitnessData_alternative <- 
       do.call(
            rbind, 
            lapply(
                 seq_len(length(MyPopulation)), # from 1 to PopSize
                 function(i) { 
                      EvaluateIndividualAndReturnAllFitnessData(
                           Individual = MyPopulation[[i]], 
                           MultiplexNetwork = Multiplex, 
                           VerticesWithNodesScores = GenesWithNodesScores
                      )
                 } 
            ))

  return (FitnessData)
}



# definition of the function to evaluate an individual 
# INPUTS: Individual - individual to evaluate
#         MultiplexNetwork - network to which the individual belongs to
#         VerticesWithNodesScores - nodes scores of the genes
# OUTPUT: Evaluation of the individual (fitness)
EvaluateIndividualAndReturnAllFitnessData <- function (Individual, MultiplexNetwork, VerticesWithNodesScores) {
  
  # verifies that there is at least one node present in the individual
  if ( length(Individual) > 0 ) {
    
    # gets the sum of the nodes scores of the subnetwork
    # I wil use fixed the layer number (1), because the nodes are ordered equally in every layer
    SumNodesScores <- 
      sum(
        VerticesWithNodesScores[VerticesWithNodesScores$gene %in% V(MultiplexNetwork[[1]])$name[Individual], 
                            "nodescore"]
      )
    
    AverageNodesScore <- SumNodesScores / length(Individual)
    
    SumDensityAllLayers <- 0

    # loop through all the layers in the multiplex
    for (layer in 1:length(MultiplexNetwork)) {
      # get the subnetwork corresponding to the individual, in the current layer
      Subnetwork <- induced_subgraph(MultiplexNetwork[[layer]], Individual)
      
      # calculate the density of the subnetwork corresponding to the individual, in the current layer
      SubnetworkDensity <- graph.density(Subnetwork)
      if (is.nan(SubnetworkDensity)) {
           SubnetworkDensity <- 0
      }

      # add the normalized subnetwork's density with respect to the density of the current layer
      SumDensityAllLayers <- SumDensityAllLayers + (SubnetworkDensity / DensityPerLayerMultiplex[layer])
    }
    
    Res <- data.frame(AverageNodesScore = AverageNodesScore, Density = SumDensityAllLayers)
  } else {
    Res <- data.frame(AverageNodesScore = 0, Density = 0)
  }
  return (Res)
}


############## FUNCTION FROM PACKAGE nsga2R ##################
# definition of the function to sort by non-domination a whole population
# INPUTS: inputData - data frame (contains the individuals' objective functions)
# OUTPUT: Rank of the population
fastNonDominatedSorting <-   function(inputData) {
  popSize = nrow(inputData)
  idxDominators = vector("list", popSize)
  idxDominatees = vector("list", popSize)
  for (i in 1:(popSize-1)) {
    for (j in (i+1):popSize) {
        xi = inputData[i, ]
        xj = inputData[j, ]
        if (all(xi >= xj) && any(xi > xj)) {  ## i dominates j
          idxDominators[[j]] = c(idxDominators[[j]], i)
          idxDominatees[[i]] = c(idxDominatees[[i]], j) 
        } else if (all(xj >= xi) && any(xj > xi)) {  ## j dominates i
          idxDominators[[i]] = c(idxDominators[[i]], j)
          idxDominatees[[j]] = c(idxDominatees[[j]], i) 
        }
    }
  }
  noDominators <- lapply(idxDominators,length);
  rnkList <- list();
  rnkList <- c(rnkList,list(which(noDominators==0)));
  solAssigned <- c();
  solAssigned <- c(solAssigned,length(which(noDominators==0)));
  while (sum(solAssigned) < popSize) {
    Q <- c();
    noSolInCurrFrnt <- solAssigned[length(solAssigned)];
    for (i in 1:noSolInCurrFrnt) {
      solIdx <- rnkList[[length(rnkList)]][i];
      hisDominatees <- idxDominatees[[solIdx]]; # A vector
      for (i in hisDominatees) {
        noDominators[[i]] <- noDominators[[i]] - 1;
        if (noDominators[[i]] == 0) {
          Q <- c(Q, i);
        }
      }
    }
    rnkList <- c(rnkList,list(sort(Q))); # sort Q before concatenating
    solAssigned <- c(solAssigned,length(Q));
  }
  
  return(rnkList)
}



# function to check for dominance of Ind1 over Ind2
# INPUTS: Ind1 - individual 1
#         Ind2 - individual 2
# OUTPUT: Dominance flag. 0: non-dominance, 1: Ind1 dominates Ind2, 2: Ind2 dominates Ind1
CheckDominance <- function(Ind1, Ind2) {
  
  # save all the fitness functions (objectives) values of the individuals
  FitnessFunctionsInd1 <- Ind1[, (colnames(Ind1) %in% MyObjectiveNames)]
  FitnessFunctionsInd2 <- Ind2[, (colnames(Ind2) %in% MyObjectiveNames)]
  
  #### In order for Ind1 to dominate Ind2, Ind1 has to have fitness functions values higher or equal
  #### than Ind2 for all the objectives, and it has to be better in at least one, i.e., have a higher
  #### value in at least one fitness function
  
  if (all(FitnessFunctionsInd1 >= FitnessFunctionsInd2) && any(FitnessFunctionsInd1 > FitnessFunctionsInd2)) {
    Dominance <- 1 # individual 1 dominates individual 2
  } else if(all(FitnessFunctionsInd2 >= FitnessFunctionsInd1) && any(FitnessFunctionsInd2 > FitnessFunctionsInd1)) {
    Dominance <- 2 # individual 2 dominates individual 1
  } else {
    Dominance <- 0 # no dominance between individuals 1 and 2
  }
  
  return(Dominance)
}


# definition of the function to perform a tournament selection
# INPUTS: TournamentSize - size of the tournament
#         PopulationForTournament - population to do the tournament with
# OUTPUT: Individual that won the tournament 
TournamentSelection <- function (TournamentSize, PopulationForTournament) {
  # randomly choose as many individuals as the tournament size indicates
  ids <- sample(1:nrow(PopulationForTournament), TournamentSize, replace=FALSE)
  
  # verify all both individuals are in the same Pareto front
  if (length(unique(PopulationForTournament$Rank[ids])) == 1) {
    # verify if these individuals have information about crowding distance (they won't if generation == 1)
    if ("CrowdingDistance" %in% colnames(PopulationForTournament)) {
      # get id of the individual with the highest crowding distance
      winner <- which(PopulationForTournament$CrowdingDistance[ids] == max(PopulationForTournament$CrowdingDistance[ids]))
    } else {
      # if there is no crowding distance and the individuals have the same rank, they are all winners
      winner <- c(1:TournamentSize) 
    }
  } else { # if the individuals are in different Pareto fronts, the winner of the tournament is the one with lowest value
    winner <- which(PopulationForTournament$Rank[ids] == min(PopulationForTournament$Rank[ids]))
  }
  
  # verify if there was a tie
  if (length(winner) > 1) {
    
    # pick a single winner randomly
    winner <- winner[ sample(1:length(winner), 1) ]
  }
  
  TournamentWinner <- PopulationForTournament[ids[winner],]
  
  return(TournamentWinner)
}


# definition of the function to get the union of neighbors from all the layers, of a list of nodes
# INPUTS: NodeList - list of nodes to look for its neighbors 
#         Multiplex - multiplex network
# OUTPUT: List of neighbors
GetNeighborsOfNodeList <- function(NodeList, Multiplex) {
  Neighbors <- NULL
  
  # loop through all the layer of the multiplex
  for (i in 1:length(Multiplex)) {
    # add to the list the neighbors of the node list in the current layer
    Neighbors <- c(Neighbors, unlist( lapply( NodeList, function(X) { neighbors(Multiplex[[i]], X) } ) ) )
  }
  
  # sort result and delete duplicates
  Neighbors <- sort(unique(Neighbors))
  
  return(Neighbors)
}


# definition of the function to perform crossover
# INPUTS: Parent1 - first parent to cross
#         Parent2 - second parent to cross
# OUTPUT: Two children (a single connected component per child) 
Crossover <- function (Parent1, Parent2) {
  Children <- NULL # intialize empty list
  
  p <- runif(1, 0, 1) # generate a random number between 0 and 1 
  
  if (p <= CrossoverRate) { # check if crossover is to be performed
    # join both parents (if they are compatible, this will be a single connected component)
    NodesInParents <- union( unlist(Parent1$Individual), unlist(Parent2$Individual)  )
    
    # first, generate the multiplex of the joint parents
    Multiplex_JointParents <- FilterMultiplex(Multiplex = Multiplex, NodesToKeep = NodesInParents)
    
    # if the size of the union of the parents is equal to the minimum number of valid nodes, then just copy the parents
    if (length(NodesInParents) == MinNumberOfNodesPerIndividual) {
      Children <- rbind(Parent1$Individual, Parent2$Individual)
    } else {
      # loop to generate two children
      for (k in 1:2) {
        ##### To choose the size of the child, we have to take into account, not only the max and min allowed,
        ##### but also the size of the joined parents, becuase if this is shorter than the maximum and 
        ##### we just randomly pick a size >min and <max, it can be out of bounds with the available nodes
        
        # randomly pick a size for the child, considering the maximum size allowed and the maximun number
        #  of nodes available
        SizeOfIndividual <- ifelse(
          length(NodesInParents) >= MaxNumberOfNodesPerIndividual,
          sample(MinNumberOfNodesPerIndividual:MaxNumberOfNodesPerIndividual, 1),
          sample(MinNumberOfNodesPerIndividual:length(NodesInParents), 1)
        )
        
        # randomly pick the root of the child
        Root <- PickRandomRoot(MyMultiplex = Multiplex_JointParents) 
        
        # randomly pick the search method: DFS (Depth First Search) or BFS (Breadth First Search)
        SearchMethod <- sample(c("DFS", "BFS"), 1)
        
        # if the selected method is DFS
        if (SearchMethod == "DFS") {
          
          # perform depth first search in the filtered multiplex
          DFS <- DFS_iterative_Multiplex(
            MyMultiplex = Multiplex_JointParents,    # use the filtered version of the multiplex for the crossover
            Root = Root,
            SizeOfIndividual = SizeOfIndividual
          )
          
          # add child to the population
          Children[k] <- list(DFS)
        } else if (SearchMethod == "BFS") {
          # perform depth first search
          BFS <- BFS_iterative_Multiplex(
            MyMultiplex = Multiplex_JointParents,    # use the filtered version of the multiplex for the crossover
            Root = Root,
            SizeOfIndividual = SizeOfIndividual
          )
          
          # add child to the population
          Children[k] <- list(BFS)
        }
      }
    }
  } else {
    # if no crossover was performed, make a copy of the parents 
    Children <- rbind(Parent1$Individual, Parent2$Individual)
  }  
  return(Children)
}


# definition of the function to perform a RANDOM breadth first search of a specified length on a MULTIPLEX network.
# This means that the nodes to be visted are always randomnly picked from all the available ones
# INPUTS:
#         Multiplex - the multiplex network
#         Root - name of the node that will be used as root for the search 
#         SizeOfIndividual - desired length of the search
# OUTPUT: Individual with the desired length
BFS_iterative_Multiplex <- function (MyMultiplex, Root, SizeOfIndividual) {
  KeepLooking <- TRUE # flag to control the execution of the current function
  Attempts <- 0
  
  while(KeepLooking == TRUE) {
    Discovered <- NULL  # initialize empty list of already visited nodes 
    Result <- NULL  # initialize empty list of the resulting individual
    Queue <- Root  # initialize the queue of the pending nodes to visit with the root of the tree
    CurrentLayer <- NULL
    
    # loop while there are pending elements to visit and the desired size of the individual hasn't been met
    while(length(Queue) > 0 && length(Result) < SizeOfIndividual) {
      # if this is the first node to be picked (CurrentLayer == NULL) or if with 50% of probability we change of layer
      if (is.null(CurrentLayer) || runif(1, 0, 1) <= 0.5) {
        AvailableLayers <- 1:length(MyMultiplex) # build a list with the layers' numbers
        
        if (length(AvailableLayers) > 1) {
          AvailableLayers <- AvailableLayers[!AvailableLayers %in% CurrentLayer] # remove current layer
        }
        
        CurrentLayer <- sample(AvailableLayers, 1) # choose a new layer
      }
      
      Node <- Queue[1]  # take a node from the queue to visit it
      Queue <- Queue[-1] # remove node from queue
      
      # verify if the node hasn't been visited yet
      if (!(Node %in% Discovered)) {
        Discovered <- c(Discovered, Node) # add node to the list of visited 
        
        # --- NOTE. Any node can be disconnected in one layer or more, therefore, if we find a node with no neighbors
        # --------  we will change of layer, until finding neighbors or visiting all layers
        KeepGoing <- TRUE
        VisitedLayers <- NULL
        AvailableLayers <- 1:length(MyMultiplex)
        
        while (KeepGoing == TRUE) {
          # get the neighbors of the current node in current layer
          MyNeighbors <- names(neighbors(MyMultiplex[[CurrentLayer]], Node)) 
          
          if (length(MyNeighbors) == 0) { 
            # add the current layer to the list of visited layers
            VisitedLayers <- c(VisitedLayers, CurrentLayer)
            
            # update the list of available layers to visit
            AvailableLayers <- AvailableLayers[!AvailableLayers %in% VisitedLayers]
            
            # verify if there are still available layers to visit
            if (length(AvailableLayers) > 0) {
              CurrentLayer <- AvailableLayers[sample(length(AvailableLayers), 1)]
            } else {
              KeepGoing <- FALSE # no neighbors could be found for this node
            }
          } else {
            KeepGoing <- FALSE
          }
        }
        
        MyNeighbors <- MyNeighbors[sample(length(MyNeighbors))] # shuffle neighbors to RANDOMLY pick them
        Queue <- c(Queue, MyNeighbors) # add shuffled neighbors to the stack
        Result <- c(Result, Node) # add node to the individual
      } 
    }
    
    # security check to make sure that the generated individual has the desired size, otherwise a new attempt will be done
    # with another root
    if (length(Result) != SizeOfIndividual) {
      # if no "good" root was found in the maximum number of attempts, generate a random individual with the original multiplex
      if (Attempts > MaxNumberOfAttempts) {
        MyMultiplex <- Multiplex # use the big original multiplex
        Root <- PickRandomRoot(MyMultiplex = MyMultiplex) # pick a random root with the big original multiplex
      } else { # if we still have attempts left
        Attempts <- Attempts + 1 # increment number of attempts
        Root <- PickRandomRoot(MyMultiplex = MyMultiplex) # pick a root with a potentially reduced version of the multiplex
      }
    } else { # if a good root was found and the individual has the desired size
      KeepLooking <- FALSE # deactivate flag to stop the search
    }
  }
  
  ##### the following conversion has to be done due to the fact that when a subnetwork is created,
  ##### the nodes get new IDs, according to the size of the new subnetwork. But the individuals
  ##### should always have IDs with respect to the global network, therefore the IDs need to be
  ##### "re-calculated"
  
  # get the names of the nodes in the local network with the local IDs
  Nodes <- names(V(MyMultiplex[[1]])[Result])
  
  # get the global IDs of the corresponding nodes
  NodesIDs <- GetNetworkIDofListOfNodes(Nodes, Multiplex[[1]])
  
  return (NodesIDs)
}



# definition of the function that filters the multiplex network, to keep only the specified list of nodes
# INPUTS: Multiplex - multiplex network to filter
#         ListOfNodes - list of nodes to keep 
# OUTPUT: Filtered version of the multiplex 
FilterMultiplex <- function(Multiplex, NodesToKeep) {
  # declare empty list to store the multiplex
  FilteredMultiplex <- list()
  
  # loop through all the layers to get the corresponding subnetwork from each of them
  for (i in 1:length(Multiplex)) {
    
    # create network with the current layer
    CurrentNetwork <- induced_subgraph(Multiplex[[i]], NodesToKeep)
    
    # add subnetwork as a layer into the multiplex network
    FilteredMultiplex[[ length(FilteredMultiplex) + 1 ]] <- CurrentNetwork
  }  
  
  return (FilteredMultiplex)
}


# definition of the function to get the list of IDs of a set of nodes
# INPUTS: ListOfNodes - list of nodes' names 
#         GlobalNetwork - the general network where we want to look for our list of nodes
# OUTPUT: List of IDs
GetNetworkIDofListOfNodes <- function(ListOfNodes, GlobalNetwork) {
  return (which(names(V(GlobalNetwork)) %in% ListOfNodes))
}


# definition of the function to perform mutation
# INPUTS: Individuals - the individuals to perform mutation with
#         Multiplex - network to where the individuals belong to 
# OUTPUT: Two mutated individuals
Mutation <- function (Individuals, Multiplex) {
  Mutants <- NULL
  
  for (i in 1:length(Individuals)) { # loop through all the individuals to be mutated
    
    p <- runif(1, 0, 1) # generate a random number between 0 and 1 
    
    if (p <= MutationRate) { # check if mutation is to be performed

      # create the corresponding subnetork in the merged network
      IndividualToMutate_Network <- induced_subgraph(Merged, Individuals[[i]])
      
      # obtain the degree per node 
      IndividualToMutate_NodesDegrees <- degree(IndividualToMutate_Network)
      
      # remove all the DE genes from the list
      IndividualToMutate_NodesDegrees <- 
        IndividualToMutate_NodesDegrees[!names(IndividualToMutate_NodesDegrees) %in% as.character(DifferentiallyExpressedGenes$gene)]
      
      # get the list of nodes with the minimum degree, i.e., the peripheral nodes
      PeripheralNodes <- which.min(IndividualToMutate_NodesDegrees)

      # get the list of node names that can be mutated
      PotentialNodesToMutate <- names(PeripheralNodes)
      
      # generates a vector with the chromosomes to mutate. Each chromosome has a percentage of getting mutated of MutationRate
      NodesToMutate <- 
        apply(
          array(rep(0, length(PotentialNodesToMutate))), 
          1, 
          function(X) { ifelse(runif(1, 0, 1) <= MutationRate, 1, 0) } 
        ) 
      
      NodesToMutate <- which(NodesToMutate == 1)  # gets the list of chromosomes to mutate
      
      # verify if at least one of the nodes will be mutated
      if (length(NodesToMutate) > 0) {
        # loop through the nodes to remove
        for (j in NodesToMutate) {
          MutatedNetwork <- delete_vertices(IndividualToMutate_Network, PotentialNodesToMutate[j])
          
          # obtain neighbors of nodes 
          Neighbors_OriginalIndividual <- 
            names( unlist( lapply( names(V(IndividualToMutate_Network)), function(X) { neighbors(Merged, X) } ) ) )
          
          Neighbors_MutatedIndividual <- 
            names( unlist( lapply( names(V(MutatedNetwork)), function(X) { neighbors(Merged, X) } ) ) )
          
          # delete from the list all the nodes that originally belonged to the individual
          Neighbors_OriginalIndividual <- 
            Neighbors_OriginalIndividual[!Neighbors_OriginalIndividual %in% names(V(IndividualToMutate_Network))]
          Neighbors_MutatedIndividual <- 
            Neighbors_MutatedIndividual[!Neighbors_MutatedIndividual %in% names(V(IndividualToMutate_Network))]
          
          # verify if the mutated network is connected
          if( is.connected(MutatedNetwork) ) {
            # save the changes in the individual
            Individuals[[i]] <- GetNetworkIDofListOfNodes(names(V(MutatedNetwork)), Multiplex[[1]])
            
            IndividualToMutate_Network <- MutatedNetwork
            
            AvailableNeighbors <- Neighbors_MutatedIndividual
          } else {
            AvailableNeighbors <- Neighbors_OriginalIndividual
          }
          
          # verify that there is at least one available neighbor to be added to the subnetwork
          if(length(AvailableNeighbors) > 0) {
            
            # get the list of DE genes that are neighbors of the current individual
            DE_Neighbors <- AvailableNeighbors[AvailableNeighbors %in% as.character(DifferentiallyExpressedGenes$gene)]
            
            # verify if there is at least one DE gene in the list of neighbors
            if (length(DE_Neighbors) > 0) {
              AvailableNeighbors <- DE_Neighbors # leave only the DE gens as available neighbors
            } 
            
            # get the number of times each node was a neighbor 
            Incidences <- table(as.factor(AvailableNeighbors))
            
            # get list of nodes that have the higheest incidence
            MaxIncidences <- which.max(Incidences)
            
            # keep as available neighbors only those with a maximum incidence
            AvailableNeighbors <- names(MaxIncidences)

            # choose a node, randomly
            RandomNode <- AvailableNeighbors[sample(1:length(AvailableNeighbors), 1)]
            
            # get the global ID of the chosen node
            NewNodeID <- GetNetworkIDofListOfNodes(RandomNode, Multiplex[[1]])
            
            # add new node to the individual
            Individuals[[i]] <- c(Individuals[[i]], NewNodeID)
          } else {
            print("Attempt to add a new neighbor, after removing one FAILED")
          }
        }
      } else { # if no nodes were removed, add a new neighbor
        
        # obtain neighbors of nodes 
        AvailableNeighbors <- 
          names( unlist( lapply( names(V(IndividualToMutate_Network)), function(X) { neighbors(Merged, X) } ) ) )
        
        # delete from the list all the nodes that originally belonged to the individual
        AvailableNeighbors <- AvailableNeighbors[!AvailableNeighbors %in% names(V(IndividualToMutate_Network))]
        
        # verify that there is at least one available neighbor to be added to the subnetwork
        if(length(AvailableNeighbors) > 0) {
          
          # get the list of DE genes that are neighbors of the current individual
          DE_Neighbors <- AvailableNeighbors[AvailableNeighbors %in% as.character(DifferentiallyExpressedGenes$gene)]
          
          # verify if there is at least one DE gene in the list of neighbors
          if (length(DE_Neighbors) > 0) {
            AvailableNeighbors <- DE_Neighbors
          } 
          
          # get the number incidences for the nodes in the neighborhood 
          Incidences <- table(as.factor(AvailableNeighbors))
          
          # get list of nodes that have a higher incidence
          MaxIncidences <- which.max(Incidences)
          
          # keep as available neighbors only those with a maximum incidence
          AvailableNeighbors <- names(MaxIncidences)
          
          # choose a node, randomly
          RandomNode <- AvailableNeighbors[sample(1:length(AvailableNeighbors), 1)]
          
          # get the global ID of the chosen node
          NewNodeID <- GetNetworkIDofListOfNodes(RandomNode, Multiplex[[1]])
          
          # add new node to the individual
          Individuals[[i]] <- c(Individuals[[i]], NewNodeID)
        } else {
          print("No nodes to be mutated. Attempt to add a new neighbor FAILED")
        }
      }
    } 
    # save the individual in the mutants' list
    Mutants[i] <- Individuals[i]
  }
  return(Mutants)
}


# definition of the function to perform replacement
# INPUTS: Parents - individuals from the previous generation 
#         Children - individuals from the new generation
# OUTPUT: The selected invididuals to keep for the next generation
Replacement <- function (Parents, Children) {
  # combine the new population (offspring) with old population (parents)
  CombinedPopulation <- rbind(Parents, Children)
  
  ## CHECK FOR DUPLICATED INDIVIDUALS HERE AND REPLACE THEM FOR NEW GENERATED ONES INSTEAD
  CombinedPopulation <- ReplaceDuplicatedIndividualsWithRandomOnes(CombinedPopulation)
  
  # order combined population by rank 
  CombinedPopulation <- CombinedPopulation[with(CombinedPopulation, order(Rank)), ]
  
  # get last rank which will be in the replacement population
  LastRank <- CombinedPopulation$Rank[PopSize]
  
  # get last rank individuals
  LastRankInds <- CombinedPopulation[CombinedPopulation$Rank == LastRank, ]

  # order by crowding distance
  LastRankInds <- LastRankInds[with(LastRankInds, order(CrowdingDistance, decreasing = TRUE)), ]
  
  # select new population for replacement
  NewPopulationForReplacement <- CombinedPopulation[CombinedPopulation$Rank < LastRank, ]
  NewPopulationForReplacement <- 
       rbind( NewPopulationForReplacement, LastRankInds[1:(PopSize - nrow(NewPopulationForReplacement)), ] )

  return(NewPopulationForReplacement)
}


# definition of the function that replaces all duplicated individuals and those above the permitted threshold of similarity 
# with newly generated random individuals
# INPUT:  CombinedPopulation - Population of size 2*N that corresponds to the union of parents and children
# OUTPUT: A garanteed diverse population of size 2*N
ReplaceDuplicatedIndividualsWithRandomOnes <- function(CombinedPopulation) {
  DiversePopulation <- CombinedPopulation
  
  #######################################################################
  ############## ---- BEGINS OPTIMIZATION OF CODE ---- ##################
  #######################################################################
  # create all the combinations of 2 numbers with the ids of individuals
  MyIndexes <- t(combn(1:nrow(DiversePopulation), 2))
  
  # create data frame to save the similarities
  Similarities <- data.frame(Index1 = MyIndexes[, 1], Index2 = MyIndexes[, 2], JS = 0)
  
  # calculate the Jaccard similarity index
  Similarities$JS <- 
    apply(
      MyIndexes, 
      1, 
      function(X) { 
        Ind1 <- as.character(unlist(DiversePopulation$Individual[X[1]]))
        Ind2 <- as.character(unlist(DiversePopulation$Individual[X[2]]))
        
        # Jaccard similarity = (intersection of A and B) / (union of A and B)
        JS <- (length(intersect(Ind1, Ind2))) / (length(union(Ind1, Ind2))) * 100
        
        return(as.double(JS))
      } 
    )
  
  # keep only the rows of "duplicated" individuals
  Similarities <- Similarities[Similarities$JS >= JaccardSimilarityThreshold, ]
  
  # verify if there is at least a pair of duplicated individuals
  if (nrow(Similarities) > 0) {
    IndividualsToRemove <- NULL
    
    # non-dominated sorting and crowding distance calculus
    PopulationWithCrowdingDistance <- 
         SortByNonDominationAndObtainCrowdingDistance(PopulationToSort = DiversePopulation)
    
    # loop through all the "duplicated" pairs of individuals
    i <- 1
    while(i < nrow(Similarities)) {
      # ID of individuals
      Ind1_ID <- Similarities[i, 1]
      Ind2_ID <- Similarities[i, 2]
      
      if ( Similarities$JS[i] < 100 &
           all(PopulationWithCrowdingDistance$Rank[c(Ind1_ID, Ind2_ID)] == 1) & 
           all(is.infinite(unique(PopulationWithCrowdingDistance$CrowdingDistance[c(Ind1_ID, Ind2_ID), ])[1]))
     ) {
        print("Very similar individuals in first Pareto front. Keeping both of them. Inf crowding distance.")
        
        print(PopulationWithCrowdingDistance[c(Ind1_ID, Ind2_ID), ])
      } else {
        # tournament between the two individuals
        IndividualToKeep <- 
          TournamentSelection(
               TournamentSize = 2, 
               PopulationForTournament = PopulationWithCrowdingDistance[c(Ind1_ID, Ind2_ID), ]
               )
        
        # get ID of the individual to remove
        IndividualsToRemove <- 
             c(IndividualsToRemove, 
               ifelse(rownames(IndividualToKeep) == row.names(PopulationWithCrowdingDistance)[Ind1_ID],
                      Ind2_ID, Ind1_ID))
        
        # get all future incidences of the individual to remove
        References <- 
             which(Similarities == IndividualsToRemove[length(IndividualsToRemove)], arr.ind = TRUE)[, "row"]
        
        # remove all future references to the new individual to remove
        Similarities <- 
             Similarities[!row.names(Similarities) %in% row.names(Similarities)[References[References > i]], ]
      }
      
      i <- i + 1
    }
    
    # remove the corresponding individuals
    DiversePopulation <- 
         DiversePopulation[!row.names(DiversePopulation) %in% row.names(DiversePopulation)[IndividualsToRemove], ]
  }
  
  # generate as many new individuals as duplicated ones   
  MyNewIndividuals <- 
    GenerateInitialPopulation(PopSize = (nrow(CombinedPopulation) - nrow(DiversePopulation)), Multiplex = Multiplex)
  
  # Evaluate individuals
  FitnessData <- EvaluatePopulation(MyNewIndividuals, Multiplex, GenesWithNodesScores)
  
  DiversePopulation <- 
       rbind(DiversePopulation, 
             data.frame("Individual" = I(MyNewIndividuals), FitnessData, Rank = 0, CrowdingDistance = 0)
             )

  DiversePopulation <- SortByNonDominationAndObtainCrowdingDistance(PopulationToSort = DiversePopulation)

  return(DiversePopulation)  
}


# definition of the function that performs the fast non dominated sorting and the calculus of the crowding distance
# INPUT:  PopulationToSort - Unsorted population 
# OUTPUT: Sorted population with ranking and crowding distance 
SortByNonDominationAndObtainCrowdingDistance <- function(PopulationToSort) {
  # sort individuals by non domination
  Ranking <- 
    fastNonDominatedSorting(PopulationToSort[, (colnames(PopulationToSort) %in% MyObjectiveNames)])
  
  # transform the output of the sorting into a matrix of 2 columns: 1.- Individual ID.  2.- Rank 
  MyResult <-
    do.call(
      rbind,
      lapply(
        1:length(Ranking),
        function(i) {
          matrix(c(unlist(Ranking[i]), rep(i, length(unlist(Ranking[i])))), ncol = 2, byrow = FALSE)
        }
      )
    )
  
  # order the matrix by individual ID
  MyResult <- MyResult[order(MyResult[, 1]), ]
  
  # add the rank to the data frame
  PopulationToSort$Rank <- MyResult[, 2]
  
  # calculate (MAX - MIN) of every objective function
  Range <- apply(PopulationToSort[, (colnames(PopulationToSort) %in% MyObjectiveNames)], 2, max) - 
    apply(PopulationToSort[, (colnames(PopulationToSort) %in% MyObjectiveNames)], 2, min)
  
  # create a matrix just removing the individuals' codes and the crowding distances
  PopulationMatrix <- as.matrix(PopulationToSort[, !colnames(PopulationToSort) %in% c("Individual", "CrowdingDistance")])
  
  PopulationToSort$CrowdingDistance <- apply(crowdingDist4frnt(pop = PopulationMatrix, rnk = Ranking, rng = Range), 1, sum)
  
  return(PopulationToSort)
}



############################################################################################################
######################### -------- FUNCTION FROM R PACKAGE nsga2R -------- #################################
############################################################################################################
# Description
# This function estimates the density of solutions surrounding a particular solution within each front.
# It calculates the crowding distances of solutions according to their objectives and those within the
# same front.
# Arguments
# pop Population matrix including decision variables, objective functions, and nondomination rank
# rnk List of solution indices for each front
# rng Vector of each objective function range, i.e. the difference between the maximum and minimum objective function value of each objective
# Value
# Return a matrix of crowding distances of all solutions
crowdingDist4frnt <- function(pop,rnk,rng) {
  popSize <- nrow(pop)
  objDim <- length(rng)
  varNo <- ncol(pop)-1-length(rng)
  cd <- matrix(Inf,nrow=popSize,ncol=objDim)
  for (i in 1:length(rnk)){
    selectRow <- pop[,ncol(pop)] == i
    len <- length(rnk[[i]])
    if (len > 2) {
      for (j in 1:objDim) {
        originalIdx <- rnk[[i]][order(pop[selectRow,varNo+j])]
        cd[originalIdx[2:(len-1)],j] = abs(pop[originalIdx[3:len],varNo+j] - pop[originalIdx[1:(len-2)],varNo+j])/rng[j]
      }
    }
  }
  return(cd)
}


# definition of the function to save in a file the best N individuals of the final population 
# INPUTS: File - file to save the individuals 
#         Population - population of individuals
#         N - number of individuals to keep
#         Network - the network to take the names from (to decode the individuals)
# OUTPUT: None
SaveTheBestIndividualsFromFinalPopulation <- function (BestIndividualsFile, Population, N, Network) {
  
  # loop through the N best individuals in the population
  for (i in 1:N) {
    Ind <- Population[[i, "Individual"]] # get the individual's code
    
    # get the names of the correponding nodes
    DecodedInd <- names(V(Network)[Ind])
    
    write(
      paste(
        c(
          DecodedInd, 
          Population[i, c("AverageNodesScore", "Density", "Rank", "CrowdingDistance")]
        ), 
        collapse=" ", sep=""), 
      file = BestIndividualsFile, 
      append = TRUE,
      sep = ",")
  }
  return(TRUE)  
}

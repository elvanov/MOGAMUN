############ MOGAMUN's functions

# definition of the function that filters the DE results, to remove duplicates
# INPUTS: DE_results: data frame with 5 columns (output of edgeR), 
#                     {gene, logFC, logCPM, PValue, FDR}
# OUTPUT: DE analysis results, with only unique values
RemoveDuplicates_DE_Results <- function(DE_results) {
    # remove NAs
    DE_results <- DE_results[!is.na(DE_results$gene), ]
    
    dup <- unique(as.character(DE_results$gene[
        which(duplicated(as.character(DE_results$gene)))]))
    
    DE_results_filtered <- DE_results[!as.character(DE_results$gene) %in% dup, ]
    
    for (d in dup) {
        DupRows <- DE_results[as.character(DE_results$gene) == d, ]
        
        DE_results_filtered <- 
            rbind(DE_results_filtered, DupRows[which.min(DupRows$FDR), ])
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
    
    # loop through all the layers to get the corresponding subnetwork 
    for (LayerFile in Files) {
        # load layer
        Layer <- data.frame(read.table(LayerFile, header = FALSE, 
                    stringsAsFactors = FALSE))
        
        # create network with the current layer
        # NOTE. All the nodes will have the same ID in every layer
        CurrentNetwork <- 
            graph_from_data_frame(d = Layer, 
                vertices = AllNodesToTake, directed = FALSE)
        
        # add subnetwork as a layer into the multiplex network
        Multiplex[[ length(Multiplex) + 1 ]] <- CurrentNetwork
    }
    
    return (Multiplex)
}


# definition of the function that generates the merged network
# INPUTS: Files - list of file names that contain the biological networks data
#         Multiplex - multiplex network
# OUTPUT: Merged network (igraph object)
GenerateMergedNetwork <- function(Files, Multiplex) {
    # declare empty data frame for the edges of the merged
    Merged <- 
        data.frame(V1 = character(), V2 = character(), Layer = character())
    
    # loop through all the layers to get all the interactions
    for (LayerFile in Files) {
        # load layer
        Layer <- data.frame(read.table(LayerFile, header = FALSE, 
            stringsAsFactors = FALSE), Layer = LayerFile)
        Merged <- rbind(Merged, Layer)
    }
    
    # create merged network, keeping the same order of the nodes as in the mx
    MergedNetwork <- 
        graph_from_data_frame(
            d = Merged, vertices = names(igraph::V(Multiplex[[1]])), 
            directed = FALSE
        )
    
    return (MergedNetwork)
}

# definition of the function that gets the sorted list of nodes 
# INPUTS: Files - list of file names that contain the biological networks data
# OUTPUT: Sorted list of nodes
GetListOfAllNodesPresentInLayers <- function(Files) {
    # declare empty variable to store the nodes
    AllNodes <- NULL
    
    # loop through all the layers to get the corresponding subnetwork 
    for (LayerFile in Files) {
        # load layer
        Layer <- data.frame(read.table(LayerFile, header = FALSE, 
            stringsAsFactors = FALSE))
        
        AllNodes <- c(AllNodes, unique(union(Layer[, 1], Layer[, 2])))
    }
    
    # remove duplicates
    AllNodes <- unique(AllNodes)
    
    # order alphabetically
    AllNodes <- sort(AllNodes)
    
    return (AllNodes)
}



# definition of the function that calculates the node score of a list of genes
# INPUTS: DE_results - differential expression analysis results 
#         ListOfGenes - list of gene names to calculate the node score
#         Measure - 'PValue' or 'FDR' for the node score calculus
# OUTPUT: nodes scores
GetNodesScoresOfListOfGenes <- function(DE_results, ListOfGenes, Measure) {
    ListOfGenes <- unique(ListOfGenes)
    
    # calculate nodes scores for the genes with formula: z = inverse CDF(1-p) 
    # if gene has an expression value, and z = -Inf otherwise
    NodesScores <- as.numeric(vapply(
        ListOfGenes,
        function(X) {
            ifelse (
                X %in% as.character(DE_results$gene), 
                qnorm( 1 - DE_results[[Measure]][DE_results$gene == X] ), 
                (-Inf)
            )
        }, 
        numeric(1)
    ))
    
    # LEAVE THE +INF AND -INF UNTIL AFTER THE NORMALIZATION!!! 
    # THESE SHOULD ALWAYS BE 1 AND 0, RESPECTIVELY
    
    # normalizes nodes scores to be in the range [0-1], ignoring +Inf and -Inf 
    NodesScores <-
        (NodesScores - min(NodesScores[which(is.finite(NodesScores))])) /
        (
            max(NodesScores[which(is.finite(NodesScores))]) - 
                min(NodesScores[which(is.finite(NodesScores))])
        )
    
    # replace +Inf with 1
    NodesScores[which(is.infinite(NodesScores) & NodesScores > 1)] <- 1
    
    # replace -Inf with 0
    NodesScores[which(is.infinite(NodesScores) & NodesScores < 1)] <- 0
    
    return(NodesScores)
}



# definition of the function to generate the initial population
# INPUTS: PopSize - population size, i.e. number of individuals to create
#         Multiplex - multiplex network (it is the search space)
#         LoadedData - general data 
# OUTPUT: Population of size PopSize with integer codification
GenerateInitialPop <- function (PopSize, Multiplex, LoadedData) {
    MinSize <- LoadedData$MinSize
    MaxSize <- LoadedData$MaxSize
    
    MyPopulation <- vector("list", PopSize) # initialize an empty vector 
    for (i in seq_len(PopSize)) { # loop from 1 to the total population size
        
        # randomly determine the size of the individual
        SizeOfIndividual <- 
            sample(
                x = MinSize:MaxSize, 
                size = 1
            )
        
        Root <- PickRoot(MyMx = Multiplex, LoadedData = LoadedData)
        
        # use DFS (Depth First Search) from the chosen root
        # NOTE. DFS will be used in the initial population because it allows 
        #       to go further from the root
        
        # makes a random DFS, i.e., the branches are shuffled before 
        # visiting them, in order to generate different individuals each time
        DFS <-
            DFS_iterative_Mx(
                MyMx = Multiplex,
                Root = Root,
                SizeOfIndividual = SizeOfIndividual, 
                LoadedData = LoadedData
            )
        
        # add individual to the population
        MyPopulation[i] <- list(DFS)
    }
    return(MyPopulation)
}


# definition of the function to perform a RANDOM depth first search of a given 
# length on a MULTIPLEX network. This means that the branch to be visted is 
# randomnly picked from all the available ones
# INPUTS: MyMx - a multiplex network
#         Root - name of the node that will be used as root for the search
#         SizeOfIndividual - desired length of the search
#         LoadedData - general data 
# OUTPUT: Individual with the desired length
DFS_iterative_Mx <- function (MyMx, Root, SizeOfIndividual, LoadedData) {
    MaxNumberOfAttempts <- LoadedData$MaxNumberOfAttempts
    Multiplex <- LoadedData$Multiplex
    KeepLooking <- TRUE # flag to control the execution of the current function
    Attempts <- 0
    
    while(KeepLooking == TRUE) {
        # make depth first search on MyMx using Root as seed
        Result <- MakeDFS(MyMx, Root, SizeOfIndividual) 
        
        # security check of individual's size
        if (length(Result) != SizeOfIndividual) {
            if (Attempts > MaxNumberOfAttempts) { # if max attemps reached
                MyMx <- Multiplex # use the big original multiplex
                Root <- PickRoot(MyMx, LoadedData) # pick a random root 
            } else { # if we still have attempts left
                Attempts <- Attempts + 1 # increment number of attempts
                Root <- PickRoot(MyMx, LoadedData) # pick a root 
            }
        } else { # if a good root was found and the ind has the desired size
            KeepLooking <- FALSE # deactivate flag to stop the search
        }
    }
    
    ##### the following conversion has to be done because when a subnetwork is 
    ##### created, the nodes get new IDs, according to the size of the new 
    ##### subnetwork, but the individuals should always have IDs with respect 
    ##### to the global network, therefore the IDs need to be "re-calculated"
    
    Nodes <- names(igraph::V(MyMx[[1]])[Result]) # get local nodes' names
    NodesIDs <- GetIDOfNodes(Nodes, Multiplex[[1]]) # get IDs
    
    return (NodesIDs)
}



# definition of the function to perform a RANDOM depth first search of a given 
# length on a MULTIPLEX network. This means that the branch to be visted is 
# randomnly picked from all the available ones
# INPUTS: MyMx - a multiplex network
#         Root - name of the node that will be used as root for the search
#         SizeOfIndividual - desired length of the search
# OUTPUT: Individual 
MakeDFS <- function(MyMx, Root, SizeOfIndividual) {
    Discovered <- Result <- CurrentLayer <- NULL  # initialization
    Stack <- Root  # initialize the stack of the pending nodes to visit 
    
    # loop while there are pending elements and the desired size hasn't been met
    while(length(Stack) > 0 && length(Result) < SizeOfIndividual) {
        # check if we will change of layer (first loop or 50% chance later)
        if (is.null(CurrentLayer) || runif(1, 0, 1) <= 0.5) {
            AvLayers <- seq_len(length(MyMx)) # list all layers
            if (length(AvLayers) > 1) { 
                AvLayers <- AvLayers[!AvLayers %in% CurrentLayer] 
            }
            
            # choose a new layer
            CurrentLayer <- 
                ifelse(length(AvLayers > 0), sample(AvLayers, 1), CurrentLayer) 
        }
        Node <- Stack[length(Stack)]  # take a node from the stack to visit it
        Stack <- Stack[-length(Stack)] # remove node from stack
        
        if (!(Node %in% Discovered)) { # verify if the node hasn't been visited
            Discovered <- c(Discovered, Node) # add node to the list of visited
            
            # when findind a disconnected node in a layer, change layer 
            KeepGoing <- TRUE
            VisitedLayers <- NULL
            AvLayers <- seq_len(length(MyMx))
            
            while (KeepGoing == TRUE) {
                MyNeighbors <- 
                    names(igraph::neighbors(MyMx[[CurrentLayer]], Node)) 
                MyNeighbors <- MyNeighbors[!MyNeighbors %in% Discovered] 
                
                if (length(MyNeighbors) == 0) {
                    VisitedLayers <- c(VisitedLayers, CurrentLayer) 
                    AvLayers <- AvLayers[!AvLayers %in% VisitedLayers]
                    if (length(AvLayers) > 0) { 
                        CurrentLayer <- AvLayers[sample(length(AvLayers), 1)]
                    } else { KeepGoing <- FALSE }
                } else { KeepGoing <- FALSE }
            }
            
            MyNeighbors <- MyNeighbors[sample(length(MyNeighbors))] # shuffle 
            Stack <- c(Stack, MyNeighbors) # add shuffled neighbors to stack
            Result <- c(Result, Node) # add node to the individual
        }
    }
    
    return(Result) # return result
}


# Definition of the function to pick a random node to be the root of a search. 
# DEG are preferred over non-DEG
# INPUT:  MyMx - network containing all the available nodes
#         LoadedData - general data 
# OUTPUT: Root node
PickRoot <- function (MyMx, LoadedData) {
    DEG <- LoadedData$DEG
    
    SearchSpaceGenes <- names(igraph::V(MyMx[[1]]))
    
    # get the list of all the DE genes (from the whole DE analysis test)
    DEGenesNames <- as.character(DEG$gene)
    
    # get the list of DE genes in the search space
    DEGAvail <- as.character(DEG[DEGenesNames %in% SearchSpaceGenes, "gene"])
    
    # verify if there is at least one DE gene
    if (length(DEGAvail) > 0) {
        # randomly pick a DEG to be the root of the tree (individual)
        Root <- DEGAvail[sample(seq_len(length(DEGAvail)), 1)]
    } else {
        # randomly pick any gene from the list of available genes
        Root <- SearchSpaceGenes[sample(seq_len(length(SearchSpaceGenes)), 1)]
    }
    
    return(Root)
}


# definition of the function to evaluate a whole population
# INPUTS: MyPopulation - population to evaluate (codes of the individuals)
#         Multiplex - network to which the population belongs to
#         LoadedData - general data 
# OUTPUT: Evaluation of the individual (fitness)
EvaluatePopulation <- function (MyPopulation, Multiplex, LoadedData) {
    FitnessData <-
        do.call(
            rbind,
            lapply(
                seq_len(length(MyPopulation)), # from 1 to PopSize
                function(i) {
                    EvaluateInd(
                        Individual = MyPopulation[[i]],
                        MultiplexNetwork = Multiplex,
                        LoadedData = LoadedData
                    )
                }
            ))
    
    return (FitnessData)
}


# definition of the function to evaluate an individual
# INPUTS: Individual - individual to evaluate
#         MultiplexNetwork - network to which the individual belongs to
#         LoadedData - general data 
# OUTPUT: Evaluation of the individual (fitness)
EvaluateInd <- function (Individual, MultiplexNetwork, LoadedData) {
    GenesNS <- LoadedData$GenesWithNodesScores
    DensityMX <- LoadedData$DensityPerLayerMultiplex
    
    # verifies that there is at least one node present in the individual
    if ( length(Individual) > 0 ) {
        
        # gets the sum of the nodes scores of the subnetwork
        SumNodesScores <- sum(GenesNS[GenesNS$gene %in% 
            igraph::V(MultiplexNetwork[[1]])$name[Individual], "nodescore"])
        
        AverageNodesScore <- SumNodesScores / length(Individual)
        
        SumDensityAllLayers <- 0
        
        # loop through all the layers in the multiplex
        for (layer in seq_len(length(MultiplexNetwork))) {
            # get the inds subnetwork corresponding in the current layer
            Subnetwork <- 
                igraph::induced_subgraph(MultiplexNetwork[[layer]], Individual)
            
            # calculate the density of the subnetwork in the current layer
            SubnetworkDensity <- igraph::graph.density(Subnetwork)
            if (is.nan(SubnetworkDensity)) {
                SubnetworkDensity <- 0
            }
            
            # add the normalized subnetwork's density with respect to 
            # the density of the current layer
            SumDensityAllLayers <- 
                SumDensityAllLayers + (SubnetworkDensity / DensityMX[layer])
        }
        
        Res <- data.frame(AverageNodesScore = AverageNodesScore, 
            Density = SumDensityAllLayers)
    } else {
        Res <- data.frame(AverageNodesScore = 0, Density = 0)
    }
    
    return (Res)
}


############## FUNCTION FROM PACKAGE nsga2R ##################
# definition of the function to sort by non-domination a whole population
# INPUTS: inputData - data frame (contains the individuals' objective functions)
# OUTPUT: Rank of the population
fastNonDomSorting <-   function(inputData) {
    popSize = nrow(inputData)
    idxDominators = vector("list", popSize)
    idxDominatees = vector("list", popSize)
    for (i in seq_len(popSize-1)) {
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
        for (i in seq_len(noSolInCurrFrnt)) {
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



# definition of the function to generate a new population from an existing one
# INPUTS: LoadedData - general data 
#         Population - full population of individuals (parents)
# OUTPUT: New population (children)
MakeNewPopulation <- function(LoadedData, Population) {
    PopSize <- LoadedData$PopSize
    MaxNumberOfAttempts <- LoadedData$MaxNumberOfAttempts
    TournSize <- LoadedData$TournamentSize
    Multiplex <- LoadedData$Multiplex
    MyNewPopulation <- vector("list", PopSize) # initialize vector 
    
    # loop to generate the new population. In each loop, 2 children are created
    for(i in seq(from = 1, to = PopSize, by = 2)) {
        Attempts <- 0 # counter for max. attemps to find parents
        KeepLooking <- TRUE  # flag to keep looking for parents
        
        while (Attempts < MaxNumberOfAttempts & KeepLooking == TRUE) {
            Parent1 <- TournamentSelection(TournSize, Population) # selection
            Parent2 <- GetParent2(Parent1, Population, LoadedData)
            
            if (is.null(Parent2)) { # verify if the search was unsuccessful
                Attempts <- Attempts + 1 # increment no. of attempts  
            } else { KeepLooking <- FALSE } # parent 2 found - stop the search
        }
        
        if (Attempts == MaxNumberOfAttempts) { # if unsuccessful search
            print("Max. no. of attemps to find compatible parents")
            
            # generate two random individuals
            Children <- GenerateInitialPop( PopSize = 2, Multiplex = Multiplex, 
                                            LoadedData = LoadedData )
        } else {
            Children <- Crossover(Parent1, Parent2, LoadedData) # mate parents
        }
        Children <- Mutation(Children, Multiplex, LoadedData) # mutate children
        MyNewPopulation[i] <- Children[1] # add child 1 to the population
        MyNewPopulation[i+1] <- Children[2] # add child 2 to the population
    }
    
    # evaluate offspring
    FitnessData <- EvaluatePopulation(MyNewPopulation, Multiplex, LoadedData)
    
    # generate data frame with the individuals and their fitness
    NewPopulation <- data.frame("Individual" = I(MyNewPopulation), FitnessData,
                                Rank = rep(0, nrow(FitnessData)), 
                                CrowdingDistance = rep(0, nrow(FitnessData)))
    NewPopulationForReplacement <- 
        Replacement(Parents = Population, Children = NewPopulation, 
                    LoadedData = LoadedData) # prepare new population

    return(NewPopulationForReplacement) # return the new population 
}



# definition of the function to choose a compatible individual with parent 1
# INPUTS: Parent1 - parent 1
#         Population - full population of individuals
#         LoadedData - general data 
# OUTPUT: Compatible parent 2 or NULL if the search was unsuccessful
GetParent2 <- function(Parent1, Population, LoadedData) {
    Multiplex <- LoadedData$Multiplex
    PopSize <- LoadedData$PopSize
    TournSize <- LoadedData$TournamentSize
    
    ### filter population to leave only inds near parent 1 
    # get nodes ids with respect to the big network
    NodesIDsOfParent1 <- sort( unlist(Parent1$Individual) )
    
    # get list of nodes IDs in parent 1 and their neighbors
    NeighborsNodesP1 <- GetNeighbors(NodesIDsOfParent1, Multiplex )
    
    # get inds that contain at least one node from the previous list
    NeighborsIndsP1 <- unlist( vapply( seq_len(PopSize), function(X) { 
        if (length(intersect(unlist(Population[X,"Individual"]), 
            NeighborsNodesP1)) > 0){ X } else { 0 }}, numeric(1) ) )
    
    # filter population and leave individuals near parent 1
    PotentialIndsParent2 <- Population[NeighborsIndsP1, ]
    
    # verify if parent 1 is in the list and if so, remove it
    if ( rownames(Parent1) %in% rownames(PotentialIndsParent2) ) {
        Parent1ID <- 
            which(rownames(PotentialIndsParent2) == rownames(Parent1))
        PotentialIndsParent2 <- PotentialIndsParent2[ -Parent1ID, ]
    }
    
    if ( nrow(PotentialIndsParent2) >= 2 ) { # if there are more than 2 inds
        
        # tournament to choose parent 2
        Parent2 <- TournamentSelection(TournSize, PotentialIndsParent2)
    } else if (nrow(PotentialIndsParent2) == 1) { 
        # if there is only 1 ind, keep it
        Parent2 <- PotentialIndsParent2
    } else { # unsuccessful search for compatible parents
        Parent2 <- NULL
    }
    
    return(Parent2)
}


# definition of the function to perform a tournament selection
# INPUTS: TournamentSize - size of the tournament
#         TournPop - population to do the tournament with
# OUTPUT: Individual that won the tournament
TournamentSelection <- function (TournamentSize, TournPop) {
    # randomly choose as many individuals as the tournament size indicates
    ids <- 
        sample(seq_len(nrow(TournPop)), TournamentSize, replace=FALSE)
    
    # verify all both individuals are in the same Pareto front
    if (length(unique(TournPop$Rank[ids])) == 1) {
        # verify if these individuals have information about crowding distance 
        # (they won't if generation == 1)
        if ("CrowdingDistance" %in% colnames(TournPop)) {
            # get id of the individual with the highest crowding distance
            winner <- which(TournPop$CrowdingDistance[ids] == 
                max(TournPop$CrowdingDistance[ids])
                )
        } else {
            # if there is no crowding distance and the individuals have the 
            # same rank, they are all considered as winners
            winner <- c(seq_len(TournamentSize))
        }
    } else { 
        # if the individuals are in different Pareto fronts, the winner 
        # of the tournament is the one with lowest value (best ranking)
        winner <- which(TournPop$Rank[ids] == min(TournPop$Rank[ids]))
    }
    
    # verify if there was a tie
    if (length(winner) > 1) {
        # pick a single winner randomly
        winner <- winner[ sample(seq_len(length(winner)), 1) ]
    }
    
    TournamentWinner <- TournPop[ids[winner],]
    
    return(TournamentWinner)
}


# definition of the function to get the neighbors of a list of nodes, 
# considering all the layers
# INPUTS: NodeList - list of nodes to look for its neighbors
#         Multiplex - multiplex network
# OUTPUT: List of neighbors
GetNeighbors <- function(NodeList, Multiplex) {
    Neighbors <- NULL
    
    # loop through all the layer of the multiplex
    for (i in seq_len(length(Multiplex))) {
        # add to the list the neighbors of the node list in the current layer
        Neighbors <- 
            c(Neighbors, unlist(lapply(NodeList, function(X) {
                igraph::neighbors(Multiplex[[i]], X)})))
    }
    
    # sort result and delete duplicates
    Neighbors <- sort(unique(Neighbors))
    
    return(Neighbors)
}


# definition of the function to perform crossover
# INPUTS: Parent1 - first parent to cross
#         Parent2 - second parent to cross
#         LoadedData - general data 
# OUTPUT: Two children (a single connected component per child)
Crossover <- function (Parent1, Parent2, LoadedData) {
    MaxNumberOfAttempts <- LoadedData$MaxNumberOfAttempts
    MinSize <- LoadedData$MinSize
    MaxSize <- LoadedData$MaxSize
    CrossoverRate <- LoadedData$CrossoverRate
    Multiplex <- LoadedData$Multiplex
    Children <- NULL # intialize empty list
    p <- runif(1, 0, 1) # generate a random number between 0 and 1
    
    if (p <= CrossoverRate) { # check if crossover is to be performed
        # join both parents
        NodesP <- union(unlist(Parent1$Individual), unlist(Parent2$Individual))
        
        # first, generate the multiplex of the joint parents
        Mx_Parents <- FilterMultiplex(Multiplex=Multiplex, NodesToKeep=NodesP)
        
        # if length is minimum, copy parents 
        if (length(NodesP) == MinSize) {
            Children <- rbind(Parent1$Individual, Parent2$Individual)
        } else {
            for (k in seq_len(2)) { # loop to generate two children
                # pick a size for the child 
                Size <- ifelse(length(NodesP) >= MaxSize, 
                    sample(MinSize:MaxSize, 1), 
                    sample(MinSize:length(NodesP), 1))
                # pick root
                Root <- PickRoot(MyMx = Mx_Parents, LoadedData = LoadedData) 
                SearchMethod <- sample(c("DFS", "BFS"), 1) # pick search method
                
                if (SearchMethod == "DFS") { # depth first search
                    DFS <- DFS_iterative_Mx(MyMx = Mx_Parents, Root = Root,
                            SizeOfIndividual = Size, LoadedData = LoadedData)
                    
                    Children[k] <- list(DFS) # add child to the population
                } else if (SearchMethod == "BFS") { # breadth first search
                    BFS <- BFS_iterative_Mx(MyMx = Mx_Parents, Root = Root, 
                            SizeOfIndividual = Size, LoadedData = LoadedData )
                    
                    Children[k] <- list(BFS) # add child to the population
                }
            }
        }
    } else { # if no crossover was performed, make a copy of the parents
        Children <- rbind(Parent1$Individual, Parent2$Individual)
    }
    
    return(Children)
}


# definition of the function to perform a RANDOM breadth first search of a 
# given length on a MULTIPLEX network. This means that the nodes to be visted  
# are always randomnly picked from all the available ones
# INPUTS: Multiplex - the multiplex network
#         Root - name of the node that will be used as root for the search
#         SizeOfIndividual - desired length of the search
#         LoadedData - general data 
# OUTPUT: Individual with the desired length
BFS_iterative_Mx <- function (MyMx, Root, SizeOfIndividual, LoadedData) {
    MaxNumberOfAttempts <- LoadedData$MaxNumberOfAttempts
    Multiplex <- LoadedData$Multiplex
    KeepLooking <- TRUE # flag to control the execution of the current function
    Attempts <- 0
    
    while(KeepLooking == TRUE) {
        Result <- MakeBFS(MyMx, Root, SizeOfIndividual)
        
        # security check of individual's size
        if (length(Result) != SizeOfIndividual) {
            if (Attempts > MaxNumberOfAttempts) { # if max attemps reached
                MyMx <- Multiplex # use the big original multiplex
                Root <- PickRoot(MyMx, LoadedData) # pick a random root 
            } else { # if we still have attempts left
                Attempts <- Attempts + 1 # increment number of attempts
                Root <- PickRoot(MyMx, LoadedData) # pick a root 
            }
        } else { # if a good root was found and the ind has the desired size
            KeepLooking <- FALSE # deactivate flag to stop the search
        }
    }
    
    ##### the following conversion has to be done because when a subnetwork is 
    ##### created, the nodes get new IDs, according to the size of the new 
    ##### subnetwork, but the individuals should always have IDs with respect 
    ##### to the global network, therefore the IDs need to be "re-calculated"
    
    # get the names of the nodes in the local network with the local IDs
    Nodes <- names(igraph::V(MyMx[[1]])[Result])
    
    # get the global IDs of the corresponding nodes
    NodesIDs <- GetIDOfNodes(Nodes, Multiplex[[1]])
    
    return (NodesIDs)
}



# definition of the function to perform a RANDOM breadth first search of a  
# given length on a MULTIPLEX network. This means that the branch to be visted 
# is randomnly picked from all the available ones
# INPUTS: MyMx - a multiplex network
#         Root - name of the node that will be used as root for the search
#         SizeOfIndividual - desired length of the search
# OUTPUT: Individual 
MakeBFS <- function(MyMx, Root, SizeOfIndividual) {
    Discovered <- Result <- CurrentLayer <- NULL  # initialization
    Queue <- Root  # initialize the queue of the pending nodes to visit 
    
    # loop while there are pending elements and the desired size hasn't been met
    while(length(Queue) > 0 && length(Result) < SizeOfIndividual) {
        # check if we will change of layer 
        if (is.null(CurrentLayer) || runif(1, 0, 1) <= 0.5) {
            AvLayers <- seq_len(length(MyMx)) # list all layers
            if (length(AvLayers) > 1) { 
                AvLayers <- AvLayers[!AvLayers %in% CurrentLayer] }
            
            # choose a new layer
            CurrentLayer <- 
                ifelse(length(AvLayers > 0), sample(AvLayers, 1), CurrentLayer) 
        }
        Node <- Queue[1]  # take a node from the queue to visit it
        Queue <- Queue[-1] # remove node from queue
        
        if (!(Node %in% Discovered)) { # verify if the node hasn't been visited 
            Discovered <- c(Discovered, Node) # add node to the list of visited
            
            # when findind a disconnected node in a layer, change layer 
            KeepGoing <- TRUE
            VisitedLayers <- NULL
            AvLayers <- seq_len(length(MyMx))
            
            while (KeepGoing == TRUE) {
                MyNeighbors <- 
                    names(igraph::neighbors(MyMx[[CurrentLayer]], Node)) 
                MyNeighbors <- MyNeighbors[!MyNeighbors %in% Discovered] 
                
                if (length(MyNeighbors) == 0) {
                    VisitedLayers <- c(VisitedLayers, CurrentLayer) 
                    AvLayers <- AvLayers[!AvLayers %in% VisitedLayers] 
                    if (length(AvLayers) > 0) { 
                        CurrentLayer <- AvLayers[sample(length(AvLayers), 1)]
                    } else { KeepGoing <- FALSE }
                } else { KeepGoing <- FALSE }
            }
            
            MyNeighbors <- MyNeighbors[sample(length(MyNeighbors))] # shuffle 
            Queue <- c(Queue, MyNeighbors) # add shuffled neighbors to the stack
            Result <- c(Result, Node) # add node to the individual
        }
    }
    return(Result) # return result    
}



# definition of the function that filters the multiplex network, to keep 
# only the specified list of nodes
# INPUTS: Multiplex - multiplex network to filter
#         ListOfNodes - list of nodes to keep
# OUTPUT: Filtered version of the multiplex
FilterMultiplex <- function(Multiplex, NodesToKeep) {
    # declare empty list to store the multiplex
    FilteredMultiplex <- list()
    
    # loop through all the layers to get the corresponding subnetwork 
    for (i in seq_len(length(Multiplex))) {
        # create network with the current layer
        CurrentNetwork <- igraph::induced_subgraph(Multiplex[[i]], NodesToKeep)
        
        # add subnetwork as a layer into the multiplex network
        FilteredMultiplex[[ length(FilteredMultiplex) + 1 ]] <- CurrentNetwork
    }
    
    return (FilteredMultiplex)
}


# definition of the function to get the list of IDs of a set of nodes
# INPUTS: ListOfNodes - list of nodes' names
#         GlobalNetwork - general network to look for our list of nodes
# OUTPUT: List of IDs
GetIDOfNodes <- function(ListOfNodes, GlobalNetwork) {
    return (which(names(igraph::V(GlobalNetwork)) %in% ListOfNodes))
}


# definition of the function to perform mutation
# INPUTS: Individuals - the individuals to perform mutation with
#         Multiplex - network to where the individuals belong to
#         LoadedData - general data 
# OUTPUT: Two mutated individuals
Mutation <- function (Individuals, Multiplex, LoadedData) {
    Merged <- LoadedData$Merged
    DEG <- LoadedData$DEG
    MutationRate <- LoadedData$MutationRate
    
    Mutants <- NULL
    
    # loop through all the individuals to be mutated
    for (i in seq_len(length(Individuals))) { 
        p <- runif(1, 0, 1) # generate a random number between 0 and 1
        
        if (p <= MutationRate) { # check if mutation is to be performed
            # make subnetwork
            IndToMutNet <- igraph::induced_subgraph(Merged, Individuals[[i]]) 
            IndToMutDeg <- igraph::degree(IndToMutNet) # get nodes' degrees
            
            # remove DEG from the list
            IndToMutDeg <- 
                IndToMutDeg[ !names(IndToMutDeg) %in% as.character(DEG$gene) ]
            
            # get the list of nodes with the min degree (peripheral nodes)
            PeripheralNodes <- which.min(IndToMutDeg)
            
            # get the list of node names that can be mutated
            PotNodesToMutate <- names(PeripheralNodes)
            
            # make vector with the nodes to mutate
            NodesToMutate <- 
                apply(array(rep(0, length(PotNodesToMutate))), 1, function(X) { 
                    ifelse(runif(1, 0, 1) <= MutationRate, 1, 0)})
            
            # gets the list of chromosomes to be mutated
            NodesToMutate <- which(NodesToMutate == 1) 
            
            # if at least one of the nodes will be mutated
            if (length(NodesToMutate) > 0) { 
                Individuals[[i]] <- 
                    MutateNodes(Individuals[[i]], IndToMutNet,  NodesToMutate, 
                                PotNodesToMutate, LoadedData) 
            } else { # if no nodes were removed, add a new neighbor
                Individuals[[i]] <- 
                    AddNode(Individuals[[i]], IndToMutNet, LoadedData) 
            }
        }
        Mutants[i] <- Individuals[i] # save the individual in the mutants' list
    }
    return(Mutants)
}



# definition of the function to mutate a given list of nodes 
# INPUTS: Ind - individual to mutate
#         NodesToMutate - list of nodes to mutate
#         PotNodesToMutate - list of potential nodes to mutate
#         LoadedData - general data 
# OUTPUT: Mutated individual
MutateNodes <- function(Ind, IndToMutNet, NodesToMutate, PotNodesToMutate, 
                        LoadedData) {
    DEG <- LoadedData$DEG # get DEG
    
    # if at least one of the nodes will be mutated
    for (j in NodesToMutate) { # loop through the nodes to remove
        MutatedNetwork <- 
            igraph::delete_vertices(IndToMutNet, PotNodesToMutate[j])
        
        # obtain neighbors of nodes
        Neighbors_OrInd <- 
            names(unlist(lapply(names( igraph::V(IndToMutNet)),
                    function(X) {igraph::neighbors(LoadedData$Merged, X)})))
        Neighbors_MutInd <- 
            names(unlist(lapply(names(igraph::V(MutatedNetwork)), 
                function(X) {igraph::neighbors(LoadedData$Merged, X)})))
        
        # delete nodes that originally belonged to the individual
        Neighbors_OrInd <- 
            Neighbors_OrInd[
                !Neighbors_OrInd %in% names(igraph::V(IndToMutNet))]
        
        Neighbors_MutInd <-
            Neighbors_MutInd[
                !Neighbors_MutInd %in% names(igraph::V(IndToMutNet))]
        
        # check if network is connected
        if( igraph::is.connected(MutatedNetwork) ) { 
            # save changes 
            Ind <- GetIDOfNodes(names(igraph::V(MutatedNetwork)), 
                                LoadedData$Multiplex[[1]])
            IndToMutNet <- MutatedNetwork
            AvNeighbors <- Neighbors_MutInd
        } else { AvNeighbors <- Neighbors_OrInd }
        
        # if there is at least 1 available neighbor to be added 
        if(length(AvNeighbors) > 0) {
            # pick a node to add
            NewNodeID <- GetNodeToAdd(AvNeighbors, LoadedData) 
            
            # add node
            Ind <- c(Ind, NewNodeID) 
        } else { print("Attempt to add a new neighbor FAILED") }
    }
    return(Ind)
}



# definition of the function to add a node to a subnetwork 
# INPUTS: IndToMutNet - subnetwork of the individual to modify
#         LoadedData - general data 
# OUTPUT: Modified individual
AddNode <- function(Ind, IndToMutNet, LoadedData) {
    DEG <- LoadedData$DEG # get DEG
    
    # obtain neighbors of nodes
    AvNeighbors <- 
        names(unlist(lapply(names(igraph::V(IndToMutNet)), 
            function(X) {igraph::neighbors(LoadedData$Merged, X)})))
    
    # delete from the list all the nodes that originally belonged to the ind
    AvNeighbors <- AvNeighbors[!AvNeighbors %in% names(igraph::V(IndToMutNet))]
    
    # if there is at least one available neighbor to be added
    if(length(AvNeighbors) > 0) {
        # pick a node to add
        NewNodeID <- GetNodeToAdd(AvNeighbors, LoadedData) 
        
        # add node
        Ind <- c(Ind, NewNodeID) 
    } else { 
        print("No nodes to be mutated. Attempt to add a new neighbor FAILED") 
    }
    
    return(Ind)
}


# definition of the function to pick a node to add 
# INPUTS: AvNeighbors - available neighbors to choose the node from
#         LoadedData - general data 
# OUTPUT: ID of the node to be added
GetNodeToAdd <- function(AvNeighbors, LoadedData) {
    DEG <- LoadedData$DEG
    Multiplex <- LoadedData$Multiplex
    
    # get the list of neighbors DEG
    DE_Neighbors <- AvNeighbors[AvNeighbors %in% as.character(DEG$gene)]
    
    # if there is at least one DE gene in the list of neighbors, keep
    if (length(DE_Neighbors) > 0) { AvNeighbors <- DE_Neighbors }
    
    # calculate the number of times each node is neighbor
    Incidences <- table(as.factor(AvNeighbors)) 
    
    # get the nodes with highest incidence
    MaxIncidences <- which.max(Incidences) 
    
    # keep nodes with a max incidence
    AvNeighbors <- names(MaxIncidences) 
    
    # pick a random node to add
    RandomNode <- AvNeighbors[sample(seq_len(length(AvNeighbors)), 1)] 
    
    # get ID of the node to add
    NewNodeID <- GetIDOfNodes(RandomNode, Multiplex[[1]]) 
    
    return (NewNodeID)
}


# definition of the function to perform replacement
# INPUTS: Parents - individuals from the previous generation
#         Children - individuals from the new generation
#         LoadedData - general data 
# OUTPUT: The selected invididuals to keep for the next generation
Replacement <- function (Parents, Children, LoadedData) {
    PopSize <- LoadedData$PopSize
    
    # combine the new population (offspring) with old population (parents)
    CombinedPopulation <- rbind(Parents, Children)
    
    ## replace duplicated individuals with random ones
    CombinedPopulation <- ReplaceDuplicatedInds(CombinedPopulation, LoadedData)
    
    # order combined population by rank
    CombinedPopulation <- CombinedPopulation[order(CombinedPopulation$Rank), ]
    
    # get last rank which will be in the replacement population
    LastRank <- CombinedPopulation$Rank[PopSize]
    
    # get last rank individuals
    LastRankInds <- CombinedPopulation[CombinedPopulation$Rank == LastRank, ]
    
    # order by crowding distance
    LastRankInds <- 
        LastRankInds[order(LastRankInds$CrowdingDistance, decreasing = TRUE), ]
    
    # select new population for replacement
    NewPopulationForReplacement <- 
        CombinedPopulation[CombinedPopulation$Rank < LastRank, ]
    NewPopulationForReplacement <- rbind(NewPopulationForReplacement, 
        LastRankInds[seq_len(PopSize - nrow(NewPopulationForReplacement)), ])
    
    return(NewPopulationForReplacement)
}


# definition of the function that replaces all duplicated individuals and 
# those above the permitted threshold of similarity with random individuals
# INPUTS:  CombinedPopulation - Population of size 2*N that corresponds to the 
#                               union of parents and children
#          LoadedData - general data 
# OUTPUT:  A garanteed diverse population of size 2*N
ReplaceDuplicatedInds <- function(CombinedPopulation, LoadedData) {
    DivPop <- CombinedPopulation
    Threshold <- LoadedData$JaccardSimilarityThreshold
    Multiplex <- LoadedData$Multiplex
    Sim <- GetDuplicatedInds(DivPop, Threshold) # get duplicated individuals
    
    if (nrow(Sim) > 0) { # verify if there are duplicated individuals
        IndsToRemove <- NULL
        
        # non-dominated sorting and crowding distance calculus
        SortedPop <- 
            NonDomSort(PopulationToSort = DivPop, LoadedData = LoadedData)
        
        i <- 1
        while(i < nrow(Sim)) { # loop through all the duplicated individuals
            # ID of individuals
            Ind1_ID <- Sim[i, 1]
            Ind2_ID <- Sim[i, 2]
            
            if (Sim$JS[i] < 100 & all(SortedPop$Rank[c(Ind1_ID, Ind2_ID)] == 1) 
                & all(is.infinite(SortedPop$CrowdingDistance[
                c(Ind1_ID, Ind2_ID)])) ) {
                print("Keeping similar individuals. Inf crowding distance.")
                print(SortedPop[c(Ind1_ID, Ind2_ID), ])
            } else { # tournament between the two individuals
                IndToKeep <- TournamentSelection(TournamentSize = 2,
                    TournPop <- SortedPop[c(Ind1_ID, Ind2_ID), ])
                IndsToRemove <- c(IndsToRemove, ifelse(rownames(IndToKeep) == 
                    row.names(SortedPop)[Ind1_ID], Ind2_ID, Ind1_ID)) 
                
                # get and remove all future incidences of the ind
                Ref <- which(Sim == IndsToRemove[length(IndsToRemove)], 
                    arr.ind = TRUE)[, "row"]
                Sim <- Sim[!row.names(Sim) %in% row.names(Sim)[Ref[Ref > i]], ]
            }
            i <- i + 1
        }
        # remove the corresponding individuals
        DivPop<-DivPop[!row.names(DivPop) %in% row.names(DivPop)[IndsToRemove],]
    }
    # generate as many new individuals as duplicated ones
    NewInds <- GenerateInitialPop(
        PopSize = (nrow(CombinedPopulation)-nrow(DivPop)), 
        Multiplex = Multiplex, LoadedData = LoadedData )
    FitnessData <- EvaluatePopulation(NewInds, Multiplex, LoadedData) # eval
    DivPop <- rbind(DivPop, data.frame("Individual" = I(NewInds), 
        FitnessData, Rank = 0, CrowdingDistance = 0))
    DivPop <- NonDomSort(PopulationToSort = DivPop, LoadedData = LoadedData)
    return(DivPop)
}


# definition of the function that finds the duplicated individuals, considering
# a given threshold
# INPUTS:  DivPop - Population
#          Threshold - Individuals over this threshold are duplicated
# OUTPUT:  Similarities and indexes of duplicated individuals
GetDuplicatedInds <- function(DivPop, Threshold) {
    # create all the combinations of 2 numbers with the ids of individuals
    Index <- t(combn(seq_len(nrow(DivPop)), 2))
    
    # calculate the Jaccard similarity index
    Sim <- data.frame(Index1 = Index[, 1], Index2 = Index[, 2], JS = 0) # sim
    Sim$JS <- apply(Index, 1, function(X) {
        Ind1 <- as.character(unlist(DivPop$Individual[X[1]]))
        Ind2 <- as.character(unlist(DivPop$Individual[X[2]]))
        # Jaccard similarity = (intersection of A and B) / (union of A and B)
        JS <- (length(intersect(Ind1, Ind2))) / (length(union(Ind1, Ind2)))*100
        return(as.double(JS)) } )
    Sim <- Sim[Sim$JS >= Threshold, ] # keep "duplicated" individuals
    
    return(Sim)
}

# definition of the function that performs the fast non dominated sorting and 
# the calculus of the crowding distance
# INPUTS:  PopulationToSort - Unsorted population
#          LoadedData - general data 
# OUTPUT:  Sorted population with ranking and crowding distance
NonDomSort <- function(PopulationToSort, LoadedData) {
    ObjectiveNames <- LoadedData$ObjectiveNames
    
    # sort individuals by non domination
    Ranking <- fastNonDomSorting(
        PopulationToSort[, (colnames(PopulationToSort) %in% ObjectiveNames)] )
    
    # transform the output of the sorting into a matrix of 2 columns: 
    # 1.- Individual ID.  2.- Rank
    MyResult <- do.call(rbind, lapply(seq_len(length(Ranking)), function(i) {
        matrix(c(unlist(Ranking[i]), rep(i, length(unlist(Ranking[i])))), 
            ncol = 2, byrow = FALSE) }))
    
    # order the matrix by individual ID
    MyResult <- MyResult[order(MyResult[, 1]), ]
    
    # add the rank to the data frame
    PopulationToSort$Rank <- MyResult[, 2]
    
    # calculate (MAX - MIN) of every objective function
    Range <- apply(PopulationToSort[, (colnames(PopulationToSort) %in% 
        ObjectiveNames)], 2, max) - apply(PopulationToSort[, 
        (colnames(PopulationToSort) %in% ObjectiveNames)], 2, min)
    
    # create a matrix removing the ind codes and the crowding distances
    PopulationMatrix <- as.matrix(PopulationToSort[ , 
        !colnames(PopulationToSort) %in% c("Individual", "CrowdingDistance")])
    
    PopulationToSort$CrowdingDistance <- 
        apply(crowdingDist4frnt(pop = PopulationMatrix, rnk = Ranking, 
            rng = Range), 1, sum)
    
    return(PopulationToSort)
}



###############################################################################
########## -------- FUNCTION FROM R PACKAGE nsga2R -------- ###################
###############################################################################
# Description
# This function estimates the density of solutions surrounding a particular 
# solution within each front. It calculates the crowding distances of solutions 
# according to their objectives and those within the same front.
# INPUTS: pop - Population matrix including decision variables, objective 
#               functions, and nondomination rank
#         rnk - List of solution indices for each front
#         rng - Vector of each objective function range, i.e. the difference
#               between the maximum and minimum objective function value of  
#               each objective
# OUTPUT: Matrix of crowding distances of all solutions
crowdingDist4frnt <- function(pop,rnk,rng) {
    popSize <- nrow(pop)
    objDim <- length(rng)
    varNo <- ncol(pop)-1-length(rng)
    cd <- matrix(Inf,nrow=popSize,ncol=objDim)
    for (i in seq_len(length(rnk))){
        selectRow <- pop[,ncol(pop)] == i
        len <- length(rnk[[i]])
        if (len > 2) {
            for (j in seq_len(objDim)) {
                originalIdx <- rnk[[i]][order(pop[selectRow,varNo+j])]
                cd[originalIdx[2:(len-1)],j] = 
                    abs(pop[originalIdx[3:len],varNo+j] - 
                            pop[originalIdx[seq_len(len-2)],varNo+j])/rng[j]
            }
        }
    }
    return(cd)
}


# definition of the function to save in a file the best N individuals of the 
# final population
# INPUTS: File - file to save the individuals
#         Population - population of individuals
#         N - number of individuals to keep
#         Network - network to take the names from (to decode the individuals)
# OUTPUT: None
SaveFinalPop <- function (BestIndividualsFile, Population, N, Network) {
    # loop through the N best individuals in the population
    for (i in seq_len(N)) {
        Ind <- Population[[i, "Individual"]] # get the individual's code
        
        # get the names of the correponding nodes
        DecodedInd <- names(igraph::V(Network)[Ind])
        
        write(paste(c(DecodedInd, Population[i, c("AverageNodesScore", 
            "Density", "Rank", "CrowdingDistance")]), collapse=" ", sep=""),
            file = BestIndividualsFile, append = TRUE, sep = ",")
    }
    
    return(TRUE)
}



# definition of a function to postprocess the results. This function: 
#    i) calculates the accumulated Pareto front, i.e. the individuals on the 
#       first Pareto front after re-ranking the results from multiple runs 
#       (NOTE. If there is a single run, the result is the set of individuals
#       in the first Pareto front),
#   ii) filters the networks to leave only the interactions between the genes  
#       that are included in the results, 
#  iii) generates some plots (scatter plots and boxplots) of interest, and 
#   iv) (optional) creates a Cytoscape file to visualize the results, merging  
#       the subnetworks with a Jaccard similarity coefficient superior to  
#       JaccardSimilarityThreshold (NOTE. Make sure to open Cytoscape if 
#       VisualizeInCytoscape is TRUE)
#
# INPUTS: ExperimentDir - folder containing the results to be processed
#         LoadedData - list containing several important variables 
#                      (see mogamun_load_data())    
#         JaccardSimilarityThreshold - subnetworks over this Jaccard similarity 
#                                      threshold are merged 
#         VisualizeInCytoscape - boolean value. If it is TRUE, a Cytoscape file 
#                                is generated
# OUTPUT: None
PostprocessResults <- function(ExperimentDir, LoadedData, 
                        JaccardSimilarityThreshold, VisualizeInCytoscape) {
    Threshold <- JaccardSimilarityThreshold
    # get list of directories (each contains the results of a single exp)
    Dirs <- list.dirs(ExperimentDir, recursive = FALSE)
    Dirs <- Dirs[grep("Experiment", Dirs)]
    
    for (Path in Dirs) {
        PathForPlots <- ExperimentsPath <- paste0(Path, "/") # path for plots
        Population <- GetIndividualsAllRuns(ExperimentsPath) # get results
        Inds1stRank_AllRuns <- Population[Population$Rank == 1, ] # filter 
        MaxRank <- max(Population$Rank) # get maximum rank in all the runs
        Nodes1stRank_AllRuns <- 
            data.frame(Nodes = character(), Run = character())
        
        for (r in seq_len(max(Population$Run))) {
            Nodes1stRank_AllRuns <- 
                rbind(Nodes1stRank_AllRuns, data.frame(
                Nodes = unique(unlist(Population[Population$Run == r & 
                Population$Rank == 1, "Individual"])), Run = paste0("Run_", r)))
        }
        
        # if there are results for several runs, make boxplot of similarities
        if (length(unique(Nodes1stRank_AllRuns$Run)) > 1) {
            SimilarityBetweenRunsBoxplot(PathForPlots, Nodes1stRank_AllRuns)
        }
        
        # calculate and plot the accumulated Pareto front 
        AccPF <- ObtainAccParetoFront(PathForPlots, Inds1stRank_AllRuns, 
            Nodes1stRank_AllRuns, LoadedData) 
        AccPF <- RemoveIdenticalInds(AccPF) # remove duplicated inds
        
        # filter networks to keep only interactions between genes in the AccPF
        FilterNetworks(ExperimentsPath, LoadedData, AccPF) 
        
        # merge "duplicated" individuals 
        DiversePopulation <- 
            JoinSimilarInds(AccPF, LoadedData, Threshold) 
        
        # save the merged individuals in a file
        SaveFilteredAccPF(DiversePopulation, ExperimentsPath, Threshold) 
    }
    
    # check if visualization in Cytoscape is to be done
    if (VisualizeInCytoscape == TRUE) {
        CytoscapeVisualization(ExperimentDir, LoadedData)
    }
}



# definition of a function to get all the individuals from all the runs
#
# INPUT:  ExperimentsPath - folder containing the results to be processed 
# OUTPUT: Data frame with the joint population from all the runs
GetIndividualsAllRuns <- function(ExperimentsPath) {
    # initialize data empty data frame
    Population <- 
        data.frame(Run = double(), Individual = list(), 
            AverageNodesScore = double(), Density = double(), 
            Rank = numeric(), CrowdingDistance = double())
    
    PatResFiles <- "MOGAMUN_Results__Run_[0-9]+.txt$"
    
    # get all files in the path that match the pattern
    ResFiles <- list.files(ExperimentsPath, pattern = PatResFiles)
    
    # loop to read all the results 
    for (counter in seq_len(length(ResFiles))) {
        # open and read file
        con <- base::file(paste0(ExperimentsPath, ResFiles[counter]), open="r")
        MyContent <-readLines(con) 
        close(con)
        Run <- as.numeric(gsub("^.*_|\\..*$", "", ResFiles[counter], 
            perl = TRUE))
        
        # loop through all the individuals
        for (i in seq_len(length(MyContent))) {
            Ind <- unlist(strsplit(MyContent[i], " "))
            
            # divide the data of the individual. Each line contains the nodes, 
            # average nodes score, density, rank and crowding distance
            MyInd <- 
                data.frame(Run = Run, 
                    Individual = I(list(sort(Ind[seq_len(length(Ind) - 4)]))), 
                    AverageNodesScore = as.numeric(Ind[(length(Ind) - 3)]),
                    Density = as.numeric(Ind[(length(Ind) - 2)]), 
                    Rank = as.numeric(Ind[(length(Ind)-1)]),
                    CrowdingDistance = as.numeric(Ind[length(Ind)])
            )
            
            # add individual to the population
            Population <- rbind(Population, MyInd)
        }
    }
    
    return(Population)
}



# definition of a function to make the boxplot of similarities between runs
#
# INPUTS: PathForPlots - folder to save the plot
#         Nodes1stRank_AllRuns - results from all the runs
# OUTPUT: None
SimilarityBetweenRunsBoxplot <- function(PathForPlots, Nodes1stRank_AllRuns) {
    AllJaccardSimilarities <- NULL
    
    # loop to calculate the JS between the results from different runs
    for (i in seq_len(length(unique(Nodes1stRank_AllRuns$Run)))) {
        # get nodes of a run
        CurrentRun <- as.character( 
            Nodes1stRank_AllRuns$Nodes[
                Nodes1stRank_AllRuns$Run == paste0("Run_", i)] )
        
        # verify that there is another run
        if ((i+1) <= length(unique(Nodes1stRank_AllRuns$Run))) {
            for (j in (i+1):length(unique(Nodes1stRank_AllRuns$Run))) {
                # get nodes of another run
                NextRun <- as.character(
                    Nodes1stRank_AllRuns$Nodes[
                        Nodes1stRank_AllRuns$Run == paste0("Run_", j)] )
                
                # calculate the intersection and union sizes. 
                # Jaccard similarity coefficient formula: intersection/union
                inter <- length(intersect(CurrentRun, NextRun))
                un <- length(union(CurrentRun, NextRun))
                JS <- inter / un
                
                # add Jaccard similarity result to the list
                AllJaccardSimilarities <- c(AllJaccardSimilarities, JS)
            }
        }
    }
    
    # make a boxplot of similarities
    svg(paste0(PathForPlots, "A_Boxplot_similarities_between_runs.svg"))
    boxplot(AllJaccardSimilarities, 
        main = "Pairwise Jaccard similarities of all the runs")
    dev.off()
}



# definition of a function to calculate and plot the accumulated Pareto front
#
# INPUTS: PathForPlots - folder to save the plot
#         Inds1stRank_AllRuns - individuals with Rank = 1 from all runs
#         Nodes1stRank_AllRuns - nodes in Inds1stRank_AllRuns
#         LoadedData - list containing several important variables 
#                      (see mogamun_load_data())    
# OUTPUT: individuals in the accumulated Pareto front
ObtainAccParetoFront <- function(PathForPlots, Inds1stRank_AllRuns, 
    Nodes1stRank_AllRuns, LoadedData) {
    # re-rank individuals in the first Pareto front
    AccPF <- NonDomSort(PopulationToSort = Inds1stRank_AllRuns,
        LoadedData = LoadedData)
    
    Colors <- c("black", rainbow(length(unique(AccPF$Rank))-1))
    
    # verify if there is more than one run
    if (length(unique(Nodes1stRank_AllRuns$Run)) > 1) {
        
        svg(paste0(PathForPlots, "A_ScatterPlot_AccPF_ALL_INDS.svg"))
        plot(AccPF$AverageNodesScore, AccPF$Density, 
            main = "Accumulated Pareto front: all individuals",
            xlab = "Average nodes score", ylab = "Density", 
            col = Colors[AccPF$Rank], pch = 16, cex = 2)
        dev.off()
        
        # make the legend of the previous plot as a separate file
        svg(paste0(PathForPlots, "A_ScatterPlot_AccPF_LEGEND_ALL_INDS.svg"))
        plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', 
            xlim = 0:1, ylim = 0:1)
        legend("center", legend = paste("Pareto front ", 
            seq_len(length(Colors))), col = Colors, pch = 16, cex = 1)
        dev.off()
    }
    
    # filter the individuals, to leave only those in the first Pareto front 
    # (i.e., the accumulated Pareto front)
    AccPF <- AccPF[AccPF$Rank == 1, ]
    
    svg(paste0(PathForPlots, "A_ScatterPlot_AccPF.svg"))
    plot(AccPF$AverageNodesScore, AccPF$Density, 
        main = "Accumulated Pareto front", xlab = "Average nodes score", 
        ylab = "Density", col = Colors[AccPF$Rank], pch = 16, cex = 2)
    dev.off()
    
    return(AccPF)
}



# definition of a function to identify and remove duplicated individuals
#
# INPUT:  AccPF - individuals in the accumulated Pareto front
# OUTPUT: filtered accumulated Pareto front
RemoveIdenticalInds <- function(AccPF) {
    # create all the combinations of 2 numbers with the ids of individuals
    MyIndexes <- t(combn(seq_len(nrow(AccPF)), 2))
    Similarities <- data.frame(Index1 = MyIndexes[, 1], Index2 = MyIndexes[, 2])
    
    # check for identical individuals
    Similarities$Identical <- apply(Similarities, 1, function(X) { 
        setequal(AccPF$Individual[[X[1]]], AccPF$Individual[[X[2]]])})
    
    # keep only the identical individuals in the comparison list
    Similarities <- Similarities[Similarities$Identical == TRUE, ]
    
    IndsToRemove <- NULL
    s <- 1
    
    # loop through the list
    while (s <= nrow(Similarities)) {
        # verify if the ID of the ind is not in the list of inds to remove
        if (!Similarities[s, 2] %in% IndsToRemove) {
            # add individual and remove references to it in the comparison list
            IndsToRemove <- c(IndsToRemove, Similarities[s, 2])
            Ref <- 
                which(Similarities == Similarities[s, 2], 
                    arr.ind = TRUE)[, "row"]
            Ref <- Ref[Ref > s]
            Similarities <- Similarities[!row.names(Similarities) %in% 
                row.names(Similarities)[Ref], ]
        }
        s <- s + 1
    }
    
    # remove identical individuals
    AccPF <- AccPF[!row.names(AccPF) %in% row.names(AccPF)[IndsToRemove], ]
    
    return(AccPF)
}



# definition of a function to filter networks to keep only the interactions 
# between the genes in the individuals from the accumulated Pareto front
#
# INPUTS: ExperimentsPath - folder to save the filtered networks
#         LoadedData - list containing several important variables 
#                      (see mogamun_load_data())    
#         AccPF - filtered accumulated Pareto front
# OUTPUT: None
FilterNetworks <- function(ExperimentsPath, LoadedData, AccPF) {
    # read the file names of the networks for the current experiment
    myLayers <- list.files(LoadedData$NetworkLayersDir, 
        pattern = paste0("^[", LoadedData$Layers, "]_"))
    
    # get the list of genes in the accumulated Pareto front
    genes <- as.character(unique(unlist(AccPF$Individual)))
    
    # filter layers to keep only interactions among the given list of genes
    for (j in seq_len(length(myLayers))) {
        # read network
        FullNetwork <- read.table(paste0(LoadedData$NetworkLayersDir, 
            myLayers[j]), sep = "\t")
        keep <- NULL # initialize null object
        
        # get all rows in the PPI network that contain these genes
        keep$V1 <- FullNetwork$V1 %in% as.character(genes) # first in column 1
        keep$V2 <- FullNetwork$V2 %in% as.character(genes) # then in column 2
        
        keep <- Reduce("&", keep) # logical AND - gets the list of rows to keep
        MyFilteredNetwork <- FullNetwork[keep, ] # filter network
        MyFilteredNetwork$TypeOfInteraction <- as.character(myLayers[j]) # type 
        
        # save filtered network in a file
        write.table(MyFilteredNetwork, file = paste0(ExperimentsPath, 
            substr(myLayers[j], 1, nchar(myLayers[j]) - 3), "_FILTERED.csv"), 
            row.names = FALSE, quote = FALSE)
    }
    
    # get files that contain the filtered networks 
    myDataFiles <- list.files(ExperimentsPath, pattern = '_FILTERED.csv')

    myFullListOfInteractions <- NULL # initialize interactions' list 
    
    # loop through all the files
    for (i in seq_len(length(myDataFiles))) {
        # get data
        data <- data.frame(read.table(paste0(ExperimentsPath, myDataFiles[i]), 
            header = TRUE))
        
        # add it to the list
        myFullListOfInteractions <- rbind(myFullListOfInteractions, data)
    }
    
    # save the full list of filtered interactions
    write.table(myFullListOfInteractions, file = paste0(ExperimentsPath, 
        "A_ALL_FILTERED_INTERACTIONS_CYTOSCAPE.csv"), 
        row.names= FALSE, quote = FALSE )
}


# definition of a function to get the IDs of pairs of similar individuals 
#
# INPUTS: DiversePop - population of individuals to compare
#         Threshold - subnetworks over this Jaccard similarity are merged 
# OUTPUT: List of pairs of similar individuals
GetSimilarInds <- function(DiversePop, Threshold) {
    # create all the combinations of 2 numbers with the ids of individuals
    MyIndexes <- t(combn(seq_len(nrow(DiversePop)), 2))
    
    # create data frame to save the similarities
    Similarities <- data.frame(I1 = MyIndexes[, 1], I2 = MyIndexes[, 2], JS = 0)
    
    # calculate the Jaccard similarity index
    Similarities$JS <- apply(MyIndexes, 1, function(X) { 
        Ind1 <- as.character(unlist(DiversePop$Individual[X[1]]))
        Ind2 <- as.character(unlist(DiversePop$Individual[X[2]]))
        
        # Jaccard similarity = (intersection of A and B) / (union of A and B)
        JS <- (length(intersect(Ind1, Ind2))) / (length(union(Ind1, Ind2)))*100
        return(as.double(JS)) } )
    
    # keep only the rows of "duplicated" individuals
    Similarities <- Similarities[Similarities$JS >= Threshold, ]
    
    return(Similarities)    
}



# definition of a function to join those individuals with a Jaccard similatiry
# coefficient superior to a given threshold
#
# INPUTS: AccPF - filtered accumulated Pareto front
#         LoadedData - list containing several important variables 
#                      (see mogamun_load_data())    
#         Threshold - subnetworks over this Jaccard similarity are merged 
# OUTPUT: None
JoinSimilarInds <- function(AccPF, LoadedData, Threshold) {
    MaxSize <- LoadedData$MaxSize
    DiversePop <- AccPF # save a copy of the accumulated Pareto front
    
    # flag to deactivate when no individuals-pair fulfill the threshold
    KeepGoing <- TRUE
    
    # loop to join similar individuals
    while(KeepGoing == TRUE) {
        if (nrow(DiversePop) > 1) { # verify that there individuals left
            # get "similar"duplicated" individuals
            Sim <-  GetSimilarInds(DiversePop, Threshold) 
            
            if (nrow(Sim) > 0) { # if there are duplicated inds
                IndsToRemove <- NULL
                i <- 1
                while(i <= nrow(Sim)) { # loop through duplicated inds
                    ID1 <- Sim[i, 1]
                    ID2 <- Sim[i, 2]
                    IndsUnion <- union(DiversePop$Individual[[ID1]], 
                        DiversePop$Individual[[ID2]])
                    
                    # verify that union of individuals respects the max size
                    if (length(IndsUnion) <= MaxSize) {
                        DiversePop$Individual[[ID1]] <- IndsUnion # 
                        IndsToRemove <- c(IndsToRemove, ID2) 
                        
                        # get all future incidences to joined individuals 
                        Ref <- 
                            union(which(Sim == ID1, arr.ind = TRUE)[, "row"], 
                                which(Sim == ID2, arr.ind = TRUE)[, "row"])
                        
                        # check and remove future incidences to joined inds
                        Ref <- Ref[Ref > i]
                        Sim <- Sim[!row.names(Sim) %in% row.names(Sim)[Ref], ]
                    } else {print("Union of individuals overpasses max size")}
                    i <- i + 1
                }
                if (is.null(IndsToRemove)) { # verify if no ind will be removed
                    KeepGoing <- FALSE
                } else {
                    # remove the individuals that were joined with another one
                    DiversePop <- DiversePop[!row.names(DiversePop) %in% 
                        row.names(DiversePop)[IndsToRemove], ]
                }
            } else { KeepGoing <- FALSE }
        } else { KeepGoing <- FALSE }
    }
    return(DiversePop)
}


# definition of a function to save in files the filtered acc Pareto front  
#
# INPUTS: DiversePopulation - population to save in the files
#         ExperimentsPath - folder to save the filtered networks
#         Threshold - threshold used to merge subnetworks 
# OUTPUT: None
SaveFilteredAccPF <- function(DiversePopulation, ExperimentsPath, Threshold) {
    # save postfiltered individuals
    for (i in seq_len(nrow(DiversePopulation))) { 
        Ind <- unlist(DiversePopulation$Individual[i]) # get the inds code
        write(paste(Ind, collapse=" ", sep=""), file = paste0(ExperimentsPath, 
            "A_AccPF_JacSimT_", Threshold, ".csv"), append = TRUE, sep = ",") 
    }
    
    # create data frame with all the nodes in the accumulated Pareto front
    AllNodesDF_Acc <- 
        data.frame(Nodes_Names = unique(unlist(DiversePopulation$Individual)))
    
    # loop through all the postfiltered individuals to create a file
    # to be processed  in cytoscape
    for (i in seq_len(nrow(DiversePopulation))) {
        NodesCurrentInd <- unlist(DiversePopulation$Individual[i])
        
        AllNodesDF_Acc[, (ncol(AllNodesDF_Acc)+1)] <- 
            ifelse(AllNodesDF_Acc[,1] %in% NodesCurrentInd, 
                paste0("Ind", i), "")
        
        colnames(AllNodesDF_Acc)[ncol(AllNodesDF_Acc)] <- paste0("Ind", i)
    }
    
    AllNodesDF_Acc$Any <- "YES"
    
    write.csv(AllNodesDF_Acc, paste0(ExperimentsPath, 
        "A_AccPF_CYTOSCAPE_JacSimT_", Threshold, ".csv"), row.names = FALSE)
}



# definition of a function to visualize the accumulated Pareto front in 
# Cytoscape. NOTE. Cytoscape needs to be open
#
# INPUTS: ExperimentDir - folder containing the results to be processed
#         LoadedData - list containing several important variables 
#                      (see mogamun_load_data())
# OUTPUT: None
CytoscapeVisualization <- function(ExperimentDir, LoadedData) {
    cytoscapePing() # verify that Cytoscape is launched
    closeSession(save.before.closing = FALSE) # close any open session 
    GeneralPath <- ExperimentDir
    
    # get list of directories. Each directory contains the results of one exp
    Dirs <- list.dirs(GeneralPath, recursive = FALSE, full.names = FALSE)
    Dirs <- Dirs[grep("Experiment", Dirs)]
    
    for (d in Dirs) { # loop through all the experiments
        ExpPath <- paste0(GeneralPath, "/", d, "/")
        
        # read network
        Network <- read.csv( paste0(ExpPath, 
            "A_ALL_FILTERED_INTERACTIONS_CYTOSCAPE.csv"), sep = " ",
            stringsAsFactors = FALSE )
        
        # change column names to match the requirements for Cytoscape
        colnames(Network) <- c("source", "target", "interaction")
        
        # get list of nodes
        Nodes <- data.frame(id = unique(c(Network$source, Network$target)))
        
        # generate network in cytoscape
        createNetworkFromDataFrames(Nodes, Network, title = d)
        
        # read expression data
        DE <- LoadedData$DE_results
        DE$DEG <- ifelse(abs(DE$logFC) > 1 & DE$FDR < 0.05, TRUE, FALSE)
        
        # remove rows without gene name and filter columns
        DE <- DE[!is.na(DE$gene), ]
        
        # load expression data in cytoscape
        loadTableData(data = DE, data.key.column = "gene", table = "node",
            table.key.column = "name", namespace = "default", network = d)
        
        # adds nodes' color and borders
        FormatNodesAndEdges(Network, d, LoadedData) 
        
        # creates the subnetworks corresponding the active modules
        CreateActiveModules(d, ExpPath)
        
        # verify if there are more exp to analyze to leave last one open
        if (which(Dirs == d) < length(Dirs)) {
            closeSession(save.before.closing = TRUE, 
                filename = paste0(ExpPath, "A_Acc_PF_", d))
        } else { saveSession(filename = paste0(ExpPath, "A_Acc_PF_", d)) }
    }
}


# definition of a function to format the nodes and edges in Cytoscape
# Cytoscape. NOTE. Cytoscape needs to be open
#
# INPUTS: Network - list of edges
#         d - Directory containing the results of one experiment
#         LoadedData - list containing several important variables 
#                      (see mogamun_load_data())
# OUTPUT: None
FormatNodesAndEdges <- function(Network, d, LoadedData) {
    # read expression data
    DE <- LoadedData$DE_results
    DE$DEG <- ifelse(abs(DE$logFC) > 1 & DE$FDR < 0.05, TRUE, FALSE)
    
    # remove rows without gene name and filter columns
    DE <- DE[!is.na(DE$gene), ]
    
    numberOfLayers <- length(LoadedData$DensityPerLayerMultiplex)
    
    # if there are 3 layers, use the same edges' colors as in the paper
    if(numberOfLayers == 3) {
        edgesColors <- c("#0033FF", "#FF6600", "#FFFF00")
    } else { # otherwise, generate a list of colors of the appropriate size
        edgesColors <- rainbow(length(LoadedData$DensityPerLayerMultiplex)) 
        
        # remove "transparency" from colors (i.e., the last 2 characters)
        edgesColors <- 
            unlist(lapply(edgesColors, function(X){ substr(X, 1, 7) }))
    }
    
    # define edges color (it is adapted to any number of layers)
    setEdgeColorMapping(table.column = "interaction", mapping.type = "d",
        table.column.values = unique(Network$interaction), colors = edgesColors)
    
    # set nodes' colors according to the logFC, from green (downregulated) 
    # to white and then to red (upregulated)
    setNodeColorMapping(table.column = "logFC",
        table.column.values = c(min(DE$logFC), 0.0, max(DE$logFC)),
        colors = c("#009933", "#FFFFFF", "#FF0000"), mapping.type = "c", 
        style.name = "default", network = d)
    
    setNodeBorderColorMapping(table.column = "name", colors = "#000000",
        mapping.type = "d", style.name = "default", default.color = "#000000",
        network = d)
    
    # set width of border to highlight DEG
    setNodeBorderWidthMapping(table.column = "DEG",
        table.column.values = c(TRUE, FALSE), widths = c(5, 0), 
        mapping.type = "d", style.name = "default", network = d)
    
}



# definition of a function that creates the subnetworks corresponding the  
# active modules from the accumulated Pareto front. 
# NOTE. Cytoscape needs to be open
#
# INPUTS: d - Directory containing the results of one experiment
#         ExperimentsPath - path to the results
# OUTPUT: None
CreateActiveModules <- function(d, ExperimentsPath) {
    
    # set layout 
    layoutNetwork(layout.name = "force-directed", network = d)
    
    # load data for active modules
    ActiveModules <- 
        data.frame(read.csv(paste0(ExperimentsPath, list.files(ExperimentsPath,
            pattern = "A_AccPF_CYTOSCAPE_JacSimT_"))))
    
    # load active modules data in cytoscape
    loadTableData(data = ActiveModules, data.key.column = "Nodes_Names",
        table = "node", table.key.column = "name", namespace = "default", 
        network = d)
    
    # create a subnetwork for each active module
    lapply( 
        colnames(ActiveModules)[grep("Ind", colnames(ActiveModules))], 
        function(am){ 
            createColumnFilter(filter.name = "ActiveModules",  column = am,  
                criterion = am, predicate = "IS", type = "node", network = d)
                                
            # create the subnetwork with the selected nodes
            createSubnetwork(subnetwork.name = paste0("ActiveModule_", am))
            
            # set layout 
            layoutNetwork(layout.name = "force-directed", 
                network = paste0("ActiveModule_", am))
        } 
    )
}


# this function is a part of run_mogamun. It is executed if the OS is not 
# Windows and it allows running MOGAMUN in parallel
RunInLinux <- function(LoadedData, Cores, NumberOfRunsToExecute, ResultsDir) {
    doParallel::registerDoParallel(cores = Cores) # in line with the phys. cores
    PopSize <- LoadedData$PopSize
    Generations <- LoadedData$Generations
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
}


# this function is a part of run_mogamun. It is executed if the OS is  
# Windows and it DOES NOT allow running MOGAMUN in parallel
RunInWindows <- function(LoadedData, NumberOfRunsToExecute, ResultsDir) {
    PopSize <- LoadedData$PopSize
    Generations <- LoadedData$Generations
    Multiplex <- LoadedData$Multiplex
    ResultsPath <- paste0(ResultsDir, "/Experiment_", Sys.Date(), "/")
    dir.create(ResultsPath, recursive = TRUE)  # create result folder 
    BestIndsPath <- paste0(ResultsPath, "MOGAMUN_Results_") # path for res
    
    RunNumber <- 1
    
    # loop to execute the algorithm many times
    for (RunNumber in seq_len(NumberOfRunsToExecute)){
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
}


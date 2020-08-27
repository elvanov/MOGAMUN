![R-CMD-check](https://github.com/elvanov/MOGAMUN/workflows/R-CMD-check/badge.svg)

# MOGAMUN - A Multi-Objective Genetic Algorithm to Find Active Modules in Multiplex Biological Networks

## Description

A package to find active modules (i.e., highly connected subnetworks with an overall deregulation) in multiplex biological networks. 

MOGAMUN is described in detail in Novoa-del-Toro, E.M., et al.: *A Multi-Objective Genetic Algorithm to Find Active Modules in Multiplex Biological Networks*, https://www.biorxiv.org // COMPLETE URL TO PAPER

### How to install

If you wish to run MOGAMUN **sequentially** (i.e., one run at a time), you can install the package by running the following commands under R console or RStudio:

```R
> library("devtools")
> install_github("elvanov/MOGAMUN")
> library(MOGAMUN)
```

To run multiple executions of MOGAMUN in **parallel**, we strongly recommend you to download the GitHub repository to your local computer, and edit four lines from *R/MOGAMUN.R*, as follows. Uncomment lines 169 and 170:

```R
library(doParallel)
registerDoParallel(cores = 1)
```
**Important note.** Make sure to specify the number of cores you want MOGAMUN to use (one core is needed per run).

Uncomment line 189 and comment out line 190, as follows:

```R
foreach(RunNumber = 1:NumberOfRunsToExecute) %dopar% {
# for (RunNumber in 1:NumberOfRunsToExecute) {
```

### Usage

MOGAMUN uses two information sources: one or more biological networks and the statistical values resulting from a differential expression analysis or any other test that gives as result *p*-values or False Discovery Rates (FDR) associated to genes. 

The first step (out of three) to run MOGAMUN is to set the values for all parameters (e.g. those related to the evolution process, such as crossover and mutation rates), using the `mogamun.init` command. To see the full list of parameters, see `?mogamun.init`. 

The second step is to provide the input data (including the path to the biological networks), using the `mogamun.load.data` command. Please note that each biological network must be in a separate file with 2-column format (separated by tabs), and we strongly recommend to replace the dashes with underscores. To see the full list of parameters, see `?mogamun.load.data`.

The last step is to run MOGAMUN, using the `mogamun.run` command, giving as parameter the path to the folder where the results will be stored. See `?mogamun.run` for further details.

Be aware that the runnig time of for 500 generations was of approximately 12 hours, in a desk computer with Intel processor i7 at 3.60GHz and 32GB of RAM.

**Important note.** You can find a full set of example files to run MOGAMUN in [MOGAMUN data github page](https://github.com/elvanov/MOGAMUN-data)


### Interpreting the results

In the results' path that you specified in `mogamun.run`, you will find two files per run (*MOGAMUN_Results_StatisticsPerGeneration_RunN.csv* and *MOGAMUN_Results__Run_N.txt*). The file *MOGAMUN_Results_StatisticsPerGeneration_RunN.csv* contains the best values of the average nodes score and density per generation (you can use them to check the convergence), and the file *MOGAMUN_Results__Run_N.txt* contains the complete final population (i.e. the subnetworks from the last generation), one per row. The number of elements in every row is variable. If *X_n* is the number of elements in the *n*-th row: the nodes of the subnetwork are the first *X_n*-4 elements. The last four elements are the average nodes score, the density, the rank and the crowding distance, respectively. The best (non-dominated) subnetworks have are those with rank = 1. 



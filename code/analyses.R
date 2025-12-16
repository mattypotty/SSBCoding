## ----setup, echo=FALSE--------------------------------------------------------
# R markdown is used to create documents (in this case, an HTML document) out of R code and output
# Code is in "chunks" that can be run, everything else is text
# This code chunk sets up how the HTML file will be created from this markdown document when we "knit" it
# By default, we want no code to run, no code to be shown, and no output to be shown
# Then for each code chunk, we will decide whether to run it, show it, and show its results

knitr::opts_chunk$set(eval=FALSE,echo=FALSE,results=FALSE,message=FALSE,warning=FALSE,tidy='styler')
# eval    - whether the code will be run (evaluated) when we create our HTML file
# echo    - whether the code will be shown (printed) when we create our HTML file
# results - whether the results of the code will be shown (printed) when we create our HTML file
# message - whether the messages output by R will be shown (printed) when we create our HTML file
# warning - whether the warnings output by R will be shown (printed) when we create our HTML file


## ----packages, eval=TRUE------------------------------------------------------
# Some functions are available in base R without installing any packages
# However, packages (which are built by outside individuals) provide additional functions
# All of the packages for this project are accessible on CRAN, a repository run by the R dev team
# This code chunk is for loading different packages so we can access their functions
# Note that we want to run this code (eval=TRUE) but we don't need to show it in our document
print("Loading packages...",quote=F)

library(this.path) # This package lets us get the filepath of the script we are running
library(ape) # Used for handling trees
library(phangorn) # Used for making MCC (maximum clade credibility) tree from posterior of trees
library(castor) # Used for checking if tree is bifurcating
library(caper) # Used for phylogenetic signal (d) analyses
library(phytools) # Used for plotting trees
library(diversitree) # Used for SSE analyses
library(MCMCglmm) # Used for MCMCglmm analyses
library(coda) # Used for analyzing MCMC output
library(MCMCvis) # Used for analyzing MCMC output
library(mcmcr) # Used for analyzing MCMC output
library(styler) # For knitting


## ----directories, eval=TRUE---------------------------------------------------
# This code chunk sets the working directory for the analyses and creates some file paths for saving output
# This can be done by hand too:
# clicking session > setting working directory > to source file location

setwd(this.path::here()) # Sets the working directory to the main SSB folder (using the "this.path" package)
print(paste0("Working directory: ",getwd()),quote=F)


## ----analysis_settings, eval=TRUE---------------------------------------------
# We will use ALL CAPS for our analysis settings
# This way, we can look at a variable and know it is a setting
print("Setting up analysis...",quote=F)

USEGLOBAL <- FALSE

args <- commandArgs(trailingOnly = TRUE)
print("Input arguments:",quote=F)
print(args)

if (length(args) != 0) {
  print("Using input arguments to set up analysis...",quote=F)
  USEGLOBAL <- TRUE
  .NUMBER <- args[1]
  .USERANDOM <- args[2]
  .ANALYSIS <- args[3]
  .SSBTYPE <- args[4]
  .HYPOTHESIS <- args[5]
  .VARIABLE <- args[6]
}

# We'll use the variable USEGLOBAL to indicate whether we want to use inputs from shell script
set_global <- function() {
  if (exists(".NUMBER")) {NUMBER <<- .NUMBER}
  if (exists(".ANALYSIS")) {ANALYSIS <<- .ANALYSIS}
  if (exists(".SSBTYPE")) {SSBTYPE <<- .SSBTYPE}
  if (exists(".HYPOTHESIS")) {HYPOTHESIS <<- .HYPOTHESIS}
  if (exists(".VARIABLE")) {VARIABLE <<- .VARIABLE}
  if (exists(".USERANDOM")) {USERANDOM <<- .USERANDOM}
}

if (USEGLOBAL) {set_global()}

if (!exists("NUMBER")) {NUMBER <- 001} # Allows us to run multiple analyses of the same type, each with a different number
set.seed(NUMBER) # For analyses requiring random numbers, we might want the same output every time

if (!exists("USERANDOM")) {USERANDOM <- FALSE} # Whether we use a random tree from the set instead of the MCC tree
# OPTIONS: FALSE, TRUE

if (!exists("ANALYSIS")) {ANALYSIS <- "NONE"}
# OPTIONS: "NONE", "MCCTREE", "SPECMATS", "MCMCGLMM", "PAGEL"

if (!exists("SSBTYPE")) {SSBTYPE <- "NONE"}
# OPTIONS: "NONE", "SSB", "FSSB", "MSSB"

if (!exists("HYPOTHESIS")) {HYPOTHESIS <- "NONE"}
# OPTIONS: "NONE", "MI", "MGS", "SDT", "RDM", "AM", "IUCN"

if (!exists("VARIABLE")) {VARIABLE <- "NONE"}
# OPTIONS: "NONE", "NAC", "SPI", "TSP", "SDT", "RDM"

# MCMC Settings
NITT <- 100000 # 100,000
BURNIN <- 10000 # 10,000
THIN <- 100 # 100
prior <- list(G=list(G1=list(V=1,nu=0.02)),
              R=list(V=1,nu=0.02))


## ----do_mcmc, echo=TRUE, eval=TRUE--------------------------------------------
do_mcmc <- function(f) {
  # File names
  analysis_name <- paste0(ANALYSIS,"_",SSBTYPE,"_",HYPOTHESIS,"_",tree_number,"_",NUMBER)
  print(paste0("Performing analysis: ",analysis_name),quote=F)
  out_file <- paste0("../output/",analysis_name,".log.csv")
  summary_file <- paste0("../output/",analysis_name,".sum.txt")
  
  #Read in data
  columns <- c("Species",all.vars(f))
  data_subset <- data[,columns]
  # Removes incomplete cases
  data_subset <- data_subset[complete.cases(data_subset),]
  
  # Running MCMCglmm
  MCMCanalysis <- MCMCglmm::MCMCglmm(f,
                  random = ~Species,
                  family = "categorical",
                  ginverse = list(Species=inv.phylo),
                  prior = prior,
                  data = data_subset,
                  nitt = NITT,
                  burnin = BURNIN,
                  thin = THIN,
                  verbose = TRUE)
  
  write.csv(MCMCanalysis$Sol,out_file,row.names=FALSE,quote=FALSE)
  out <- capture.output(cat(paste0("ESS: ",ess(MCMCanalysis$Sol),"\n")),
                        print(summary(MCMCanalysis$Sol)))
  write(out,summary_file)
}


## ----do_pagel, echo=TRUE, eval=TRUE-------------------------------------------
do_pagel <- function() {
  # File names
  analysis_name <- paste0(ANALYSIS,"_",SSBTYPE,"_",VARIABLE,"_",tree_number,"_",NUMBER)
  print(paste0("Performing analysis: ",analysis_name),quote=F)
  summary_file <- paste0("../output/",analysis_name,".sum.txt")
  
  #Read in data
  columns <- c("Species",SSBTYPE,VARIABLE)
  data_subset <- data[,columns]
  # Removes incomplete cases
  data_subset <- data_subset[complete.cases(data_subset),]
  phy_subset <- ape::keep.tip(phy, data_subset$Species)
  
  #Assign species names
  SSB <- setNames(data_subset[[SSBTYPE]],data_subset$Species)
  VAR <- setNames(data_subset[[VARIABLE]],data_subset$Species)
  
  write(paste0("Pagel's Directional Test for ",SSBTYPE," and ",VARIABLE),summary_file)
  write("--------------------",summary_file,append=TRUE)

  # Fit models where SSB and NAC are totally independent or dependent
  fit_SSB_VAR <- phytools::fitPagel(phy_subset,SSB,VAR)
  out1 <- capture.output(print(fit_SSB_VAR))
  write(paste0("Both Independent or Dependent"),summary_file,append=TRUE)
  write(out1,summary_file,append=TRUE)
  write("--------------------",summary_file,append=TRUE)
  
  # Fit model where SSB depends on NAC
  fit_SSB <- phytools::fitPagel(phy_subset,SSB,VAR,dep.var="x")
  out2 <- capture.output(print(fit_SSB))
  write(paste0("Dependent ",SSBTYPE),summary_file,append=TRUE)
  write(out2,summary_file,append=TRUE)
  write("--------------------",summary_file,append=TRUE)
  
  # Fit model where NAC depends on SSB
  fit_VAR <- phytools::fitPagel(phy_subset,SSB,VAR,dep.var="y")
  out3 <- capture.output(print(fit_VAR))
  write(paste0("Dependent ",VARIABLE),summary_file,append=TRUE)
  write(out3,summary_file,append=TRUE)
  write("--------------------",summary_file,append=TRUE)
  
  # Comparing the goodness of all 4 models using AIC
  aic <- setNames(c(fit_SSB_VAR$independent.AIC,
                    fit_SSB$dependent.AIC,
                    fit_VAR$dependent.AIC,
                    fit_SSB_VAR$dependent.AIC),
                  c("independent","dependent SSB",paste0("dependent ",VARIABLE),paste0("dependent SSB & ",VARIABLE)))
  out4 <- capture.output(print(aic))
  out5 <- capture.output(print(aic.w(aic)))
  write(paste0("AIC"),summary_file,append=TRUE)
  write(out4,summary_file,append=TRUE)
  write(paste0("Weights"),summary_file,append=TRUE)
  write(out5,summary_file,append=TRUE)
}


## ----dataset, eval=TRUE, results="SHOW"---------------------------------------
# This code chunk reads in our data files
# Note that we want to run this code (eval=TRUE) and show the output (results="SHOW"), but we don't need to see it
print("Reading in data...",quote=F)

data_file <- read.csv("../data/FINALdata_20251124.csv", header = TRUE, sep = ",") # Reads the data, which is a csv file (comma separated values), and saves to "data" variable
#print(data_file$Column.Name.Key[data_file$Column.Name.Key != ""]) # Prints the column name key
#print(data_file$Mating.Group.Structures.Key[data_file$Mating.Group.Structures.Key != ""]) # prints the group mating structures key

data_full <- data.frame(data_file[,-which(colnames(data_file)=="Column.Name.Key" | colnames(data_file)=="Mating.Group.Structures.Key")]) # Removing keys from dataset

# Fixing any problems with the dataset
data_full$Species[which(data_full$Species=="Otaria flavescens")] <- "Otaria_flavescens"

# Removing species which are not in the tree
phylogeny_changes <- read.csv("../data/Dataset_to_Phylogeny_Changes.csv", header = TRUE, sep = ",")[,1:2]
for (i in 1:nrow(phylogeny_changes)) {
  if (phylogeny_changes$Name.used.for.Phylacine[i] == "NLV") {
    name <- phylogeny_changes$X1700.Species.Name[i]
    data_full <- data_full[data_full$Species != name,]
  }
}

SOC <- rep(NA,nrow(data_full))
for (i in 1:nrow(data_full)) {
  if (data_full$GL[i] == "N") {SOC[i] <- "SOL"}
  if (data_full$GL[i] == "Y") {SOC[i] <- "GL"}
  if (data_full$SOG[i] == "Y") {SOC[i] <- "SOG"}
}
data_full$SOC <- SOC

quantitative <- c("GP","BM","AFM","MR","PC","FWA","FAM","MAM") # quantitative columns
binary <- c("RDM","TSP","NB","SDT","GL","SPI","NAC","GIM","GFE","MMC","BG","DHE","DCA","DPMI","DIN","DFR","DGR","DCE","DMOC","SOG","FSSB","MSSB","SSB") # binary columns
yn <- c("RDM","TSP","NB","SDT","GL","SPI","NAC","GIM","GFE","MMC","BG","DHE","DCA","DPMI","DIN","DFR","DGR","DCE","DMOC","SOG") # binary columns with y/n
multiple <- c("MGS","IUCN","Modifiers", "SOC") # categorical columns with multiple options
categorical <- c(binary,multiple) # all categorical columns

for (i in 1:ncol(data_full)) { # looping through all the columns to format each one correctly
  col <- colnames(data_full)[i] # getting the name of the column we're working on
  data_full[[col]] <- gsub("#N/A","",data_full[[col]]) # replacing any of the "#N/A" generated by excel with ""
  data_full[[col]] <- gsub("^$",NA,data_full[[col]]) # replacing the blank "" with NA
  if (col %in% quantitative) {data_full[[col]] <- as.numeric(data_full[[col]])} # formatting quantitative columns as numeric
  if (col %in% yn) {data_full[[col]] <- gsub("Y","1",data_full[[col]]); data_full[[col]] <- gsub("N","0",data_full[[col]])} # swapping y/n for 1/0
  if (col %in% categorical) {data_full[[col]] <- as.factor(data_full[[col]])} # formatting all categorical columns as factors
}


## ----make_tree, echo=TRUE-----------------------------------------------------
# This code chunk reads our giant list of trees and combines them to make one "best" tree
# We want to see this code (echo=TRUE), but we don't want to run it every time we do our analyses (it takes a while)
ANALYSIS <- "MCMCTREE"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "MCCTREE") {
  print("Generating MCC tree...",quote=F)
  trees <- ape::read.nexus("../data/phylogeny.nex")
  phy_mcc <- phangorn::maxCladeCred(trees) # Generates a maximum clade credibility (best) tree from the 1000 trees
  cat(paste0("Is Ultrametric? ", ape::is.ultrametric(phy_mcc)))
  cat(paste0("Is Bifurcating? ", castor::is_bifurcating(phy_mcc)))
  n_taxa <- length(phy_mcc$tip.label) # Gets the number of species from looking at the tip labels
  phy_mcc$node.label <- c((n_taxa + 1):(n_taxa + phy_mcc$Nnode)) # Gives number labels to nodes
  ape::write.tree(phy_mcc,"../data/mcc_tree.txt") # Saves this tree to a file called "mcc_tree.txt"
}


## ----read_tree, eval=TRUE-----------------------------------------------------
# This code chunk reads in the full set of trees
# It also reads in our MCC tree (one best tree)

print("Reading trees...",quote=F)
phy_mcc <- ape::read.tree("../data/mcc_tree.txt") # Reads our MCC tree from a file
n_taxa <- length(phy_mcc$tip.label)
n_node <- phy_mcc$Nnode

phy_all <- ape::read.nexus("../data/phylogeny.nex") # Reads our posterior sample of trees from a file
# Giving number labels to nodes in all 1000 phylogenies
for (i in 1:length(phy_all)) {
  phy_all[[i]]$node_label <- c((n_taxa + 1):(n_taxa + n_node))
}


## ----phylogeny, eval=TRUE-----------------------------------------------------
# This code chunk chooses a phylogeny to use with our current analysis
# There are two options: use the MCC tree, or use a random tree

print("Choosing tree...",quote=F)
phy_full <- phy_mcc # This is the default behavior
tree_number <- "TREEMCC"
# If we want a random tree instead, we do this
if (USERANDOM == TRUE) {
  num <- sample(1:1000,1)
  phy_full <- phy_all[[num]]
  tree_number <- paste0("TREE",num)
}

# Renaming or removing Phylacine species that do not match our dataset
phylogeny_changes <- read.csv("../data/Dataset_to_Phylogeny_Changes.csv", header = TRUE, sep = ",")[,1:2]
for (i in 1:nrow(phylogeny_changes)) {
  if (phylogeny_changes$Name.used.for.Phylacine[i] != "NLV") {
    new_name <- phylogeny_changes$X1700.Species.Name[i]
    old_name <- phylogeny_changes$Name.used.for.Phylacine[i]
    phy_idx <- which(phy_full$tip.label==old_name)
    phy_full$tip.label[phy_idx] <- new_name
  }
}
print(paste0("Chosen tree: ",tree_number),quote=F)


## ----get_mismatches-----------------------------------------------------------
# This code checks for all the species that are in the data but not in the tree
print("Checking for mismatched species...",quote=F)
ANALYSIS <- "SPECMATS"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "SPECMATS") {
  LON <- data_full$Species # Reduces species from my data to LON ("List of Names")
  check <- LON[which(!(LON %in% phy_full$tip.label))] # Lists species which are in mine that are not in theirs (species to check)
  print(paste0("Species to check: ",paste(check,collapse=", ")),quote=F) # Prints out the list of species that need to be checked
}


## ----specmats, eval=TRUE, results="SHOW"--------------------------------------
# This code chunk makes sure our species data and our phylogeny have the same taxa in them
# We will use this smaller tree (phy) for most of our analyses

print("Getting tree subset...",quote=F)
LON <- data_full$Species # Reduces species from my data to LON ("List of Names")
SPECMATS <- LON[which((LON %in% phy_full$tip.label))] # Creates a variable SPECMATS for all the species of mine which match those in the full trees
phy <- ape::keep.tip(phy_full, SPECMATS) # Creates a new tree reducing phy to just SPECMATS, can then plot this with plot.phylo(SPECMATSTREE)
print(paste0("Total species: ",length(phy$tip.label)),quote=F) # Prints out the number of species in the tree

# Read tree and verify ultrametric, form matrix for evolution time
inv.phylo <- inverseA(phy, "TIPS")$Ainv

# THIS SHOULD NOT BE REQUIRED ONCE ALL OF THE NAMES MATCH!!!
data <- data_full[which((LON %in% phy$tip.label)),] # Creates dataset of SPECMATS (rows that match with the tree).


## ----plot_tree, eval=TRUE, results="SHOW", fig.align="center", fig.width=12, fig.height=12----
phytools::plotTree(phy, type="fan", fsize=.5, asp=1) # This plots the tree that we are using (with only matching taxa) as a circle phylogeny


## ----MI_mcmcglmm, echo=TRUE---------------------------------------------------
#Runs MCMCglmm for all our Maternal investment variables. Individual binary variables each get their own suite of pagels below.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "MI"; SSBTYPE <- "SSB"; if (USEGLOBAL) {set_global()}
f <- SSB ~ NAC + PC + SPI + GP + FWA + MR + TSP + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "MI" & SSBTYPE == "SSB") {do_mcmc(f)}




## ----MIf_mcmcglmm, echo=TRUE--------------------------------------------------
#Runs MCMCglmm for all our Maternal investment variables. Individual binary variables each get their own suite of pagels below.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "MI"; SSBTYPE <- "FSSB"; if (USEGLOBAL) {set_global()}
f <- FSSB ~ NAC + PC + SPI + GP + FWA + MR + TSP + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "MI" & SSBTYPE == "FSSB") {do_mcmc(f)}




## ----MIm_mcmcglmm, echo=TRUE--------------------------------------------------
#Runs MCMCglmm for all our Maternal investment variables. Individual binary variables each get their own suite of pagels below.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "MI"; SSBTYPE <- "MSSB"; if (USEGLOBAL) {set_global()}
f <- MSSB ~ NAC + PC + SPI + GP + FWA + MR + TSP + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "MI" & SSBTYPE == "MSSB") {do_mcmc(f)}




## ----NAC_pagel, echo=TRUE-----------------------------------------------------
ANALYSIS <- "PAGEL"; VARIABLE <- "NAC"; SSBTYPE <- "SSB"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "PAGEL" & VARIABLE == "NAC" & SSBTYPE == "SSB") {do_pagel()}




## ----NACf_pagel, echo=TRUE----------------------------------------------------
ANALYSIS <- "PAGEL"; VARIABLE <- "NAC"; SSBTYPE <- "FSSB"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "PAGEL" & VARIABLE == "NAC" & SSBTYPE == "FSSB") {do_pagel()}




## ----NACm_pagel, echo=TRUE----------------------------------------------------
ANALYSIS <- "PAGEL"; VARIABLE <- "NAC"; SSBTYPE <- "MSSB"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "PAGEL" & VARIABLE == "NAC" & SSBTYPE == "MSSB") {do_pagel()}




## ----SPI_pagel, echo=TRUE-----------------------------------------------------
ANALYSIS <- "PAGEL"; VARIABLE <- "SPI"; SSBTYPE <- "SSB"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "PAGEL" & VARIABLE == "SPI" & SSBTYPE == "SSB") {do_pagel()}




## ----SPIf_pagel, echo=TRUE----------------------------------------------------
ANALYSIS <- "PAGEL"; VARIABLE <- "SPI"; SSBTYPE <- "FSSB"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "PAGEL" & VARIABLE == "SPI" & SSBTYPE == "FSSB") {do_pagel()}




## ----SPIm_pagel, echo=TRUE----------------------------------------------------
ANALYSIS <- "PAGEL"; VARIABLE <- "SPI"; SSBTYPE <- "MSSB"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "PAGEL" & VARIABLE == "SPI" & SSBTYPE == "MSSB") {do_pagel()}




## ----TSP_pagel, echo=TRUE-----------------------------------------------------
ANALYSIS <- "PAGEL"; VARIABLE <- "TSP"; SSBTYPE <- "SSB"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "PAGEL" & VARIABLE == "TSP" & SSBTYPE == "SSB") {do_pagel()}




## ----TSPf_pagel, echo=TRUE----------------------------------------------------
ANALYSIS <- "PAGEL"; VARIABLE <- "TSP"; SSBTYPE <- "FSSB"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "PAGEL" & VARIABLE == "TSP" & SSBTYPE == "FSSB") {do_pagel()}




## ----TSPm_pagel, echo=TRUE----------------------------------------------------
ANALYSIS <- "PAGEL"; VARIABLE <- "TSP"; SSBTYPE <- "MSSB"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "PAGEL" & VARIABLE == "TSP" & SSBTYPE == "MSSB") {do_pagel()}




## ----MGS_mcmcglmm, echo=TRUE--------------------------------------------------
#Runs MCMCglmm for all our Sociality variables.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "MGS"; SSBTYPE <- "SSB"; if (USEGLOBAL) {set_global()}
f <- SSB ~ MGS + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "MGS" & SSBTYPE == "SSB") {do_mcmc(f)}




## ----MGSf_mcmcglmm, echo=TRUE-------------------------------------------------
#Runs MCMCglmm for all our Sociality variables.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "MGS"; SSBTYPE <- "FSSB"; if (USEGLOBAL) {set_global()}
f <- FSSB ~ MGS + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "MGS" & SSBTYPE == "FSSB") {do_mcmc(f)}




## ----MGSm_mcmcglmm, echo=TRUE-------------------------------------------------
#Runs MCMCglmm for all our Sociality variables.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "MGS"; SSBTYPE <- "MSSB"; if (USEGLOBAL) {set_global()}
f <- MSSB ~ MGS + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "MGS" & SSBTYPE == "MSSB") {do_mcmc(f)}




## ----SDT_mcmcglmm, echo=TRUE--------------------------------------------------
#Runs MCMCglmm for all our Territoriality variables. Individual binary variables each get their own suite of pagels below.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "SDT"; SSBTYPE <- "SSB"; if (USEGLOBAL) {set_global()}
f <- SSB ~ SDT + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "SDT" & SSBTYPE == "SSB") {do_mcmc(f)}




## ----SDTf_mcmcglmm, echo=TRUE-------------------------------------------------
#Runs MCMCglmm for all our Territoriality variables. Individual binary variables each get their own suite of pagels below.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "SDT"; SSBTYPE <- "FSSB"; if (USEGLOBAL) {set_global()}
f <- FSSB ~ SDT + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "SDT" & SSBTYPE == "FSSB") {do_mcmc(f)}




## ----SDTm_mcmcglmm, echo=TRUE-------------------------------------------------
#Runs MCMCglmm for all our Territoriality variables. Individual binary variables each get their own suite of pagels below.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "SDT"; SSBTYPE <- "MSSB"; if (USEGLOBAL) {set_global()}
f <- MSSB ~ SDT + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "SDT" & SSBTYPE == "MSSB") {do_mcmc(f)}




## ----SDT_pagel, echo=TRUE-----------------------------------------------------
ANALYSIS <- "PAGEL"; VARIABLE <- "SDT"; SSBTYPE <- "SSB"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "PAGEL" & VARIABLE == "SDT" & SSBTYPE == "SSB") {do_pagel()}




## ----SDTf_pagel, echo=TRUE----------------------------------------------------
ANALYSIS <- "PAGEL"; VARIABLE <- "SDT"; SSBTYPE <- "FSSB"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "PAGEL" & VARIABLE == "SDT" & SSBTYPE == "FSSB") {do_pagel()}




## ----SDTm_pagel, echo=TRUE----------------------------------------------------
ANALYSIS <- "PAGEL"; VARIABLE <- "SDT"; SSBTYPE <- "MSSB"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "PAGEL" & VARIABLE == "SDT" & SSBTYPE == "MSSB") {do_pagel()}




## ----SOC_mcmcglmm, echo=TRUE--------------------------------------------------
#Runs MCMCglmm for all our Territoriality variables. Individual binary variables each get their own suite of pagels below.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "SOC"; SSBTYPE <- "SSB"; if (USEGLOBAL) {set_global()}
f <- SSB ~ SOC + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "SOC" & SSBTYPE == "SSB") {do_mcmc(f)}




## ----SOCf_mcmcglmm, echo=TRUE-------------------------------------------------
#Runs MCMCglmm for all our Territoriality variables. Individual binary variables each get their own suite of pagels below.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "SOC"; SSBTYPE <- "FSSB"; if (USEGLOBAL) {set_global()}
f <- FSSB ~ SOC + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "SOC" & SSBTYPE == "FSSB") {do_mcmc(f)}




## ----SOCm_mcmcglmm, echo=TRUE-------------------------------------------------
#Runs MCMCglmm for all our Territoriality variables. Individual binary variables each get their own suite of pagels below.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "SOC"; SSBTYPE <- "MSSB"; if (USEGLOBAL) {set_global()}
f <- MSSB ~ SOC + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "SOC" & SSBTYPE == "MSSB") {do_mcmc(f)}




## ----RDM_mcmcglmm, echo=TRUE--------------------------------------------------
#Runs MCMCglmm for all our Confusion variables. Individual binary variables each get their own suite of pagels below.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "RDM"; SSBTYPE <- "SSB"; if (USEGLOBAL) {set_global()}
f <- SSB ~ RDM + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "RDM" & SSBTYPE == "SSB") {do_mcmc(f)}




## ----RDMf_mcmcglmm, echo=TRUE-------------------------------------------------
#Runs MCMCglmm for all our Confusion variables. Individual binary variables each get their own suite of pagels below.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "RDM"; SSBTYPE <- "FSSB"; if (USEGLOBAL) {set_global()}
f <- FSSB ~ RDM + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "RDM" & SSBTYPE == "FSSB") {do_mcmc(f)}




## ----RDMm_mcmcglmm, echo=TRUE-------------------------------------------------
#Runs MCMCglmm for all our Confusion variables. Individual binary variables each get their own suite of pagels below.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "RDM"; SSBTYPE <- "MSSB"; if (USEGLOBAL) {set_global()}
f <- MSSB ~ RDM + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "RDM" & SSBTYPE == "MSSB") {do_mcmc(f)}




## ----RDM_pagel, echo=TRUE-----------------------------------------------------
ANALYSIS <- "PAGEL"; VARIABLE <- "RDM"; SSBTYPE <- "SSB"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "PAGEL" & VARIABLE == "RDM" & SSBTYPE == "SSB") {do_pagel()}




## ----RDMf_pagel, echo=TRUE----------------------------------------------------
ANALYSIS <- "PAGEL"; VARIABLE <- "RDM"; SSBTYPE <- "FSSB"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "PAGEL" & VARIABLE == "RDM" & SSBTYPE == "FSSB") {do_pagel()}




## ----RDMm_pagel, echo=TRUE----------------------------------------------------
ANALYSIS <- "PAGEL"; VARIABLE <- "RDM"; SSBTYPE <- "MSSB"; if (USEGLOBAL) {set_global()}
if (ANALYSIS == "PAGEL" & VARIABLE == "RDM" & SSBTYPE == "MSSB") {do_pagel()}




## ----AM_mcmcglmm, echo=TRUE---------------------------------------------------
#Runs MCMCglmm for all our Maturity variables.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "AM"; SSBTYPE <- "SSB"; if (USEGLOBAL) {set_global()}
f <- SSB ~ FAM + MAM + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "AM" & SSBTYPE == "SSB") {do_mcmc(f)}




## ----AMf_mcmcglmm, echo=TRUE--------------------------------------------------
#Runs MCMCglmm for all our Maturity variables.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "AM"; SSBTYPE <- "FSSB"; if (USEGLOBAL) {set_global()}
f <- FSSB ~ FAM + MAM + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "AM" & SSBTYPE == "FSSB") {do_mcmc(f)}




## ----AMm_mcmcglmm, echo=TRUE--------------------------------------------------
#Runs MCMCglmm for all our Maturity variables.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "AM"; SSBTYPE <- "MSSB"; if (USEGLOBAL) {set_global()}
f <- MSSB ~ FAM + MAM + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "AM" & SSBTYPE == "MSSB") {do_mcmc(f)}




## ----IUCN_mcmcglmm, echo=TRUE-------------------------------------------------
#Runs MCMCglmm for all our IUCN variables.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "IUCN"; SSBTYPE <- "SSB"; if (USEGLOBAL) {set_global()}
f <- SSB ~ IUCN + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "IUCN" & SSBTYPE == "SSB") {do_mcmc(f)}




## ----IUCNf_mcmcglmm, echo=TRUE------------------------------------------------
#Runs MCMCglmm for all our IUCN variables.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "IUCN"; SSBTYPE <- "FSSB"; if (USEGLOBAL) {set_global()}
f <- FSSB ~ IUCN + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "IUCN" & SSBTYPE == "FSSB") {do_mcmc(f)}




## ----IUCNm_mcmcglmm, echo=TRUE------------------------------------------------
#Runs MCMCglmm for all our IUCN variables.
ANALYSIS <- "MCMCGLMM"; HYPOTHESIS <- "IUCN"; SSBTYPE <- "MSSB"; if (USEGLOBAL) {set_global()}
f <- MSSB ~ IUCN + 1
if (ANALYSIS == "MCMCGLMM" & HYPOTHESIS == "IUCN" & SSBTYPE == "MSSB") {do_mcmc(f)}


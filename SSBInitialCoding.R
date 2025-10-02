#Start by clicking session > setting working directory > to source file location
data <- read.csv("FirstRCodeTest.csv", header = TRUE, sep = ",") #reads my csv collection
library(ape) #keeps ape package functions accessible
library(phytools) #keeps phytools package functions accessible
library(diversitree) #keeps diversitree package functions accessible
library(caper) #keeps caper package functions accessible
#library(phylolm) #keeps phylolm package functions accessible (this one is glitched rn, not installing properly)
library(MCMCglmm) #keeps MCMCglmm package functions accessible


#read in tree for comparison:
phy <- read.tree("Complete_phylogeny.nex") #reduces reading full 1000 trees to just phy command
phytree1 <- phy[[1]] #reduces to just one of the 1000 trees

LON <- data$Species #reduces species from my data to LON ("List of Names")
phytree1$tip.label <-gsub("_"," ",phytree1$tip.label) #Converts underscores to space in the full trees data to match with mine
LON[which(!(LON %in% phytree1$tip.label))] #tells you which are in mine that are not in theirs

SPECMATS <- LON[which((LON %in% phytree1$tip.label))] #creates a variable SPECMATS for all the species of mine which match those in the full trees
SPECMATSTREE <- keep.tip(phytree1, SPECMATS) #creates a new tree reducing phytree1 to just SPECMATS, can then plot this with plot.phylo(SPECMATSTREE)
dataSPECMATS <- data[which((LON %in% phytree1$tip.label)),] #creates dataset of SPECMATS (rows that match with the tree).

#ape analysis
#We use a CTMC to determine correlation between different variable in a matrix.
dataAFSSBmtrx <- matrix(c(0, 1, 2, 0),2) #forms matrix for future use by AFSSB
dataAFSSB <- ace(dataSPECMATS$Female_Homosexuality,SPECMATSTREE, type = "discrete", model = dataAFSSBmtrx) #reconstructed ancestral states of female SSB

for (i in 1:nrow(dataSPECMATS)) { #create for loop for finding four states of dual relation between FSSB and TSP
  dataifFSSB <- dataSPECMATS$Female_Homosexuality[i]
  dataifTSP <- dataSPECMATS$Typically.single.progeny.[i]
  if (dataifFSSB == 0) {
    if (dataifTSP == "N") {
      state <- 1 #no FSSB, no TSP
    } else {
      state <- 2 #no FSSB, yes TSP
    }
  }else {
    if (dataifTSP == "N") {
      state <- 3 #yes FSSB, no TSP
    } else {
      state <- 4 #yes FSSB, yes TSP
    }
  }
  dataSPECMATS$FSSB_TSP[i] <- state #takes new state info and adds it as a row to data named FSSB_TSP
}

dataAFSSB_TSPmtrx <- matrix(c(0, 1, 2, 0, 3, 0, 0, 4, 5, 0, 0, 6, 0, 7, 8, 0),4, byrow=TRUE) #forms matrix for use in FSSB_TSP plotting, manually adjusted 0 terms for things we determine possible
dataAFSSB_TSP <- ace(dataSPECMATS$FSSB_TSP,SPECMATSTREE, type = "discrete", model = dataAFSSB_TSPmtrx) #plots matric of AFSSB_TSP


#caper analysis
#going to do a phylogenetic signal test here
#D=1 complete randomness of trait (no phylogenetic signal), all rest have phylogenetic signal (D=0 trait follows brownian, D positve under 1 more different than brownian, d negative more similar than brownian)
phylo.d(dataSPECMATS,SPECMATSTREE, Species, Typically.single.progeny.) #tests for a D value for a single trait. If D positive, more divergence than brownian, if negative less so.)
phylo.d(dataSPECMATS,SPECMATSTREE, Species, Female_Homosexuality) #same analysis for female SSB. Positive D shows more differentiation between species than expected.
phylo.d(dataSPECMATS,SPECMATSTREE, Species, Male_Homosexuality) #same analysis for male SSB. Positive D (even higher) shows more differentiation between species than expected.

#MCMCglmm analysis
plot(SPECMATSTREE) #just brings this tree back in so R doesn't get mad at me
inv.phylo <- inverseA(force.ultrametric(SPECMATSTREE), "TIPS")$Ainv #reads a tree and creates a matrix that tells us how much time two species have spent evolving together vs separately.
prior <- list(G=list(G1=list(V=1,nu=0.02), #makes these things called priors, which are prior assumptions on variance distribution for the following model
              G1=list(V=1,nu=0.02)),
          R=list(V=1,nu=0.02))
MCanalysis <- MCMCglmm(Female_Homosexuality~1, #runs MCMCglmm fancy linear regression analysis
         random = ~Typically.single.progeny.+Progeny.Count,
         ginverse = list(Species=inv.phylo),
         prior = prior,
         data = dataSPECMATS)
summary(MCanalysis)




##################################################################################################################



data <- read.csv("temp_data.csv", header = TRUE, sep = ",") #reads my csv collection

names <- colnames(data)
type <- c(0,NA,1,1,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

ind_names <- names[c(3:13,15:length(names))]
ind_type <- type[c(3:13,15:length(names))]

cat_names <- ind_names[which(ind_type==0)]
#cat_names <- cat_names[1]
quant_names <- ind_names[which(ind_type==1)]

SSB <- as.factor(data$SSB)

tests <- data.frame(matrix(NA,ncol=3,nrow=0))
for (name in names) {
  if (name %in% cat_names) {
    ind <- as.factor(data[,which(colnames(data)==name)])
    test <- chisq.test(x=SSB,y=ind)$p.value
    if (test <= 0.05) {
      row <- c(name=name,p=test,keep="Y")
    } else {
      row <- c(name=name,p=test,keep="N")
    }
    tests <- rbind(tests,row)
    print(row)
    print(chisq.test(x=SSB,y=ind)$observed)
  }
  if (name %in% quant_names) {
    ind <- as.numeric(data[,which(colnames(data)==name)])
    test <- t.test(ind ~ SSB, data=cbind(SSB,ind),na.rm=TRUE)$p.value
    if (test <= 0.05) {
      row <- c(name=name,p=test,keep="Y")
    } else {
      row <- c(name=name,p=test,keep="N")
    }
    tests <- rbind(tests,row)
    print(row)
    print(t.test(ind ~ SSB, data=cbind(SSB,ind),na.rm=TRUE)$estimate)
  }
}

colnames(tests) <- c("name","p","keep")
print(tests)



dim <- as.factor(data$Radically.Dimorphic.Morphology.)

chisq.test(x=SSB,y=dim)

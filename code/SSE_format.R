library(ape)

tree_filepath <- "data/mcc_tree.txt"
print(paste0("Reading in tree from ",tree_filepath),quote=F)
phy_full <- ape::read.tree(tree_filepath)

strong_files <- list.files("SSE_Clades/Strong")
weak_files <- list.files("SSE_Clades/Weak")

names <- c()
locations <- c()

for (i in 1:length(strong_files)) {
  locations <- c(locations,paste0("SSE_Clades/Strong/",strong_files[i]))
  names <- c(names, tools::file_path_sans_ext(strong_files[i]))
}
for (i in 1:length(weak_files)) {
  locations <- c(locations,paste0("SSE_Clades/Weak/",weak_files[i]))
  names <- c(names, tools::file_path_sans_ext(weak_files[i]))
}

for (i in 1:length(names)) {

  data_filepath <- locations[i]
  out_filepath <- paste0("SSE_Clades/formatted/",names[i],".csv")
  out_treepath <- paste0("SSE_Clades/formatted/",names[i],".txt")

  print(paste0("Reading in data from ",data_filepath),quote=F)
  data <- read.csv(data_filepath, header = TRUE, sep = ",")[c("Binomial.1.2","SSB")]
  colnames(data) <- c("Taxon","SSB")
  data[is.na(data)] <- "?"
  write.csv(data,out_filepath,row.names=FALSE,quote=FALSE)

  print("Getting tree subset...",quote=F)
  LON <- data$Taxon
  SPECMATS <- LON[which((LON %in% phy_full$tip.label))]
  if (length(LON) > length(SPECMATS)) {print("Warning: there are taxon in the dataset which are not in the tree!")}
  phy <- ape::keep.tip(phy_full, SPECMATS)
  print(paste0("Total species: ",length(phy$tip.label)),quote=F)
  write.tree(phy,out_treepath)

}

phy <- read.tree("SSE_Clades/formatted/Rodentia_Supramyomorpha_Myomorphi.txt")
file.copy("SSE_Clades/formatted/Rodentia_Supramyomorpha_Myomorphi.txt",
          "SSE_Clades/formatted/old_Rodentia_Supramyomorpha_Myomorphi.txt")

node_ind <- which(phy$edge.length == min(phy$edge.length))
node <- phy$edge[node_ind,]
parent_ind <- which(phy$edge[,2]==node[1])
parent <- phy$edge[parent_ind,]
node_length <- phy$edge.length[node_ind]
parent_length <- phy$edge.length[parent_ind]
new_node_length <- node_length + 0.0001
new_parent_length <- parent_length - 0.0001

phy$edge.length[node_ind] <- new_node_length
phy$edge.length[parent_ind] <- new_parent_length

write.tree(phy,"SSE_Clades/formatted/Rodentia_Supramyomorpha_Myomorphi.txt")

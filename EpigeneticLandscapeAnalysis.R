  library(igraph)
  library(BoolNet)
  source("Functions.R")
  
  #Name of the network file in BoolNet format
  name = "EMTgenericModel.bnet"
  #Output file name
  nameF <- "EMTgeneric"
  #number of repetitions for each network state
  repet <- 10000
  #probability of error
  noise <- 0.01
  #gen names that will be mutate
  # "WT" means that the WT model will be simulated
  genM <- list(c("WT"),c("p53"),c("p16"))
  #vector to specify the types of mutations, "0" for loss of function and "1" for gain of function. The length of the vector must be the same as that of the genM vector.
  mut <- list(c(1),c(0),c(0))
  #To get the matrix of the net transition rates among the attractors of the network for each specified mutation
  for(i in 1:length(genM)){
    if(genM[i]== "WT"){
      network<- loadNetwork(name)
      getMatrixTrans(network, repet, noise, nameF, genM[i], mut[i])
      next
    }
    network<- loadNetwork(name)
    network <- fixGenes(network, unlist(genM[i]), unlist(mut[i]))
    getMatrixTrans(network, repet, noise, nameF, unlist(genM[i]), unlist(mut[i]))
  }
  
  
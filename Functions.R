library(igraph)
library(BoolNet)
library(markovchain)

#To get the basins of attraction and the transition table of the network
#Print the attractors of the network
dataExtract <- function(network, nameFile){
  attrs <- getAttractors(network)
  pdf(file= nameFile)
  plotAttractors(attrs, drawLegend = FALSE)
  dev.off()
  nF <- paste("Table", nameFile)
  s <- paste(toString(attrs$stateInfo[5]),toString(attrs$stateInfo[6]),toString(attrs$attractors) )
  write.table(s, file=nF, append=FALSE, row.names= FALSE, sep= " " , col.names = FALSE, eol = "\n")
  noAttrs <- length(attrs$attractors)
  TransTable <- getTransitionTable(attrs)
  vecBasinSize <- vector("integer", noAttrs)
  for(i in 1:noAttrs){
    vecBasinSize[i] <- attrs$attractors[[i]]$basinSize
  }
  rm(attrs)
  gc()
  Inits <- TransTable[[1]]
  noStates <- 2^length(network$genes)
  for(i in 2:length(network$genes)) Inits <- cbind(Inits, TransTable[[i]])
  CharStateSpace <- apply(Inits, 1, function(i) paste(i, collapse=""))
  rm(Inits)
  gc()
  AttractsVector <- TransTable$attractorAssignment
  inicio <- length(network$genes)+1
  nextSt <- TransTable[[inicio]]
  final <- length(network$genes)-1
  for(i in inicio+1:final) nextSt <- cbind(nextSt, TransTable[[i]])
  CharStateSpaceNext <- apply(nextSt, 1, function(i) paste(i, collapse=""))
  rm(TransTable)
  rm(nextSt)
  gc()
  l <- matrix(c(CharStateSpaceNext, AttractsVector), ncol = 2)
  dic <- vector("integer", noStates)
  for(i in 1:length(CharStateSpace)){
    index <- strtoi(CharStateSpace[i], base = 2) + 1
    dic[index] = AttractsVector[i]
  }
  return (list(l, vecBasinSize, dic))
}

#To get the Markov matrix and the matrix of the net transition rates. Print both matrices.
getMatrixTrans <- function(network, repet, noise, name, genM, mut){
  n <-""
  for (i in 1:length(genM)) {
    if(genM[i] != "WT"){
      if(mut[i] == "0"){
        n <- paste(n, "KO", genM[i])
      }else{
        n <- paste(n,"OE", genM[i])
      }
    } else {
      n <- "WT"
    }
  }
  r <- toString(noise)
  r <- sub("\\.","_",r)
  nameF <- paste("Attractors", name, repet, "rep", n, r, "prob.pdf")
  listStatesBasinsDic <- dataExtract(network,nameF)
  noAttrs <- length(listStatesBasinsDic[[2]])
  noStates <- length(listStatesBasinsDic[[1]][,1])
  noNodes <- nchar(listStatesBasinsDic[[1]][1])
  mutNode <-c()
  for(i in 1:length(genM)){
    mutNode <- c(mutNode, which(network$genes == genM[i]))
  }
  if(length(mutNode) == 0){
    mutNode <- c(-2)
  }

  matrixTrans <- matrix(0,noAttrs, noAttrs)
  prob <- rbinom(10000*repet*noNodes,1,noise)
  count <- 1
  for(i in 1:noStates){
    if(i%%10000 == 0){
      prob <- rbinom(10000*repet*noNodes,1,noise)
      count <- 1
      cat(i, " ")
    }
    for(j in 1:repet){
      b = FALSE
      est <- listStatesBasinsDic[[1]][i,1]
      originalAtt <- as.integer(listStatesBasinsDic[[1]][i,2])
      for(k in 1:nchar(est)){
        if(prob[count] == 1){
          dif <- TRUE
          for(m in 1:length(mutNode)){
            if(mutNode[m] != k){
              next
            } else{
              dif <- FALSE
              break
            }
          }
          if(dif){
            b=TRUE
            if(substr(est, k,k) == '1'){
              substr(est,k,k) <- '0'
            }else{
              substr(est,k,k) <- '1'
            }
          }
        }
        count <- count + 1
      }
      if(b){
        index <- strtoi(est, base=2)+1
        newAtt <- listStatesBasinsDic[[3]][index]
      } else {
        newAtt <- originalAtt
      }
      matrixTrans[originalAtt, newAtt] = matrixTrans[originalAtt, newAtt]+1
    }
  }

  for(i in 1:noAttrs){
    for(j in 1:noAttrs){
      matrixTrans[i,j] <- matrixTrans[i,j]/(listStatesBasinsDic[[2]][i]*repet)
    }
  }
  
  s <- paste("MarkovMatrix",name,n,repet,"rep",r,"prob","R.txt")
  write.table(matrixTrans, file=s, append=FALSE, row.names= FALSE, sep= " " , col.names = FALSE, eol = "\n")
  
  rm(listStatesBasinsDic)
  rm(prob)
  gc()
  
  mc <- new("markovchain",states=letters[1:noAttrs], transitionMatrix=matrixTrans)
  MFPT <- meanFirstPassageTime(mc)
  netTransMat <- (1/MFPT) - 1/t(MFPT)
  diag(netTransMat) <- 0
  s <- paste("transitionRatesMatrix",name,n,repet,"rep",r,"prob","R.txt")
  write.table(netTransMat, file=s, append=FALSE, row.names= FALSE, sep= " " , col.names = FALSE, eol = "\n")
}


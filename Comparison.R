library(BoolNet)

#Nodes must be in the same order in both .bnet files.
n1 = loadNetwork("reducedhEMTModelUsingGinsim.bnet") 
n2 <- loadNetwork("reducedhEMTModelUsingBooleanAlgebra.bnet")


#to compare regulators and truth tables for each node
for(i in n1$genes){
  print(i)
  r1 <- n1$interactions[[i]][1]
  r2 <- n2$interactions[[i]][1]
  if(!identical(r1, r2)){
    s <- paste("Regulators are different ", " node1: ", r1, " node2: ", r2)
    print(s)
  }
  t1 <- n1$interactions[[i]][2]
  t2 <- n2$interactions[[i]][2]
  if(!identical(t1, t2)){
    s <- paste("Truth tables are different ", " node1: ", t1, " node2: ", t2)
    print(s)
  }
}

#To obtain the attractors of both networks
att1 <- getAttractors(n1)
plotAttractors(att1)
att2 <- getAttractors(n2)
plotNetworkWiring(n2)
plotAttractors(att2)


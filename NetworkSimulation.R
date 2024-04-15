library(BoolNet)
#simulate the extended and reduced models using synchronous and asynchronous schemes.


#Loading the extended model
r <- loadNetwork("CompletehEMTmodel.bnet")


startStateNumber <- 5000000
startStateNumberA <- 100000
pdf(file="Simulation of networks.pdf")

attr0 <- getAttractors(r, type= "synchronous", method="random",  startStates=startStateNumber)
tit <- paste("Attractors of the extended model, Synchronous scheme")
plotAttractors(attr0, title= tit, drawLegend = FALSE)

attr0 <- getAttractors(r, type= "asynchronous", method="random",  startStates=startStateNumberA)
tit <- paste("Attractors of the extended model, Asynchronous scheme")
plotAttractors(attr0, title= tit, drawLegend = FALSE)


#Loading the reduced network
r2 <- loadNetwork("reducedhEMTModelUsingBooleanAlgebra.bnet")

#In this case, the simulation is exhaustive
attr0 <- getAttractors(r2, type= "synchronous")
tit <- paste("Attractors of the reduced model, Synchronous scheme")
plotAttractors(attr0, title= tit, drawLegend = FALSE)

attr0 <- getAttractors(r2, type= "asynchronous", method="random",  startStates=startStateNumberA)
tit <- paste("Attractors of the reduced model, Asynchronous scheme")
plotAttractors(attr0, title= tit, drawLegend = FALSE)

dev.off()

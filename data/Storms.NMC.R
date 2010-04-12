load("Storms.NMC.Rdata")
Storms.NMC <- results
class(Storms.NMC) <- "CCI.object"
rm(results); gc(reset=TRUE)


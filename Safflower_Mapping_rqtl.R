#Analyse the Safflower data sets

#Mapping library
library(qtl)
library(readr)

saff_F2 <- read.cross("csvr", file="safflower_x395_SNP_1row.csv", estimate.map = FALSE, genotypes = c(0,2,1), na.strings = c("-","X"))


summary(saff_F2)

#visualisation of crosses - marks indicate missing data
par(mfrow=c(1,1), las=1)
plotMissing(saff_F2)

par(mfrow=c(1,2), las=1)
plot(ntyped(saff_F2), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(saff_F2, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")
par(mfrow=c(1,1), las=1) #(restore single plot)


#visualisation of present markers and crosses without poor samples
#NB: number at the end determines the MINIMUM CUTOFF for the number of markers in a genotype
saff_F2_no_bad_data<- subset(saff_F2, ind=(ntyped(saff_F2)>1800)) #NB: NB: This is approximately 90% complete markers for individuals

#NB: number at the end determines the MINIMUM CUTOFF for the number of individuals for a marker
nt.bymar <- ntyped(saff_F2_no_bad_data, "mar")
todrop <- names(nt.bymar[nt.bymar < 75]) #NB: This is so each marker covers approximately 80% of the samples
saff_F2_no_bad_data <- drop.markers(saff_F2_no_bad_data, todrop)

summary(saff_F2_no_bad_data)
par(mfrow=c(1,2), las=1)
plot(ntyped(saff_F2_no_bad_data), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(saff_F2_no_bad_data, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")
par(mfrow=c(1,1), las=1) #(restore single plot)


#Identify duplicate genotypes
cg <- comparegeno(saff_F2_no_bad_data)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. genotypes matching another")
rug(cg[lower.tri(cg)])


#provide a matrix of all identical genotypes (100%)
#NB: for % similar (not identical), replace 'cg >= 1' with 'cg > .##'
wh <- which(cg > .95, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh


#identify duplicate markers and drop them
dup <- findDupMarkers(saff_F2_no_bad_data, exact.only=TRUE)
saff_F2_nodup_markers <- drop.markers(saff_F2_no_bad_data, unlist(dup))


#Study individuals' genotype frequencies
g <- pull.geno(saff_F2_nodup_markers)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))


#Show the frequency of genotype distribution
#NB: For a RIL population, these should practically be a mirror of one another
par(mfrow=c(1,3), las=1)
for(i in 1:3){
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA","AB", "BB")[i],ylim=c(0,1))
}
par(mfrow=c(1,1), las=1)

saff_F2_map <- est.rf(saff_F2_nodup_markers)
checkAlleles(saff_F2_map, threshold = 5)

plotMap(saff_F2_nodup_markers, saff_F2_map)

#This plots the recombination factor against 
rf <- pull.rf(saff_F2_map)
lod <- pull.rf(saff_F2_map, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")


#Form Linkage Groups - This will be based on the recombination factor and LOD score
lg <- formLinkageGroups(saff_F2_map, max.rf=0.35, min.lod=10)
table(lg[,2])
length(table(lg[,2]))

#Examining changing the minimum LOD scores, 14.5 creates a distinct linkage group at groups 1 and 2.
saff_F2_map <- formLinkageGroups(saff_F2_map, max.rf=0.35, min.lod=14.5, reorgMarkers = TRUE)
summary(saff_F2_map)

par(mfrow=c(1,1), las=1)
plotRF(saff_F2_map, alternate.chrid=TRUE)


#Observe marker names 
#mn<x> <- markernames(mapthis, chr=<x>)

mn2 <- markernames(saff_F2_map, chr=2)

par(mfrow=c(2,1), las=1)
plot(rf, mn2[1], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
abline(h=0.5, lty=2)
plot(lod, mn2[1], bandcol="gray70", alternate.chrid=TRUE)





#Phenotyping
saff_F2_map <- calc.errorlod(saff_F2_map, error.prob=0.01)
saff_F2_map <- calc.genoprob(saff_F2_map, step=2)
out.2p <- scanone(saff_F2_map, pheno.col="Flowering", model="2part",
                  upper=TRUE)
summary(out.2p)
plot(out.2p)



#This allows multithreading support
library("mvabund", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")

#Mapping library
library(ASMap)
library(qtl)
library(readr)

data("mapF2", package = "ASMap")


#data import
test_SBS <- read.cross("csvr", file="/data/cloudStor/Shared/IWYP60/Data/Genomics/SBS_Mapping_Gx_Data/SBS_RIL_marker_matrix_CIMMYT2017.csv", 
                       estimate.map = FALSE, crosstype = "riself", na.strings = c("-","X"))

#test data for the QTL packages
data("mapthis")

#examine details of crosses
summary(test_SBS)
summary(mapthis)

#visualisation of crosses - marks indicate missing data
par(mfrow=c(1,1), las=1)
plotMissing(test_SBS)
plotMissing(mapthis)

par(mfrow=c(1,2), las=1)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")
par(mfrow=c(1,1), las=1)


#visualisation of present markers and crosses
par(mfrow=c(1,2), las=1)
plot(ntyped(test_SBS), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(test_SBS, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")
par(mfrow=c(1,1), las=1) #(restore single plot)

#visualisation of present markers and crosses without poor samples
#NB: number at the end determines the MINIMUM CUTOFF for the number of markers in a genotype
test_SBS_nomissing <- subset(test_SBS, ind=(ntyped(test_SBS)>1100))

#NB: number at the end determines the MINIMUM CUTOFF for the number of individuals for a marker
nt.bymar <- ntyped(test_SBS_nomissing, "mar")
todrop <- names(nt.bymar[nt.bymar < 150])
test_SBS_nomissing <- drop.markers(test_SBS_nomissing, todrop)


summary(test_SBS_nomissing)
#Display the same as above but with the poor data removed
par(mfrow=c(1,2), las=1)
plot(ntyped(test_SBS_nomissing), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(test_SBS_nomissing, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")
par(mfrow=c(1,1), las=1)



#Identify duplicate genotypes
cg <- comparegeno(test_SBS_nomissing)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. genotypes matching another")
rug(cg[lower.tri(cg)])

#provide a matrix of all identical genotypes (100%)
#NB: for % similar (not identical), replace 'cg >= 1' with 'cg > .##'
wh <- which(cg >= 1, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh

#show number of matching markers
#NB: rows/columns correspond to alleles
g <- pull.geno(test_SBS_nomissing)
table(g[19,], g[24,])

#remove mismatching genotypes
for(i in 1:nrow(wh)) {
  tozero <- !is.na(g[wh[i,1],]) & !is.na(g[wh[i,2],]) & g[wh[i,1],] != g[wh[i,2],]
  test_SBS_nomissing$geno[[1]]$data[wh[i,1],tozero] <- NA
}

#remove one of the multiple genotypes
test_SBS_nomissing <- subset(test_SBS_nomissing, ind=-wh[,2])
summary(test_SBS_nomissing)

#identify duplicate markers and drop them
dup <- findDupMarkers(test_SBS_nomissing, exact.only=TRUE)
test_SBS_nodupmrks <- drop.markers(test_SBS_nomissing, unlist(dup))



#identify markers with distorted segregation patterns
#NB: 0.05 = confidence interval
gt <- geno.table(test_SBS_nodupmrks)
gt[gt$P.value < 0.05/totmar(test_SBS_nodupmrks),]
#Q: Should we be omitting the worst of these? Should we expect an even split?


#Study individuals' genotype frequencies
g <- pull.geno(test_SBS_nodupmrks)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:2)))
gfreq <- t(t(gfreq) / colSums(gfreq))


#Show the frequency of genotype distribution
#NB: For a RIL population, these should practically be a mirror of one another
par(mfrow=c(1,1), las=1)
for(i in 1:1){
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA","BB")[i],ylim=c(0,1))
}

test_SBS_map <- est.rf(test_SBS_nodupmrks)
checkAlleles(test_SBS_map, threshold = 3)

#This plots the recombination factor against 
rf <- pull.rf(test_SBS_map)
lod <- pull.rf(test_SBS_map, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")


#Form Linkage Groups - This will be based on the recombination factor and LOD score
lg <- formLinkageGroups(test_SBS_map, max.rf=0.5, min.lod=6)
table(lg[,2])
length(table(lg[,2]))


test_SBS_map <- formLinkageGroups(test_SBS_map, max.rf=0.5, min.lod=6, reorgMarkers=TRUE)


plotRF(test_SBS_map, alternate.chrid=TRUE)
rf <- pull.rf(test_SBS_map)
lod <- pull.rf(test_SBS_map, what="lod")
mn4 <- markernames(test_SBS_map, chr=4)
par(mfrow=c(2,1))
plot(rf, mn4[24], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
abline(h=0.5, lty=2)
plot(lod, mn4[3], bandcol="gray70", alternate.chrid=TRUE)



test_SBS_map <- orderMarkers(test_SBS_map, chr=1)

test_SBS_map_5ish <- test_SBS_map

test_SBS_map_5ish <- orderMarkers(test_SBS_map, chr=5)
plotRF(test_SBS_map_5ish, alternate.chrid=TRUE)


lg <- formLinkageGroups(test_SBS_map_5ish, max.rf=0.5, min.lod=6, reorgMarkers = TRUE)
table(lg[,2])

#Show marker order (in "Chromosome 5")
pull.map(test_SBS_map_5ish, chr=5)
rip5 <- ripple(test_SBS_map_5ish, chr=5, window=7)

rip5lik <- ripple(test_SBS_map_5ish, chr=5, window=3, method="likelihood", error.prob = 0.005)

compareorder(test_SBS_map_5ish, chr=5, c(1:66,68,67,69:82), error.prob=0)
test_SBS_map_5ish <- switch.order(test_SBS_map_5ish, chr=5, c(1:66,68,67,69:82), error.prob=0.005)
pull.map(test_SBS_map_5ish, chr=5)


plotMap(test_SBS_map_5ish, show.marker.names=TRUE)
plotRF(test_SBS_map_5ish)

dropone <- droponemarker(test_SBS_map_5ish, error.prob=0.005)


#dist.fun - kosambi = moderate interference, haldane = no interference
#objective.fun - COUNT = minimal sum, ML = maximum likelihood

cat 
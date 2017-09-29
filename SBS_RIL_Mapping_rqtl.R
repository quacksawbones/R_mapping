#Mapping library
library(qtl)
library(readr)


#data import
setwd("/data/cloudStor/Shared/IWYP60/Data/Genomics/SBS_Mapping_Gx_Data")


SBS_RIL <- read.cross("csvsr", genfile="SBS_RIL_marker_matrix_CIMMYT2017.csv", phefile = "SBS_TILL_17_Phenotyping_mod_R1_t_DC.csv", estimate.map = TRUE, crosstype = "riself", na.strings = c("-","X"))




#examine details of crosses
summary(SBS_RIL)

#visualisation of crosses - marks indicate missing data
par(mfrow=c(1,1), las=1)
plotMissing(SBS_RIL)

#Number of 
plotPheno(SBS_RIL, pheno.col = "j")

colnames(SBS_RIL$pheno)

#visualisation of present markers and crosses
par(mfrow=c(1,2), las=1)
plot(ntyped(SBS_RIL), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(SBS_RIL, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")
par(mfrow=c(1,1), las=1) #(restore single plot)


#visualisation of present markers and crosses without poor samples
#NB: number at the end determines the MINIMUM CUTOFF for the number of markers in a genotype
SBS_RIL_no_bad_data<- subset(SBS_RIL, ind=(ntyped(SBS_RIL)>1100))

#NB: number at the end determines the MINIMUM CUTOFF for the number of individuals for a marker
nt.bymar <- ntyped(SBS_RIL_no_bad_data, "mar")
todrop <- names(nt.bymar[nt.bymar < 150])
SBS_RIL_no_bad_data <- drop.markers(SBS_RIL_no_bad_data, todrop)

#At this point, the SBS data set should not have any markers and 
summary(SBS_RIL_no_bad_data)


#Display the same as above but with the poor data removed
par(mfrow=c(1,2), las=1)
plot(ntyped(SBS_RIL_no_bad_data), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(SBS_RIL_no_bad_data, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")
par(mfrow=c(1,1), las=1)



#Identify duplicate genotypes
cg <- comparegeno(SBS_RIL_no_bad_data)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. genotypes matching another")
rug(cg[lower.tri(cg)])

#provide a matrix of all identical genotypes (100%)
#NB: for % similar (not identical), replace 'cg >= 1' with 'cg > .##'
wh <- which(cg >= 1, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh

#show number of matching markers
#NB: rows/columns correspond to alleles
g <- pull.geno(SBS_RIL_no_bad_data)

#NB: Use this to examine the similarities between the pairs
#table(g[19,], g[24,])

#remove mismatching genotypes
for(i in 1:nrow(wh)) {
  tozero <- !is.na(g[wh[i,1],]) & !is.na(g[wh[i,2],]) & g[wh[i,1],] != g[wh[i,2],]
  SBS_RIL_no_bad_data$geno[[1]]$data[wh[i,1],tozero] <- NA
}

#remove one of the multiple genotypes
SBS_RIL_no_bad_data <- subset(SBS_RIL_no_bad_data, ind=-wh[,2])
summary(SBS_RIL_no_bad_data)

#identify duplicate markers and drop them
dup <- findDupMarkers(SBS_RIL_no_bad_data, exact.only=TRUE)
SBS_RIL_nodup_markers <- drop.markers(SBS_RIL_no_bad_data, unlist(dup))



#identify markers with distorted segregation patterns
#NB: 0.05 = confidence interval
gt <- geno.table(SBS_RIL_nodup_markers)
gt[gt$P.value < 0.05/totmar(SBS_RIL_nodup_markers),]
#Q: Should we be omitting the worst of these? Should we expect an even split?


#Study individuals' genotype frequencies
g <- pull.geno(SBS_RIL_nodup_markers)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:2)))
gfreq <- t(t(gfreq) / colSums(gfreq))


#Show the frequency of genotype distribution
#NB: For a RIL population, these should practically be a mirror of one another
par(mfrow=c(1,1), las=1)
for(i in 1:1){
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA","BB")[i],ylim=c(0,1))
}

SBS_RIL_map <- est.rf(SBS_RIL_nodup_markers)
checkAlleles(SBS_RIL_map, threshold = 3)

plotMap(SBS_RIL_nodup_markers, SBS_RIL_map)

#This plots the recombination factor against 
rf <- pull.rf(SBS_RIL_map)
lod <- pull.rf(SBS_RIL_map, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")


#Form Linkage Groups - This will be based on the recombination factor and LOD score
lg <- formLinkageGroups(SBS_RIL_map, max.rf=0.5, min.lod=6)
table(lg[,2])
length(table(lg[,2]))



SBS_RIL_map <- formLinkageGroups(SBS_RIL_map, max.rf=0.5, min.lod=6, reorgMarkers=TRUE)

par(mfrow=c(1,1), las=1)
plotRF(SBS_RIL_map, alternate.chrid=TRUE)
rf <- pull.rf(SBS_RIL_map)
lod <- pull.rf(SBS_RIL_map, what="lod")
mn4 <- markernames(SBS_RIL_map, chr=4)

#For this 
par(mfrow=c(2,1), las=1)
plot(rf, mn4[24], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
abline(h=0.5, lty=2)
plot(lod, mn4[24], bandcol="gray70", alternate.chrid=TRUE)


#Change chr=# to whichever linkage group needs to be changed
#NB: This will take a LONG time!
SBS_RIL_map <- orderMarkers(SBS_RIL_map, chr=1)

plotRF(SBS_RIL_map, alternate.chrid=TRUE)


lg <- formLinkageGroups(SBS_RIL_map, max.rf=0.5, min.lod=6, reorgMarkers = TRUE)
table(lg[,2])

#Show marker order (in "Chromosome 5")
pull.map(SBS_RIL_map, chr=5)

#Examine the number of different permutations of markers
#NB: Takes time... but use this to estimate the window size for below
rip5 <- ripple(SBS_RIL_map, chr=5, window=7)

#Examine whether there are better marker orders available
#NB: Takes A LONG time... increasing the window size (exponentially?) increases the processing time (i.e. more permutations)
rip5lik <- ripple(SBS_RIL_map, chr=5, window=3, method="likelihood", error.prob = 0.005)
summary(rip5lik)

#Examine how changing the order of a marker or two will change the chromosome size, using c(1:x-1, x, y, y+1)
compareorder(SBS_RIL_map, chr=5, c(1:66,68,67,69:82), error.prob=0)

#switch the order of markers, using c(1:x-1, x, y, y+1), error.prob will change the estimations
SBS_RIL_map <- switch.order(SBS_RIL_map, chr=5, c(1:66,68,67,69:82), error.prob=0.005)
pull.map(SBS_RIL_map, chr=5)


plotMap(SBS_RIL_map, show.marker.names=TRUE)
plotRF(SBS_RIL_map)


#Grind through and see the effect of dropping each marker and see what the effect is on the map
dropone <- droponemarker(SBS_RIL_map, error.prob=0.005)


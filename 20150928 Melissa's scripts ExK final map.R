
########################## Modified from Julian's quick script

library(qtl)
library(ASMap)

#map <- read.cross(file = "ExK.csv",  format = "csvr", genotypes = c("A","B"),
#                  na.strings = c("-","U"), estimate.map = FALSE)

map <- read.cross(file = "/data/cloudStor/Shared/IWYP60/Data/Genomics/ExK RIL Genotypic data/ExK_Darren_Temp.csv",  format = "csvr", genotypes = c("AA","BB"),
                  na.strings = c("-","U"), estimate.map = FALSE)


#272 ind and 3538 markers

## checking marker set vitals .. missing value plot

plot.missing(map)

## some genotypes have large amounts of missing values .. are they reliable?
## omit genotypes with more than 25% of the markers missing

map1 <- subset(map, ind = ntyped(map) > (totmar(map)*75/100))
nind(map1)

#dropped 50 lines, 222 lines kept

## check marker stats

#clones

gc <- genClones(map1, tol = 0.90)

# The statistics for the pairs of genotypes is stored in cgd

gc$cgd


#14 pairs of clones identified, but pair 4 have a high proportion of na
#remove pair 4 from cgd

cgd <- gc$cgd[-4, ]
cgd

#Consensus for remaining clones

map2 <- fixClones(map1, cgd, consensus = TRUE)

levels(map2$pheno[[1]])[grep("_", levels(map2$pheno[[1]]))]

#Marker Profile

pm <- profileMark(map2, stat.type = c("seg.dist","prop","miss"), layout = c(1, 4), type = "l")

## removing markers with proportion of missing values greater than 20%

dm <- markernames(map2)[pm$marker$miss > 0.20]
map3 <- drop.markers(map2, dm)
summary(map3)

#3526 markers remaining

## check vitals again

pm2 <- profileMark(map3, stat.type = c("seg.dist","prop","miss"), layout = c(1, 4), type = "l")


#Removing distorted at p=0.05.. THIS IS CONSERVATIVE

print("Removing distorted at p=0.05.. THIS IS CONSERVATIVE")
gt1 <- geno.table(map3)
nrow(gt1[gt1$P.value < 0.05,])   ##at p=0.05 
map4 <- drop.markers(map3, rownames(gt1[gt1$P.value < 0.05,]))
rm(gt1)
print(summary(map4))
print(summary(map3))

#3526-3210= 316 markers dropped

#Removing distorted markers - Expecting a 0.5 genotype frequency, removing markers with >0.6 or <0.4 of one genotype

#mm <- statMark(map3, stat.type = "marker")$marker$AB
#map4 <- drop.markers(map3, c(markernames(map3)[mm > 0.6],
#markernames(map3)[mm < 0.4]))


## check vitals again

pm3 <- profileMark(map4, stat.type = c("seg.dist","prop","miss"), layout = c(1, 4), type = "l", crit.val="bonf")

## turn into "riself" population

map5 <- convert2riself(map4)
names(map5$pheno)
class(map5)

## construct map

map6 <- mstmap(map5, bychr = FALSE, dist.fun = "kosambi", p.value = 1e-8)

#32 chromosomes

nmar(map6)[nmar(map6) < 5]

# four chromosomes with less than 5 markers

heatMap(map6)
chrlen(map6)
heatMap(map6, chr = "L.3")
summary(map6)

#209 inds 3210 markers

pg <- profileGen(map6, bychr = FALSE, layout = c(1, 3), lty = 2)

#couting xo

plot(countXO(map6)); abline(h=75)

pg <- profileGen(map6, bychr = FALSE, stat.type = c("xo", "dxo",
                                                    "miss"), id = "Genotype", xo.lambda = 60, layout = c(1, 3), lty = 2, cex =
                   0.7)

#Remove lines with too many crossovers

map7 <- subsetCross(map6, ind = !pg$xo.lambda)
map8 <- mstmap(map7, bychr = TRUE, dist.fun = "kosambi", trace = TRUE,
               p.value = 1e-8)
summary(map8)

#208 lines #3210 markers

chrlen(map8)
heatMap(map8)


## remove linkage groups with one marker

map9 <- map8
map9$geno <- map9$geno[nmar(map9) > 2]
nmar(map9)
summary(map9)

#208 lines 3206 markers


write.cross(map9, filestem="20150922_ExK_map_MG", format="csvr")
######################### 

#To move markers slightly so that they are not at the same position in the map

map10<-jittermap(map9, amount = 1e-6)

write.cross(map10, filestem="20150922_ExK_map_MG_jitter", format="csvr")

length(map10)
chrlen(map10)
sum(chrlen(map10))

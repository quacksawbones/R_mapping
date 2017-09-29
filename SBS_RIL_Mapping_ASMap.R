#Mapping library

library(ASMap)
library(readr)

#Import a cross object as a 
SBS_RIL <- read.cross("csvr", file="/data/cloudStor/Shared/IWYP60/Data/Genomics/SBS_Mapping_Gx_Data/SBS_RIL_marker_matrix_CIMMYT2017.csv", 
                      estimate.map = FALSE, crosstype = "riself", na.strings = c("-","X"))

SBS_RIL_existing <- read.cross("csvr", file="/data/cloudStor/Shared/IWYP60/Data/Genomics/SBS_Mapping_Gx_Data/SBS_RIL_marker_matrix_CIMMYT2017_existing_map.csv", 
                      estimate.map = FALSE, crosstype = "riself", na.strings = c("-","X"))


SBS_RIL_df <- t(read.csv("/data/cloudStor/Shared/IWYP60/Data/Genomics/SBS_Mapping_Gx_Data/SBS_RIL_marker_matrix_CIMMYT2017.csv",
                       header = FALSE, stringsAsFactors = FALSE))



#Construct a linkage map from the existing marker map r/qtl object
SBS_RIL_ASMap <- mstmap.cross(SBS_RIL, id="Genotypes", bychr = FALSE, dist.fun = "haldane", trace = TRUE)

SBS_RIL_ASMap_CIMMYT <- mstmap.cross(SBS_RIL, id="Genotypes", bychr = TRUE, dist.fun = "haldane", trace = TRUE)

length(nmar(SBS_RIL_ASMap))

heatMap(SBS_RIL_ASMap, lmax = 50)

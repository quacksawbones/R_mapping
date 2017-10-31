#Mapping library
library(qtl)
library(readr)

#stage="TILL"
stage="TILL"
#load the cross into the old version of rqtl
# SBS_RIL_BOOT_R1 <- read.cross("csvsr", genfile="SBS_RIL_marker_matrix_CIMMYT2017_existing_map.csv", phefile = "SBS_BOOT_17_Phenotyping_mod_R1_t_DC.csv",
#                       estimate.map = FALSE, crosstype = "riself", na.strings = c("-","X"))
# 
# 
# SBS_RIL_BOOT_R2 <- read.cross("csvsr", genfile="SBS_RIL_marker_matrix_CIMMYT2017_existing_map.csv", phefile = "SBS_BOOT_17_Phenotyping_mod_R2_t_DC.csv",
#                               estimate.map = FALSE, crosstype = "riself", na.strings = c("-","X"))


assign(paste0("SBS_RIL_",stage,"_R1"), read.cross("csvsr", genfile="SBS_RIL_marker_matrix_CIMMYT2017_existing_map.csv", 
                                                  phefile = paste0("SBS_",stage,"_17_Phenotyping_mod_R1_t_DC.csv"),
                                                  estimate.map = FALSE, crosstype = "riself", na.strings = c("-","X")))

assign(paste0("SBS_RIL_",stage,"_R2"), read.cross("csvsr", genfile="SBS_RIL_marker_matrix_CIMMYT2017_existing_map.csv", 
                                                  phefile = paste0("SBS_",stage,"_17_Phenotyping_mod_R2_t_DC.csv"),
                                                  estimate.map = FALSE, crosstype = "riself", na.strings = c("-","X")))


library(qtl2)

#convert the old qtl cross into the new rqtl2 version
# SBS_RIL_BOOT_R1 <- convert2cross2(SBS_RIL_BOOT_R1)
# SBS_RIL_BOOT_R2 <- convert2cross2(SBS_RIL_BOOT_R2)

assign(paste0("SBS_RIL_",stage,"_R1"), convert2cross2(get(paste0("SBS_RIL_",stage,"_R1"))))
assign(paste0("SBS_RIL_",stage,"_R2"), convert2cross2(get(paste0("SBS_RIL_",stage,"_R2"))))


# SBS_RIL_BOOT_R1_map <- insert_pseudomarkers(SBS_RIL_BOOT_R1$gmap, step=1)
# SBS_RIL_BOOT_R2_map <- insert_pseudomarkers(SBS_RIL_BOOT_R2$gmap, step=1)

assign(paste0("SBS_RIL_",stage,"_R1_map"), insert_pseudomarkers(get(paste0("SBS_RIL_",stage,"_R1"))[["gmap"]],step=1))
assign(paste0("SBS_RIL_",stage,"_R2_map"), insert_pseudomarkers(get(paste0("SBS_RIL_",stage,"_R2"))[["gmap"]],step=1))


# pr_R1 <- calc_genoprob(SBS_RIL_BOOT_R1, SBS_RIL_BOOT_R1_map, err=0.002, cores=0)
# pr_R2 <- calc_genoprob(SBS_RIL_BOOT_R2, SBS_RIL_BOOT_R2_map, err=0.002, cores=0)

pr_R1 <- calc_genoprob(get(paste0("SBS_RIL_",stage,"_R1")), get(paste0("SBS_RIL_",stage,"_R1_map")), err=0.002, cores=0)
pr_R2 <- calc_genoprob(get(paste0("SBS_RIL_",stage,"_R2")), get(paste0("SBS_RIL_",stage,"_R2_map")), err=0.002, cores=0)


apr_R1 <- genoprob_to_alleleprob(pr_R1)
apr_R2 <- genoprob_to_alleleprob(pr_R2)

#Use Create Kinship Matrix or the second one, not both...
kinship_R1 <- calc_kinship(pr_R1)
kinship_R2 <- calc_kinship(pr_R2)


# grid_R1 <- calc_grid(SBS_RIL_BOOT_R1_map$gmap, step=1)
# grid_R2 <- calc_grid(SBS_RIL_BOOT_R2_map$gmap, step=1)

grid_R1 <- calc_grid(get(paste0("SBS_RIL_",stage,"_R1_map"))[["gmap"]], step=1)
grid_R2 <- calc_grid(get(paste0("SBS_RIL_",stage,"_R2_map"))[["gmap"]], step=1)


pr_grid_R1 <- probs_to_grid(pr_R1, grid_R1)
pr_grid_R2 <- probs_to_grid(pr_R2, grid_R2)

kinship_grid_R1 <- calc_kinship(pr_grid_R1)
kinship_grid_R2 <- calc_kinship(pr_grid_R2)

kinship_loco_R1 <- calc_kinship(pr_R1, "loco", cores=0)
kinship_loco_R2 <- calc_kinship(pr_R2, "loco", cores=0)


library(qtl2scan)
 
# out_R1 <- scan1(pr_R1, SBS_RIL_BOOT_R1$pheno, cores=0)
# out_R2 <- scan1(pr_R2, SBS_RIL_BOOT_R2$pheno, cores=0)

out_R1 <- scan1(pr_R1, get(paste0("SBS_RIL_",stage,"_R1"))[["pheno"]], cores=0)
out_R2 <- scan1(pr_R2, get(paste0("SBS_RIL_",stage,"_R2"))[["pheno"]], cores=0)

# out_kinship_R1 <- scan1(pr_R1, SBS_RIL_BOOT_R1$pheno, kinship_loco_R1, cores=0)
# out_kinship_R2 <- scan1(pr_R2, SBS_RIL_BOOT_R2$pheno, kinship_loco_R2, cores=0)

out_kinship_R1 <- scan1(pr_R1, get(paste0("SBS_RIL_",stage,"_R1"))[["pheno"]], kinship_loco_R1, cores=0)
out_kinship_R2 <- scan1(pr_R2, get(paste0("SBS_RIL_",stage,"_R2"))[["pheno"]], kinship_loco_R2, cores=0)



library(qtl2plot)

par(mfrow=c(1,1), las=1)

ymx <- max(maxlod(out_R1),maxlod(out_R2)) # overall maximum LOD score

phenotypes <- c("q2rda","q2rddm","vcmax25","q2rdfm","vcmax25.Narea")
# 
# plot(out_R1, SBS_RIL_BOOT_R1_map, lodcolumn="q2rda", col="blue", ylim=c(0, ymx*1.02), lty=3, lwd=2, main="Non-Kinship")
# plot(out_kinship_R2, SBS_RIL_BOOT_R2_map, lodcolumn="q2rda", col="blue", ylim=c(0, ymx*1.02), lty=3, lwd=2, main="Kinship")
# 
# plot(out_R2, SBS_RIL_BOOT_R2_map, lodcolumn="q2rda", col="red", ylim=c(0, ymx*1.02), lty=3, lwd=2, add=TRUE)


for(i in 1:length(phenotypes)) {

  jpeg(paste0("./rqtl2_outputs/QTL_",stage,"_KIN_",phenotypes[i],".jpg"), width=800, height=600)
  
  plot(out_kinship_R1, get(paste0("SBS_RIL_",stage,"_R1_map")), lodcolumn=i, col="blue", ylim=c(0, ymx*1.02), main=paste("QTLs of",stage,"stage (using Kinship Matrix) -",phenotypes[i]), lty=3, lwd=2)
  plot(out_kinship_R2, get(paste0("SBS_RIL_",stage,"_R2_map")), lodcolumn=i, col="red", ylim=c(0, ymx*1.02), add=TRUE,lwd=2, lty=3)
  
  dev.off()
  
  jpeg(paste0("./rqtl2_outputs/QTL_",stage,"_NOKIN_",phenotypes[i],".jpg"), width=800, height=600)
  
  plot(out_R1, get(paste0("SBS_RIL_",stage,"_R1_map")), lodcolumn=i, col="blue", ylim=c(0, ymx*1.02), main=paste("QTLs of",stage,"stage (without Kinship Matrix) -",phenotypes[i]), lty=3, lwd=2)
  plot(out_R2, get(paste0("SBS_RIL_",stage,"_R2_map")), lodcolumn=i, col="red", ylim=c(0, ymx*1.02), add=TRUE,lwd=2, lty=3)

  dev.off()  
}
  


# legend("topright", lwd=2, col=c("#e41a1c", "#377eb8","#4daf4a","#984ea3","#ff7f00"),
       # c("q2rda","q2rddm","vcmax25","q2rdfm","vcmax25.Narea"), bg="gray90",cex = .75)














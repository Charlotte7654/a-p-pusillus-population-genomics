# ==============================================================
#   SETUP
# ==============================================================

cat("Setting working directory...\n")
setwd("~/ch2/05_stats")

cat("Loading packages...\n")
library(devtools)
library(SNPRelate)
library(adegenet)
library(poppr)
library(hierfstat)
library(StAMPP)
library(BiocManager)
library(gdsfmt)
library(vcfR)
library(tidyr)
library(dplyr)
library(dunn.test, lib.loc="~/R/libraries")

# ==============================================================
#   READ AND PREPARE DATA
# ==============================================================

cat("Reading VCF file...\n")
vcf_neutral <- read.vcfR("~/ch2/04_outlierdetection/03.65_HG_neutral.vcf")

cat("Converting VCF to genind...\n")
genind_neutral <- vcfR2genind(vcf_neutral)

cat("Assigning populations...\n")
populations <- c(
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "Cape Cross",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "False Bay",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Kleinzee",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Lambert's Bay",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point",
  "Pelican Point"
)


genind_neutral@pop <- as.factor(populations)

cat("Is length of populations and nInd the same?...\n")
length(populations) == nInd(genind_neutral)

cat("Converting genind to hierfstat format...\n")
hf_neutral <- genind2hierfstat(genind_neutral)

cat("Is genind_neutral genind?...\n")
is.genind(genind_neutral)
#TRUE

# ==============================================================
#   Generate bootstraps
# ==============================================================

cat("Generate bootstraps for fis...\n")
#Generate bootstraps - fine to use with basic.stats outputs. Halfway CI is not FIS. 
populations.structure_boot_fis<-boot.ppfis(dat=genind_neutral,nboot=1000,quant=c(0.025,0.975),diploid=TRUE,dig=4)
populations.structure_boot_fis


# ==============================================================
#   Generate FST
# ==============================================================

cat("Load .gl object...\n")
#on computer: 
#populations.structure.neutral.gl<-gi2gl(genind_neutral)
#save(populations.structure.neutral.gl, file = "populations.structure.neutral.gl.RData")

load("populations.structure.neutral.gl.RData")

cat("Generate bootstraps for FST...\n")
stamppFst_values <- stamppFst(populations.structure.neutral.gl, nboots = 1000, percent = 95, nclusters = 1) #NB takes a long time - saved to computer

cat("saving stampp...\n")

# Save the object to a file on your computer
save(stamppFst_values, file = "stamppFst_values_neutral.RData")

# Load the saved object from the file
#load("stamppFst_values.RData")

cat("printing STAMPP...\n")
options(max.print = 100)
stamppFst_values

cat("Printing FST...\n")
FST_pvalues<-stamppFst_values[["Pvalues"]]
FST_pvalues




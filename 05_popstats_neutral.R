# House keeping ####
#Clear memory
rm(list = ls())

#Set working directory
setwd("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/05_r_scripts")

# Generate citations ####
#For R
citation()

#For R Studio
if (!requireNamespace("rstudioapi", quietly = TRUE))
  install.packages("rstudioapi")
library(rstudioapi)

versionInfo()$citation

#For packages
#citation("pkgname")

# Install and load packages ####
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
library(devtools)

BiocManager::install("SNPRelate") #doesn't work normal way
library(SNPRelate)

if (!require("adegenet", quietly = TRUE))
  install.packages("adegenet")
library(adegenet)

if (!require("poppr", quietly = TRUE))
  install.packages("poppr")
library(poppr)

if (!require("hierfstat", quietly = TRUE))
  install.packages("hierfstat")
library(hierfstat)

if (!require("StAMPP", quietly = TRUE))
  install.packages("StAMPP")
library(StAMPP)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

if (!require("gdsfmt", quietly = TRUE))
  install.packages("gdsfmt")
library(gdsfmt)

if (!require("dartR", quietly = TRUE))
  install.packages("dartR")
library(dartR)

if (!require("dunn.test", quietly = TRUE))
  install.packages("dunn.test")
library(dunn.test)

if (!require("vcfR", quietly = TRUE))
  install.packages("vcfR")
library(vcfR)

# Basic stats ####

# read vcf 
vcf_neutral <- read.vcfR("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/03.65_HG_neutral.vcf")
genind_neutral <- vcfR2genind(vcf_neutral)

#used excel to add commas 
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

length(populations) == nInd(genind_neutral) # should be true 

genind_neutral@pop <- as.factor(populations)

# convert genind obj to hierfstat object 
hf_neutral <- genind2hierfstat(genind_neutral)
is.genind(genind_neutral)
#TRUE

#Takes long
basic_stats_neutral <- basic.stats(hf_neutral)

# Save the object to a file on your computer
save(basic_stats_neutral, file = "basic_stats_neutral_data.RData")

# Load the saved object from the file
#load("basic_stats_neutral_data.RData")

#write a standard error function with NA omitted 
#se <- function(x) sd(x)/sqrt(length(x)) 
#produces same outputs in outlier file as funciton above 
se <- function(x) {
  x <- na.omit(x)            # remove NAs
  sd(x) / sqrt(length(x))    # now safe
}

#print overall stats 
basic_stats_neutral
#$overall
#Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
#0.3143  0.2812  0.2813  0.0000  0.2813  0.0000  0.0001  0.0002 -0.1176  0.0001 

#Ho ####
#mean observed heterozygosities

Ho(hf_neutral) #takes long
#Cape Cross     False Bay      Kleinzee Lambert's Bay Pelican Point 
#    0.3149659     0.3141191     0.3154913     0.3110079     0.3160419 

Ho<-basic_stats_neutral[["Ho"]]
#View(Ho)

#Prep data for se calulation 
Ho_CC<-Ho[,1,drop=FALSE]
Ho_PP<-Ho[,5,drop=FALSE]
Ho_KZ<-Ho[,3,drop=FALSE]
Ho_LB<-Ho[,4,drop=FALSE]
Ho_FB<-Ho[,2,drop=FALSE]

#se calculation
se(Ho_CC) #0.0009474728
se(Ho_PP) #0.0009468468
se(Ho_KZ) #0.0009476008
se(Ho_LB) #0.00095045
se(Ho_FB) #0.0009735249

#Testing normality (Plots) 
#Skipped shapiro.test b/c data too long
hist(Ho_CC, main="Histogram of Ho_CC", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_CC)
qqline(Ho_CC, col = "red")
#CC not normal

hist(Ho_PP, main="Histogram of Ho_PP", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_PP)
qqline(Ho_PP, col = "red")
#PP not normal

hist(Ho_KZ, main="Histogram of Ho_KZ", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_KZ)
qqline(Ho_KZ, col = "red")
#KZ not normal

hist(Ho_LB, main="Histogram of Ho_LB", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_LB)
qqline(Ho_LB, col = "red")
#LB not normal

hist(Ho_FB, main="Histogram of Ho_FB", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Ho_FB)
qqline(Ho_FB, col = "red")
#FB not normal

#Non-parametric test for two multiple independent groups (Kruskal-Wallis Multiple Comparisons test):
#reformat data to long format 
# Assuming your data frame is called 'df' and population names are the column headers
# Convert from wide to long format
Ho_df <- as.data.frame(Ho)

Ho_long <- Ho_df %>%
  pivot_longer(cols = everything(), # Select all columns
               names_to = "population", # Name of the new 'population' column
               values_to = "Ho") # Name of the new 'Ho' column
# View the reshaped data
View(Ho_long)

kruskal.test(Ho ~ factor(population), data = Ho_long)
#Kruskal-Wallis rank sum test
#data:  Ho by factor(population)
#Kruskal-Wallis chi-squared = 22.352, df = 4, p-value = 0.0001705

#Posthoc test for Kruskal-Wallis (Dunn's test)
# Perform Dunn's test
dunn_results_Ho <- dunn.test(Ho_long$Ho, Ho_long$population, method = "bh")  # 'bh' for Benjamini-Hochberg adjustment

Kruskal-Wallis rank sum test

data: x and group
Kruskal-Wallis chi-squared = 22.3523, df = 4, p-value = 0


Comparison of x by group                            
(Benjamini-Hochberg)                              
Col Mean-|
  Row Mean |   Cape Cro   False Ba   Kleinzee   Lambert'
---------+--------------------------------------------
False Ba |   1.857253
         |     0.0527
         |
Kleinzee |  -0.443615  -2.300869
         |     0.3287    0.0214*
         |
Lambert' |   3.038561   1.181307   3.482177
|    0.0040*     0.1696    0.0012*
  |
  Pelican  |  -0.903328  -2.760582  -0.459712  -3.941889
|     0.2290    0.0072*     0.3587    0.0004*
  
  alpha = 0.05
Reject Ho if p <= alpha/2

# Create a data frame from the results
dunn_results_Ho_df <- data.frame(
  Comparison = dunn_results_Ho$comparisons,  # Pairwise comparisons
  Z = dunn_results_Ho$Z,                     # Z-statistics
  P.adjusted = dunn_results_Ho$P.adjusted    # Adjusted p-values
)

# Print the data frame
print(dunn_results_Ho_df)

#Comparison          Z   P.adjusted
#1         Cape Cross - False Bay  1.8572538 0.0527292109
#2          Cape Cross - Kleinzee -0.4436157 0.3286602252
#3           False Bay - Kleinzee -2.3008695 0.0213990080 #
#4     Cape Cross - Lambert's Bay  3.0385617 0.0039618409 #
#5      False Bay - Lambert's Bay  1.1813079 0.1696288915
#6       Kleinzee - Lambert's Bay  3.4821774 0.0012433854 #
#7     Cape Cross - Pelican Point -0.9033283 0.2289698166
#8      False Bay - Pelican Point -2.7605821 0.0072123076 #
#9       Kleinzee - Pelican Point -0.4597125 0.3587347563
#10 Lambert's Bay - Pelican Point -3.9418899 0.0004042104 #

# Significant differences: FB-KZ, CC-LB, KL_LB, FB-PP, LB-PP

#Hs ####
#mean gene diversities within population Hs
Hs(hf_neutral) #Takes long
#   Cape Cross     False Bay      Kleinzee Lambert's Bay Pelican Point 
# 0.2809507     0.2817042     0.2819653     0.2795921     0.2819953 

Hs<-basic_stats_neutral[["Hs"]]

#Prep data for se calulation 
Hs_CC<-Hs[,1,drop=FALSE]
Hs_PP<-Hs[,5,drop=FALSE]
Hs_KZ<-Hs[,3,drop=FALSE]
Hs_LB<-Hs[,4,drop=FALSE]
Hs_FB<-Hs[,2,drop=FALSE]

se(Hs_CC) 
se(Hs_PP) 
se(Hs_KZ) 
se(Hs_LB) 
se(Hs_FB) 

> se(Hs_CC) 
[1] 0.0006454861
> se(Hs_PP) 
[1] 0.0006419165
> se(Hs_KZ) 
[1] 0.0006449155
> se(Hs_LB) 
[1] 0.000653761
> se(Hs_FB) 
[1] 0.0006650558

#Testing normality (Plots) 
#Skipped shapiro.test b/c data too long
hist(Hs_CC, main="Histogram of Hs_CC", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_CC)
qqline(Hs_CC, col = "red")
#CC not normal

hist(Hs_PP, main="Histogram of Hs_PP", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_PP)
qqline(Hs_PP, col = "red")
#PP not normal

hist(Hs_KZ, main="Histogram of Hs_KZ", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_KZ)
qqline(Hs_KZ, col = "red")
#KZ almost normal

hist(Hs_LB, main="Histogram of Hs_LB", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_LB)
qqline(Hs_LB, col = "red")
#LB not normal

hist(Hs_FB, main="Histogram of Hs_FB", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(Hs_FB)
qqline(Hs_FB, col = "red")
#FB not normal

#Non-parametric test for two multiple independent groups (Kruskal-Wallis Multiple Comparisons test):
#reformat data to long format 
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
library(tidyr)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

# Assuming your data frame is called 'df' and population names are the column headers
# Convert from wide to long format
Hs_df <- as.data.frame(Hs)

Hs_long <- Hs_df %>%
  pivot_longer(cols = everything(), # Select all columns
               names_to = "population", # Name of the new 'population' column
               values_to = "Hs") # Name of the new 'Hs' column

# View the reshaped data
View(Hs)
View(Hs_long)

kruskal.test(Hs ~ factor(population), data = Hs_long)
#Kruskal-Wallis rank sum test
#data:  Hs by factor(population)
#Kruskal-Wallis chi-squared = 10.681, df = 4, p-value = 0.03039

#PostHoc test for Kruskal-Wallis (Dunn's test)
if (!require("dunn.test", quietly = TRUE))
  install.packages("dunn.test")
library(dunn.test)

#Dunn's test
dunn_results_Hs <- dunn.test(Hs_long$Hs, Hs_long$population, method = "bh")  # 'bh' for Benjamini-Hochberg adjustment
#The procedure of using Dunn's (z-)tests with Bonferroni-correction of the p-values is called "Dunn-Bonferroni"
# c("none", "bonferroni", "sidak", "holm", "hs", "hochberg", "bh", "by")

#Kruskal-Wallis rank sum test
#data: x and group
#Kruskal-Wallis chi-squared = 10.6811, df = 4, p-value = 0.03

#Comparison of x by group                            
#(Benjamini-Hochberg)                              
#Col Mean-|
#  Row Mean |   Cape Cro   False Ba   Kleinzee   Lambert'
---------+--------------------------------------------
#False Ba |  -2.301085
#         |     0.0535
#         |
#Kleinzee |  -1.159397   1.141710
#         |     0.2052     0.1811
#         |
#Lambert' |   0.692849   2.993921   1.852246
#|     0.2713    0.0138*     0.1067
#|
#  Pelican  |  -1.104713   1.196392   0.054683  -1.797562
#|     0.1683     0.2315     0.4782     0.0903

#alpha = 0.05
#Reject Ho if p <= alpha/2

# Create a data frame from the results
dunn_results_Hs_df <- data.frame(
  Comparison = dunn_results_Hs$comparisons,  # Pairwise comparisons
  Z = dunn_results_Hs$Z,                     # Z-statistics
  P.adjusted = dunn_results_Hs$P.adjusted    # Adjusted p-values
)

# Print the data frame
print(dunn_results_Hs_df)
Comparison           Z P.adjusted
#1         Cape Cross - False Bay -2.30108544 0.05346700 
#2          Cape Cross - Kleinzee -1.15939714 0.20524528
#3           False Bay - Kleinzee  1.14171045 0.18112457
#4     Cape Cross - Lambert's Bay  0.69284907 0.27133571
#5      False Bay - Lambert's Bay  2.99392128 0.01377085 #
#6       Kleinzee - Lambert's Bay  1.85224622 0.10665079 #
#7     Cape Cross - Pelican Point -1.10471392 0.16830223
#8      False Bay - Pelican Point  1.19639264 0.23154337
#9       Kleinzee - Pelican Point  0.05468323 0.47819542
#10 Lambert's Bay - Pelican Point -1.79756299 0.09030786 #

#Comparisons that are significant: FB-LB, KZ-LB, LB-PP

#FIS ####
#Inbreeding Coeficcient 
#subset data to find FIS using basic.stats
#poppr for popsub
if (!require("poppr", quietly = TRUE))
  install.packages("poppr")
library(poppr)

FIS_CC<-
  popsub(
    genind_neutral,
    sublist = "Cape Cross",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

FIS_PP<-
  popsub(
    genind_neutral,
    sublist = "Pelican Point",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

FIS_KZ<-
  popsub(
    genind_neutral,
    sublist = "Kleinzee",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

FIS_LB<-
  popsub(
    genind_neutral,
    sublist = "Lambert's Bay",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

FIS_FB<-
  popsub(
    genind_neutral,
    sublist = "False Bay",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

#Find FIS using basic.stats
basic.stats_CC <-basic.stats(FIS_CC,diploid=TRUE,digits=5)
basic.stats_CC
#$overall
#Ho       Hs       Ht      Dst      Htp     Dstp      Fst     Fstp      Fis     Dest 
#0.31497  0.28095  0.28095  0.00000      NaN      NaN  0.00000      NaN -0.12107      NaN 

FIS_CC #-0.12107

basic.stats_PP<-basic.stats(FIS_PP,diploid=TRUE,digits=5)
basic.stats_PP
#$overall
#Ho       Hs       Ht      Dst      Htp     Dstp      Fst     Fstp      Fis     Dest 
#0.31604  0.28200  0.28200  0.00000      NaN      NaN  0.00000      NaN -0.12073      NaN 

FIS_PP #-0.12073

basic.stats_KZ <-basic.stats(FIS_KZ,diploid=TRUE,digits=5)
basic.stats_KZ
#$overall
#Ho       Hs       Ht      Dst      Htp     Dstp      Fst     Fstp      Fis     Dest 
#0.31549  0.28197  0.28197  0.00000      NaN      NaN  0.00000      NaN -0.11890      NaN 

FIS_KZ #-0.11890

basic.stats_LB <-basic.stats(FIS_LB,diploid=TRUE,digits=5)
basic.stats_LB
#$overall
#Ho       Hs       Ht      Dst      Htp     Dstp      Fst     Fstp      Fis     Dest 
#0.31101  0.27959  0.27959  0.00000      NaN      NaN  0.00000      NaN -0.11236      NaN 

FIS_LB #-0.11236

basic.stats_FB <-basic.stats(FIS_FB,diploid=TRUE,digits=5)
basic.stats_FB
#$overall
#Ho       Hs       Ht      Dst      Htp     Dstp      Fst     Fstp      Fis     Dest 
#0.31412  0.28171  0.28171  0.00000      NaN      NaN  0.00000      NaN -0.11507      NaN 

FIS_FB #-0.11507

# I did this in the HPC 
#Generate bootstraps - fine to use with basic.stats outputs. Halfway CI is not FIS. 
populations.structure_boot_fis<-boot.ppfis(dat=genind_neutral,nboot=10000,quant=c(0.025,0.975),diploid=TRUE,dig=4)
View(populations.structure_boot_fis)
populations.structure_boot_fis

View(populations.structure_boot_fis)
pops <- levels(genind_outlier@pop)
pops

#$call
#boot.ppfis(dat = genind_neutral, nboot = 1000, quant = c(0.025,
                                                         0.975), diploid = TRUE, dig = 4)
#$fis.ci
#ll      hl
#1 -0.1249 -0.1170 CC 
#2 -0.1189 -0.1109 FB 
#3 -0.1229 -0.1151 KZ 
#4 -0.1165 -0.1084 LB 
#5 -0.1248 -0.1169 PP 

#NB check order of pops ! 
# HPC done 

#Subset datasets to prepare to test nomality 
#Can probably also use basic.stats_PP$FIS
FIS<-basic_stats_neutral[["Fis"]]
View(FIS)
View(basic_stats_neutral)

#Prep data for se calulation 
FIS_CC<-FIS[,1,drop=FALSE]
FIS_PP<-FIS[,5,drop=FALSE]
FIS_KZ<-FIS[,3,drop=FALSE]
FIS_LB<-FIS[,4,drop=FALSE]
FIS_FB<-FIS[,2,drop=FALSE]

#se calculation
se(FIS_CC) #NA
se(FIS_PP) #NA
se(FIS_KZ) #NA
se(FIS_LB) #NA
se(FIS_FB) #NA

#Testing normality (Plots) 
#Skipped shapiro.test b/c data too long
hist(FIS_CC, main="Histogram of FIS_CC", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(FIS_CC)
qqline(FIS_CC, col = "red")
#CC not normal

hist(FIS_PP, main="Histogram of FIS_PP", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(FIS_PP)
qqline(FIS_PP, col = "red")
#PP not normal

hist(FIS_KZ, main="Histogram of FIS_KZ", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(FIS_KZ)
qqline(FIS_KZ, col = "red")
#KZ not normal

hist(FIS_LB, main="Histogram of FIS_LB", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(FIS_LB)
qqline(FIS_LB, col = "red")
#LB not normal

hist(FIS_FB, main="Histogram of FIS_FB", xlab="Value", ylab="Frequency") #very right-skewed
qqnorm(FIS_FB)
qqline(FIS_FB, col = "red")
#FB not normal

#Non-parametric test for two multiple independent groups (Kruskal-Wallis Multiple Comparisons test):
#reformat data to long format 
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
library(tidyr)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

# Assuming your data frame is called 'df' and population names are the column headers
# Convert from wide to long format
FIS_df <- as.data.frame(FIS)

FIS_long <- FIS_df %>%
  pivot_longer(cols = everything(), # Select all columns
               names_to = "population", # Name of the new 'population' column
               values_to = "FIS") # Name of the new 'FIS' column
# View the reshaped data
View(FIS_long)

kruskal.test(FIS ~ factor(population), data = FIS_long)
#Kruskal-Wallis rank sum test
#data:  FIS by factor(population)
#Kruskal-Wallis chi-squared = 5.9031, df = 4, p-value = 0.2065

#Posthoc test for Kruskal-Wallis (Dunn's test)
if (!require("dunn.test", quietly = TRUE))
  install.packages("dunn.test")
library(dunn.test)

#DON'T NEED POST-HOC (p>0.05)
# Perform Dunn's test
dunn_results_FIS <- dunn.test(FIS_long$FIS, FIS_long$population, method = "bh")  # 'bh' for Benjamini-Hochberg adjustment
#The procedure of using Dunn's (z-)tests with Bonferroni-correction of the p-values is called "Dunn-Bonferroni"
# c("none", "bonferroni", "sidak", "holm", "hs", "hochberg", "bh", "by")
#Kruskal-Wallis rank sum test
#data: x and group
#Kruskal-Wallis chi-squared = 5.9031, df = 4, p-value = 0.21
#Comparison of x by group                            
#(Benjamini-Hochberg)                              
#Col Mean-|
#  Row Mean |   Cape Cro   False Ba   Kleinzee   Lambert'
---------+--------------------------------------------
  #False Ba |  -1.273456
  #         |     0.2536
#         |
#Kleinzee |   1.039496   2.282292
#         |     0.2986     0.1124
#         |
#Lambert' |  -0.624919   0.651022  -1.650900
#|     0.2956     0.3219     0.2469
#|
#  Pelican  |   0.151207   1.396004  -0.866118   0.761725
#|     0.4399     0.2712     0.3220     0.3187
#
#alpha = 0.05
#Reject Ho if p <= alpha/2

# Create a data frame from the results
dunn_results_FIS_df <- data.frame(
  Comparison = dunn_results_FIS$comparisons,  # Pairwise comparisons
  Z = dunn_results_FIS$Z,                     # Z-statistics
  P.adjusted = dunn_results_FIS$P.adjusted    # Adjusted p-values
)

# Print the data frame
print(dunn_results_FIS_df)

#No comparisons are significant. No differences! 

#FST ####
#FST from Stamp 
if (!require("StAMPP", quietly = TRUE))
  install.packages("StAMPP")
library(StAMPP)

# Dependancies for dartR (not installing without them)
BiocManager::install("SNPRelate") #doesn't work normal way
library(SNPRelate)

if (!require("gdsfmt", quietly = TRUE))
  install.packages("gdsfmt")
library(gdsfmt)

if (!require("dartR", quietly = TRUE))
  install.packages("dartR")
library(dartR)


#neutral 
#NB ran bootstraps in the HPC - see appropriate scripts 
#populations.structure.neutral.gl<-gi2gl(genind_neutral)
#save(populations.structure.neutral.gl, file = "populations.structure.neutral.gl.RData")

#stamppFst_values <- stamppFst(populations.structure.gl, nboots = 1000, percent = 95, nclusters = 1) #NB takes a long time - saved to computer

# Save the object to a file on your computer
#save(stamppFst_values, file = "stamppFst_values_neutral.RData")

# Load the saved object from the file
load("stamppFst_values_neutral.RData")

options(max.print = 100)
stamppFst_values
#View(stamppFst_values)

# looking at CIs
# Check the components.
names(stamppFst_values)

# Extract the bootstrap results.
bootstrap_results <- stamppFst_values[["Bootstraps"]]
#View(bootstrap_results)

# Show only the useful summary columns, excluding the 1,000 replicate columns.
FST_confidence_intervals <- bootstrap_results[, c(
  "Population1",
  "Population2",
  "Fst",
  "Lower bound CI limit",
  "Upper bound CI limit",
  "p-value"
)]

FST_confidence_intervals

# ==============================================================
#   Correct for multiple comparisons
# ==============================================================
# Extract the raw pairwise P-values.
FST_pvalues_raw <- stamppFst_values$Pvalues
#FST_pvalues_raw <- stamppFst_values[["Pvalues"]]
#FST_pvalues_raw
#View(FST_pvalues_raw)

# Identify the ten cells containing comparisons (exclude NAs)
tested_pairs <- !is.na(FST_pvalues_raw)
tested_pairs

# Copy the matrix and correct those ten P-values.
FST_pvalues_holm <- FST_pvalues_raw
FST_pvalues_holm #still orginal
FST_pvalues_holm[tested_pairs] <- p.adjust(
  FST_pvalues_raw[tested_pairs],
  method = "holm"
)

cat("Raw P-values:\n")
FST_pvalues_raw

cat("Holm-corrected P-values:\n")
FST_pvalues_holm

# Add the FDR corrected p-values back to your data frame 
stamppFst_values$Holm_P_Value <- FST_pvalues_holm

# Identify significant results
# You can set your desired significance level (e.g., 0.05) for the FDR-corrected p-values
alpha_level <- 0.05
stamppFst_values$Significant_Holm <- stamppFst_values$Holm_P_Value <= alpha_level

## Remove the incorrectly named FDR fields.
#stamppFst_values$Significant_FDR <- NULL
#stamppFst_values$FDR_P_Value <- NULL

# Check the remaining fields.
#names(stamppFst_values)
#print(stamppFst_values)

# View the results 
# Display the data frame with the new columns
print(stamppFst_values)

#$Fsts
#              Cape Cross     False Bay        Kleinzee Lambert's Bay Pelican Point
#Cape Cross                NA            NA              NA            NA            NA
#False Bay     -0.00001077735            NA              NA            NA            NA
#Kleinzee      -0.00007377278  0.0001872277              NA            NA            NA
#Lambert's Bay  0.00041168711 -0.0001392380  0.000360941768            NA            NA
#Pelican Point  0.00009720389  0.0003591899 -0.000003402574  0.0003114581            NA
#
#$Pvalues
#Cape Cross False Bay Kleinzee Lambert's Bay Pelican Point
#Cape Cross            NA        NA       NA            NA            NA
#False Bay          0.532        NA       NA            NA            NA
#Kleinzee           0.611     0.240       NA            NA            NA
#Lambert's Bay      0.036     0.712    0.074            NA            NA
#Pelican Point      0.327     0.077    0.461         0.097            NA
#
#$Bootstraps
#Population1 Population2 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48
#49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98
#[ reached 'max' / getOption("max.print") -- omitted 906 columns ]
#[ reached 'max' / getOption("max.print") -- omitted 10 rows ]
#
#$Holm_P_Value
#Cape Cross False Bay Kleinzee Lambert's Bay Pelican Point
#Cape Cross            NA        NA       NA            NA            NA
#False Bay           1.00        NA       NA            NA            NA
#Kleinzee            1.00     1.000       NA            NA            NA
#Lambert's Bay       0.36     1.000    0.666            NA            NA
#Pelican Point       1.00     0.666    1.000         0.679            NA
#
#$Significant_Holm
#Cape Cross False Bay Kleinzee Lambert's Bay Pelican Point
#Cape Cross            NA        NA       NA            NA            NA
#False Bay          FALSE        NA       NA            NA            NA
#Kleinzee           FALSE     FALSE       NA            NA            NA
#Lambert's Bay      FALSE     FALSE    FALSE            NA            NA
#Pelican Point      FALSE     FALSE    FALSE         FALSE            NA

#no FST comparisons are statistically significant

# Export the results to a new CSV file 
#this wont work becaue stamppFst_values is a list containing several differently shaped tables, so it cannot be written directly as one CSV. 
#write.csv(stamppFst_values, "stampp_fst_results_with_holm_correction_full.csv", row.names = FALSE)
#Create one row per population comparison instaed 

# Avoid scientific notation.
options(scipen = 999)

# Locate the ten pairwise comparisons.
tested_pairs <- !is.na(stamppFst_values$Pvalues)
positions <- which(tested_pairs, arr.ind = TRUE)

# Combine the relevant results into one table.
fst_results <- data.frame(
  Population_1 = rownames(stamppFst_values$Pvalues)[positions[, 1]],
  Population_2 = colnames(stamppFst_values$Pvalues)[positions[, 2]],
  FST = stamppFst_values$Fsts[positions],
  Raw_P_value = stamppFst_values$Pvalues[positions],
  Holm_P_value = stamppFst_values$Holm_P_Value[positions],
  Significant_Holm = stamppFst_values$Significant_Holm[positions]
)

print(fst_results)

#Population_1  Population_2             FST Raw_P_value Holm_P_value Significant_Holm
#1      False Bay    Cape Cross -0.000010777349       0.532        1.000            FALSE
#2       Kleinzee    Cape Cross -0.000073772779       0.611        1.000            FALSE
#3  Lambert's Bay    Cape Cross  0.000411687109       0.036        0.360            FALSE
#4  Pelican Point    Cape Cross  0.000097203886       0.327        1.000            FALSE
#5       Kleinzee     False Bay  0.000187227713       0.240        1.000            FALSE
#6  Lambert's Bay     False Bay -0.000139238043       0.712        1.000            FALSE
#7  Pelican Point     False Bay  0.000359189851       0.077        0.666            FALSE
#8  Lambert's Bay      Kleinzee  0.000360941768       0.074        0.666            FALSE
#9  Pelican Point      Kleinzee -0.000003402574       0.461        1.000            FALSE
#10 Pelican Point Lambert's Bay  0.000311458141       0.097        0.679            FALSE

write.csv(
  fst_results,
  "stampp_fst_results_with_holm_correction_neutral.csv",
  row.names = FALSE
)



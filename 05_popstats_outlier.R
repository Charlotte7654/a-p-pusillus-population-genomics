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
vcf_outlier <- read.vcfR("/Users/Char/Documents/Stellenbosch_University/Masters/ch2/bioinformatics/00_data/03.65_HG_outlier.vcf")
genind_outlier <- vcfR2genind(vcf_outlier)

# Create a population vector
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

length(populations) == nInd(genind_outlier) # should be true 

genind_outlier@pop <- as.factor(populations)

# convert genind obj to hierfstat object 
hf_outlier <- genind2hierfstat(genind_outlier)
is.genind(genind_outlier)
#TRUE

basic_stats_outlier <- basic.stats(hf_outlier)

# Save the object to a file on your computer
save(basic_stats_outlier, file = "basic_stats_outlier_data.RData")

#load("basic_stats_outlier_data.RData")

#write a standard error function with NA omitted 
#se <- function(x) sd(x)/sqrt(length(x)) 
#produces same outputs in outlier file as funciton above 
se <- function(x) {
  x <- na.omit(x)            # remove NAs
  sd(x) / sqrt(length(x))    # now safe
}

basic_stats_outlier
#$overall
#Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
#0.1944  0.1829  0.1862  0.0033  0.1870  0.0041  0.0177  0.0220 -0.0627  0.0050 

#Ho ####
#mean observed heterozygosities

Ho(hf_outlier) #takes long
#Cape Cross     False Bay      Kleinzee     Lambert's Bay Pelican Point 
#    0.1993944     0.1820191     0.2229904     0.1845985     0.1829491 

Ho<-basic_stats_outlier[["Ho"]]
View(Ho)

#Prep data for se calulation 
Ho_CC<-Ho[,1,drop=FALSE]
Ho_PP<-Ho[,5,drop=FALSE]
Ho_KZ<-Ho[,3,drop=FALSE]
Ho_LB<-Ho[,4,drop=FALSE]
Ho_FB<-Ho[,2,drop=FALSE]

#se calculation
se(Ho_CC) #0.009651682
se(Ho_PP) #0.01008328
se(Ho_KZ) #0.0102673
se(Ho_LB) #0.009875877
se(Ho_FB) #0.01062011

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
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
library(tidyr)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

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
#	Kruskal-Wallis rank sum test
#data:  Ho by factor(population)
#Kruskal-Wallis chi-squared = 18.729, df = 4, p-value = 0.0008885

#Interpretation: Since the p-value is much less than 0.05, 
#you reject the null hypothesis that the $\text{H}_O$ distributions are the same across all five populations. 
#This means at least one population is significantly different from the others.
#So, now we need to see where the diffs are 

#Posthoc test for Kruskal-Wallis (Dunn's test)
if (!require("dunn.test", quietly = TRUE))
  install.packages("dunn.test")
library(dunn.test)

# Perform Dunn's test
dunn_results_Ho <- dunn.test(Ho_long$Ho, Ho_long$population, method = "bh")  # 'bh' for Benjamini-Hochberg adjustment
#The procedure of using Dunn's (z-)tests with Bonferroni-correction of the p-values is called "Dunn-Bonferroni"
# c("none", "bonferroni", "sidak", "holm", "hs", "hochberg", "bh", "by")

#Comparison of x by group                            
#(Benjamini-Hochberg)                              
#Col Mean-|
#  Row Mean |   Cape Cro   False Ba   Kleinzee   Lambert'
---------+--------------------------------------------
#False Ba |   2.195442
#         |     0.0352
#         |
#Kleinzee |  -1.568167  -3.763610
#         |     0.0974    0.0008*
#         |
#Lambert' |   1.381542  -0.813900   2.949709
#|     0.1194     0.2598    0.0053*
#  |
#  Pelican  |   1.719468  -0.475974   3.287635   0.337925
#|     0.0855     0.3523    0.0025*     0.3677
#
#alpha = 0.05
#Reject Ho if p <= alpha/2

# Create a data frame from the results
dunn_results_Ho_df <- data.frame(
  Comparison = dunn_results_Ho$comparisons,  # Pairwise comparisons
  Z = dunn_results_Ho$Z,                     # Z-statistics
  P.adjusted = dunn_results_Ho$P.adjusted    # Adjusted p-values
)

# Print the data frame
print(dunn_results_Ho)

# Create a data frame from the desired components
dunn_table_Ho <- data.frame(
  Comparison = dunn_results_Ho$comparisons,
  Z_Statistic = dunn_results_Ho$Z,
  P_Value_Unadjusted = dunn_results_Ho$P,
  P_Value_Adjusted = dunn_results_Ho$P.adjusted
)

# Print the resulting table
print(dunn_table_Ho)

# The following comparisons are significant: CC-FB, FB-KZ, KZ-LB, KZ-PP
#     Comparison                    Z_Statistic  P_Value_Unadjusted  P_Value_Adjusted
#1         Cape Cross - False Bay   2.1954428       1.406593e-02     0.0351648149 #
#2          Cape Cross - Kleinzee  -1.5681674       5.842104e-02     0.0973683962
#3           False Bay - Kleinzee  -3.7636101       8.373895e-05     0.0008373895 #
#4     Cape Cross - Lambert's Bay   1.3815421       8.355617e-02     0.1193659630
#5      False Bay - Lambert's Bay  -0.8139007       2.078509e-01     0.2598136539 
#6       Kleinzee - Lambert's Bay   2.9497095       1.590364e-03     0.0053012148 #
#7     Cape Cross - Pelican Point   1.7194681       4.276459e-02     0.0855291757 
#8      False Bay - Pelican Point  -0.4759747       3.170462e-01     0.3522735594
#9       Kleinzee - Pelican Point   3.2876355       5.051629e-04     0.0025258143 #
#10 Lambert's Bay - Pelican Point   0.3379260       3.677095e-01     0.3677094788

#Interpretation (Final Conclusion): By comparing the Adjusted P-values to $\alpha = 0.05$, 
#you correctly identified four significant differences (marked with a # in your final table)
#This indicates that the $\text{H}_O$ in Kleinzee is significantly different from three other populations, 
#making it the most genetically distinct population in your sample based on this metric.
                                                                                                                                                 
#Hs ####
#mean gene diversities within population Hs
Hs(hf_outlier) #Takes long
#Cape Cross     False Bay      Kleinzee Lambert's Bay Pelican Point 
#    0.1872139     0.1729988     0.2070389     0.1763577     0.1710966 

Hs<-basic_stats_outlier[["Hs"]]

#Prep data for se calulation 
Hs_CC<-Hs[,1,drop=FALSE]
Hs_PP<-Hs[,5,drop=FALSE]
Hs_KZ<-Hs[,3,drop=FALSE]
Hs_LB<-Hs[,4,drop=FALSE]
Hs_FB<-Hs[,2,drop=FALSE]

#se calculation
se(Hs_CC) #0.008010025
se(Hs_PP) #0.00827669
se(Hs_KZ) #0.008068789
se(Hs_LB) #0.008382698
se(Hs_FB) #0.009035088

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
#Kruskal-Wallis chi-squared = 17.596, df = 4, p-value = 0.00148

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
#Kruskal-Wallis chi-squared = 17.5961, df = 4, p-value = 0

#Comparison of x by group                            
#(Benjamini-Hochberg)                              
#Col Mean-|
#  Row Mean |   Cape Cro   False Ba   Kleinzee   Lambert'
---------+--------------------------------------------
#False Ba |   1.884722
#         |     0.0743
#         |
#Kleinzee |  -1.687129  -3.571852
#         |     0.0763    0.0018*
#         |
#Lambert' |   1.201454  -0.683268   2.888583
#|     0.1640     0.3090    0.0064*
#  |
#  Pelican  |   1.693876  -0.190846   3.381005   0.492422
#|     0.0903     0.4243    0.0018*     0.3458

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
#                     Comparison       Z      P.adjusted
#1         Cape Cross - False Bay  1.8847225 0.074334140 
#2          Cape Cross - Kleinzee -1.6871296 0.076315361
#3           False Bay - Kleinzee -3.5718521 0.001772328 #
#4     Cape Cross - Lambert's Bay  1.2014542 0.163982188
#5      False Bay - Lambert's Bay -0.6832683 0.309023323
#6       Kleinzee - Lambert's Bay  2.8885838 0.006449682 #
#7     Cape Cross - Pelican Point  1.6938763 0.090288797 
#8      False Bay - Pelican Point -0.1908462 0.424323033
#9       Kleinzee - Pelican Point  3.3810059 0.001805525 #
#10 Lambert's Bay - Pelican Point  0.4924221 0.345789443

#Comparisons that are significant: FB-KZ, KZ-LB, KZ-PP

#FIS ####
#Inbreeding Coeficcient 
#subset data to find FIS using basic.stats
#poppr for popsub
if (!require("poppr", quietly = TRUE))
  install.packages("poppr")
library(poppr)

FIS_CC<-
  popsub(
    genind_outlier,
    sublist = "Cape Cross",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

FIS_PP<-
  popsub(
    genind_outlier,
    sublist = "Pelican Point",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

FIS_KZ<-
  popsub(
    genind_outlier,
    sublist = "Kleinzee",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

FIS_LB<-
  popsub(
    genind_outlier,
    sublist = "Lambert's Bay",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

FIS_FB<-
  popsub(
    genind_outlier,
    sublist = "False Bay",
    exclude = "NULL",
    blacklist = NULL,
    mat = NULL,
    drop = TRUE
  )

#Find FIS using basic.stats
basic.stats_PP<-basic.stats(FIS_PP,diploid=TRUE,digits=5)
basic.stats_PP
$overall
Ho       Hs       Ht      Dst      Htp     Dstp      Fst     Fstp      Fis     Dest 
0.18295  0.17110  0.17110  0.00000      NaN      NaN  0.00000      NaN -0.06926      NaN 

FIS_PP #-0.06926

basic.stats_CC <-basic.stats(FIS_CC,diploid=TRUE,digits=5)
basic.stats_CC
$overall
Ho       Hs       Ht      Dst      Htp     Dstp      Fst     Fstp      Fis     Dest 
0.19940  0.18722  0.18722  0.00000      NaN      NaN  0.00000      NaN -0.06505      NaN 

FIS_CC #-0.06505

basic.stats_KZ <-basic.stats(FIS_KZ,diploid=TRUE,digits=5)
basic.stats_KZ
$overall
Ho       Hs       Ht      Dst      Htp     Dstp      Fst     Fstp      Fis     Dest 
0.22299  0.20704  0.20704  0.00000      NaN      NaN  0.00000      NaN -0.07705      NaN 

FIS_KZ #-0.07705

basic.stats_LB <-basic.stats(FIS_LB,diploid=TRUE,digits=5)
basic.stats_LB
$overall
Ho       Hs       Ht      Dst      Htp     Dstp      Fst     Fstp      Fis     Dest 
0.18460  0.17636  0.17636  0.00000      NaN      NaN  0.00000      NaN -0.04673      NaN 

FIS_LB #-0.04673

basic.stats_FB <-basic.stats(FIS_FB,diploid=TRUE,digits=5)
basic.stats_FB
$overall
Ho       Hs       Ht      Dst      Htp     Dstp      Fst     Fstp      Fis     Dest 
0.18202  0.17300  0.17300  0.00000      NaN      NaN  0.00000      NaN -0.05214      NaN 

FIS_FB #-0.05214 

#Generate bootstraps - fine to use with basic.stats outputs. Halfway CI is not FIS. 
populations.structure_boot_fis<-boot.ppfis(dat=genind_outlier,nboot=1000,quant=c(0.025,0.975),diploid=TRUE,dig=4)
populations.structure_boot_fis

View(populations.structure_boot_fis)
pops <- levels(genind_outlier@pop)
pops

#$call
#boot.ppfis(dat = genind_outlier, nboot = 1000, quant = c(0.025, 
#                                                         0.975), diploid = TRUE, dig = 4)
#$fis.ci
#ll      hl
#1 -0.1120 -0.0198 CC
#2 -0.1024  0.0038 FB 
#3 -0.1232 -0.0237 KZ 
#4 -0.0977  0.0046 LB 
#5 -0.1187 -0.0152 PP

#Subset datasets to prepare to test nomality 
#Can probably also use basic.stats_PP$FIS
FIS<-basic_stats_outlier[["Fis"]]
View(FIS)
View(basic_stats_outlier)

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


#OUTLIER
populations.structure.gl<-gi2gl(genind_outlier)
stamppFst_values <- stamppFst(populations.structure.gl, nboots = 1000, percent = 95, nclusters = 1) #NB takes a long time - saved to computer

# Save the object to a file on your computer
save(stamppFst_values, file = "stamppFst_values_outlier.RData")

# Load the saved object from the file
load("stamppFst_values_outlier.RData")

#optional: print stampp values 
options(max.print = 30)
stamppFst_values
#View(stamppFst_values)

FST_pvalues<-stamppFst_values[["Pvalues"]]
FST_pvalues
#I think the decimals are shown as the least necessary 

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
#Cape Cross  False Bay   Kleinzee Lambert's Bay Pelican Point
#Cape Cross             NA         NA         NA            NA            NA
#False Bay     0.014752075         NA         NA            NA            NA
#Kleinzee      0.009573261 0.02044062         NA            NA            NA
#Lambert's Bay 0.023828236 0.03210104 0.01596217            NA            NA
#Pelican Point 0.028625617 0.04191970 0.01983426    0.01527553            NA
#
#$Pvalues
#Cape Cross False Bay Kleinzee Lambert's Bay Pelican Point
#Cape Cross            NA        NA       NA            NA            NA
#False Bay          0.009        NA       NA            NA            NA
#Kleinzee           0.004         0       NA            NA            NA
#Lambert's Bay      0.000         0        0            NA            NA
#Pelican Point      0.000         0        0             0            NA
#
#$Bootstraps
#Population1 Population2 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28
#[ reached 'max' / getOption("max.print") -- omitted 976 columns ]
#[ reached 'max' / getOption("max.print") -- omitted 10 rows ]
#
#$Holm_P_Value
#Cape Cross False Bay Kleinzee Lambert's Bay Pelican Point
#Cape Cross            NA        NA       NA            NA            NA
#False Bay          0.009        NA       NA            NA            NA
#Kleinzee           0.008         0       NA            NA            NA
#Lambert's Bay      0.000         0        0            NA            NA
#Pelican Point      0.000         0        0             0            NA

#$Significant_Holm
#Cape Cross False Bay Kleinzee Lambert's Bay Pelican Point
#Cape Cross            NA        NA       NA            NA            NA
#False Bay           TRUE        NA       NA            NA            NA
#Kleinzee            TRUE      TRUE       NA            NA            NA
#Lambert's Bay       TRUE      TRUE     TRUE            NA            NA
#Pelican Point       TRUE      TRUE     TRUE          TRUE            NA

#all FST comparisons are statistically significant

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

options(max.print = 100)
print(fst_results)

#    Population_1  Population_2         FST Raw_P_value Holm_P_value Significant_Holm
#1      False Bay    Cape Cross 0.014752075       0.009        0.009             TRUE
#2       Kleinzee    Cape Cross 0.009573261       0.004        0.008             TRUE
#3  Lambert's Bay    Cape Cross 0.023828236       0.000        0.000             TRUE
#4  Pelican Point    Cape Cross 0.028625617       0.000        0.000             TRUE
#5       Kleinzee     False Bay 0.020440620       0.000        0.000             TRUE
#6  Lambert's Bay     False Bay 0.032101037       0.000        0.000             TRUE
#7  Pelican Point     False Bay 0.041919703       0.000        0.000             TRUE
#8  Lambert's Bay      Kleinzee 0.015962170       0.000        0.000             TRUE
#9  Pelican Point      Kleinzee 0.019834263       0.000        0.000             TRUE
#10 Pelican Point Lambert's Bay 0.015275526       0.000        0.000             TRUE

write.csv(
  fst_results,
  "stampp_fst_results_with_holm_correction_outlier.csv",
  row.names = FALSE
)





 

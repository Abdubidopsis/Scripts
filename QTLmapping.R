# <Abdul Saboor Khan - QTL Mapping - Arabis F3 Drought experiment >
# <last updated on 07/25>


#useful links
#https://rqtl.org/tutorials/new_multiqtl.pdf
#https://rqtl.org/manual/qtl-manual.pdf
#https://rqtl.org/tutorials/rqtltour.pdf



#load libraries
library(ASMap)
library(qtl)
library(here)
library(dplyr)
library(tidyr)





#Prepare data----
#for successful matching during the read.cross process, it is necessary that there are the same IDs present in both files
here()

#load data frames
pheno <- read.csv("data/pheno2.csv")
geno <- read.csv2("data/map4_imputed_223_1086_gen.csv", header = T, check.names = FALSE) 
colnames(geno)[3] <- ""

#save IDs of both files
geno.id <- colnames(geno)[4:ncol(geno.id)] #"id", "", "", the first three columnnames are not needed here
pheno.id <- colnames(pheno)[2:ncol(pheno.id)]

#check common IDs
common_ids <- intersect(geno.id, pheno.id)

#missing in geno
pheno_not_in_geno <- setdiff(pheno.id, geno.id) #465_1

#missing in pheno
geno_not_in_pheno <- setdiff(geno.id, pheno.id) #"105" "107" "110" "115" "118" "19"  "330" "340" "347" "360" "361" "363" "367" "383" "-"   "52"  "62" "72"  "8"   "94"  "98" 

#filter both IDs
pheno_filtered <- pheno[pheno$id %in% common_ids, ]
geno_filtered <- geno[, c(1:3, which(colnames(geno) %in% common_ids))]

#confirm filtering
geno.id <- colnames(geno_filtered)[4:length(geno.id)] #"id", "", "", the first three columnnames are not needed here
pheno.id <- pheno_filtered$id
common_ids <- intersect(geno.id, pheno.id)
pheno_not_in_geno <- setdiff(pheno.id, geno.id) #empty
geno_not_in_pheno <- setdiff(geno.id, pheno.id) #empty

write.csv(pheno_filtered, "./data/pheno_filtered.csv", row.names = F, quote = FALSE)
write.csv(geno_filtered, "./data/geno_filtered.csv", row.names = F, quote = FALSE)



#read cross object----
cross <- read.cross(format = "csvsr", 
                   dir = "./data",
                   genfile = "geno_filtered.csv",
                   phefile = "pheno_filtered.csv",
                   sep =  ";",
                   na.strings = "NA",
                   check.names = FALSE,
                   genotypes = c("NN", "NS", "SS"))


write.cross(cross, format = "csvsr", 
            filestem = "./data/Cross_20250704",
            digits = 3)

write.cross(cross, format = "csv", 
            filestem = "./data/Cross_20250704",
            digits = 5)

#check cross----
cross <- jittermap(cross)
summary(cross)
plot(cross)
plotMissing(cross_2022j)


plotGeno(cross)
cross <- drop.nullmarkers(cross)
geno.image(cross, col = c("#66c2a5","#fc8d62","#8da0cb"))

#calculate conditional Qtl genotype probabilities
cross <- calc.genoprob(cross, step = 1)
plotGeno(cross)
plotInfo(cross)



##Survival----

####Single QTL scan: scanone----
#permutation tests with different numbers of permutation replicates
Surv.operm.5k <- scanone(cross, pheno.col = as.numeric(cross$pheno$survival=="alive"), method = "hk", model = "binary", n.perm = 5000)
summary(Surv.operm.5k, alpha = 0.05) # 5% 3.73
plot(Surv.operm.5k)


#genome scan with scanone
Surv.hk <- scanone(cross, pheno.col = as.numeric(cross$pheno$survival=="alive"), model = "binary")
summary(Surv.hk)

summary(Surv.hk, perms = Surv.operm.5k, alpha = 0.05, pvalues = T) 
max(Surv.hk) #1_513_11962149   1  12 3.02

# There were no LOD peaks above the threshold.



require(iRF)
require(dplyr)

setwd("/exeh_3/kai/data/heterogenous_treatment_effects")

dataset <- read.csv('TCGA_survival_CDR.csv', header = TRUE)

# some clearning
dataset <- dataset[dataset$PFI.time!=0,]
dataset[dataset == "#N/A"] <- NA
dataset <- subset(dataset, !is.na(type))

# specify cancer type here
outcome <- as.factor(as.numeric(dataset$type == "BRCA"))
patients <- cbind(outcome, dataset)

# **********************************************************************************
#                                mutation data
# **********************************************************************************
load("TCGA_PanCanAtlas_No_Of_Mutations_NoSilentMut.rda")
bcr_patient_barcode <- substr(mut_adj.rev$Tumor_Sample_Barcode, start=1, stop=12)
mut_adj.rev <- cbind(bcr_patient_barcode, mut_adj.rev)

# merge patients with mutitation dataset by "bcr_patient_barcode"
merged.dataset <- inner_join(patients, mut_adj.rev, by="bcr_patient_barcode")

# find mutations with mutation count greater than the defined mutation.threhold
pos <- which(colnames(merged.dataset) == "5S_rRNA")
mutation.mat <- as.matrix(merged.dataset[, pos:ncol(merged.dataset)])
# mutation.count <- apply(mutation.mat, 2, function(x) sum(x != 0))

save(fitte.obj.RData')

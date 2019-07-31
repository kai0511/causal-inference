# require(RSarules)
require(dplyr)
require(data.table)
require(methods)

# source('/exeh_3/kai/R.code/HTE/RSarules.R')
source('/exeh_3/kai/R.code/HTE/RSarules_1.0.R')

setwd("/exeh_3/kai/data/heterogenous_treatment_effects")

patients.info <- read.csv('TCGA_survival_CDR.csv', header = TRUE)

# some clearning
patients.info <- patients.info[patients.info$PFI.time!=0,]
patients.info[patients.info == "#N/A"] <- NA
patients.info <- subset(patients.info, !is.na(PFI.time) & !is.na(PFI) & !is.na(age_at_initial_pathologic_diagnosis))

# load mutation data
load("TCGA_PanCanAtlas_No_Of_Mutations_NoSilentMut.rda")
bcr_patient_barcode <- substr(mut_adj.rev$Tumor_Sample_Barcode, start=1, stop=12)
mut_adj.rev <- cbind(bcr_patient_barcode, mut_adj.rev)

# merge dataset with mutitation dataset by "bcr_patient_barcode"
merged <- inner_join(patients.info, mut_adj.rev, by="bcr_patient_barcode")

outcomes <- as.numeric(merged$type ==  'BRCA')
merged <- cbind(merged, outcomes)

pos <- which(colnames(merged) == "5S_rRNA")

# feature selection before running model
cases <- subset(merged, merged$outcomes == 1)
mutation.count <- apply(cases[, pos:ncol(cases)], 2, function(x) sum(!(x == 0)))
selected_features <- mutation.count >= 30 

mutation.data <- as.matrix(merged[, pos: ncol(merged)])[, selected_features]
mutation.data[!mutation.data == 0] = 1

# M stands for the numb of association rule to be sampled.
# ig is the value for tuning parameter
M = 10000 
ig = 500000 
rsar.obj <- RSarules(mutation.data, M = M, ig = ig, rhs = ncol(mutation.data))

save(rsar.obj, file = 'rsar.obj.RData')
print('mission done!')

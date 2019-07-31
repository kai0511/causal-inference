require(dplyr)

setwd('/Users/zhaokai/Documents/projects/HTE/dataset/Metabric_BRAC')

patients <- read.table('patientData.txt', sep = '\t', header = T)
patients$daysFollowup <- as.Date(patients$dateLastFollowup, format="%Y-%m-%d") - as.Date(patients$dateOfDiagnosis, format="%Y-%m-%d")


# mutation <- read.table('somaticMutations.txt', sep = '\t', header = T)

sm <- read.table('somaticMutations_incNC.txt', sep = '\t', header = T)
sampMutation <- cbind(sm$sample, sm$gene)
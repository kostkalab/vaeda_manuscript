#! /usr/bin/Rscript

suppressMessages(library('SingleCellExperiment'))
suppressMessages(library('scds'))
suppressMessages(library('reticulate'))

args = commandArgs(trailingOnly=TRUE)

scs <- import("scipy.sparse")

pth <- args[1]
mtx <- t(scs$load_npz(pth))

sce = SingleCellExperiment(assays=list(counts=mtx))
sce <- scds::cxds(sce)
sce <- scds::bcds(sce)
sce <- scds::cxds_bcds_hybrid(sce)

CD  = colData(sce)

cxds_scores = CD$cxds_score
bcds_scores = CD$bcds_score
hybrid_scores = CD$hybrid_score

save_pth <- substr(pth,1,nchar(pth)-4)
print('HERE')
print(save_pth)

write.table(cxds_scores, paste0(save_pth, '_cxds.csv'), quote = F, sep = ",", row.names = F, col.names = F)
write.table(bcds_scores, paste0(save_pth, '_bcds.csv'), quote = F, sep = ",", row.names = F, col.names = F)
write.table(hybrid_scores, paste0(save_pth, '_hybrid.csv'), quote = F, sep = ",", row.names = F, col.names = F)

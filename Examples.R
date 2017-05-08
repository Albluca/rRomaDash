library(rRomaDash)

rRomaDash(Interactive = FALSE)

Data <- readRDS("/Users/newmac-luca/Google Drive/Datasets/Mosaic/rRoma-Apoptosis.rds")
Data$Groups <- factor(as.character(Data$Groups), levels = c("7", "7+2", "7+3", "7+4", "7+7", "7+10", "7+15"))
names(Data$Groups) <- colnames(Data$ExpMat)
saveRDS(Data, "/Users/newmac-luca/Google Drive/Datasets/Mosaic/rRoma-Apoptosis.rds")


Data <- readRDS("/Users/newmac-luca/Google Drive/Datasets/Mosaic/rRoma-CellCycleMap.rds")
Data$Groups <- factor(as.character(Data$Groups), levels = c("7", "7+2", "7+3", "7+4", "7+7", "7+10", "7+15"))
names(Data$Groups) <- colnames(Data$ExpMat)
saveRDS(Data, "/Users/newmac-luca/Google Drive/Datasets/Mosaic/rRoma-CellCycleMap.rds")


Data <- readRDS("/Users/newmac-luca/Google Drive/Datasets/Mosaic/rRoma-DNARepair.rds")
Data$Groups <- factor(as.character(Data$Groups), levels = c("7", "7+2", "7+3", "7+4", "7+7", "7+10", "7+15"))
names(Data$Groups) <- colnames(Data$ExpMat)
saveRDS(Data, "/Users/newmac-luca/Google Drive/Datasets/Mosaic/rRoma-DNARepair.rds")


Data <- readRDS("/Users/newmac-luca/Google Drive/Datasets/Mosaic/rRoma-Hallmarks.rds")
Data$Groups <- factor(as.character(Data$Groups), levels = c("7", "7+2", "7+3", "7+4", "7+7", "7+10", "7+15"))
names(Data$Groups) <- colnames(Data$ExpMat)
saveRDS(Data, "/Users/newmac-luca/Google Drive/Datasets/Mosaic/rRoma-Hallmarks.rds")


Data <- readRDS("/Users/newmac-luca/Google Drive/Datasets/Mosaic/rRoma-Biocarta.rds")
Data$Groups <- factor(as.character(Data$Groups), levels = c("7", "7+2", "7+3", "7+4", "7+7", "7+10", "7+15"))
names(Data$Groups) <- colnames(Data$ExpMat)
saveRDS(Data, "/Users/newmac-luca/Google Drive/Datasets/Mosaic/rRoma-Biocarta.rds")



Data <- readRDS("/Users/newmac-luca/Google Drive/Datasets/Mosaic/rRoma-InfoSigMap.rds")
Data$Groups <- factor(as.character(Data$Groups), levels = c("7", "7+2", "7+3", "7+4", "7+7", "7+10", "7+15"))
names(Data$Groups) <- colnames(Data$ExpMat)
saveRDS(Data, "/Users/newmac-luca/Google Drive/Datasets/Mosaic/rRoma-InfoSigMap.rds")

SubSet_rRomaRDS(SourceFile = "/Users/newmac-luca/Google Drive/Datasets/Mosaic/rRoma-InfoSigMap.rds",
                TargetFile = "/Users/newmac-luca/Google Drive/Datasets/Mosaic/rRoma-InfoSigMap_.rds")



library(readr)
GSE72056_melanoma_single_cell_revised_v2 <- read_delim("~/Google Drive/Datasets/Tirosh et al - Human oligodendroglioma immune cells/GSE72056_melanoma_single_cell_revised_v2.txt",
                                                       "\t", escape_double = FALSE, trim_ws = TRUE)


GSE72056_melanoma_single_cell_revised_v2 <- GSE72056_melanoma_single_cell_revised_v2[-c(1:3),]

GSE72056_melanoma_single_cell_revised_v2 <- GSE72056_melanoma_single_cell_revised_v2[(rowSums(GSE72056_melanoma_single_cell_revised_v2[,-1]>0) > 10), ]

GSE72056_melanoma_single_cell_revised_v2 <-
  GSE72056_melanoma_single_cell_revised_v2[!duplicated(GSE72056_melanoma_single_cell_revised_v2[,1]),]

GSE72056_melanoma_single_cell_revised_v2[,-1] <- log10(GSE72056_melanoma_single_cell_revised_v2[,-1] + 1)

write_delim(GSE72056_melanoma_single_cell_revised_v2, delim = "\t",
            path = "~/Google Drive/Datasets/Tirosh et al - Human oligodendroglioma immune cells/CleanedData.txt")



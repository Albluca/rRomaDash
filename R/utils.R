#' Subset an rRomaDash save file
#'
#' @param SourceFile the source RDS file
#' @param TargetFile the target RDS file
#' @param KeyWord the keywork to use to subset the original save file
#'
#' @return
#' @export
#'
#' @examples
SubSet_rRomaRDS <- function(SourceFile, TargetFile, KeyWord) {

  SourceRDS <- readRDS(SourceFile)

  Selected <- grep(pattern = KeyWord, x = rownames(SourceRDS$RomaData$ModuleMatrix), ignore.case = TRUE)

  if(length(Selected) == 0){
    print("No geneset found")
  }

  print(paste(length(Selected), "genesets found"))

  TargetRoma <- SourceRDS$RomaData

  TargetRoma$ModuleMatrix <- TargetRoma$ModuleMatrix[Selected, ]
  TargetRoma$ProjMatrix <- TargetRoma$ProjMatrix[Selected, ]
  TargetRoma$ModuleSummary <- TargetRoma$ModuleSummary[Selected]
  TargetRoma$WeigthList <- TargetRoma$WeigthList[Selected]
  TargetRoma$PVVectMat <- TargetRoma$PVVectMat[Selected, ]
  TargetRoma$OutLiersList <- TargetRoma$OutLiersList[Selected]

  TargetRDS <- SourceRDS
  TargetRDS$RomaData <- TargetRoma

  saveRDS(TargetRDS, TargetFile)

}

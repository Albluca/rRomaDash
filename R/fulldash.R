#' Title
#'
#' @param RomaData
#' @param ExpMat
#' @param Groups
#'
#' @return
#'
#' @examples
Initialize <- function(RomaData, ExpMat, Groups){

  if(is.null(RomaData)){
    return(list(FoundSamp = NULL, Groups = Groups, AddInfo = NULL,
                SelList = list(), SelListAF = list(), GSList = list(),
                PCAProj = NULL, ProcessedSamples = NULL))
  }

  ProcessedSamples <- colnames(RomaData$SampleMatrix)

  if(!is.null(Groups)){

    FoundSamp <- intersect(colnames(RomaData$SampleMatrix), names(Groups))
    Groups <- Groups[FoundSamp]

    if(length(FoundSamp) > 2){
      if(!is.factor(Groups)){
        Groups <- as.factor(Groups)
      }
      AddInfo <- data.frame(Groups = Groups)
    } else {
      AddInfo <- NULL
    }
    SelList <- list("By sample" = "sample", "By group" = "group")
    SelListAF <- list("Mean" = "mean", "Median" = "median", "Std. dev." = "sd", "IQR" = "IQR", "mad" = "mad")

  } else {

    Groups <- rep("N/A", length(ProcessedSamples))
    names(Groups) <- ProcessedSamples
    FoundSamp = NULL

    AddInfo <- NULL
    SelList <- list("By sample" = "sample")
    SelListAF <- list()

  }

  if(length(unique(Groups))<2){
    SelList <- list("By sample" = "sample")
    SelListAF <- list()
  }

  if(!is.null(RomaData)){

    GSList <- as.list(1:nrow(RomaData$PVVectMat))
    names(GSList) <- as.list(rownames(RomaData$SampleMatrix))
    GSList <- GSList[order(names(GSList))]

    PCAProj <- irlba::prcomp_irlba(t(ExpMat[,ProcessedSamples]), 2, retx = TRUE)$x
    rownames(PCAProj) <- ProcessedSamples
    colnames(PCAProj) <- c("PC1", "PC2")

  } else {
    GSList <- list()
    PCAProj <- NULL
  }

  return(list(FoundSamp = FoundSamp, Groups = Groups, AddInfo = AddInfo,
              SelList = SelList, SelListAF = SelListAF, GSList = GSList,
              PCAProj = PCAProj, ProcessedSamples = ProcessedSamples))
}












#' Launch the rRoma dashboard
#'
#' @param RomaData
#' @param ExpMat
#' @param Groups
#' @param Interactive boolean, scalar. Should interactivity (via Plotly) be enabled?
#'
#' @return
#' @export
#'
#' @examples
rRomaDash <- function(RomaData = NULL,
                      ExpMat = NULL,
                      Groups = NULL,
                      Interactive = FALSE) {

  # library(rRoma)
  # library(shiny)
  # library(shinyFiles)

   # preprocess data ---------------------------------------------------------

  if(R.utils::isPackageLoaded("plotly")){
    print("Detaching plotly.")
    detach("package:plotly", unload=TRUE)
  }

  if(Interactive){
    library(plotly)
  }



  # if(Interactive){
  #   print("Using plotly. This can cause problems on some systems. Try setting 'Interactive = FALSE' if errors are encountered")
  # } else {
  #
  # }

  StartTime <- Sys.time()

  TimeVect <- rep(StartTime, 3)
  names(TimeVect) <- c("Upload", "Local", "Analysis")

  InitReturn <- Initialize(RomaData, ExpMat, Groups)

  FoundSamp <- InitReturn$FoundSamp
  Groups <- InitReturn$Groups
  AddInfo <- InitReturn$AddInfo
  SelList <- InitReturn$SelList
  SelListAF <- InitReturn$SelListAF
  GSList <- InitReturn$GSList
  PCAProj <- InitReturn$PCAProj
  ProcessedSamples <- InitReturn$ProcessedSamples

  GeneList <- list()

  tSNEProj <- PCAProj

  InternalGMTList <- list(
    "Molecular signature DB (v6.0)" = "MsigDB",
    "ACSN Globlal Map (v1.1)" = "ACSN_Global",
    "ACSN Apoptosis Map (v1.1)" = "ACSN_Apoptosis",
    "ACSN Cell Cycle Map (v1.1)" = "ACSN_CellCycle",
    "ACSN DNA Repair Map (v1.1)" = "ACSN_DNARepair",
    "ACSN EMT Map (v1.1)" = "ACSN_EMT",
    "ACSN Survival Map (v1.1)" = "ACSN_Survival",
    "InfoSigMap Conseved" = "InfoSig_Conserved",
    "InfoSigMap" = "InfoSig_Informative"
  )

  MapList <- list(
    "Atlas of Cancer Signalling Network global map" = "https://acsn.curie.fr/navicell/maps/acsn/master/index.php",
    "Apoptosis and mitochondria metabolism map" = "https://acsn.curie.fr/navicell/maps/apoptosis/master/index.php",
    "Cell survival map" = "https://acsn.curie.fr/navicell/maps/survival/master/index.php",
    "EMT and cell motility map" = "https://acsn.curie.fr/navicell/maps/emtcellmotility/master/index.php",
    "Cell cycle map" = "https://acsn.curie.fr/navicell/maps/cellcycle/master/index.php",
    "DNA repair map" = "https://acsn.curie.fr/navicell/maps/dnarepair/master/index.php",
    "InfoSigMap" = "https://navicell.curie.fr/navicell/newtest/maps/infosigmap/master/index.php"
  )

  nProcList <- as.list(as.character(1:32))
  names(nProcList) <- as.character(1:32)

  InitialSave <- paste(getwd(), '/', 'rRoma-', Sys.Date(), '.rds', sep = "")

  # define ui ---------------------------------------------------------

  ui <- navbarPage("rRoma dashboard",

             # Perform analysis (Top Tab 1) ---------------------------------------------------------
             tabPanel("Analyze Data",

                      pageWithSidebar(

                        # Application title
                        headerPanel(""),

                        # Sidebar with a slider input
                        sidebarPanel(
                          actionButton("doROMA", "Execute rROMA"),
                          hr(),
                          htmlOutput("ROMAOut"),
                          htmlOutput("ROMAOut2")
                        ),

                        # Show a plot of the generated distribution
                        mainPanel(
                          tabsetPanel(

                            # Data input ---------------------------------------------------------

                            tabPanel("Input",
                                     fluidPage(

                                       fluidRow(
                                         titlePanel("Expression matrix"),
                                         column(8,
                                                fileInput("expmatfile", "Choose an expression matrix (TSV file)", accept = c(".tsv", ".txt"))
                                                )
                                       ),

                                       fluidRow(
                                         titlePanel("Sample groups"),
                                         column(8,
                                                fileInput("groupfile", "Choose a group matrix (TSV file)", accept = c(".tsv", ".txt"))
                                         ),
                                         column(4,
                                                checkboxInput("usegroups", "Use groups", TRUE)
                                         )
                                       ),

                                       fluidRow(
                                         titlePanel("Geneset list"),
                                         column(6,

                                                selectInput("gmtsrc", "Geneset source:",
                                                            list("Local file" = "File",
                                                                 "Internal DB" = "Internal")),

                                                conditionalPanel(
                                                  condition="input.gmtsrc == 'File'",
                                                  fileInput("gmtfile", "Choose a GMT file", accept = c(".gmt", ".txt"))
                                                ),

                                                conditionalPanel(
                                                  condition="input.gmtsrc == 'Internal'",
                                                  selectInput("gmtlist", "Available geneset list:", InternalGMTList)
                                                )

                                                ),

                                         column(6,
                                                textInput("msigkw", "Keywords", ""),
                                                checkboxInput("msigkwall", "search all keywords", FALSE),
                                                checkboxInput("loadwei", "load weights", FALSE)
                                         ),

                                         column(12,
                                                actionButton("searchDB", "Apply")
                                         )

                                       ),

                                       fluidRow(
                                         titlePanel("Available Genesets:"),
                                         column(12,
                                                dataTableOutput("PrintGeneSets"))
                                         )
                                       )

                                     ),

                            # Parameters ---------------------------------------------------------

                            tabPanel("Parameters",

                                     fluidPage(title = "Base parameters",
                                               titlePanel("Base parameters"),
                                                column(6, selectInput("par_FixedCenter", "FixedCenter",
                                                                      list("TRUE" = "TRUE", "FALSE" = "FALSE"), selected = "FALSE"),
                                                       selectInput("par_PCSignMode", "PCSignMode",
                                                                   list("none" = "none",
                                                                        "PreferActivation" = "PreferActivation",
                                                                        "UseAllWeights" = "UseAllWeights",
                                                                        "UseKnownWeights" = "UseKnownWeights",
                                                                        "CorrelateAllWeightsByGene" = "CorrelateAllWeightsByGene",
                                                                        "CorrelateKnownWeightsByGene" = "CorrelateKnownWeightsByGene"),
                                                                   selected = "CorrelateAllWeightsByGene"),
                                                       textInput("par_nSamples", "nSamples", "100"),
                                                       textInput("par_GeneOutThr", "GeneOutThr", "5")
                                                       ),
                                                column(6,
                                                       selectInput("par_UseParallel", "UseParallel",
                                                                      list("TRUE" = "TRUE", "FALSE" = "FALSE"), selected = "TRUE"),
                                                       conditionalPanel(condition="input.par_UseParallel == 'TRUE'",
                                                         selectInput("par_nCores", "nCores", nProcList,
                                                                     selected = as.character(parallel::detectCores() - 1)),
                                                         selectInput("par_ClusType", "ClusType",
                                                                     list("PSOCK", "FORK"), selected = "PSOCK")
                                                       )

                                                       )
                                     ),

                                     fluidPage(title = "Advanced parameters",
                                               titlePanel("Advanced parameters"),
                                                column(4, selectInput("par_UseWeigths", "UseWeigths",
                                                                      list("TRUE" = "TRUE", "FALSE" = "FALSE"), selected = "FALSE"),
                                                       selectInput("par_FullSampleInfo", "FullSampleInfo",
                                                                   list("TRUE" = "TRUE", "FALSE" = "FALSE"), selected = "FALSE"),
                                                       selectInput("par_centerData", "centerData",
                                                                   list("TRUE" = "TRUE", "FALSE" = "FALSE"), selected = "TRUE"),
                                                       selectInput("par_GeneSelMode", "GeneSelMode",
                                                                   list("All" = "All", "Others" = "Others"), selected = "TRUE"),
                                                       textInput("par_Ncomp", "Ncomp", "5"),
                                                       textInput("par_DefaultWeight", "DefaultWeight", "1")
                                                       ),

                                                column(4, selectInput("par_SampleFilter", "SampleFilter",
                                                                      list("TRUE" = "TRUE", "FALSE" = "FALSE"), selected = "TRUE"),
                                                       selectInput("par_ExpFilter", "ExpFilter",
                                                                   list("TRUE" = "TRUE", "FALSE" = "FALSE"), selected = "FALSE"),
                                                       selectInput("par_MoreInfo", "MoreInfo",
                                                                   list("TRUE" = "TRUE", "FALSE" = "FALSE"), selected = "FALSE"),
                                                       selectInput("par_GeneOutDetection", "GeneOutDetection",
                                                                   list("L1OutVarPerc" = "L1OutVarPerc",
                                                                        "L1OutVarDC" = "L1OutVarDC",
                                                                        "L1OutExpOut" = "L1OutExpOut",
                                                                        "L1OutSdMean" = "L1OutSdMean"), selected = "L1OutExpOut"),
                                                       textInput("par_OutGeneNumber", "OutGeneNumber", "5"),
                                                       textInput("par_OutGeneSpace", "OutGeneSpace", "NULL")
                                                       ),

                                                column(4, textInput("par_MinGenes", "MinGenes", "10"),
                                                       textInput("par_MaxGenes", "MaxGenes", "500"),
                                                       textInput("par_ApproxSamples", "ApproxSamples", "5"),
                                                       textInput("par_PCSignThr", "PCSignThr", "NULL"),
                                                       selectInput("par_CorMethod", "CorMethod",
                                                                   list("pearson" = "pearson",
                                                                        "kendall" = "kendall",
                                                                        "spearman" = "spearman"), selected = "pearson"),
                                                       selectInput("par_PCAType", "PCAType",
                                                                   list("DimensionsAreGenes" = "DimensionsAreGenes",
                                                                        "DimensionsAreSamples" = "DimensionsAreSamples"),
                                                                   selected = "DimensionsAreGenes")
                                                       )
                                     )


                                     )

                          )
                          )
                        )
                      ),

             # Data summary (Top Tab 2) ---------------------------------------------------------
             tabPanel("Summarize Info",

                      fluidPage(
                        titlePanel(""),

                        fluidRow(
                          tabsetPanel(
                            tabPanel("Parameters", dataTableOutput("PrintPar")),
                            tabPanel("Genesets", dataTableOutput("PrintGMT")),
                            tabPanel("Groups", dataTableOutput("PrintGroup"))
                          )
                        )
                      )

                    ),

             # Results  (Top Tab 3) ---------------------------------------------------------
             tabPanel("Visualize Results",

                      pageWithSidebar(

                        # Application title
                        headerPanel(NULL),

                        # sidebar ---------------------------------------------------------

                        sidebarPanel(
                          conditionalPanel(
                            condition="input.ResTabs == 'Modules'",
                            selectInput("prjt", "Projection type:",
                                        list("PCA" = "PCA", "tSNE" = "tSNE"))
                          ),

                          conditionalPanel(
                            condition="input.ResTabs == 'Modules' && input.prjt == 'tSNE'",
                            sliderInput("perp", "tSNE perplexity:",
                                        min = 0,  max = 100,  value = 0, step = 1),
                            numericInput("initial_dims", "Initial dimensions:", min = 2, max = NA, step = 1, value = 50),
                            actionButton("dotSNE", "Compute tSNE"),
                            hr()
                          ),

                          conditionalPanel(
                            condition="input.ResTabs == 'Modules'",
                            selectInput("boxcomp", "Group comparison:",
                                        c("Show all", "Show none", "Show significant")),
                            hr()
                          ),

                          conditionalPanel(
                            condition="input.ResTabs == 'Modules' || input.ResTabs == 'Gene contribution'",
                            selectInput("gs", "Geneset:",
                                        GSList),
                            htmlOutput("info"),
                            hr()
                          ),

                          conditionalPanel(
                            condition="input.ResTabs == 'Heatmap' || input.ResTabs == 'Correlation'",
                            selectInput("htype", "Heatmap type:", SelList)
                          ),

                          conditionalPanel(
                            condition="(input.ResTabs == 'Heatmap' || input.ResTabs == 'Correlation') && input.htype == 'group'",
                            selectInput("aggfun", "Aggregating function:", SelListAF)
                          ),


                          conditionalPanel(
                            condition="input.ResTabs == 'Correlation' || input.ResTabs == 'Gene contribution'",
                            selectInput("cortype", "Correlation method:",
                                        list("Pearson" = "pearson",
                                             "Kendall" = "kendall",
                                             "Spearman" = "spearman"))
                          ),

                          conditionalPanel(
                            condition="input.ResTabs == 'Gene contribution'",
                            selectInput("corlelvel", "Confidence level:",
                                        list(".95", ".99", ".999", ".9999")),
                            actionButton("doCorr", "Compute correlations"),
                            hr()
                          ),

                          conditionalPanel(
                            condition="input.ResTabs == 'Heatmap'",
                            checkboxInput("gscol", "Samples on columns", FALSE),
                            checkboxInput("saclus", "Cluster samples / groups", FALSE)
                          ),

                          conditionalPanel(
                            condition="input.ResTabs == 'Heatmap' || input.ResTabs == 'Correlation'",
                            checkboxInput("gsclus", "Cluster genesets", FALSE),
                            selectInput("GSOrdeMode", "Geneset information:",
                                        choices = c("None", "Gene number",
                                                    "Overdispersion pv", "Underdispersion pv",
                                                    "Overcoordination pv", "Undercoordination pv",
                                                    "Overexpression pv", "Underexpression pv")),
                            hr()
                          ),

                          selectInput("disp", "Dispersion filter:",
                                      list("Overdispersed" = "Over",
                                           "Underdispersed" = "Under",
                                           "None" = "None")),
                          sliderInput("pdisp", "Log 10 p-value threshold:",
                                      max = -1.3,  min = -5,  value = -1.3, step = .1),
                          hr(),

                          selectInput("coord", "Coordination filter:",
                                      list("Overcoordinated" = "Over",
                                           "Undercoordinated" = "Under",
                                           "None" = "None"), selected = "None"),
                          sliderInput("pcoord", "Log 10 p-value threshold:",
                                      max = -1.3,  min = -5,  value = -1.3, step = .1),
                          hr(),

                          selectInput("exp", "Expression filter:",
                                      list("Overexpressed" = "Over",
                                           "Underexpressed" = "Under",
                                           "None" = "None"), selected = "None"),
                          sliderInput("pexp", "Log 10 p-value threshold:",
                                      max = -1.3,  min = -5,  value = -1.3, step = .1),
                          hr(),



                          conditionalPanel(
                            condition="input.ResTabs == 'Heatmap'",
                            sliderInput("llim", "Lower limit",
                                        max = 0,  min = -100,  value = -100, step = .1),
                            sliderInput("ulim", "Upper limit",
                                        max = 100,  min = 0,  value = 100, step = .1)
                          ),

                          conditionalPanel(
                            condition="input.ResTabs == 'Correlation'",
                            selectInput("gs_x", "Geneset (x axis):",
                                        GSList),
                            selectInput("gs_y", "Geneset (y axis):",
                                        GSList),
                            sliderInput("lcor", "Lower limit",
                                        max = 0,  min = -1,  value = -1, step = .1),
                            sliderInput("ucor", "Upper limit",
                                        max = 1,  min = 0,  value = 1, step = .1)
                          )


                        ),

                        # main plots ---------------------------------------------------------

                        mainPanel(
                          tabsetPanel(id = "ResTabs",
                                      # SubTab 1 (PCA/tSNE) ---------------------------------------------------------
                                      tabPanel(title = "Modules",
                                               if(Interactive){
                                                 tabPanel("Plot",
                                                          plotlyOutput("scatPlot", height = "900px"),
                                                          plotOutput("boxPlot", height = "500px"),
                                                          plotOutput("SamplesBoxPlot", height = "400px")
                                                 )
                                               } else {
                                                 tabPanel("Plot",
                                                          plotOutput("scatPlot", height = "900px"),
                                                          plotOutput("boxPlot", height = "500px"),
                                                          plotOutput("SamplesBoxPlot", height = "400px")
                                                 )
                                               }
                                      ),

                                      # SubTab 2 (Module heatmap) ---------------------------------------------------------
                                      tabPanel("Heatmap", id = "tab2",
                                               plotOutput("hmPlot", height = "1600px")
                                               ),

                                      # SubTab 3 ---------------------------------------------------------
                                      tabPanel("Gene contribution", id = "tab3",
                                               selectInput("contribType", "Contribution mode:",
                                                           list("Positive", "Negative")),
                                               plotOutput("CorrCI", height = "1000px"),
                                               # dataTableOutput("CorrTable"),
                                               selectInput("availGenes", "", GeneList),
                                               plotOutput("ExpProj", height = "1000px")
                                      ),

                                      # SubTab 4 ---------------------------------------------------------
                                      tabPanel("Correlation", id = "tab4",
                                               if(Interactive){
                                                 tabPanel("Plot",
                                                          plotOutput("CorHmPlot", height = "1500px"),
                                                          plotlyOutput("CorScatPlot", height = "1000px")
                                                 )
                                               } else {
                                                 tabPanel("Plot",
                                                          plotOutput("CorHmPlot", height = "1500px"),
                                                          plotOutput("CorScatPlot", height = "1000px")
                                                 )
                                               }

                                               ),

                                      # SubTab 5 ---------------------------------------------------------
                                      tabPanel("ACSN (Selection)", id = "tab5",
                                               fluidPage(

                                                 fluidRow(
                                                   headerPanel("Map to use:"),
                                                   column(12,
                                                          selectInput("mapurl", "", MapList, width = "75%")
                                                   ),

                                                   headerPanel("Selected genesets:"),
                                                   column(12,
                                                          dataTableOutput("SelGSTable")
                                                   ),

                                                   headerPanel("Sample(s) to use:"),
                                                   column(12,
                                                          checkboxGroupInput("selSamples", "", inline = TRUE,
                                                                             choices = character(0))
                                                   )
                                                 ),

                                                 fluidRow(
                                                   headerPanel("Group Selection:"),
                                                   column(6,
                                                          selectInput("selgroup", NULL, as.list(unique(Groups)))
                                                   ),
                                                   column(6,
                                                          actionButton("selByGroup", "Select group"),
                                                          actionButton("selNone", "Clear selection")
                                                   )
                                                 ),

                                                 fluidRow(
                                                   headerPanel("Projection options:"),
                                                   column(6,
                                                          selectInput("projtype", "Display mode:",
                                                                      list("Module" = "Module", "Gene" = "Gene")),
                                                          textInput("ACSNWeiFil", "Filtering threshold:", "20")
                                                   ),
                                                   column(6,
                                                          selectInput("scoreaggfun", "Score aggregation:",
                                                                      list("mean" = "mean", "median" = "median",
                                                                           "sd" = "sd", "IQR" = "IQR", "mad" = "mad")),
                                                          selectInput("geneaggfun", "Gene aggregation:",
                                                                      list("mean" = "mean", "median" = "median",
                                                                           "sd" = "sd", "IQR" = "IQR", "mad" = "mad"))
                                                   )
                                                 ),

                                                 hr(),

                                                 fluidRow(
                                                   column(6,
                                                          actionButton("doACSN", "Plot on ACSN map")
                                                   ),
                                                   column(6,
                                                          htmlOutput("ACSNStatus")
                                                   )
                                                 ),

                                                 hr()
                                               )


                                      ),

                                      # SubTab 6 ---------------------------------------------------------
                                      tabPanel("ACSN (Info)", id = "tab5",
                                               fluidPage(

                                                 fluidRow(
                                                   column(6,
                                                          plotOutput("WeiVarBP")
                                                          ),
                                                   column(6,
                                                          plotOutput("ScoreDist")
                                                          )
                                                 ),

                                                 fluidRow(
                                                   column(6,
                                                          plotOutput("GeneMult")
                                                   ),
                                                   column(6,
                                                          plotOutput("ScoreVar")
                                                   )
                                                 )

                                               )


                                      )
                          )

                        )

                      )

                      ),

             # Save / Load  (Top Tab 4) ---------------------------------------------------------
             tabPanel("Save/Load",

                      fluidPage(
                        titlePanel(""),

                        fluidRow(

                          column(6,
                              wellPanel(
                                helpText("Upload a previously performed rRoma analysis"),
                                fileInput("prev.rRoma", "Choose an rRoma file", accept = c(".rds"))
                              ),
                              wellPanel(
                                helpText("Download the rRoma analysis"),
                                downloadButton("downloadData", "Download")
                              )
                            ),

                          column(6,
                            wellPanel(
                              helpText("Load rRoma analysis locally"),
                              shinyFilesButton("load", "Load data", "Please select a file", FALSE)
                            ),
                            wellPanel(
                              helpText("Save rRoma analysis locally"),
                              shinySaveButton("save", "Save data", "Save file as",
                                              filetype=list(rds="rds"))
                            )
                          )
                        )

                      )
                    )

  )




  # define server ---------------------------------------------------------

  server <- function(input, output, session) {

    options(shiny.maxRequestSize=1000*1024^2)

    Volumes = list("Working directory" = getwd(), "Volumes"= getVolumes())

    shinyFileChoose(input, "load",
                    roots=Volumes,
                    session=session, filetypes=c('rds'))

    shinyFileSave(input, "save",
                  roots=Volumes,
                  session=session)

    # Load GMT file ---------------------------------------------------------

    GetGMTFile <- reactive({

      inFile <- input$gmtfile

      if(is.null(inFile)){
        return(NULL)
      } else {
        print(paste("Loading", inFile$datapath))

        if(trimws(input$msigkw) == ""){
          KWStrings <- ""
        } else {
          KWStrings <- unlist(strsplit(trimws(input$msigkw), split = " ", fixed = TRUE))
        }

        LoadedData <- ReadGMTFile(FileLocation = inFile$datapath,
                                  SearchString = KWStrings,
                                  Mode = ifelse(input$msigkwall, "ALL", "ANY"))

        if(!input$loadwei){
          LoadedData <- lapply(LoadedData, function(GS) {
            GS$Weigths[!is.na(GS$Weigths)] <- NA
            return(GS)
          })
        }

        return(LoadedData)
      }

    })


    # Get GMT list ---------------------------------------------------------

    GetModuleList <- eventReactive(input$searchDB, {

      if(input$gmtsrc == "File"){

        ModuleList <- GetGMTFile()
        return(ModuleList)

      }

      if(input$gmtsrc == "Internal"){

        if(trimws(input$msigkw) == ""){
          KWStrings <- ""
        } else {
          KWStrings <- unlist(strsplit(trimws(input$msigkw), split = " ", fixed = TRUE))
        }

        FoundGS <- SelectFromInternalDB(SearchString = KWStrings,
                                        Mode = ifelse(input$msigkwall, "ALL", "ANY"),
                                        BDName = input$gmtlist, Version = NULL)

        if(!input$loadwei){
          FoundGS <- lapply(FoundGS, function(GS) {
            GS$Weigths[!is.na(GS$Weigths)] <- NA
            return(GS)
          })
        }

        return(FoundGS)

      }

    }, ignoreInit = FALSE, ignoreNULL = FALSE)

    # Print selected genesets ---------------------------------------------------------

    output$PrintGeneSets <- renderDataTable({

      ModuleList <- GetModuleList()

      ModuleDF <- data.frame(Names = unlist(lapply(ModuleList, "[[", "Name")),
                             Genes = unlist(lapply(lapply(ModuleList, "[[", "Genes"), length)),
                             Weighted = unlist(lapply(lapply(ModuleList, "[[", "Weigths"), function(x){sum(!is.na(x))})))

      ModuleDF

    })



    # Load expression matrix ---------------------------------------------------------

    GetExpMat <- reactive({

      inFile <- input$expmatfile

      if(is.null(inFile)){
        return(NULL)
      } else {
        print(paste("Loading", inFile$datapath))

        PlainFile.Head <- readr::read_tsv(file = inFile$datapath, col_names = TRUE, n_max = 1)

        ListSpec <- list('c', 'd')
        names(ListSpec) <- c(colnames(PlainFile.Head)[1], '.default')

        PlainFile <- readr::read_tsv(file = inFile$datapath,
                                     col_types = do.call(what = readr::cols, args = ListSpec),
                                     col_names = TRUE)

        EmptyColumns <- colSums(is.na(PlainFile)) == nrow(PlainFile)
        EmptyRows <- rowSums(is.na(PlainFile)) >= ncol(PlainFile) - 1

        if(any(EmptyColumns) | any(EmptyRows)){
          print("Filtering empty rows and columns")
          PlainFile <- PlainFile[!EmptyRows, !EmptyColumns]
        }

        ExpMat <- data.matrix(PlainFile[,-1])
        rownames(ExpMat) <- unlist(PlainFile[,1])

        return(ExpMat)
      }

    })



    # Load Group file ---------------------------------------------------------

    GetGroups <- reactive({

      if(!input$usegroups){
        return(NULL)
      }

      inFile <- input$groupfile

      if(is.null(inFile)){
        return(NULL)
      } else {
        print(paste("Loading", inFile$datapath))

        PlainFile <- readr::read_tsv(file = inFile$datapath, col_names = FALSE)

        GroupVect <- unlist(PlainFile[,2])
        names(GroupVect) <- unlist(PlainFile[,1])

        return(GroupVect)
      }

    })

    # Load previos ROMA file link ---------------------------------------------------------

    LoadPastRoma <- reactive({

      print("Loading data")

      FileName <- input$prev.rRoma

      if(!is.null(FileName)){
        print("Setting upload timestamp")
        TimeVect["Upload"] <<- Sys.time()
        return(FileName)
      }

      return(NULL)

    })


    # Do rROMA ---------------------------------------------------------

    RunROMA <- eventReactive(input$doROMA, {

      print("Running ROMA")

      # Get expression matrix
      MatData <- GetExpMat()

      if(is.null(MatData)){
        return(NULL)
      }

      # Get groups
      GroupVect <- GetGroups()

      # Get Module list
      ModuleList <- GetModuleList()

      if(is.null(ModuleList)){
        return(NULL)
      }

      Cleaned_OutGeneSpace <- NULL
      if(input$par_OutGeneSpace != "NULL"){
        Cleaned_OutGeneSpace <- as.numeric(input$par_OutGeneSpace)
      }

      Cleaned_PCSignThr <- NULL
      if(input$par_PCSignThr != "NULL"){
        Cleaned_PCSignThr <- as.numeric(input$par_PCSignThr)
      }


      RomaData <- rRoma.R(
        ExpressionMatrix = MatData, ModuleList = ModuleList,
        FixedCenter = eval(parse(text=input$par_FixedCenter)),
        PCSignMode = input$par_PCSignMode,
        nSamples = as.integer(input$par_nSamples),
        GeneOutThr = as.integer(input$par_GeneOutThr),
        UseParallel = eval(parse(text=input$par_UseParallel)),
        nCores = as.integer(input$par_nCores),
        ClusType = input$par_ClusType,
        UseWeigths = eval(parse(text=input$par_UseWeigths)),
        FullSampleInfo = eval(parse(text=input$par_FullSampleInfo)),
        centerData = eval(parse(text=input$par_centerData)),
        GeneSelMode = input$par_GeneSelMode,
        Ncomp = as.integer(input$par_Ncomp),
        DefaultWeight = as.numeric(input$par_DefaultWeight),
        SampleFilter = eval(parse(text=input$par_SampleFilter)),
        ExpFilter = eval(parse(text=input$par_ExpFilter)),
        MoreInfo = eval(parse(text=input$par_MoreInfo)),
        GeneOutDetection = input$par_GeneOutDetection,
        OutGeneNumber = as.numeric(input$par_OutGeneNumber),
        OutGeneSpace = Cleaned_OutGeneSpace,
        MinGenes = as.integer(input$par_MinGenes),
        MaxGenes = as.integer(input$par_MaxGenes),
        ApproxSamples = as.numeric(input$par_ApproxSamples),
        CorMethod = input$par_CorMethod,
        PCSignThr = Cleaned_PCSignThr,
        PlotData = FALSE, PCADims = 2, SamplingGeneWeights = NULL, FillNAMethod = list(),
        Grouping = NULL, GroupPCSign = FALSE,
        PCAType = input$par_PCAType
      )

      TimeVect["Analysis"] <<- Sys.time()
      return(RomaData)

    }, ignoreInit = FALSE, ignoreNULL = FALSE)



    # Do print ROMA results ---------------------------------------------------------

    output$ROMAOut <- renderUI({

      ExpMat <- GetExpMat()

      GroupMat <- GetGroups()

      # Get expression matrix
      if(!is.null(ExpMat)){
        str1 <- "Expression matrix loaded"
      } else {
        str1 <- "Expression matrix missing"
      }

      if(!is.null(GroupMat)){
        if(input$usegroups){
          str2 <- "Group data loaded"
        } else {
          str2 <- "Group data loaded, but ignored"
        }
      } else {
        str2 <- "Group information missing"
      }

      return(HTML(paste(str1, str2, sep = '<br/>')))

    })

    output$ROMAOut2 <- renderUI({

      ModuleList <- GetModuleList()

      # Get expression matrix
      if(!is.null(ModuleList)){
        return(HTML("Geneset list loaded"))
      }

    })



    # data table for GMT ---------------------------------------------------------

    output$PrintGMT <- renderDataTable({

      RomaData <- GetData()$RomaData

      if(is.null(RomaData)){
        return(NULL)
      }

      ModuleDF <- data.frame(Names = unlist(lapply(RomaData$ModuleSummary, "[[", "ModuleName")),
                             Genes = unlist(lapply(lapply(RomaData$ModuleSummary, "[[", "UsedGenes"), length)))

      ModuleDF[order(ModuleDF$Names),]

    })


    # parameter list ---------------------------------------------------------

    output$PrintPar <- renderDataTable({

      RomaData <- GetData()$RomaData

      if(is.null(RomaData)){
        return(NULL)
      }

      ToUse <- !(names(RomaData$InputPars) %in% c("ModuleList", "SamplingGeneWeights", "Grouping"))

      ModuleDF <- data.frame(Parameter = names(RomaData$InputPars[ToUse]),
                             Value = as.character(RomaData$InputPars[ToUse])
                             )

      ModuleDF[order(ModuleDF$Parameter),]

    })

    # Group information ---------------------------------------------------------

    output$PrintGroup <- renderDataTable({

      RomaData <- GetData()

      if(is.null(RomaData$Groups)){
        return(NULL)
      }

      ModuleDF <- data.frame(Sample = names(RomaData$Groups),
                             Group = RomaData$Groups)
      rownames(ModuleDF) <- NULL

      ModuleDF

    })



    # download data ---------------------------------------------------------

    output$downloadData <- downloadHandler(

      filename = function() {
        paste('rRoma-', Sys.Date(), '.rds', sep='')
      },

      content = function(con) {
        rRomaDashData <- list(RomaData = GetData()$RomaData,
                              ModuleList = GetModuleList(),
                              ExpMat = GetData()$ExpMat,
                              Groups = GetData()$Groups)
        saveRDS(rRomaDashData, con)
      }

    )

    # load data ---------------------------------------------------------

    GetData <- reactive({

      print("Getting Data")

      print("Trying to load uploaded ROMA data")
      LoadDataStatus <- LoadPastRoma()

      print("Trying to load ROMA input")
      LoadInputStatus <- RunROMA()

      print("Trying to load local ROMA data")
      LoadServerStatus <- LoadFromServer()

      if(all(TimeVect == StartTime)){

        return(list(
          RomaData = RomaData,
          ExpMat = ExpMat,
          FoundSamp = FoundSamp,
          Groups = Groups,
          AddInfo = AddInfo,
          GSList = GSList,
          PCAProj = PCAProj,
          ProcessedSamples = ProcessedSamples
        ))

      }

      if(max(TimeVect) == TimeVect["Upload"]){

        print(paste("Loading", LoadDataStatus$datapath))

        LoadedData <- readRDS(LoadDataStatus$datapath)

        InitReturn <- Initialize(LoadedData$RomaData, LoadedData$ExpMat, LoadedData$Groups)

        SelList <<- InitReturn$SelList
        SelListAF <<- InitReturn$SelListAF
        updateSelectInput(session, inputId = "htype", choices = SelList)
        updateSelectInput(session, inputId = "aggfun", choices = SelListAF)

        print("Passing on the information from the uploaded rds")

        return(list(
          RomaData = LoadedData$RomaData,
          ExpMat = LoadedData$ExpMat,
          FoundSamp = InitReturn$FoundSamp,
          Groups = InitReturn$Groups,
          AddInfo = InitReturn$AddInfo,
          GSList = InitReturn$GSList,
          PCAProj = InitReturn$PCAProj,
          ProcessedSamples = InitReturn$ProcessedSamples
        ))

      }

      if(max(TimeVect) == TimeVect["Analysis"]){

        print("Passing on the information from the analysis")

        InitReturn <- Initialize(LoadInputStatus, GetExpMat(), GetGroups())

        SelList <<- InitReturn$SelList
        SelListAF <<- InitReturn$SelListAF
        updateSelectInput(session, inputId = "htype", choices = SelList)
        updateSelectInput(session, inputId = "aggfun", choices = SelListAF)

        return(list(
          RomaData = LoadInputStatus,
          ExpMat = GetExpMat(),
          FoundSamp = InitReturn$FoundSamp,
          Groups = InitReturn$Groups,
          AddInfo = InitReturn$AddInfo,
          GSList = InitReturn$GSList,
          PCAProj = InitReturn$PCAProj,
          ProcessedSamples = InitReturn$ProcessedSamples
        ))

      }

      if(max(TimeVect) == TimeVect["Local"]){

        print(paste("Loading", LoadServerStatus))

        LoadedData <- readRDS(LoadServerStatus)

        InitReturn <- Initialize(LoadedData$RomaData, LoadedData$ExpMat, LoadedData$Groups)

        SelList <<- InitReturn$SelList
        SelListAF <<- InitReturn$SelListAF
        updateSelectInput(session, inputId = "htype", choices = SelList)
        updateSelectInput(session, inputId = "aggfun", choices = SelListAF)

        print("Passing on the information from the local rds")

        return(list(
          RomaData = LoadedData$RomaData,
          ExpMat = LoadedData$ExpMat,
          FoundSamp = InitReturn$FoundSamp,
          Groups = InitReturn$Groups,
          AddInfo = InitReturn$AddInfo,
          GSList = InitReturn$GSList,
          PCAProj = InitReturn$PCAProj,
          ProcessedSamples = InitReturn$ProcessedSamples
        ))

      }

    })

    # select module ---------------------------------------------------------

    SelectedGS <- reactive({
      as.integer(input$gs)
    })

    # Do tSNE ---------------------------------------------------------

    GettSNEGS <- eventReactive(input$dotSNE, {

      if(input$prjt == "tSNE"){
        ExpMat <- GetData()$ExpMat
        ProcessedSamples <- GetData()$ProcessedSamples

        initial_dims <- as.integer(input$initial_dims)

        if(initial_dims < 2){
          initial_dims = 2
        }

        updateNumericInput(session, "initial_dims", value = initial_dims)

        tSNEProj <- Rtsne::Rtsne(X = t(ExpMat[,ProcessedSamples]), perplexity = input$perp, initial_dims = initial_dims)$Y
        rownames(tSNEProj) <- ProcessedSamples
        colnames(tSNEProj) <- c("PC1", "PC2")
        print("tSNE computed")
        return(tSNEProj)
      }

    }, ignoreInit = TRUE)

    # Select modules using conditions ---------------------------------------------------------

    SelectedIdx <- reactive({

      RomaData <- GetData()$RomaData

      if(is.null(RomaData)){
        return(NULL)
      }

      SelectGeneSets(RomaData = RomaData,
                     VarThr = 10^input$pdisp, VarMode = "Wil", VarType = input$disp,
                     MedThr = 10^input$pexp, MedMode = "Wil", MedType = input$exp,
                     RatThr = 10^input$pcoord, RatMode = "Wil", RatType = input$coord)
    })

    # Update single module input ---------------------------------------------------------

    observe({

      RomaData <- GetData()$RomaData
      Groups <- GetData()$Groups
      Idx <- SelectedIdx()

      if(length(Idx)>0){
        GSList <- Idx
        names(GSList) <- as.list(rownames(RomaData$SampleMatrix)[Idx])
        GSList <- GSList[order(names(GSList))]
      } else {
        GSList <- list(" " = "")
      }

      updateSelectInput(session, "gs", choices = GSList, selected = GSList[[1]])
      updateSelectInput(session, "gs_x", choices = GSList, selected = GSList[[1]])
      updateSelectInput(session, "gs_y", choices = GSList, selected = GSList[[1]])
      updateSelectInput(session, "selgroup", choices = as.list(unique(Groups)))

    })

    # print module information ---------------------------------------------------------

    output$info <- renderUI({

      RomaData <- GetData()$RomaData

      if(is.null(RomaData)){
        return(NULL)
      }

      str1 <- paste("Number of genes = ", length(RomaData$ModuleSummary[[SelectedGS()]]$UsedGenes))
      str2 <- paste("Overdispersed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 1]))
      str3 <- paste("Underdispersed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 2]))
      str4 <- paste("Overcoordinated PV = ", signif(RomaData$PVVectMat[SelectedGS(), 3]))
      str5 <- paste("Undercoordinated PV = ", signif(RomaData$PVVectMat[SelectedGS(), 4]))
      str6 <- paste("Overexpressed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 5]))
      str7 <- paste("Underexpressed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 6]))
      HTML(paste(str1, str2, str3, str4, str5, str6, str7, sep = '<br/>'))
    })


    # sample diag boxplots ---------------------------------------------------------

    output$SamplesBoxPlot <- renderPlot({

      RomaData <- GetData()$RomaData

      if(is.null(RomaData)){
        return(NULL)
      }

      SeID <- SelectedGS()

      Sampled.DF <- t(rbind(sapply(RomaData$ModuleSummary[[SeID]]$SampledExp, "[[", "ExpVar")[1:2,],
                            sapply(RomaData$ModuleSummary[[SeID]]$SampledExp, "[[", "MedianExp")))

      Sampled.DF <- cbind(Sampled.DF, Sampled.DF[,1]/Sampled.DF[,2])

      colnames(Sampled.DF) <- c("L1", "L2", "Exp", "L1/L2")

      Sampled.DF <- reshape::melt(Sampled.DF)
      Sampled.DF <- cbind(Sampled.DF, rep("Samples", nrow(Sampled.DF)))

      colnames(Sampled.DF)[4] <- "Org"

      Sampled.DF$Org <- factor(as.character(Sampled.DF$Org), levels = c("Samples", "Data"))

      Sampled.DF <- rbind(Sampled.DF,
                          c(max(Sampled.DF$X1)+1,
                            "L1",
                            RomaData$ModuleSummary[[SeID]]$ExpVarBase[1],
                            "Data"),
                          c(max(Sampled.DF$X1)+2,
                            "L2",
                            RomaData$ModuleSummary[[SeID]]$ExpVarBase[2],
                            "Data"),
                          c(max(Sampled.DF$X1)+3,
                            "L1/L2",
                            RomaData$ModuleSummary[[SeID]]$ExpVarBase[1]/RomaData$ModuleSummary[[SeID]]$ExpVarBase[2],
                            "Data"),
                          c(max(Sampled.DF$X1)+4,
                            "Exp",
                            RomaData$ModuleMatrix[SeID,"Median Exp"],
                            "Data")
                          )

      Sampled.DF <- data.frame(Sampled.DF)

      Sampled.DF$X2 <- factor(as.character(Sampled.DF$X2), levels = c("L1", "L2", "L1/L2", "Exp"))
      Sampled.DF$value <- as.numeric(as.character(Sampled.DF$value))

      p1 <- ggplot2::ggplot(data = Sampled.DF[Sampled.DF$Org == "Samples" &
                                               Sampled.DF$X2 %in% c("L1", "L2"), ],
                           ggplot2::aes(y = value, x = X2, color = Org)) +
        ggplot2::geom_boxplot() + ggplot2::facet_wrap(~X2, scales = "free_x", ncol = 4) +
        ggplot2::scale_color_manual(values = c(Data="red", Samples="black")) +
        ggplot2::geom_point(data = Sampled.DF[Sampled.DF$Org == "Data" &
                                                Sampled.DF$X2 %in% c("L1", "L2"), ],
                   ggplot2::aes(y = value, x = X2, color = Org), inherit.aes = FALSE, size = 3) +
        ggplot2::labs(x="", y="Explained variance") +
        ggplot2::guides(color=FALSE) +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank())


        p2 <- ggplot2::ggplot(data = Sampled.DF[Sampled.DF$Org == "Samples" &
                                                 Sampled.DF$X2 %in% c("L1/L2"), ],
                             ggplot2::aes(y = value, x = X2, color = Org)) +
          ggplot2::geom_boxplot() + ggplot2::facet_wrap(~X2, scales = "free", ncol = 4) +
          ggplot2::scale_color_manual(values = c(Data="red", Samples="black")) +
          ggplot2::geom_point(data = Sampled.DF[Sampled.DF$Org == "Data" &
                                                  Sampled.DF$X2 %in% c("L1/L2"), ],
                              ggplot2::aes(y = value, x = X2, color = Org), inherit.aes = FALSE, size = 3) +
          ggplot2::labs(x="", y="Ratio") +
          ggplot2::guides(color=FALSE) +
          ggplot2::scale_y_log10() +
          ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                         axis.text.x=ggplot2::element_blank(),
                         axis.ticks.x=ggplot2::element_blank())


        p3 <- ggplot2::ggplot(data = Sampled.DF[Sampled.DF$Org == "Samples" &
                                                 Sampled.DF$X2 %in% c("Exp"), ],
                             ggplot2::aes(y = value, x = X2, color = "Samples")) +
          ggplot2::geom_boxplot() + ggplot2::facet_wrap(~X2, scales = "free", ncol = 4) +
          ggplot2::scale_color_manual("", values = c(Data="red", Samples="black")) +
          ggplot2::geom_point(data = Sampled.DF[Sampled.DF$Org == "Data" &
                                                  Sampled.DF$X2 %in% c("Exp"), ],
                              ggplot2::aes(y = value, x = X2, color = Org),
                              inherit.aes = FALSE, size = 3) +
          ggplot2::labs(x="", y="Median Expression") +
          ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                         axis.text.x=ggplot2::element_blank(),
                         axis.ticks.x=ggplot2::element_blank())

        gridExtra::grid.arrange(p1, p2, p3, ncol = 3, widths=c(2,1,1.5))

    })


    # single sample boxplot ---------------------------------------------------------

    output$boxPlot <- renderPlot({

      if(is.na(SelectedGS())){
        return(NULL)
      }

      RomaData <- GetData()$RomaData
      Groups <- GetData()$Groups
      ProcessedSamples <- GetData()$ProcessedSamples

      DataToPlot <- data.frame(Score = RomaData$SampleMatrix[SelectedGS(), ProcessedSamples],
                               Group = Groups[ProcessedSamples])

      p <- ggplot2::ggplot(data = DataToPlot, ggplot2::aes(x = Group, y = Score)) +
        ggplot2::geom_boxplot()

      if(length(unique(Groups[ProcessedSamples]))>1 & input$boxcomp != "Show none"){

        if(input$boxcomp == "Show all"){
          ToDisplay <- GetComb(unique(Groups[ProcessedSamples]))
        }

        if(input$boxcomp == "Show significant"){
          PWComp <- pairwise.wilcox.test(DataToPlot$Score, DataToPlot$Group)$p.value
          ExtPWComp <- matrix(rep(NA, length(unique(DataToPlot$Group))^2), nrow = length(unique(DataToPlot$Group)))
          colnames(ExtPWComp) <- sort(unique(DataToPlot$Group))
          rownames(ExtPWComp) <- colnames(ExtPWComp)

          ExtPWComp[rownames(PWComp), colnames(PWComp)] <- PWComp

          ToDisplay <- apply(which(ExtPWComp <= .05, arr.ind = TRUE), 1, list)
          ToDisplay <- lapply(ToDisplay, function(x){as.integer(unlist(x))})
          # names(ToDisplay) <- NULL
        }

        p <- p + ggsignif::geom_signif(comparisons = ToDisplay, map_signif_level=TRUE, test = "wilcox.test", step_increase = .1)

      }

      print(p)

    })

    # PCA / tSNE plot ---------------------------------------------------------

    if(Interactive){
      output$scatPlot <- renderPlotly({

        RomaData <- GetData()$RomaData
        Groups <- GetData()$Groups
        ProcessedSamples <- GetData()$ProcessedSamples

        if(is.null(RomaData)){
          return(NULL)
        }

        input$dotSNE

        if(input$prjt == "PCA"){
          print("using PCA")
          Projs <- GetData()$PCAProj
        }

        if(input$prjt == "tSNE"){
          print("using tSNE")
          Projs <- GettSNEGS()
        }

        p <- ggplot2::ggplot(data = data.frame(Comp1 = Projs[,1],
                                                       Comp2 = Projs[,2],
                                                       Score = RomaData$SampleMatrix[SelectedGS(), ProcessedSamples],
                                                       Group = Groups[ProcessedSamples]),
                                     ggplot2::aes(x = Comp1, y = Comp2, shape = Group, color = Score)) +
          ggplot2::scale_color_gradient2(low = "blue", high = "red", mid = "white") +
          ggplot2::labs(x = "Component 1", y = "Component 2", shape = "",
                        title = rownames(RomaData$SampleMatrix)[SelectedGS()]) +
          ggplot2::geom_point(ggplot2::aes(text = rownames(Projs)))

        print(ggplotly(p))
      })
    } else {
      output$scatPlot <- renderPlot({

        RomaData <- GetData()$RomaData
        Groups <- GetData()$Groups
        ProcessedSamples <- GetData()$ProcessedSamples

        if(is.null(RomaData)){
          return(NULL)
        }

        input$dotSNE

        if(input$prjt == "PCA"){
          print("using PCA")
          Projs <- GetData()$PCAProj
        }

        if(input$prjt == "tSNE"){
          print("using tSNE")
          Projs <- GettSNEGS()
        }

        p <- ggplot2::ggplot(data = data.frame(Comp1 = Projs[,1],
                                               Comp2 = Projs[,2],
                                               Score = RomaData$SampleMatrix[SelectedGS(), ProcessedSamples],
                                               Group = Groups[ProcessedSamples]),
                             ggplot2::aes(x = Comp1, y = Comp2, shape = Group, color = Score)) +
          ggplot2::scale_color_gradient2(low = "blue", high = "red", mid = "white") +
          ggplot2::labs(x = "Component 1", y = "Component 2", shape = "",
                        title = rownames(RomaData$SampleMatrix)[SelectedGS()]) +
          ggplot2::geom_point(size = 3)
        print(p)
      })
    }

    # heatmap (modules) plot ---------------------------------------------------------

    output$hmPlot <- renderPlot({

      RomaData <- GetData()$RomaData
      Groups <- GetData()$Groups
      ProcessedSamples <- GetData()$ProcessedSamples
      AddInfo <- GetData()$AddInfo
      FoundSamp <- GetData()$FoundSamp

      if(is.null(RomaData)){
        return(NULL)
      }

      BaseCol <- colorRamps::blue2red(54)

      Idx <- SelectedIdx()

      PlotMat <- RomaData$SampleMatrix[Idx,]

      GSCat <- rep(NA, nrow(RomaData$SampleMatrix))
      names(GSCat) <- rownames(RomaData$SampleMatrix)

      if(input$GSOrdeMode == "None"){
        GSOrdering <- order(rownames(RomaData$SampleMatrix[Idx,]))
        GSCat <- NA
      }

      if(input$GSOrdeMode == "Gene number"){
        nGenes <- unlist(lapply(lapply(RomaData$ModuleSummary[Idx], "[[", "UsedGenes"), length))
        GSOrdering <- order(nGenes)
        GSCat[] = nGenes
        GSCat <- data.frame(Genes = signif(GSCat, 2))
      }

      if(input$GSOrdeMode == "Overdispersion pv"){
        GSOrdering <- order(RomaData$PVVectMat[Idx, "L1 WT less pv"])
        GSCat[] = log10(RomaData$PVVectMat[Idx, "L1 WT less pv"])
        GSCat <- data.frame("OD lpv" = signif(GSCat, 2))
      }

      if(input$GSOrdeMode == "Underdispersion pv"){
        GSOrdering <- order(RomaData$PVVectMat[Idx, "L1 WT greater pv"])
        GSCat[] = log10(RomaData$PVVectMat[Idx, "L1 WT greater pv"])
        GSCat <- data.frame("UD lpv" = signif(GSCat, 2))
      }

      if(input$GSOrdeMode == "Overcoordination pv"){
        GSOrdering <- order(RomaData$PVVectMat[Idx, "L1/L2 WT less pv"])
        GSCat[] = log10(RomaData$PVVectMat[Idx, "L1/L2 WT less pv"])
        GSCat <- data.frame("OC lpv" = signif(GSCat, 2))
      }

      if(input$GSOrdeMode == "Undercoordination pv"){
        GSOrdering <- order(RomaData$PVVectMat[Idx, "L1/2 WT greater pv"])
        GSCat[] = log10(RomaData$PVVectMat[Idx, "L1/2 WT greater pv"])
        GSCat <- data.frame("UC lpv" = signif(GSCat, 2))
      }

      if(input$GSOrdeMode == "Overexpression pv"){
        GSOrdering <- order(RomaData$PVVectMat[Idx, "Median Exp WT less pv"])
        GSCat[] = log10(RomaData$PVVectMat[Idx, "Median Exp WT less pv"])
        GSCat <- data.frame("OE lpv" = signif(GSCat, 2))
      }

      if(input$GSOrdeMode == "Underexpression pv"){
        GSOrdering <- order(RomaData$PVVectMat[Idx, "Median Exp WT greater pv"])
        GSCat[] = log10(RomaData$PVVectMat[Idx, "Median Exp WT greater pv"])
        GSCat <- data.frame("UE lpv" = signif(GSCat, 2))
      }


      if(input$htype == "sample"){

        MinMax <- range(PlotMat)

        # print(MinMax)

        if(input$llim < MinMax[1]){
          updateSliderInput(session, "llim", value = ifelse(MinMax[1] < 0, floor(10*MinMax[1])/10, 0))
        }

        updateSliderInput(session, "llim", min = ifelse(MinMax[1] < 0, floor(10*MinMax[1])/10, 0))



        if(input$ulim > MinMax[2]){
          updateSliderInput(session, "ulim", value = ifelse(MinMax[2] > 0, ceiling(10*MinMax[2])/10, 0))
        }

        updateSliderInput(session, "ulim", max = ifelse(MinMax[2] > 0, ceiling(10*MinMax[2])/10, 0))


        SatLL <- input$llim

        if(MinMax[1] < 0){
          DoLow <- TRUE
          LowBrk <- c(MinMax[1]+SatLL/28, seq(from = SatLL*(26/27), to = 0, by = -SatLL/27))
        } else {
          DoLow <- FALSE
        }

        SatUL <- input$ulim

        if(MinMax[2] > 0){
          DoHigh <- TRUE
          HighBrk <- c(seq(from = SatUL/27, to = (26/27)*SatUL, by = SatUL/27), MinMax[2] + SatUL/28)
        } else {
          DoHigh <- FALSE
        }

        if(DoLow){
          MyBreaks <- LowBrk
          if(SatLL < 0){
            UseCol <- c(1:27)
          } else {
            UseCol <- c(1)
          }
        } else {
          MyBreaks <- 0
          UseCol <- NULL
        }

        if(DoHigh){
          if(SatUL > 0){
            MyBreaks <- c(MyBreaks, HighBrk)
            UseCol <- c(UseCol, 28:54)
          } else {
            MyBreaks <- c(MyBreaks, HighBrk[2])
            UseCol <- c(UseCol, 54)
          }

        }

        if(length(Idx)>1){

          if(input$gscol){
            pheatmap::pheatmap(PlotMat[GSOrdering,], color = BaseCol[UseCol], breaks = MyBreaks,
                               cluster_rows = input$gsclus, cluster_cols = input$saclus,
                               annotation_col = AddInfo,
                               annotation_row = GSCat,
                               main = "Module scores across samples")
          } else {
            pheatmap::pheatmap(t(PlotMat[GSOrdering,]), color = BaseCol[UseCol], breaks = MyBreaks,
                               cluster_rows = input$saclus, cluster_cols = input$gsclus,
                               annotation_row = AddInfo,
                               annotation_col = GSCat,
                               main = "Module scores across samples")
          }


        }

        if(length(Idx) == 1){

          PlotMat <- matrix(PlotMat, nrow = 1)

          colnames(PlotMat) <- colnames(RomaData$SampleMatrix)

          if(input$gscol){
            pheatmap::pheatmap(PlotMat, color = BaseCol[UseCol], breaks = MyBreaks,
                               cluster_cols = input$saclus, cluster_rows = FALSE,
                               annotation_col = AddInfo,
                               main = paste("Score of", rownames(RomaData$SampleMatrix)[Idx]))
          } else {
            pheatmap::pheatmap(t(PlotMat), color = BaseCol[UseCol], breaks = MyBreaks,
                               cluster_rows = input$saclus, cluster_cols = FALSE,
                               annotation_row = AddInfo,
                               main = paste("Score of", rownames(RomaData$SampleMatrix)[Idx]))
          }

          # pheatmap::pheatmap(PlotMat, BaseCol[UseCol], breaks = MyBreaks,
          #                    cluster_rows = FALSE, cluster_cols = FALSE,
          #                    main = paste("Score of", rownames(RomaData$SampleMatrix)[Idx]))

        }
      }

      if(input$htype == "group"){

        if(length(Idx) > 1){

          SplitData <- split(data.frame(t(PlotMat[,FoundSamp])), f=AddInfo$Groups)

          Aggmat <- sapply(SplitData, function(x) {
            apply(x, 2, get(input$aggfun))
          })

          MinMax <- range(Aggmat)

          if(input$llim < MinMax[1]){
            updateSliderInput(session, "llim", value = ifelse(MinMax[1] < 0, floor(10*MinMax[1])/10, 0))
          }

          updateSliderInput(session, "llim", min = ifelse(MinMax[1] < 0, floor(10*MinMax[1])/10, 0))



          if(input$ulim > MinMax[2]){
            updateSliderInput(session, "ulim", value = ifelse(MinMax[2] > 0, ceiling(10*MinMax[2])/10, 0))
          }

          updateSliderInput(session, "ulim", max = ifelse(MinMax[2] > 0, ceiling(10*MinMax[2])/10, 0))



          SatLL <- input$llim

          if(MinMax[1] < 0){
            DoLow <- TRUE
            LowBrk <- c(MinMax[1]+SatLL/28, seq(from = SatLL*(26/27), to = 0, by = -SatLL/27))
          } else {
            DoLow <- FALSE
          }

          SatUL <- input$ulim

          if(MinMax[2] > 0){
            DoHigh <- TRUE
            HighBrk <- c(seq(from = SatUL/27, to = (26/27)*SatUL, by = SatUL/27), MinMax[2] + SatUL/28)
          } else {
            DoHigh <- FALSE
          }

          if(DoLow){
            MyBreaks <- LowBrk
            if(SatLL < 0){
              UseCol <- c(1:27)
            } else {
              UseCol <- c(1)
            }
          } else {
            MyBreaks <- 0
            UseCol <- NULL
          }

          if(DoHigh){
            if(input$ulim > 0){
              MyBreaks <- c(MyBreaks, HighBrk)
              UseCol <- c(UseCol, 28:54)
            } else {
              MyBreaks <- c(MyBreaks, HighBrk[2])
              UseCol <- c(UseCol, 54)
            }

          }


          if(input$gscol){
            pheatmap::pheatmap(Aggmat[GSOrdering,], color = BaseCol[UseCol], breaks = MyBreaks,
                               cluster_rows = input$gsclus, cluster_cols = input$saclus,
                               annotation_row = GSCat, main = "Module scores across groups")
          } else {
            pheatmap::pheatmap(t(Aggmat[GSOrdering,]), color = BaseCol[UseCol], breaks = MyBreaks,
                               cluster_rows = input$saclus, cluster_cols = input$gsclus,
                               annotation_col = GSCat, main = "Module scores across groups")
          }

        }

        if(length(Idx) == 1){

          # names(PlotMat) <- colnames(RomaData$SampleMatrix)
          SplitData <- split(data.frame(PlotMat[FoundSamp]), f=AddInfo$Groups)

          Aggmat <- sapply(SplitData, function(x) {
            do.call(input$aggfun, list(unlist(x)))
          })

          MinMax <- range(Aggmat)

          if(input$llim < MinMax[1]){
            updateSliderInput(session, "llim", value = ifelse(MinMax[1] < 0, floor(10*MinMax[1])/10, 0))
          }

          updateSliderInput(session, "llim", min = ifelse(MinMax[1] < 0, floor(10*MinMax[1])/10, 0))



          if(input$ulim > MinMax[2]){
            updateSliderInput(session, "ulim", value = ifelse(MinMax[2] > 0, ceiling(10*MinMax[2])/10, 0))
          }

          updateSliderInput(session, "ulim", max = ifelse(MinMax[2] > 0, ceiling(10*MinMax[2])/10, 0))


          if(MinMax[1] < 0){
            DoLow <- TRUE
            LowBrk <- c(MinMax[1], seq(from = input$llim*(26/27), to = 0, by = -input$llim/27))
            if(LowBrk[1] != min(LowBrk)){
              LowBrk[1] <- min(LowBrk) + input$llim/27
            }
          } else {
            DoLow <- FALSE
          }

          if(MinMax[2] > 0){
            DoHigh <- TRUE
            HighBrk <- c(seq(from = input$ulim/27, to = (26/27)*input$ulim, by = input$ulim/27), MinMax[2])
            if(HighBrk[length(HighBrk)] != max(HighBrk)){
              HighBrk[length(HighBrk)] <- max(HighBrk) + input$ulim/27
            }
            if(HighBrk[1] == 0){
              HighBrk <- HighBrk[-1]
            }
          } else {
            DoHigh <- FALSE
          }

          if(DoLow){
            MyBreaks <- LowBrk
            if(input$llim < 0){
              UseCol <- c(1:27)
            } else {
              UseCol <- c(1,27)
            }
          } else {
            MyBreaks <- 0
            UseCol <- 27
          }

          if(DoHigh){
            MyBreaks <- c(MyBreaks, HighBrk)
            if(input$ulim > 0){
              UseCol <- c(UseCol, 28:54)
            } else {
              UseCol <- c(UseCol, 28, 54)
            }

          }

          pheatmap::pheatmap(t(Aggmat), color = BaseCol[UseCol], breaks = MyBreaks,
                             cluster_rows = FALSE, cluster_cols = FALSE,
                             main = paste("Score of", rownames(RomaData$SampleMatrix)[Idx]))
        }

      }

    })






    # corr heatmap ---------------------------------------------------------

    output$CorHmPlot <- renderPlot({

      RomaData <- GetData()$RomaData
      Groups <- GetData()$Groups
      ProcessedSamples <- GetData()$ProcessedSamples
      AddInfo <- GetData()$AddInfo
      FoundSamp <- GetData()$FoundSamp

      if(is.null(RomaData)){
        return(NULL)
      }

      BaseCol <- colorRamps::blue2red(54)

      Idx <- SelectedIdx()

      GSCat <- rep(NA, nrow(RomaData$SampleMatrix))
      names(GSCat) <- rownames(RomaData$SampleMatrix)

      if(input$GSOrdeMode == "None"){
        GSOrdering <- order(rownames(RomaData$SampleMatrix[Idx,]))
        GSCat <- NA
      }

      if(input$GSOrdeMode == "Gene number"){
        nGenes <- unlist(lapply(lapply(RomaData$ModuleSummary[Idx], "[[", "UsedGenes"), length))
        GSOrdering <- order(nGenes)
        GSCat[] = nGenes
        GSCat <- data.frame(Genes = signif(GSCat, 2))
      }

      if(input$GSOrdeMode == "Overdispersion pv"){
        GSOrdering <- order(RomaData$PVVectMat[Idx, "L1 WT less pv"])
        GSCat[] = log10(RomaData$PVVectMat[Idx, "L1 WT less pv"])
        GSCat <- data.frame("OD lpv" = signif(GSCat, 2))
      }

      if(input$GSOrdeMode == "Underdispersion pv"){
        GSOrdering <- order(RomaData$PVVectMat[Idx, "L1 WT greater pv"])
        GSCat[] = log10(RomaData$PVVectMat[Idx, "L1 WT greater pv"])
        GSCat <- data.frame("UD lpv" = signif(GSCat, 2))
      }

      if(input$GSOrdeMode == "Overcoordination pv"){
        GSOrdering <- order(RomaData$PVVectMat[Idx, "L1/L2 WT less pv"])
        GSCat[] = log10(RomaData$PVVectMat[Idx, "L1/L2 WT less pv"])
        GSCat <- data.frame("OC lpv" = signif(GSCat, 2))
      }

      if(input$GSOrdeMode == "Undercoordination pv"){
        GSOrdering <- order(RomaData$PVVectMat[Idx, "L1/2 WT greater pv"])
        GSCat[] = log10(RomaData$PVVectMat[Idx, "L1/2 WT greater pv"])
        GSCat <- data.frame("UC lpv" = signif(GSCat, 2))
      }

      if(input$GSOrdeMode == "Overexpression pv"){
        GSOrdering <- order(RomaData$PVVectMat[Idx, "Median Exp WT less pv"])
        GSCat[] = log10(RomaData$PVVectMat[Idx, "Median Exp WT less pv"])
        GSCat <- data.frame("OE lpv" = signif(GSCat, 2))
      }

      if(input$GSOrdeMode == "Underexpression pv"){
        GSOrdering <- order(RomaData$PVVectMat[Idx, "Median Exp WT greater pv"])
        GSCat[] = log10(RomaData$PVVectMat[Idx, "Median Exp WT greater pv"])
        GSCat <- data.frame("UE lpv" = signif(GSCat, 2))
      }

      LowBrk <- c(-1.01, seq(from = input$lcor*(26/27), to = 0, by = -input$lcor/27))
      HighBrk <- c(seq(from = input$ucor/27, to = (26/27)*input$ucor, by = input$ucor/27), 1.01)

      MyBreaks <- LowBrk
      if(input$lcor < 0){
        UseCol <- c(1:27)
      } else {
        UseCol <- c(1)
      }

      if(input$ucor > 0){
        MyBreaks <- c(MyBreaks, HighBrk)
        UseCol <- c(UseCol, 28:54)
      } else {
        MyBreaks <- c(MyBreaks, 1.01)
        UseCol <- c(UseCol, 54)
      }

      if(input$htype == "sample"){

        if(length(Idx)>1){
          pheatmap::pheatmap(cor(t(RomaData$SampleMatrix[Idx,]), method = input$cortype)[GSOrdering, GSOrdering],
                             color = BaseCol[UseCol], breaks = MyBreaks,
                             annotation_row = GSCat, annotation_col = GSCat,
                             cluster_rows = input$gsclus, cluster_cols = input$gsclus,
                             main = "Correlation across samples")
        }

        if(length(Idx) == 1){
          NULL
        }
      }

      if(input$htype == "group"){

        if(length(Idx) > 1){

          SplitData <- split(data.frame(t(RomaData$SampleMatrix[Idx,FoundSamp])), f=AddInfo$Groups)

          Aggmat <- sapply(SplitData, function(x) {
            apply(x, 2, get(input$aggfun))
          })

          pheatmap::pheatmap(cor(t(Aggmat), method = input$cortype)[GSOrdering, GSOrdering],
                             color = BaseCol[UseCol], breaks = MyBreaks,
                             cluster_rows = input$gsclus, cluster_cols = input$gsclus,
                             annotation_row = GSCat, annotation_col = GSCat,
                             main = "Correlation across groups")

        }

        if(length(Idx) == 1){
          NULL
        }

      }

    })



    # corr scatplot ---------------------------------------------------------

    if(Interactive){
      output$CorScatPlot <- renderPlotly({

        RomaData <- GetData()$RomaData
        Groups <- GetData()$Groups
        ProcessedSamples <- GetData()$ProcessedSamples
        AddInfo <- GetData()$AddInfo
        FoundSamp <- GetData()$FoundSamp

        if(is.null(RomaData)){
          return(NULL)
        }

        XLab <- paste("Scores - ", RomaData$ModuleSummary[[as.integer(input$gs_x)]]$ModuleName, " - ",
                      length(RomaData$ModuleSummary[[as.integer(input$gs_x)]]$UsedGenes), " genes (",
                      length(intersect(RomaData$ModuleSummary[[as.integer(input$gs_x)]]$UsedGenes,
                                       RomaData$ModuleSummary[[as.integer(input$gs_y)]]$UsedGenes)), " shared)",
                      sep = "")

        YLab <- paste("Scores - ", RomaData$ModuleSummary[[as.integer(input$gs_y)]]$ModuleName, " - ",
                      length(RomaData$ModuleSummary[[as.integer(input$gs_y)]]$UsedGenes), " genes (",
                      length(intersect(RomaData$ModuleSummary[[as.integer(input$gs_x)]]$UsedGenes,
                                       RomaData$ModuleSummary[[as.integer(input$gs_y)]]$UsedGenes)), " shared)",
                      sep = "")

        CTitle <- cor(RomaData$SampleMatrix[as.integer(input$gs_x), ProcessedSamples],
                     RomaData$SampleMatrix[as.integer(input$gs_y), ProcessedSamples],
                     method = input$cortype)


        if(input$htype == "sample"){

          p <- ggplot2::ggplot(data = data.frame(XVal = RomaData$SampleMatrix[as.integer(input$gs_x), ProcessedSamples],
                                                 YVal = RomaData$SampleMatrix[as.integer(input$gs_y), ProcessedSamples],
                                                 Group = Groups[ProcessedSamples]),
                               ggplot2::aes(x = XVal, y = YVal, color = Group)) +
            ggplot2::labs(x = XLab, y = YLab, shape = "", title = paste("Corr =", signif(CTitle, 5))) +
            ggplot2::geom_point(ggplot2::aes(text = ProcessedSamples))

        }


        if(input$htype == "group"){

          AggX <- aggregate(RomaData$SampleMatrix[as.integer(input$gs_x), ProcessedSamples], list(AddInfo$Groups), get(input$aggfun))
          AggY <- aggregate(RomaData$SampleMatrix[as.integer(input$gs_y), ProcessedSamples], list(AddInfo$Groups), get(input$aggfun))

          p <- ggplot2::ggplot(data = data.frame(XVal = AggX[,2],
                                                 YVal = AggY[,2],
                                                 Group = AggX[,1]),
                               ggplot2::aes(x = XVal, y = YVal, color = Group)) +
            ggplot2::labs(x = XLab, y = YLab, shape = "", title = paste("Corr =", signif(CTitle, 5))) +
            ggplot2::geom_point(ggplot2::aes(text = AggX[,1]))

        }

        print(ggplotly(p))

      })
    } else {

      output$CorScatPlot <- renderPlot({

        RomaData <- GetData()$RomaData
        Groups <- GetData()$Groups
        ProcessedSamples <- GetData()$ProcessedSamples
        AddInfo <- GetData()$AddInfo
        FoundSamp <- GetData()$FoundSamp

        if(is.null(RomaData)){
          return(NULL)
        }

        XLab <- paste("Scores - ", RomaData$ModuleSummary[[as.integer(input$gs_x)]]$ModuleName, " - ",
                      length(RomaData$ModuleSummary[[as.integer(input$gs_x)]]$UsedGenes), " genes (",
                      length(intersect(RomaData$ModuleSummary[[as.integer(input$gs_x)]]$UsedGenes,
                                       RomaData$ModuleSummary[[as.integer(input$gs_y)]]$UsedGenes)), " shared)",
                      sep = "")

        YLab <- paste("Scores - ", RomaData$ModuleSummary[[as.integer(input$gs_y)]]$ModuleName, " - ",
                      length(RomaData$ModuleSummary[[as.integer(input$gs_y)]]$UsedGenes), " genes (",
                      length(intersect(RomaData$ModuleSummary[[as.integer(input$gs_x)]]$UsedGenes,
                                       RomaData$ModuleSummary[[as.integer(input$gs_y)]]$UsedGenes)), " shared)",
                      sep = "")

        CTitle <- cor(RomaData$SampleMatrix[as.integer(input$gs_x), ProcessedSamples],
                      RomaData$SampleMatrix[as.integer(input$gs_y), ProcessedSamples],
                      method = input$cortype)

        if(input$htype == "sample"){

          p <- ggplot2::ggplot(data = data.frame(XVal = RomaData$SampleMatrix[as.integer(input$gs_x), ProcessedSamples],
                                                 YVal = RomaData$SampleMatrix[as.integer(input$gs_y), ProcessedSamples],
                                                 Group = Groups[ProcessedSamples]),
                               ggplot2::aes(x = XVal, y = YVal, color = Group)) +
            ggplot2::labs(x = XLab, y = YLab, shape = "", title = paste("Corr =", signif(CTitle, 5))) +
            ggplot2::geom_point()

        }


        if(input$htype == "group"){

          AggX <- aggregate(RomaData$SampleMatrix[as.integer(input$gs_x), ProcessedSamples], list(AddInfo$Groups), get(input$aggfun))
          AggY <- aggregate(RomaData$SampleMatrix[as.integer(input$gs_y), ProcessedSamples], list(AddInfo$Groups), get(input$aggfun))


           p <- ggplot2::ggplot(data = data.frame(XVal = AggX[,2],
                                                 YVal = AggY[,2],
                                                 Group = AggX[,1]),
                               ggplot2::aes(x = XVal, y = YVal, color = Group)) +
            ggplot2::labs(x = XLab, y = YLab, shape = "", title = paste("Corr =", signif(CTitle, 5))) +
            ggplot2::geom_point()

        }

        print(p)

      })
    }

    # Selected GS datatable ---------------------------------------------------------

    output$SelGSTable <- renderDataTable({

      RomaData <- GetData()$RomaData

      if(is.null(RomaData)){
        return(NULL)
      } else {
        SelIdxs <- SelectedIdx()

        RetDF <- data.frame(Names = unlist(lapply(RomaData$ModuleSummary[SelIdxs], "[[", "ModuleName"), use.names = FALSE))

        return(RetDF)
      }

    })

    # Update available samples ---------------------------------------------------------

    observe({

      RomaData <- GetData()$RomaData

      if(!is.null(RomaData)){
        Samples <- colnames(RomaData$SampleMatrix)
        updateCheckboxGroupInput(session, "selSamples", choices = Samples, inline = TRUE)
      }

    })

    # Plot on ACSN ---------------------------------------------------------

    GetACSN <- eventReactive(input$doACSN, {

      print("Initiating ACSN plotting")

      SelIdxs <- SelectedIdx()
      RomaData <- GetData()$RomaData

      if(!is.null(input$selSamples) & !is.null(RomaData) & !is.null(SelIdxs)){
        RetData <- PlotOnACSN(
          RomaData = RomaData, SampleName = input$selSamples, AggScoreFun = input$scoreaggfun,
          MapURL = input$mapurl, Selected = SelIdxs, FilterByWei = as.numeric(input$ACSNWeiFil),
          AggGeneFun = input$geneaggfun, DispMode = input$projtype, DefDispVal = 0,
          PlotInfo = FALSE, ReturnInfo = TRUE, LocalRange = FALSE)
      } else {
        RetData  <- NULL
      }

      return(RetData)

    }, ignoreInit = TRUE, ignoreNULL = FALSE)


    # Include Status (ACSN) ---------------------------------------------------------

    output$ACSNStatus <- renderUI({
      Test <- GetACSN()

      if(is.null(Test)){
        HTML("<i>Data not exported on ACSN</i>")
      } else {
        HTML("<b>Data exported to ACSN</b>")
      }
    })

    # Plot weigth variation (ACSN) ---------------------------------------------------------

    output$WeiVarBP <- renderPlot({

      Data <- GetACSN()

      if(!is.null(Data)){

        p <- ggplot2::ggplot(data = data.frame(y = Data$GenesVar,
                                               x = rep(1, length(Data$GenesVar))),
                             mapping = ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_boxplot(mapping = ggplot2::aes(color = "Data")) +
          ggplot2::labs(y = "Variance of gene weight") + ggplot2::scale_y_log10() +
          ggplot2::scale_color_manual("", values = c(Data="black", Outliers="red")) +
          ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                         axis.text.x=ggplot2::element_blank(),
                         axis.ticks.x=ggplot2::element_blank())

        if(sum(Data$GeneOut) > 1){

          p <- p + ggplot2::geom_point(data = data.frame(x = rep(1, length(which(Data$GeneOut))),
                                                y = Data$GenesVar[names(which(Data$GeneOut))]),
                              mapping = ggplot2::aes(x = x, y = y, color="Outliers"), size = 3,
                              inherit.aes = FALSE)

        }

        print(p)
      }

    })

    # Plot gene multiplicity (ACSN) ---------------------------------------------------------

    output$GeneMult <- renderPlot({

      Data <- GetACSN()

      if(!is.null(Data)){
        if(!is.null(Data$GeneMult)){
          p <- ggplot2::ggplot(data = data.frame(Data$GeneMult),
                               mapping = ggplot2::aes(x = Freq)) +
            ggplot2::geom_histogram(binwidth = 1) +
            ggplot2::labs(y = "Frequency", x = "Gene multiplicity")

          print(p)
        }
      }

    })

    # Plot score distribution (ACSN) ---------------------------------------------------------

    output$ScoreDist <- renderPlot({

      Data <- GetACSN()

      if(!is.null(Data)){

        if((input$projtype == "Module") & !is.null(Data$GeneMult)){
          p <- ggplot2::ggplot(data = data.frame(Data$ScoreDist,
                                                 x = rep(1, length(Data$ScoreDist))),
                               mapping = ggplot2::aes(x = x, y = data)) +
            ggplot2::geom_boxplot() +
            ggplot2::labs(y = "Gene score (Module)", x = "") +
            ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                           axis.text.x=ggplot2::element_blank(),
                           axis.ticks.x=ggplot2::element_blank())

          print(p)
        }

        if((input$projtype == "Gene") & !is.null(Data$WeiDist)){
          p <- ggplot2::ggplot(data = data.frame(Data$WeiDist,
                                                 x = rep(1, length(Data$WeiDist))),
                               mapping = ggplot2::aes(x = x, y = data)) +
            ggplot2::geom_boxplot() +
            ggplot2::labs(y = "Gene score (Weight)", x = "") +
            ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                           axis.text.x=ggplot2::element_blank(),
                           axis.ticks.x=ggplot2::element_blank())

          print(p)
        }
      }

    })


    # Plot score variance (ACSN) ---------------------------------------------------------

    output$ScoreVar <- renderPlot({

      Data <- GetACSN()

      if(!is.null(Data)){

        if((input$projtype == "Module") & !is.null(Data$ScoreVar)){
          p <- ggplot2::ggplot(data = data.frame(y = Data$ScoreVar,
                                                 x = rep(1, length(Data$ScoreVar))),
                               mapping = ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_boxplot() +
            ggplot2::labs(y = "Gene score variance (per gene with multiplicity > 1)", x = "") +
            ggplot2::scale_y_log10() +
            ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                           axis.text.x=ggplot2::element_blank(),
                           axis.ticks.x=ggplot2::element_blank())

          print(p)
        }

      }

    })

    # Select samples by group (ACSN) ---------------------------------------------------------

    observeEvent(input$selByGroup, {

      Groups <- GetData()$Groups

      updateCheckboxGroupInput(session, "selSamples", selected = unique(c(names(Groups[Groups == input$selgroup]),
                                                                          input$selSamples)))

    }, ignoreNULL = FALSE, ignoreInit = TRUE)


    # Clear selection (ACSN) ---------------------------------------------------------

    observeEvent(input$selNone, {

      updateCheckboxGroupInput(session, "selSamples", selected = "")

    }, ignoreNULL = FALSE, ignoreInit = TRUE)


    # Compute Correlations (Gene contribution) ---------------------------------------------------------

    GetCorrs <- eventReactive(input$doCorr, {

      ExpMat <- GetData()$ExpMat
      RomaData <- GetData()$RomaData
      ModuleID <- SelectedGS()

      CorInfo <- GetCorrelations(RomaData = RomaData, Selected = ModuleID, MatData = ExpMat,
                                 Methods = input$cortype, ConfLevel = as.numeric(input$corlelvel))

      return(list(CorInfo = CorInfo, Method = input$cortype, ModuleID = ModuleID))

    }, ignoreInit = TRUE)

    output$CorrCI <- renderPlot({

      AllGenesCorr <- GetCorrs()$CorInfo$Genes

      if(GetCorrs()$Method != input$cortype | GetCorrs()$ModuleID != SelectedGS()){
        updateSelectInput(session, "availGenes", choices = list())
        return(NULL)
      }

      # Filter genes

      AllGenesCorr <- AllGenesCorr[as.numeric(AllGenesCorr[,"p.val"]) <= (1 - as.numeric(input$corlelvel)), ]

      # Get direction

      SelGenes <- NULL

      if(input$contribType == "Positive"){
        SelGenes <- which(as.numeric(AllGenesCorr[,"cor"]) > 0)
      }

      if(input$contribType == "Negative"){
        SelGenes <- which(as.numeric(AllGenesCorr[,"cor"]) < 0)
      }

      if(length(SelGenes) == 0){
        return(NULL)
        updateSelectInput(session, "availGenes", choices = list())
      }

      GeneNames <- AllGenesCorr[SelGenes, "gene"]
      GeneLabels <- paste(GeneNames, signif(as.numeric(AllGenesCorr[SelGenes, "cor"]), 4),
                          sep = ' | cor = ')

      RetList <- as.list(GeneNames)
      names(RetList) <- GeneLabels

      updateSelectInput(session, "availGenes",
                        choices = RetList[order(as.numeric(AllGenesCorr[SelGenes, "cor"]))])


      AllGenesCorr.DF <- data.frame(AllGenesCorr[SelGenes,])
      AllGenesCorr.DF$cor <- as.numeric(as.character(AllGenesCorr.DF$cor))
      AllGenesCorr.DF$p.val <- as.numeric(as.character(AllGenesCorr.DF$p.val))
      AllGenesCorr.DF$gene <- factor(as.character(AllGenesCorr.DF$gene),
                                     levels = as.character(AllGenesCorr.DF$gene)[order(AllGenesCorr.DF$cor)])

      if(input$cortype == "pearson"){
        AllGenesCorr.DF$ci.low <- as.numeric(as.character(AllGenesCorr.DF$ci.low))
        AllGenesCorr.DF$ci.high <- as.numeric(as.character(AllGenesCorr.DF$ci.high))

        p <- ggplot2::ggplot(data = AllGenesCorr.DF,
                             ggplot2::aes(x = gene, y = cor, ymin = ci.low, ymax = ci.high)) +
          ggplot2::geom_pointrange() + ggplot2::geom_point() + ggplot2::coord_flip() +
          ggplot2::labs(y = "Correlation")

      } else {

        p <- ggplot2::ggplot(data = AllGenesCorr.DF,
                             ggplot2::aes(x = gene, y = cor)) +
          ggplot2::geom_point() + ggplot2::coord_flip() +
          ggplot2::labs(y = "Correlation")

      }


      if(input$contribType == "Positive"){
        p <- p + ggplot2::scale_y_continuous(limits = c(0, 1))
      }

      if(input$contribType == "Negative"){
        p <- p + ggplot2::scale_y_continuous(limits = c(-1, 0))
      }

      print(p)

    })


    output$ExpProj <- renderPlot({

      if(input$availGenes != ""){

        GeneExp <- GetData()$ExpMat[input$availGenes, ]
        ModScore <- GetData()$RomaData$SampleMatrix[SelectedGS(),]

        SampleNames <- intersect(names(GeneExp), names(ModScore))

        NewDF <- data.frame(Samples = SampleNames,
                            Expression = GeneExp[SampleNames],
                            Score = ModScore[SampleNames],
                            Groups = GetData()$Groups)

        p <- ggplot2::ggplot(data = NewDF, ggplot2::aes(x = Score, y = Expression)) +
          ggplot2::geom_smooth() +
          ggplot2::geom_point(mapping = ggplot2::aes(color = Groups)) +
          ggplot2::labs(x = "Module score", y = "Gene expression")

        print(p)
      }


    })


    # Save data locally ---------------------------------------------------------

    observe({

      fileinfo <- parseSavePath(Volumes, input$save)

      if (nrow(fileinfo) > 0) {

        RomaData <- GetData()$RomaData
        ModuleList <- GetModuleList()
        ExpMat <- GetData()$ExpMat
        Groups <- GetData()$Groups

        print(is.null(RomaData))
        print(is.null(ModuleList))
        print(is.null(ExpMat))
        print(is.null(Groups))

        rRomaDashData <- list(RomaData = RomaData,
                              ModuleList = ModuleList,
                              ExpMat = ExpMat,
                              Groups = Groups)

        print(paste("Saving to", as.character(fileinfo$datapath)))

        saveRDS(rRomaDashData, as.character(fileinfo$datapath))
      }

    })




    # Load data locally ---------------------------------------------------------

    LoadFromServer <- reactive({

      fileinfo <- parseFilePaths(Volumes, input$load)

      if (nrow(fileinfo) > 0) {
        print("Setting local timestamp")
        TimeVect["Local"] <<- Sys.time()
        return(as.character(fileinfo$datapath))
      } else {
        return(NULL)
      }

    })

  }









  # Run the application
  shinyApp(ui = ui, server = server)

}




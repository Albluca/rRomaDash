

Initialize <- function(RomaData, ExpMat, Groups){


  if(is.null(RomaData)){
    return(list(FoundSamp = NULL, Groups = Groups, AddInfo = NULL,
                SelList = list(), SelListAF = list(), GSList = list(),
                PCAProj = NULL, ProcessedSamples = NULL))
  }

  ProcessedSamples <- colnames(RomaData$ProjMatrix)

  if(!is.null(Groups)){

    FoundSamp <- intersect(colnames(RomaData$ProjMatrix), names(Groups))
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
    SelListAF <- list("Mean" = "mean", "Median" = "median", "Std. dev." = "sd", "IQR" = "IQR")

  } else {

    Groups <- rep("N/A", length(ProcessedSamples))
    names(Groups) <- ProcessedSamples

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
    names(GSList) <- as.list(rownames(RomaData$ProjMatrix))
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












rRomaDash <- function(RomaData = NULL,
                      ExpMat = NULL,
                      Groups = NULL,
                      Interactive = FALSE) {

  library(rRoma)
  library(shiny)

   # preprocess data ---------------------------------------------------------

  library(plotly)
  if(Interactive){
    print("Using plotly. This can cause problems on some systems. Try setting 'Interactive = FALSE' if errors are encountered")
  } else {
    if(R.utils::isPackageLoaded("plotly")){
      print("Detaching plotly.")
      detach("package:plotly", unload=TRUE)
    }
  }

  InitReturn <- Initialize(RomaData, ExpMat, Groups)

  FoundSamp <- InitReturn$FoundSamp
  Groups <- InitReturn$Groups
  AddInfo <- InitReturn$AddInfo
  SelList <- InitReturn$SelList
  SelListAF <- InitReturn$SelListAF
  GSList <- InitReturn$GSList
  PCAProj <- InitReturn$PCAProj
  ProcessedSamples <- InitReturn$ProcessedSamples

  tSNEProj <- PCAProj

  # define ui ---------------------------------------------------------

  ui <- navbarPage("rRoma dashboard",

             # Perform analysis ---------------------------------------------------------
             tabPanel("Analyze Data",

                      pageWithSidebar(

                        # Application title
                        headerPanel(""),

                        # Sidebar with a slider input
                        sidebarPanel(
                          actionButton("doROMA", "Execute rROMA")
                        ),

                        # Show a plot of the generated distribution
                        mainPanel(
                          tabsetPanel(

                            tabPanel("Genesets",

                                     selectInput("gmtsrc", "Geneset source:",
                                                 list("File" = "File",
                                                      "MSig" = "MSig",
                                                      "Predefined", "Predefined")),

                                     conditionalPanel(
                                       condition="input.gmtsrc == 'File'",
                                       fileInput("gmtfile", "Choose a GMT file", accept = c(".gmt"))
                                     ),

                                     conditionalPanel(
                                       condition="input.gmtsrc == 'Predefined'",
                                       textInput("msigkw", "Keywords", "hallmark"),
                                       selectInput("gmtlist", "Available genesets:",
                                                   list("A" = "A",
                                                        "B" = "B",
                                                        "C", "C"))
                                     ),
                                     textInput("msigkw", "Keywords", "hallmark"),
                                     checkboxInput("msigkwall", "search all keywords", FALSE)


                                     ),

                            tabPanel("Data"


                                     ),

                            tabPanel("Parameters"


                                     ),



                            tabPanel("Output"


                            )

                          )
                        )
                      )

                      ),

             # Data summary ---------------------------------------------------------
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

             # Results ---------------------------------------------------------
             tabPanel("Visualize Results",

                      pageWithSidebar(

                        # Application title
                        headerPanel(NULL),

                        # sidebar ---------------------------------------------------------

                        sidebarPanel(
                          conditionalPanel(
                            condition="input.ResTabs == 'Modules'",
                            selectInput("prjt", "Projectin type:",
                                        list("PCA" = "PCA", "tSNE" = "tSNE"))
                          ),

                          conditionalPanel(
                            condition="input.ResTabs == 'Modules' && input.prjt == 'tSNE'",
                            sliderInput("perp", "tSNE perplexity:",
                                        min = 0,  max = 50,  value = 0, step = .1),
                            actionButton("dotSNE", "Compute tSNE")
                          ),

                          conditionalPanel(
                            condition="input.ResTabs == 'Modules'",
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
                            condition="input.ResTabs == 'Correlation'",
                            selectInput("cortype", "Correlation method:",
                                        list("Pearson" = "pearson",
                                             "Kendall" = "kendall",
                                             "Spearman" = "spearman"))
                          ),

                          conditionalPanel(
                            condition="input.ResTabs == 'Heatmap'",
                            checkboxInput("gscol", "Samples on columns", FALSE),
                            checkboxInput("saclus", "Cluster samples / groups", FALSE)
                          ),

                          conditionalPanel(
                            condition="input.ResTabs == 'Heatmap' || input.ResTabs == 'Correlation'",
                            checkboxInput("gsclus", "Cluster genesets", FALSE),
                            hr()
                          ),

                          selectInput("disp", "Dispersion filter:",
                                      list("Overdispersed" = "Over",
                                           "Underdispersed" = "Under",
                                           "None" = "None")),
                          sliderInput("pdisp", "Log 10 p-value threshold:",
                                      max = 0,  min = -5,  value = -1, step = .1),
                          hr(),

                          selectInput("coord", "Coordination filter:",
                                      list("Overcoordinated" = "Over",
                                           "Undercoordinated" = "Under",
                                           "None" = "None"), selected = "None"),
                          sliderInput("pcoord", "Log 10 p-value threshold:",
                                      max = 0,  min = -5,  value = 0, step = .1),
                          hr(),

                          selectInput("exp", "Expression filter:",
                                      list("Overexpressed" = "Over",
                                           "Underexpressed" = "Under",
                                           "None" = "None"), selected = "None"),
                          sliderInput("pexp", "Log 10 p-value threshold:",
                                      max = 0,  min = -5,  value = 0, step = .1),
                          hr(),

                          conditionalPanel(
                            condition="input.ResTabs == 'Modules'",
                            selectInput("gs", "GeneSet:",
                                        GSList),
                            htmlOutput("info")
                          ),

                          conditionalPanel(
                            condition="input.ResTabs == 'Heatmap'",
                            sliderInput("llim", "Lower limit",
                                        max = 0,  min = -5,  value = -5, step = .1),
                            sliderInput("ulim", "Upper limit",
                                        max = 5,  min = 0,  value = 5, step = .1)
                          ),

                          conditionalPanel(
                            condition="input.ResTabs == 'Correlation'",
                            selectInput("gs_x", "GeneSet (x axis):",
                                        GSList),
                            selectInput("gs_y", "GeneSet (y axis):",
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
                                      tabPanel(title = "Modules",
                                               if(Interactive){
                                                 tabPanel("Plot",
                                                          plotlyOutput("scatPlot"),
                                                          plotOutput("boxPlot"),
                                                          plotOutput("SamplesBoxPlot")
                                                 )
                                               } else {
                                                 tabPanel("Plot",
                                                          plotOutput("scatPlot"),
                                                          plotOutput("boxPlot"),
                                                          plotOutput("SamplesBoxPlot")
                                                 )
                                               }
                                      ),

                                      tabPanel("Heatmap", id = "tab2",
                                               plotOutput("hmPlot", height = "800px")
                                               ),

                                      tabPanel("Correlation", id = "tab3",
                                               if(Interactive){
                                                 tabPanel("Plot",
                                                          plotOutput("CorHmPlot", height = "800px"),
                                                          plotlyOutput("CorScatPlot")
                                                 )
                                               } else {
                                                 tabPanel("Plot",
                                                          plotOutput("CorHmPlot", height = "800px"),
                                                          plotOutput("CorScatPlot")
                                                 )
                                               }

                                               )
                          )

                        )

                      )

                      ),

             # Save / Load ---------------------------------------------------------
             tabPanel("Save/Load",

                      fluidPage(
                        titlePanel(""),

                        fluidRow(
                          column(6, wellPanel(
                              helpText("Uploada previously performed rRoma analysis"),
                              fileInput("prev.rRoma", "Choose an rRoma file", accept = c(".rds"))
                            )
                          ),

                          column(6,wellPanel(
                            helpText("Download the rRoma analysis"),
                            downloadButton("downloadData", "Download")
                          )
                          )
                        )

                      )
                    )

  )




  # define server ---------------------------------------------------------

  server <- function(input, output, session) {

    options(shiny.maxRequestSize=250*1024^2)


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

      RomaData <- GetData()$RomaData

      if(is.null(RomaData)){
        return(NULL)
      }

      ModuleDF <- data.frame(Sample = names(RomaData$InputPars$Grouping),
                             Group = RomaData$InputPars$Grouping)
      rownames(ModuleDF) <- NULL

      ModuleDF[order(ModuleDF$Sample),]

    })



    # download data ---------------------------------------------------------

    output$downloadData <- downloadHandler(

      filename = function() {
        paste('rRoma-', Sys.Date(), '.rds', sep='')
      },

      content = function(con) {
        rRomaDashData <- list(RomaData = GetData()$RomaData,
                              ExpMat = GetData()$ExpMat,
                              Groups = GetData()$Groups)
        saveRDS(rRomaDashData, con)
      }

    )

    # load data ---------------------------------------------------------

    GetData <- reactive({

      inFile <- input$prev.rRoma

      if(is.null(inFile)){

        return(list(
          RomaData = RomaData,
          Groups = Groups,
          ExpMat = ExpMat,
          FoundSamp = FoundSamp,
          Groups = Groups,
          AddInfo = AddInfo,
          GSList = GSList,
          PCAProj = PCAProj,
          ProcessedSamples = ProcessedSamples
        ))

      } else {

        print(paste("Loading", inFile$datapath))

        LoadedData <- readRDS(inFile$datapath)

        InitReturn <- Initialize(LoadedData$RomaData, LoadedData$ExpMat, LoadedData$Groups)

        SelList <<- InitReturn$SelList
        SelListAF <<- InitReturn$SelListAF
        updateSelectInput(session, inputId = "htype", choices = SelList)
        updateSelectInput(session, inputId = "aggfun", choices = SelListAF)

        return(list(
          RomaData = LoadedData$RomaData,
          Groups = LoadedData$Groups,
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

        tSNEProj <- Rtsne::Rtsne(X = t(ExpMat[,ProcessedSamples]), perplexity = input$perp)$Y
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
      Idx <- SelectedIdx()

      if(length(Idx)>0){
        GSList <- Idx
        names(GSList) <- as.list(rownames(RomaData$ProjMatrix)[Idx])
        GSList <- GSList[order(names(GSList))]
      } else {
        GSList <- list(" " = "")
      }

      updateSelectInput(session, "gs", choices = GSList, selected = GSList[[1]])
      updateSelectInput(session, "gs_x", choices = GSList, selected = GSList[[1]])
      updateSelectInput(session, "gs_y", choices = GSList, selected = GSList[[1]])

    })

    # print module information ---------------------------------------------------------

    output$info <- renderUI({

      RomaData <- GetData()$RomaData

      if(is.null(RomaData)){
        return(NULL)
      }

      str1 <- paste("Overdispersed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 1]))
      str2 <- paste("Underdispersed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 2]))
      str3 <- paste("Overcoordinated PV = ", signif(RomaData$PVVectMat[SelectedGS(), 3]))
      str4 <- paste("Undercoordinated PV = ", signif(RomaData$PVVectMat[SelectedGS(), 4]))
      str5 <- paste("Overexpressed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 5]))
      str6 <- paste("Underexpressed PV = ", signif(RomaData$PVVectMat[SelectedGS(), 6]))
      HTML(paste(str1, str2, str3, str4, str5, str6, sep = '<br/>'))
    })


    # sample diag boxplots ---------------------------------------------------------

    output$SamplesBoxPlot <- renderPlot({

      RomaData <- GetData()$RomaData

      if(is.null(RomaData)){
        return(NULL)
      }

      SeID <- SelectedGS()

      Sampled.DF <- t(rbind(sapply(RomaData$ModuleSummary[[SeID]]$SampledExp, "[[", "ExpVar"),
                            sapply(RomaData$ModuleSummary[[SeID]]$SampledExp, "[[", "MedianExp")))

      Sampled.DF <- cbind(Sampled.DF, Sampled.DF[,2]/Sampled.DF[,1])

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

      p <- ggplot2::ggplot(data = Sampled.DF[Sampled.DF$Org == "Samples", ], ggplot2::aes(y = value, x = X2, color = Org)) +
        ggplot2::geom_boxplot() + ggplot2::facet_wrap(~X2, scales = "free", ncol = 4) +
        ggplot2::geom_point(data = Sampled.DF[Sampled.DF$Org == "Data", ],
                   ggplot2::aes(y = value, x = X2, color = Org), inherit.aes = FALSE, size = 3)

        print(p)

    })


    # single sample boxplot ---------------------------------------------------------

    output$boxPlot <- renderPlot({

      RomaData <- GetData()$RomaData
      Groups <- GetData()$Groups
      ProcessedSamples <- GetData()$ProcessedSamples

      GetComb <- function(GrpLevs) {
        RetList <- list()
        for(i in 1:length(GrpLevs)){
          for(j in 1:length(GrpLevs)){
            if(i<j){
              RetList[[length(RetList)+1]] <- c(i, j)
            }
          }
        }
        return(RetList)
      }

      p <- ggplot2::ggplot(data = data.frame(Score = RomaData$ProjMatrix[SelectedGS(), ProcessedSamples],
                                             Group = Groups[ProcessedSamples]),
                           ggplot2::aes(x = Group, y = Score)) +
        ggplot2::geom_boxplot()

      if(length(unique(Groups[ProcessedSamples]))>1){
        p <- p + ggsignif::geom_signif(comparisons = GetComb(unique(Groups[ProcessedSamples])),
                                       map_signif_level=TRUE, test = "wilcox.test", step_increase = .1)
      }

      if(!is.na(SelectedGS())){
        print(p)
      }

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

        p <- ggplot2::ggplot(data = data.frame(Comp1 = Projs[,2],
                                                       Comp2 = Projs[,2],
                                                       Score = RomaData$ProjMatrix[SelectedGS(), ProcessedSamples],
                                                       Group = Groups[ProcessedSamples]),
                                     ggplot2::aes(x = Comp1, y = Comp2, shape = Group, color = Score)) +
          ggplot2::scale_color_gradient2(low = "blue", high = "red", mid = "white") +
          ggplot2::labs(x = "Component 1", y = "Component 2", shape = "",
                        title = rownames(RomaData$ProjMatrix)[SelectedGS()]) +
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
                                               Score = RomaData$ProjMatrix[SelectedGS(), ProcessedSamples],
                                               Group = Groups[ProcessedSamples]),
                             ggplot2::aes(x = Comp1, y = Comp2, shape = Group, color = Score)) +
          ggplot2::scale_color_gradient2(low = "blue", high = "red", mid = "white") +
          ggplot2::labs(x = "Component 1", y = "Component 2", shape = "",
                        title = rownames(RomaData$ProjMatrix)[SelectedGS()]) +
          ggplot2::geom_point(size = 3)
        print(p)
      })
    }



    # heatmap ---------------------------------------------------------

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

      PlotMat <- RomaData$ProjMatrix[Idx,]

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
            pheatmap::pheatmap(PlotMat, color = BaseCol[UseCol], breaks = MyBreaks,
                               cluster_rows = input$gsclus, cluster_cols = input$saclus,
                               annotation_col = AddInfo, main = "Module scores across samples")
          } else {
            pheatmap::pheatmap(t(PlotMat), color = BaseCol[UseCol], breaks = MyBreaks,
                               cluster_rows = input$saclus, cluster_cols = input$gsclus,
                               annotation_row = AddInfo, main = "Module scores across samples")
          }


        }

        if(length(Idx) == 1){

          names(PlotMat) <- colnames(RomaData$ProjMatrix)

          pheatmap::pheatmap(t(PlotMat), BaseCol[UseCol], breaks = MyBreaks,
                             cluster_rows = FALSE, cluster_cols = FALSE,
                             main = paste("Score of", rownames(RomaData$ProjMatrix)[Idx]))

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
            pheatmap::pheatmap(Aggmat, color = BaseCol[UseCol], breaks = MyBreaks,
                               cluster_rows = input$gsclus, cluster_cols = input$saclus,
                               main = "Module scores across groups")
          } else {
            pheatmap::pheatmap(t(Aggmat), color = BaseCol[UseCol], breaks = MyBreaks,
                               cluster_rows = input$saclus, cluster_cols = input$gsclus,
                               main = "Module scores across groups")
          }

        }

        if(length(Idx) == 1){

          # names(PlotMat) <- colnames(RomaData$ProjMatrix)
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
                             main = paste("Score of", rownames(RomaData$ProjMatrix)[Idx]))
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
          pheatmap::pheatmap(cor(t(RomaData$ProjMatrix[Idx,]), method = input$cortype),
                             color = BaseCol[UseCol], breaks = MyBreaks,
                             cluster_rows = input$gsclus, cluster_cols = input$gsclus,
                             main = "Correlation across samples")
        }

        if(length(Idx) == 1){
          NULL
        }
      }

      if(input$htype == "group"){

        if(length(Idx) > 1){

          SplitData <- split(data.frame(t(RomaData$ProjMatrix[Idx,FoundSamp])), f=AddInfo$Groups)

          Aggmat <- sapply(SplitData, function(x) {
            apply(x, 2, get(input$aggfun))
          })

          pheatmap::pheatmap(cor(t(Aggmat), method = input$cortype), color = BaseCol[UseCol], breaks = MyBreaks,
                             cluster_rows = input$gsclus, cluster_cols = input$gsclus,
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

        CTitle <- cor(RomaData$ProjMatrix[as.integer(input$gs_x), ProcessedSamples],
                     RomaData$ProjMatrix[as.integer(input$gs_y), ProcessedSamples],
                     method = input$cortype)


        if(input$htype == "sample"){

          p <- ggplot2::ggplot(data = data.frame(XVal = RomaData$ProjMatrix[as.integer(input$gs_x), ProcessedSamples],
                                                 YVal = RomaData$ProjMatrix[as.integer(input$gs_y), ProcessedSamples],
                                                 Group = Groups[ProcessedSamples]),
                               ggplot2::aes(x = XVal, y = YVal, shape = Group)) +
            ggplot2::labs(x = XLab, y = YLab, shape = "", title = paste("Corr =", signif(CTitle, 5))) +
            ggplot2::geom_point(ggplot2::aes(text = ProcessedSamples))

        }


        if(input$htype == "group"){

          AggX <- aggregate(RomaData$ProjMatrix[as.integer(input$gs_x), ProcessedSamples], list(AddInfo$Groups), get(input$aggfun))
          AggY <- aggregate(RomaData$ProjMatrix[as.integer(input$gs_y), ProcessedSamples], list(AddInfo$Groups), get(input$aggfun))

          p <- ggplot2::ggplot(data = data.frame(XVal = AggX[,2],
                                                 YVal = AggY[,2],
                                                 Group = AggX[,1]),
                               ggplot2::aes(x = XVal, y = YVal, shape = Group)) +
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

        CTitle <- cor(RomaData$ProjMatrix[as.integer(input$gs_x), ProcessedSamples],
                      RomaData$ProjMatrix[as.integer(input$gs_y), ProcessedSamples],
                      method = input$cortype)

        if(input$htype == "sample"){

          p <- ggplot2::ggplot(data = data.frame(XVal = RomaData$ProjMatrix[as.integer(input$gs_x), ProcessedSamples],
                                                 YVal = RomaData$ProjMatrix[as.integer(input$gs_y), ProcessedSamples],
                                                 Group = Groups[ProcessedSamples]),
                               ggplot2::aes(x = XVal, y = YVal, shape = Group)) +
            ggplot2::labs(x = XLab, y = YLab, shape = "", title = paste("Corr =", signif(CTitle, 5))) +
            ggplot2::geom_point()

        }


        if(input$htype == "group"){

          AggX <- aggregate(RomaData$ProjMatrix[as.integer(input$gs_x), ProcessedSamples], list(AddInfo$Groups), get(input$aggfun))
          AggY <- aggregate(RomaData$ProjMatrix[as.integer(input$gs_y), ProcessedSamples], list(AddInfo$Groups), get(input$aggfun))


           p <- ggplot2::ggplot(data = data.frame(XVal = AggX[,2],
                                                 YVal = AggY[,2],
                                                 Group = AggX[,1]),
                               ggplot2::aes(x = XVal, y = YVal, shape = Group)) +
            ggplot2::labs(x = XLab, y = YLab, shape = "", title = paste("Corr =", signif(CTitle, 5))) +
            ggplot2::geom_point()

        }

        print(p)

      })
    }
  }

  # Run the application
  shinyApp(ui = ui, server = server)

}




#' @import dplyr plyr bnutil pgUpstream ggplot2 shinyjs
#' @export
operatorProperties <- function() {
  list(
    list("Kinase_family", list("PTK", "STK")),
    list("Lock_kinase_family", list("No", "Yes"))
  )
}

#' @export
shinyServerRun <- function(input, output, session, context) {
  getFolderReactive <- context$getRunFolder()
  getDataReactive <- context$getData()
  getPropertiesReactive <- context$getPropertiesAsMap()
  output$body <- renderUI({
    sidebarLayout(
      sidebarPanel(
        tags$div(HTML("<strong><font color = #6895d1>Upstream Kinase Analysis - January 2022</font></strong>")),
        tags$hr(),
        actionButton("start", "Start"),
        tags$hr(),
        textOutput("status")
      ),
      mainPanel(
        shinyjs::useShinyjs(),
        tabsetPanel(
          tabPanel(
            "Basic Settings",
            tags$hr(),
            selectInput("kinasefamily", label = "Kinase family", choices = c("PTK", "STK")),
            sliderInput("scan", "Scan Rank From-To", min = 1, max = 20, value = c(4, 12), round = TRUE),
            sliderInput("nperms", "Number of permutations", min = 100, max = 1000, value = 500, step = 100, round = TRUE)
          ),
          tabPanel(
            "Advanced Settings",
            tags$hr(),
            helpText("Set weights for database types:"),
            sliderInput("wIviv", "In Vitro / In Vivo", min = 0, max = 1, value = 1),
            sliderInput("wPhosphoNET", "In Silico (PhosphoNET)", min = 0, max = 1, value = 1),
            tags$hr(),
            helpText("Set minimal sequence similarity required to link a peptide to a phosphosite"),
            sliderInput("seqHom", "Minimal sequence similarity", min = 0, max = 1, value = 0.9),
            tags$hr(),
            numericInput("minPScore", "Minimal PhosphoNET prediction score", value = 300),
            tags$hr(),
            checkboxInput("seed", "Set seed for random permutations")
          )
        )
      )
    )
  })

  getGroupingLabel <- function(aData) {
    if (length(aData$colorLabels) > 1) {
      stop("Need at most one Color to define the grouping")
    } else if (length(aData$colorLabels) == 1) {
      groupingLabel <- aData$colorColumnNames
      grp <- as.factor((aData$data[[groupingLabel]]))
      if (length(levels(grp)) != 2) {
        stop(paste("Wrong number of groups, found:", levels(grp)))
      }

      df <- subset(aData$data, rowSeq = 1)
      grp <- as.factor((df[[groupingLabel]]))
      if (min(summary(grp)) < 2) {
        stop("Need at least 2 observations per group.")
      }
      return(groupingLabel)
    } else {
      if (max(aData$data$colSeq) > 1) {
        stop(paste("Can not run this data without a grouping factor.", aData$colorLabels))
      } else {
        return(NULL)
      }
    }
  }

  # DB = UpstreamAppTest2021::UpstreamDatabase # note: hard code reference to the package name!
  DB <- load(system.file("extdata", "220106-86402-87102_UpstreamDb.RData", package = "UpstreamApp2"))
  DB <- get(DB)
  nid <- showNotification("Press Start to start the analysis.", duration = NULL, type = "message", closeButton = FALSE)
  updateSliderInput(session, "seqHom", min = min(DB$PepProtein_SeqSimilarity))

  observe({
    getData <- getDataReactive$value
    getFolder <- getFolderReactive$value
    getProperties <- getPropertiesReactive$value
    if (is.null(getProperties)) {
      return()
    }
    prop <- getProperties()
    updateSelectInput(session, "kinasefamily", selected = prop$Kinase_family)

    if (is.null(getData)) {
      return()
    }
    adf <- getData()
    output$status <- renderText({
      grp <- getGroupingLabel(adf)
      if (input$start == 0) {
        if (prop$Lock_kinase_family == "Yes") {
          shinyjs::disable("kinasefamily")
        }
        if (!is.null(grp)) {
          return(paste("Grouping factor:", grp))
        } else {
          return("Grouping factor: none")
        }
      }

      shinyjs::disable("kinasefamily")
      shinyjs::disable("scan")
      shinyjs::disable("nperms")
      shinyjs::disable("wIviv")
      shinyjs::disable("wPhosphoNET")
      shinyjs::disable("seqHom")
      shinyjs::disable("minPScore")

      if (!adf$hasUniqueDataMapping) stop("Mapping error: not all input data is unique")

      if (adf$hasMissingCells) stop("Missing values are not allowed.")

      if (adf$getMaxNPerCell() > 1) stop("More than one value per cell in the cross-tab view is not allowed.")
      df <- adf$getData()

      if (adf$hasZeroScaleRows) {
        zIdx <- adf$getZeroScaleRows()
        df <- df %>% filter(!(rowSeq %in% zIdx))
        msg <- paste("Warning:", length(zIdx), "peptides with zero scale have been removed from the data")
        showNotification(ui = msg, duration = NULL, type = "warning")
      }

      if (input$kinasefamily == "PTK") {
        DB <- DB %>% filter(PepProtein_Residue == "Y")
      } else if (input$kinasefamily == "STK") {
        DB <- DB %>% filter(PepProtein_Residue == "S" | PepProtein_Residue == "T")
      } else {
        stop("Unknown value for kinase family")
      }

      DB <- DB %>% filter(PepProtein_SeqSimilarity >= input$seqHom)
      DB <- DB %>% filter(Kinase_PKinase_PredictorVersion2Score >= input$minPScore | Database == "iviv")
      DB <- DB %>% filter(Kinase_Rank <= input$scan[2])

      nCores <- detectCores()
      msg <- paste("Please wait ... running analysis. Using", nCores, "cores.")
      showNotification(ui = msg, id = nid, type = "message", closeButton = FALSE, duration = NULL)

      if (input$seed) {
        set.seed(42)
      }

      if (!is.null(grp)) {
        df$grp <- as.factor(df[[grp]])
        aResult <- pgScanAnalysis2g(df,
          dbFrame = DB,
          scanRank = input$scan[1]:input$scan[2],
          nPermutations = input$nperms,
          dbWeights = c(
            iviv = input$wIviv,
            PhosphoNET = input$wPhosphoNET
          )
        )
      } else {
        aResult <- pgScanAnalysis0(df,
          dbFrame = DB,
          scanRank = input$scan[1]:input$scan[2],
          nPermutations = input$nperms,
          dbWeights = c(
            iviv = input$wIviv,
            PhosphoNET = input$wPhosphoNET
          )
        )
      }
      showNotification(ui = "Done", id = nid, type = "message", closeButton = FALSE)
      aFull <- ldply(aResult, .fun = function(.) {
        return(data.frame(.$aResult, mxRank = .$mxRank))
      })


      settings <- data.frame(
        setting = c("Kinase family", "ScanRank Min", "ScanRank Max", "Number of Permutations", "In Vitro In Vitro weight", "PhosphoNET weight", "Min PhosphoNet score", "Min Sequence Similarity"),
        value = c(input$kinasefamily, input$scan[1], input$scan[2], input$nperms, input$wIviv, input$wPhosphoNET, input$minPScore, input$seqHom)
      )

      spath <- file.path(getFolder(), "runData.RData")
      save(file = spath, df, aResult, aFull, settings)
      dpath <- file.path(getFolder(), "runDb.RData")
      save(file = dpath, DB)
      out <- data.frame(rowSeq = 1, colSeq = 1, dummy = NaN)
      meta <- data.frame(labelDescription = c("rowSeq", "colSeq", "dummy"), groupingType = c("rowSeq", "colSeq", "QuantitationType"))
      result <- AnnotatedData$new(data = out, metadata = meta)
      context$setResult(result)
      return("Done")
    })
  })
}

#' @export
#' @import colourpicker kinaseTreeParser
shinyServerShowResults <- function(input, output, session, context) {
  getFolderReactive <- context$getRunFolder()
  output$body <- renderUI(
    sidebarLayout(
      sidebarPanel(
        tags$div(HTML("<strong><font color = #6895d1>Upstream Kinase Analysis - January 2022</font></strong>")),
        tags$hr(),
        textOutput("grpText"),
        tags$hr(),
        sliderInput("minsetsize", "Include results based on peptide sets with minimal size of:",
          min = 1,
          max = 5,
          step = 1,
          value = 3
        ),
        width = 3
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Upstream Kinase Score",
            helpText("This plot shows putative upstream kinases ranked by their Final Score (median) or the value of the Kinase Statistic (median)."),
            helpText("Each point is the result of an individual analysis with a different rank cut-off for adding upstream kinases for peptides."),
            helpText("The size of the points indicates the size of the peptide set used for a kinase in the corresponding analysis."),
            helpText("The color of the points indicates the specificity score resulting from the corresponding analysis."),
            helpText("Red colored kinase labels on the y-axis mean the kinases could be dual-specificity (both PTK and STK). Caution is advised, check literature to see whether the result is valid."),
            selectInput("spsort", "Sort score plot on", choices = list(Score = "score", Statistic = "stat")),
            actionLink("saveScorePlot", "Save score plot"),
            plotOutput("scorePlot", height = "1400px")
          ),
          tabPanel(
            "Kinase Volcano",
            helpText("This plot shows the median Final Score (y-axis) versus the mean Kinase Statistic (x-axis) of putative upstream kinases."),
            helpText("The size of the kinase names indicates the size of the peptide set used for a kinase in the corresponding analysis."),
            helpText("The color of the kinase names indicates the specificity score resulting from the corresponding analysis."),
            actionLink("saveVolcanoPlot", "Save volcano plot"),
            plotOutput("volcanoPlot", height = "800px")
          ),
          tabPanel(
            "Kinase Details",
            helpText("These are kinase details"),
            selectInput("showKinase", "Select kinase", choices = ""),
            tags$hr(),
            textOutput("kinGroup"),
            actionLink("uniprot", "..."),
            tags$hr(),
            tableOutput("kinaseSummary"),
            tabsetPanel(
              tabPanel(
                "Details Table",
                helpText(""),
                actionLink("saveDetailsTable", "Save details table"),
                tableOutput("kinaseDetails")
              ),
              tabPanel(
                "Per peptide plot",
                helpText(""),
                actionLink("savePerpeptidePlot", "Save per peptide plot"),
                plotOutput("perPeptidePlot", height = "800px")
              )
            )
          ),
          tabPanel(
            "Report",
            helpText("The table below shows the settings that were used for this analysis."),
            tableOutput("InfoSettings"),
            helpText("The below shows the summary results of the analysis"),
            actionLink("saveSummaryResults", "Save summary results"),
            tableOutput("SummaryTable")
          )
          # ,tabPanel("Kinase tree",
          #          tags$a(href = "http://kinhub.org/kinmap/","The summary data can be saved as a file that can be used to map the data to a phylogenetic tree using the external websit Kinmap."),
          #          helpText(""),
          #          actionButton("saveKinMap", "Save data as KinMap file"),
          #          tags$hr(),
          #          helpText("The specificty score will be mapped to the symbol size"),
          #          textInput("sclow", "Scale score from:", 0),
          #          textInput("schigh","Scale score to:", 2),
          #          textInput("scmax", "max symbol size", 50),
          #          tags$hr(),
          #          helpText("The kinase statistic will be mapped to the symbol color"),
          #          textInput("stlow", "Scale kinase statistic from:", -1),
          #          textInput("stmid", "Scale kinase statistic midpoint:", 0),
          #          textInput("sthigh","Scale kinase statistic to:", 1),
          #          colourpicker::colourInput("cllow", label = "Low color", value = "green"),
          #          colourpicker::colourInput("clmid", label = "Mid color", value = "black"),
          #          colourpicker::colourInput("clhigh", label = "High color", value = "red")
          # )
        )
      )
    )
  )

  observe({
    getFolder <- getFolderReactive$value
    if (is.null(getFolder)) {
      return()
    }
    spath <- file.path(getFolder(), "runData.RData")
    load(file = spath)
    # load the run-db "DB" that was used for the analysis
    dpath <- file.path(getFolder(), "runDb.RData")
    load(file = dpath)
    kinase2uniprot <- DB %>%
      group_by(Kinase_UniprotID) %>%
      dplyr::summarise(Kinase_Name = Kinase_Name[1])

    output$grpText <- renderText({
      grp <- df[["grp"]]

      if (!is.null(grp)) {
        txt <- paste("Grouping factor with levels", levels(grp)[1], "and", levels(grp)[2])
      } else {
        return("Grouping factor: none")
      }
    })

    xaxText <- function() {
      grp <- df[["grp"]]
      if (!is.null(grp)) {
        txt <- paste(">0 indicates higher activity in the", as.character(levels(grp)[2]), "group")
      } else {
        return(NULL)
      }
    }

    scorePlot <- reactive({
      aSub <- subset(aFull, nFeatures >= input$minsetsize)
      group_names <- DB %>% distinct(Kinase_Name, Kinase_group)
      aSub <- left_join(aSub, group_names, by = c("ClassName" = "Kinase_Name"))

      cs <- makeScorePlot(aSub, input$spsort)
      xax <- paste("Normalized kinase statistic (", xaxText(), ")", sep = "")

      axis_order <- ggplot_build(cs)$layout$panel_ranges[[1]]$y.labels
      kinLabelGroup <- aSub %>%
        distinct(ClassName, Kinase_group) %>%
        arrange(match(ClassName, axis_order)) %>%
        pull(Kinase_group)

      if (sum(grepl("Y", DB$PepProtein_Residue)) > 0) {
        kinLabelColors <- ifelse(kinLabelGroup == "TK", "black", "red")
      } else if (sum(grepl("Y", DB$PepProtein_Residue)) == 0) {
        kinLabelColors <- ifelse(kinLabelGroup == "TK", "red", "black")
      }

      cs <- cs + ylab(xax) + theme(axis.text.y = element_text(colour = kinLabelColors))
      # cs = cs + ylab(xax)
    })

    perPeptidePlot <- reactive({
      aPlot <- makePerPeptidePlot(df, DB %>% filter(Kinase_Name == input$showKinase))
      return(aPlot)
    })

    output$scorePlot <- renderPlot({
      print(scorePlot())
    })

    aSummary <- reactive({
      aSub <- subset(aFull, nFeatures >= input$minsetsize)
      aSum <- makeSummary(aSub)
      updateSelectInput(session, "showKinase", choices = aSum$ClassName)
      aSum
    })

    volcanoPlot <- reactive({
      vp <- makeVolcanoPlot(aSummary())
    })

    output$volcanoPlot <- renderPlot({
      print(volcanoPlot())
    })

    output$perPeptidePlot <- renderPlot({
      print(perPeptidePlot())
    })

    output$kinaseSummary <- renderTable({
      aSum <- aSummary()
      aKin <- subset(aSum, ClassName == input$showKinase)
      aTable <- data.frame(
        name = c("Mean Kinase Statistic", "Mean Specificity Score", "Mean Significance Score", "Median Combined Score"),
        value = c(aKin$meanStat, aKin$meanFeatScore, aKin$meanPhenoScore, aKin$medianScore)
      )
      colnames(aTable) <- c(input$showKinase, "value")
      return(aTable)
    })

    observeEvent(input$showKinase, {
      upid <- kinase2uniprot %>%
        filter(Kinase_Name == input$showKinase) %>%
        .$Kinase_UniprotID
      aLabel <- paste(input$showKinase, " (", upid, ") on Uniprot.org", sep = "")
      updateActionButton(session, "uniprot", label = aLabel)
    })

    output$kinGroup <- renderText({
      upid <- kinase2uniprot %>%
        filter(Kinase_Name == input$showKinase) %>%
        .$Kinase_UniprotID

      kinGroup <- DB %>%
        filter(Kinase_UniprotID == upid) %>%
        select(Kinase_group) %>%
        distinct() %>%
        pull()
      kinGroupStr <- paste("Kinase", input$showKinase, "belongs to the", kinGroup, "kinase group.")
      return(kinGroupStr)
    })


    detailsTable <- reactive({
      aTable <- makeDetailsTable(df, DB %>% filter(Kinase_Name == input$showKinase))
    })

    output$kinaseDetails <- renderTable({
      aTable <- detailsTable()
    })

    output$InfoSettings <- renderTable({
      getSettingsInfo(settings)
    })

    observeEvent(input$saveDetailsTable, {
      aTable <- detailsTable()
      filename <- file.path(getFolder(), paste(gsub("/", "-", input$showKinase), format(Sys.time(), "%Y%m%d-%H%M.txt")))
      write.table(aTable, filename, sep = "\t", quote = FALSE, row.names = FALSE)
      shell.exec(getFolder())
    })

    observeEvent(input$saveScorePlot, {
      filename <- file.path(getFolder(), paste("UpstreamScorePlot", format(Sys.time(), "%Y%m%d-%H%M.png")))
      ggsave(filename, scorePlot(), device = "png", units = "cm", height = 40, width = 28)
      shell.exec(getFolder())
    })

    observeEvent(input$saveVolcanoPlot, {
      filename <- file.path(getFolder(), paste("UpstreamVolcanoPlot", format(Sys.time(), "%Y%m%d-%H%M.png")))
      ggsave(filename, volcanoPlot(), device = "png", units = "cm", height = 30, width = 30)
      shell.exec(getFolder())
    })

    observeEvent(input$savePerpeptidePlot, {
      filename <- file.path(getFolder(), paste(gsub("/", "-", input$showKinase), "_peptides", format(Sys.time(), "%Y%m%d-%H%M.png")))
      ggsave(filename, perPeptidePlot(), device = "png", units = "cm", height = 40, width = 28)
      shell.exec(getFolder())
    })

    observeEvent(input$uniprot, {
      upid <- kinase2uniprot %>%
        filter(Kinase_Name == input$showKinase) %>%
        .$Kinase_UniprotID
      browseURL(paste("http://www.uniprot.org/uniprot/", upid, sep = ""))
    })



    summaryResultTable <- reactive({
      df <- aSummary()
      df <- left_join(kinase2uniprot, df, by = c("Kinase_Name" = "ClassName"))
      group_names <- DB %>% distinct(Kinase_Name, Kinase_group, Kinase_family)
      df <- left_join(df, group_names, by = c("Kinase_Name" = "Kinase_Name"))
      # relocate doesn't exist in dplyr 0.7.1 so manual arranging
      col_order <- c("Kinase_UniprotID", "Kinase_Name", "Kinase_group", "Kinase_family", "meanFeatScore", "meanPhenoScore", "meanScore", "medianScore", "meanStat", "medianStat", "sdStat", "meanSetSize")
      df <- df[, col_order]
      df <- df %>%
        filter(!is.na(medianScore)) %>%
        arrange(-medianScore)
      df <- df %>% dplyr::rename(
        "Kinase Uniprot ID" = Kinase_UniprotID,
        "Kinase Name" = Kinase_Name,
        "Kinase Group" = Kinase_group,
        "Kinase Family" = Kinase_family,
        "Mean Specificity Score" = meanFeatScore,
        "Mean Significance Score" = meanPhenoScore,
        "Mean Final Score" = meanScore,
        "Median Final score" = medianScore,
        "Mean Kinase Statistic" = meanStat,
        "Median Kinase Statistic" = medianStat,
        "SD Kinase Statitistic" = sdStat,
        "Mean peptide set size" = meanSetSize
      )
    })

    observeEvent(input$saveSummaryResults, {
      aTable <- summaryResultTable()
      filename <- file.path(getFolder(), paste("Summaryresults", format(Sys.time(), "%Y%m%d-%H%M.txt")))
      write.table(aTable, filename, sep = "\t", quote = FALSE, row.names = FALSE)
      shell.exec(getFolder())
    })

    output$SummaryTable <- renderTable({
      summaryResultTable()
    })

    observeEvent(input$saveKinMap, {
      df <- summaryResultTable()
      colnames(df) <- make.names(colnames(df)) # convert formatted column names to valid variable names
      filename <- file.path(getFolder(), paste("KinMap file", format(Sys.time(), "%Y%m%d-%H%M.txt")))
      szFrom <- c(as.numeric(input$sclow), as.numeric(input$schigh))
      szTo <- c(0, as.numeric(input$scmax))
      szScale <- list(from = szFrom, to = szTo)
      clScale <- list(low = as.numeric(input$stlow), mid = as.numeric(input$stmid), high = as.numeric(input$sthigh))
      clr <- c(input$cllow, input$clmid, input$clhigh)
      print(colnames(df))
      treeFile(fname = filename, mappings = Kinase.Name ~ Mean.Kinase.Statistic + Mean.Specificity.Score, data = df, szScale = szScale, clScale = clScale, clValues = clr)
      shell.exec(getFolder())
    })
  })
}
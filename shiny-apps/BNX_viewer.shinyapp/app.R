# BNX_viewer.shinyapp
# A R/shiny tool to create scatterplots or histograms from BNX data
# designed to work with BNX 1.2 format
#
# Stephane Plaisance, VIB Nucleomics Core
# visit our Git: https://github.com/Nucleomics-VIB
# version: 2017-03-10_v1.0
# version 1.1, added density histogram
# version 1.2, added L50 and N50 data and lines
# version 1.3, data.table speedup and progress indicators
# version 1.4, add filtering for labels and remove useless table rows, add table to plot
# © by using this tool, you accept the licence saved under ./www/licence.pdf

library("shiny")
library("shinyBS")
library("readr")
library("stringr")
library("ggplot2")
library("data.table")
library("grid")
library("gridExtra")

# limit upload to 2GB
# the limit was kept low to make the online shiniapps.io version work with a sample BNX file
# such sample can be downloaded from the github page (look into the BNX_viewer.shinyapp Data folder)

# you may uncomment the next line to allow large input files
options(shiny.maxRequestSize=2*1000*1024^2)
# the following test checks if we are running on shinnyapps.io to limit file size dynamically
# ref: https://stackoverflow.com/questions/31423144/how-to-know-if-the-app-is-running-at-local-or-on-server-r-shiny/31425801#31425801
#if ( Sys.getenv('SHINY_PORT') == "" ) { options(shiny.maxRequestSize=1000*1024^2) }

script.version="1.4"

# defaults for controllers
def.selin1 <- 'molLength'
def.selin2 <- 'molAvgIntensity'
def.min1 <- 0
def.max1 <- 2500
def.val1 <- c(100, 2500)
def.min2 <- 0
def.max2 <- 2
def.val2 <- c(0, 1)
def.min3 <- 0
def.max3 <- 100
def.val3 <- c(0, 100)
def.min4 <- 0
def.max4 <- 1000
def.val4 <- c(0, 100)

# custom functions
## convert data.table column types
convert.magic <- function(obj, types){
  for (i in 1:length(obj)){
    FUN <- switch(types[i],
                  character = as.character,
                  numeric = as.numeric,
                  factor = as.factor)
    obj[[i]] <- FUN(obj[[i]])
  }
  obj
}

# round data.frame numbers
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

## return L50 and N50 from a vector
calc.n50 <- function (x){
  if(is.numeric(x)){
    sorted <- sort(x, decreasing = TRUE)
    avg.len <- sum(x)/2
    csum <- cumsum(sorted)
    gt.avg <- as.vector(csum >= avg.len)
    L50 <- min(which(gt.avg == TRUE))
    N50 <- round(sorted[L50], 1)
  } else {
    L50=NA
    N50=NA
  }
  return(list(L50=L50, N50=N50))
}

# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application header
  headerPanel(
  "BNX molecule distribution viewer"
  ),

  # Application title
  titlePanel(
    windowTitle = "BNX molecule distribution and filtering",
    tags$a(href="https://corefacilities.vib.be/nc", target="_blank", img(src='logo.png', align = "right", width="150", height="58.5", alt="VIB Nucleomics Core"))
  ),

  # Choose input BNX file
  sidebarLayout(
    # show file import and molecule filters
    sidebarPanel(
      tags$h4(paste("code version: ", script.version, sep="")),
      downloadButton("downloadData", label = "Download test data"),
      tags$br(),
      tags$a(href="license.pdf", target="_blank", "usage licence"),
      tags$hr(),
      fileInput('file1', 'Choose BNX File', accept='.bnx'),

      tags$h4("choose different axes for the scatterplot"),

      selectInput('Xaxis', 'X-axis:',
                  c("Molecule length" = "molLength",
                    "Molecule Average Intensity" = "molAvgIntensity",
                    "label SNR" = "labelSNR",
                    "number of labels" = "NumberofLabels",
                    "labels per 100k" = "label100kdens"),
                  selected = def.selin1
      ),

      selectInput('Yaxis', 'Y-axis:',
                  c("Molecule length" = "molLength",
                    "Molecule Average Intensity" = "molAvgIntensity",
                    "label SNR" = "labelSNR",
                    "number of labels" = "NumberofLabels",
                    "labels per 100k" = "label100kdens",
                    "NA: histogram density for the X-axis" = "histo"),
                  selected = def.selin2
      ),

      hr(),
      tags$h4("modify ranges below and click ", tags$em("Filter Data")),

      sliderInput('length',
                  "Molecule length (kb):",
                  min = def.min1,
                  max = def.max1,
                  step = 10,
                  value = def.val1
      ),

      sliderInput('mAI',
                  "Molecule AvgIntensity (au):",
                  min = def.min2,
                  max = def.max2,
                  step = 0.05,
                  value = def.val2
      ),

      sliderInput('mSNR',
                  "Molecule SNR (au):",
                  min = def.min3,
                  max = def.max3,
                  step = 1,
                  value = def.val3
      ),
      
      sliderInput('mlabels',
                  "Label count:",
                  min = def.min4,
                  max = def.max4,
                  step = 1,
                  value = def.val4
      ),

      actionButton(inputId='FilterButton', "Filter Data", style='padding:4px; font-weight: bold; font-size:150%'),
      
      actionButton(inputId='PlotButton', "Refresh Plot", style='padding:4px; font-weight: bold; font-size:150%'),

      hr(),
      
      textInput('outfile', 'name for output File:', value="densityPlot"),
      
      selectInput('format', 'Output format (png or pdf):', c("png", "pdf"), selected="pdf"),
      
      downloadButton('downloadPlot', 'Download Plot')
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput('plot1'),
      br(),
      textOutput('cntData'),
      textOutput('n50x'),
      textOutput('n50y'),
      br(),
      div(DT::dataTableOutput('summaryzero'), style = "font-size: 75%; width: 75%")
    )
  )
)

# load data and create scatterplot
server <- function(input, output) {
  output$downloadData <- downloadHandler(
    filename <- function() {
      paste("sample", "bnx", sep=".")
    },

  content <- function(file) {
      file.copy("Data/sample.bnx", file)
    },
    contentType = "application/zip"
  )

  load.data <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    # show progress
    n <- 4
    progress <- shiny::Progress$new()
    progress$set(message = "Importing data", value = 0)
    progress$set(value = 1/n, detail = "reading file in")
    on.exit(progress$close())
    dat <- readLines(inFile$datapath, skipNul=TRUE)
    Sys.sleep(1)
    progress$set(value = 2/n, detail = "filtering rows")
    dat <- dat[grepl("^0\t*", dat, perl = TRUE)]
    Sys.sleep(1)
    progress$set(value = 3/n, detail = "creating table")
    dat <- data.table(do.call('rbind', strsplit(as.character(dat),'\t', fixed=TRUE)),
                      stringsAsFactors=FALSE)
    # adjust data type per column
    dat <- convert.magic(dat, c(rep('numeric', 9), 'character', rep('numeric',3)))
    colnames(dat) <- c("LabelChannel", "MoleculeId", "molLength", "molAvgIntensity",
                       "labelSNR", "NumberofLabels", "OriginalMoleculeId", "ScanNumber", "ScanDirection",
                       "ChipId", "Flowcell", "RunId", "GlobalScanNumber")
    # add label dentity per 100k
    dat$label100kdens <- 100000*dat$NumberofLabels/dat$molLength
    # remove useless rows
    useless <- c("LabelChannel", "MoleculeId", "OriginalMoleculeId", "ScanNumber",
                 "ScanDirection", "ChipId", "Flowcell", "RunId", "GlobalScanNumber")
    dat <- dat[, .SD, .SDcols=-useless]
    Sys.sleep(1)
    progress$set(value = 4/n, detail = "creating indexes")
    # set keys for fast access
    key(dat)
    keycols = c("molLength", "molAvgIntensity", "labelSNR", "NumberofLabels")
    setkeyv(dat, keycols)   # rather than key(DT)<-keycols (which copies entire table)
    Sys.sleep(1)
    return(dat)
  })

  filter.data <- eventReactive(input$FilterButton, {
    if (is.null(load.data())) return(NULL)

    # set filters for length & AI & SNR from UI
    minlen <- min(input$length)
    maxlen <- max(input$length)
    minAI <- min(input$mAI)
    maxAI <- max(input$mAI)
    minSNR <- min(input$mSNR)
    maxSNR <- max(input$mSNR)
    minlab <- min(input$mlabels)
    maxlab <- max(input$mlabels)
    
    # 4x subset
    n <- 4
    progress <- shiny::Progress$new()
    progress$set(message = "Filtering by", value = 0)
    progress$set(value = 1/n, detail = "by molLength")
    on.exit(progress$close())
    subset <- load.data()[molLength %between% list(minlen*1000 ,maxlen*1000)]
    Sys.sleep(1)
    progress$set(value = 2/n, detail = "by molAvgIntensity")
    subset <- subset[molAvgIntensity %between% list(minAI, maxAI)]
    Sys.sleep(1)
    progress$set(value = 3/n, detail = "by labelSNR")
    subset <- subset[labelSNR %between% list(minSNR, maxSNR)]
    Sys.sleep(1)
    progress$set(value = 4/n, detail = "by NumberofLabels")
    subset <- subset[NumberofLabels %between% list(minlab, maxlab)]
    Sys.sleep(1)
    return(subset)
  })

  output$cntData <- renderText({
    if (is.null(filter.data())) return(NULL)
    paste("molecules in plot: ", nrow(filter.data()))
  })
  
  # output N50 and L50 values as separate text lines
  output$n50x <- renderText({
    if (is.null(sumData()[input$Xaxis,8])) return(NULL)
    paste("N50 for '",input$Xaxis, "': ", sumData()[input$Xaxis,8], " (L50=", sumData()[input$Xaxis,7], ") is shown with the grey dashed line.", sep="")
  })
  
  output$n50y <- renderText({
    if (input$Yaxis != 'histo') {
      paste("N50 for '",input$Yaxis, "': ", sumData()[input$Yaxis,8], " (L50=", sumData()[input$Yaxis,7], ") is shown with the grey dashed line.", sep="")
    }
  })
  sumData <- reactive({
    if (is.null(filter.data())) return(NULL)
    sum <- as.data.frame(
      do.call(cbind,
              lapply(filter.data(),
                     summary)))
    col.names <- c(row.names(sum), "L50", "N50")
    # add L50 and N50
    more <- as.data.frame(sapply(filter.data(), function(x) calc.n50(x)), stringsAsFactors=FALSE)
    more <- convert.magic(more, c(rep('numeric', 9), 'character', rep('numeric',4)))
    sum <- rbindlist(list(sum, more), fill = TRUE)
    sum <- t(sum)
    colnames(sum) <- col.names
    # round average values
    sum <- round_df(sum, 2)
    return(sum)
  })

  output$summaryzero = DT::renderDataTable({
    if (is.null(sumData())) return(NULL)
    sumData()
    }, options = list(dom = 't',
                      columnDefs=list(list(targets=1:8, class="dt-right")))
  )
  
  plotInput <- eventReactive(input$PlotButton, {
    if (is.null(filter.data())) return(NULL)

    # table of applied filters
    filterSet <- data.frame( min=c(min(input$length), min(input$mAI), min(input$mSNR), min(input$mlabels), "NA"), 
                         max=c(max(input$length), max(input$mAI), max(input$mSNR), max(input$mlabels), round(nrow(filter.data()),0))
                         )
    row.names(filterSet) <- c("molLength", "molAvgIntensity" , "labelSNR", "NumberofLabels", "TotalMolecules")
    
    tt <- gridExtra::ttheme_default(
      core = list(fg_params=list(cex = 0.75)),
      colhead = list(fg_params=list(cex = 0.75)),
      rowhead = list(fg_params=list(cex = 0.75)))
    
    # test plot type from X and Y
    if (input$Yaxis != 'histo') {
      # scatterplot
      p <- ggplot(data=filter.data(), aes_string(x = input$Xaxis, y = input$Yaxis))
      p <- p + geom_point(aes(colour = labelSNR, size=label100kdens)) +
        scale_colour_gradient(low="grey95", high = "grey15") +
        scale_size(range = c(0, 3))
      p <- p + geom_vline(xintercept = as.numeric(sumData()[input$Xaxis,8]), colour="gray50", linetype="dashed", size=0.5)
      p <- p + geom_hline(yintercept = as.numeric(sumData()[input$Yaxis,8]), colour="gray50", linetype="dashed", size=0.5)
      p
    } else {
      # density histogram from Xaxis data
      p <- ggplot(data=filter.data(), aes_string(x = input$Xaxis))
      p <- p + geom_density()
      p <- p + geom_vline(xintercept = as.numeric(sumData()[input$Xaxis,8]), colour="gray50", linetype="dashed", size=0.5)    
    }
    # add filters next to plot
    gs <- list(tableGrob(filterSet, theme=tt), p)
    hlay <- rbind(c(1,2,2,2),
                  c(NA,2,2,2))
    select_grobs <- function(lay) {
      id <- unique(c(t(lay))) 
      id[!is.na(id)]
    } 
    g <- grid.arrange(grobs=gs[select_grobs(hlay)], 
                      layout_matrix=hlay)
    return(g)
  })
  
  output$plot1 <- renderPlot({
    plot(plotInput(), width="640px", height="480px")
  })
  
  plotFile <- reactive({
    if (is.null(sumData())) return(NULL)
    if (is.null(plotInput())) return(NULL)
    plt <- plotInput()
    tbl <- as.data.frame(sumData())
    
    # print plot and table to file
    tt <- gridExtra::ttheme_default(
        core = list(fg_params=list(cex = 0.75)),
        colhead = list(fg_params=list(cex = 0.75)),
        rowhead = list(fg_params=list(cex = 0.75)))
    tbl <- tableGrob(tbl, theme=tt)
    empty <- textGrob("")

    # Plot chart and table into one object
    grid.arrange(plt, tbl,
                 nrow=2,
                 as.table=TRUE,
                 heights=c(4,1),
                 vp=viewport(width=0.95, height=0.95))
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste(input$outfile, input$format, sep=".") },
    content = function(file) {
      if(input$format == "png")
        png(file, width = 600, height = 480, units = "px") # open the png device
      else
        pdf(file, width = 10, height = 8) # open the pdf device
      plot(plotFile())
      dev.off()  # turn the device off
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)

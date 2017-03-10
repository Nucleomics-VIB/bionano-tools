# BNX_viewer.shinyapp
# A R/shiny tool to create a scatterplot from BNX data
# designed to work with BNX 1.2 format
#
# Stephane Plaisance, VIB Nucleomics Core
# visit our Git: https://github.com/Nucleomics-VIB 
# version: 2017-03-10_v1.0
# © by using this tool, you accept the licence saved under ./www/licence.pdf

library(shiny)
library("readr")
library("stringr")
library("ggplot2")
library("data.table")

# limit upload to 1GB
options(shiny.maxRequestSize=1000*1024^2) 
script.version="1.0"

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

convert.magic <- function(obj, types){
  for (i in 1:length(obj)){
    FUN <- switch(types[i],
                  character = as.character, 
                  numeric = as.numeric, 
                  factor = as.factor)
    obj[,i] <- FUN(obj[,i])
  }
  obj
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
      tags$a(href="license.pdf", target="_blank", "usage licence"),
      hr(),
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
                    "labels per 100k" = "label100kdens"),
                  selected = def.selin2
      ),
      
      hr(),
      tags$h4("modify ranges below and click ", tags$em("Calculate & Plot")),
      
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
      
      hr(),

      actionButton(inputId='goButton', "Calculate & Plot", style='padding:4px; font-weight: bold; font-size:150%')
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput('plot1'),
      textOutput('min1'),
      textOutput('max1'),
      textOutput('min2'),
      textOutput('max2'),
      textOutput('min3'),
      textOutput('max3'),
      textOutput('cntData'),
      tableOutput('summaryzero')
    )
  )
)

# load data and create scatterplot
server <- function(input, output) {

  output$min1 <- renderText({ 
    paste("min length set to: ", min(input$length))
  })
  
  output$max1 <- renderText({ 
    paste("max length set to: ", max(input$length))
  })
  
  output$min2 <- renderText({ 
    paste("min AvgIntensity set to: ", min(input$mAI))
  })
  
  output$max2 <- renderText({ 
    paste("max AvgIntensity set to: ", max(input$mAI))
  })
  
  output$min3 <- renderText({ 
    paste("min SNR set to: ", min(input$mSNR))
  })
  
  output$max3 <- renderText({ 
    paste("max SNR set to: ", max(input$mSNR))
  })
  
  load.data <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    dat <- readLines(inFile$datapath, skipNul=TRUE)
    dat <- dat[grepl("^0\t*", dat, perl = TRUE)]
    dat <- data.frame(do.call('rbind', strsplit(as.character(dat),'\t', fixed=TRUE)), 
                      stringsAsFactors=FALSE)
    # adjust data type per column
    dat <- convert.magic(dat, c(rep('numeric', 9), 'character', rep('numeric',3)))
    colnames(dat) <- c("LabelChannel", "MoleculeId", "molLength", "molAvgIntensity", 
                       "labelSNR", "NumberofLabels", "OriginalMoleculeId", "ScanNumber", "ScanDirection", 
                       "ChipId", "Flowcell", "RunId", "GlobalScanNumber")
    # add label dentity per 100k
    dat$label100kdens <- 100000*dat$NumberofLabels/dat$molLength
    dat
  })
  
  filter.data <- eventReactive(input$goButton, {
    if (is.null(load.data())) return(NULL)
    
    # set filters for length & AI & SNR from UI
    minlen <- min(input$length)
    maxlen <- max(input$length)
    minAI <- min(input$mAI)
    maxAI <- max(input$mAI)
    minSNR <- min(input$mSNR)
    maxSNR <- max(input$mSNR)
    
    # subset
    subset(load.data(), ( load.data()$molLength>=minlen*1000 & 
                            load.data()$molLength<=maxlen*1000 & 
                            load.data()$molAvgIntensity>=minAI & 
                            load.data()$molAvgIntensity<=maxAI &
                            load.data()$labelSNR>=minSNR & 
                            load.data()$labelSNR<=maxSNR) )
  })

  output$cntData <- reactive({
    if (is.null(filter.data())) return(NULL)
    paste("molecules in plot: ", nrow(filter.data()))
  })
  
  sumData <- reactive({
    if (is.null(filter.data())) return(NULL)
    t(as.data.frame(
      do.call(cbind, 
              lapply(filter.data(), 
                     summary))
    )
    )
  })
  
  output$summaryzero <- renderTable({ 
    sumData() 
  }, rownames=TRUE, colnames=TRUE, digits=2)
  
  output$plot1 <- renderPlot({
    if (is.null(filter.data())) return(NULL)
    
    # scatterplot
    p <- ggplot(data=filter.data(), aes_string(x = input$Xaxis, y = input$Yaxis))
    p <- p + geom_point(aes(colour = labelSNR, size=label100kdens)) +
      scale_colour_gradient(low="grey95", high = "grey15") +
      scale_size(range = c(0, 3))
    p
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

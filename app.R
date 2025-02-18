source("SingleChrom.R")
suppressMessages(library(xcms)) #will make chromatograms 
library(SummarizedExperiment)
suppressMessages(library(msPurity))
suppressMessages(library(Spectra)) #will make spectra
library(dplyr)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Mass Spec Analysis"),
  p("Before you load your MS data there will be an error message. This will go away once you load the MS data. It will take a little bit of time to load."),
  p("If you see 'error in evaluating the argument 'object' in selecting a method for function 'msLevel': subscript contains out-of-bounds indices' in the spectra section this means that there are no MS2 for the peak"),
  p("If you scroll down further there will be a table of the peaks in the chromatogram file. There are also tables of the chromatogram and spectra that can be downloaded."),  
  fluidRow(
      column(6,
      h2("MS Chromatogram"),
      plotOutput("Chrom"),
      downloadLink("downloadPlot", "Download Chromatogram"),
      ),
      column(6,
      h2("Spectra Plot"),
      plotOutput("SpectraPlot"),
      downloadLink("downloadSpPlot", "Download Spectra plot"),
      ),
    ),
    wellPanel(
    fluidRow(
      column(3,
      h4("Chromatogram Setup"),             
      helpText("Create chromatogram from LCMS data"),
      fileInput("file1", h5("LCMS File input"), accept = ".mzML"),
#      fileInput("file2", h3("Data File input"), accept = ".csv"),
      numericInput("num", 
                   label = "Target mass",
                   value = 315.2319), 
      sliderInput("range", 
                  label = "Retention time range:",
                  min = 0, max = 1000, value = c(0, 1000)),
      numericInput("yrange", 
                 label = "Y Max Value",
                 value = 100000),
      ),
      column(3,
      h4("Chromatogram Design"),
      helpText("Change the dimensions of the graph and font size"),
      numericInput("width", 
                   h5("Width pdf file"),
                   value = 8),
      numericInput("height", 
                   h5("Height pdf file"),
                   value = 4),
      numericInput("fontsize", 
                   h5("Font size of chromatogram"),
                   value = 15),
      ),
      column(2,
      h4("Spectra Input"),
       numericInput("peaknm", 
                   h5("Pick row in peak file"),
                   value = 1),
      numericInput("Spec", 
                   h5("Pick spectra in peak"),
                   value = 1),
      numericInput("label", 
                   h5("Peak label setting"),
                   value = 1000000)
      ),
      column(4,
      br(),
      br(),
      br(),
      helpText("Using the peak file displayed on the screen pick which row you would like to display the spectra of"),
      br(),
      helpText("Each row that was picked has multiple sepctra associated with it. Change the value below to change which specra is displayed."),
      br(),
      helpText("You can change the labeling of the specta file by changing the value below. Note: only intensities above the value will be displayed."),
      
    ),
    ),
    ),
    fluidRow(
      h3("Tables for download"),
    ),
    wellPanel(
      fluidRow(
      h4("Table of Peak Values"),
      downloadLink("downloadPeak", "Download Peak File"),
      tableOutput("Peaks"),
      br(),
    ),
    ),
    fluidRow(
      column(6,
      h4("Chromatogram Table"),
      p("Table of Chromatogram values intensity and retention time values for download."),             
      downloadLink("downloadFile", "Download Chromatogram csv"),
      tableOutput("ChromCSV")
      ),
      column(6,
      h4("Spectra Table"),
      p("Table spectra intensity values for Download."),             
      downloadLink("downloadSpCSV", "Download Spectra csv"),
      tableOutput("SpectraCsv")
      )
    )
  )



# Server logic ----
server <- function(input, output) {
  #Make it so we can upload larger files. Default is ~30 I think 
  options(shiny.maxRequestSize=99*1024^2)
  #Reactive elements that process all the big data one time rather than each time the inputs are changed
  dataInput <- reactive({
    file1 <- input$file1
    Exp<-readMSData(files = file1$datapath, mode = "onDisk")
  })
  dataChrom <- reactive({
    setup(dataInput(), input$num, input$range)
  })
  
  #make the files you want to be able to download as reactive elements 
  makeChrom <- reactive({
    peakfile(dataChrom())
  })
  makeFile <- reactive({
    chromfile(dataChrom())
  })
  makechromplot <- reactive({
    chromplot(makeFile(), input$fontsize, input$yrange)
  })
  
  setspectra <- reactive({
    file1<-input$file1
    getspectra(file1$datapath, makeChrom(), input$num, input$peaknm)
  })
 # spval<-reactiveValues(setspectra())
  getspectracsv <- reactive({
    spectrachromcsv(setspectra(), input$Spec)
  })
  getspectraplot <- reactive({
   spectrachromplot(setspectra(), input$Spec, input$label)
  })
  
  #Outputs that match with the UI section
  output$Chrom <- renderPlot({
    makechromplot()
  })
  output$downloadPlot <- downloadHandler(
    filename = function() { "Chromatogram.pdf" },
    content = function(file) {
      pdf(file, height = input$height, width = input$width)
      plot(makechromplot())
      dev.off()
    }
  )
  output$ChromCSV <- renderTable({
    makeFile()
  })
  output$downloadPeak <- downloadHandler(
    filename = function() { "ChromFile.csv" },
    content = function(file) {
      write.csv(makeFile(), file)
    }
  )
  output$Peaks <- renderTable({
    makeChrom()
  })
  output$downloadFile <- downloadHandler(
    filename = function() { "Peaks.csv" },
    content = function(file) {
      write.csv(makeChrom(), file)
    }
  )
  output$SpectraCsv <- renderTable({
    getspectracsv()
  })
  output$downloadSpCSV <- downloadHandler(
    filename = function() { "Spectra.csv" },
    content = function(file) {
      write.csv(getspectracsv(), file)
    }
  )
  output$SpectraPlot <- renderPlot({
    getspectraplot()
  })
  output$downloadSpPlot <- downloadHandler(
    filename = function() { "SpectraPlot.pdf" },
    content = function(file) {
      pdf(file)
      spectrachromplot(setspectra(), input$Spec, input$label)
      dev.off()

    }
  )
}

# Run app ----
shinyApp(ui, server)
library(shiny)
library(dplyr)
library(circlize)

# function to circlize
dp_circlize_blast <- function(test_blast1) {
  blast_header <- c("queryId", "subjectId", "percIdentity", "alnLength", "mismatchCount", "gapOpenCount", "queryStart", "queryEnd", "subjectStart", "subjectEnd", "eVal", "bitScore")
  colnames(test_blast1) <- blast_header
  length_query <- test_blast1 %>% group_by(queryId) %>% summarise(chrlength = max(queryEnd)) %>% rename(chr = queryId) %>% mutate(type = "query")
  length_subject <- test_blast1 %>% group_by(subjectId) %>% summarise(chrlength = max(subjectEnd)) %>% rename(chr = subjectId) %>% mutate(type = "subject")
  
  allchr_length <- rbind(length_query, length_subject)
  bed2 <- test_blast1 %>% select(queryId, queryStart, queryEnd) %>% rename(chr = queryId, start = queryStart, end = queryEnd) %>% mutate(color1 = as.numeric(as.factor(chr)))
  bed1 <- test_blast1 %>% select(subjectId, subjectStart, subjectEnd) %>% rename(chr = subjectId, start = subjectStart, end = subjectEnd) %>% mutate(color1 = as.numeric(as.factor(chr)))
  
  circosname = paste0("circos_", format(Sys.time(), "%Y-%b-%d-%H.%M.%S"), ".png")
  png(filename = circosname, width = 7, height = 7, units = "in", res = 1200)
  circos.clear()
  circos.par(cell.padding = c(0.02, 0, 0.02, 0))
  circos.initialize(factors = unique(allchr_length$chr), xlim = matrix(c(rep(0, length(unique(allchr_length$chr))), allchr_length$chrlength), ncol = 2))
  
  col_text <- "grey40"
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.text(mean(xlim), mean(ylim), chr, cex = 0.3, col = col_text, facing = "clockwise", niceFacing = TRUE)
  },
  bg.col = (as.numeric(as.factor(allchr_length$type)) + 2), bg.border = F, track.height = 0.06)
  
  circos.genomicLink(bed1, bed2, col = add_transparency(bed2$color1, 0.8), border = NA)
  dev.off()
  print(circosname)
  return(circosname)  # Return the name of the generated file
}

# Define UI
ui <- fluidPage(
  # Add CSS to change the background color of the page
  tags$style(HTML("
    body {
      background-color: #faf6fd;  
    }
    .well {
      background-color: #e1dde3;  
    }
  ")),
  titlePanel("BLAST Output Circos Plot Generator"),
  h4("(c) Dev Paudel, dpaudel@outlook.com"),
  p("This application allows you to upload BLAST output file in tabular format and generate a Circos plot."),
    p("Please ensure your file is in the correct format: [blastn -outfmt 6]"),
    p("Only top 200 rows will be evaluated."),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose BLAST Output File", accept = c(".txt", ".csv")),
      actionButton("loadSample", "Load Sample BLAST Output"),
      actionButton("generate", "Generate Circos Plot")
    ),
    mainPanel(
      textOutput("message"),
      imageOutput("circosPlot")  
    )
  )
)



# Define server logic
server <- function(input, output, session) {
  # Function to generate sample BLAST data directly in the app
  generate_sample_blast <- function() {
    sample_data <- data.frame(
      queryId = c("seq1", "seq1", "seq2", "seq3"),
      subjectId = c("seqA", "seqB", "seqC", "seqA"),
      percIdentity = c(99, 95, 85, 90),
      alnLength = c(100, 100, 150, 200),
      mismatchCount = c(1, 5, 15, 10),
      gapOpenCount = c(0, 1, 2, 0),
      queryStart = c(1, 1, 1, 1),
      queryEnd = c(100, 100, 150, 200),
      subjectStart = c(1, 10, 5, 1),
      subjectEnd = c(100, 110, 155, 200),
      eVal = c("1.00E-50", "1.00E-30", "1.00E-20", "1.00E-10"),
      bitScore = c(200, 150, 100, 80)
    )
    return(sample_data)
  }
  
  blast_data <- reactiveVal(NULL)  # Create a reactive value to store BLAST data
  
  observeEvent(input$loadSample, {
    # Directly load the sample BLAST output programmatically
    blast_data(generate_sample_blast())
    output$message <- renderText("Sample BLAST output loaded. Click Generate Circos Plot on sidebar")
  })
  
  observeEvent(input$generate, {
    # Use sample data if loaded, or require a file upload
    if (is.null(blast_data())) {
      req(input$file)  # If no sample is loaded, require file input
      # Read the uploaded file
      blast_data(read.table(input$file$datapath, header = TRUE, sep = "\t", nrows = 200))
    }
    
    # Generate the circos plot
    circos_file <- dp_circlize_blast(blast_data())
    
    # Construct the complete path to the output file
    output_path <- normalizePath(circos_file)
    
    output$message <- renderText({
      paste("Circos plot generated", "") # remove name for online
    })
    
    output$circosPlot <- renderImage({
      list(
        src = circos_file,
        contentType = 'image/png',
        alt = "Circos plot",
        width = 0.2 * 1600,  # Adjust based on original dimensions
        height = 0.2 * 1600   # Adjust based on original dimensions
      )
    }, deleteFile = FALSE)
  })
}



# Run the application
shinyApp(ui = ui, server = server)

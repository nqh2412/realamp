library(shiny)
library(shinyBS)
library(imager)
library(colocr)
library(dplyr)
library(qpcR)
source('modified_functions.R')
source('modified_functions_qpcr.R')
# Define UI for application that draws a histogram
ui <- navbarPage(
  
  # Application page Title
  title = 'ReaLAMP',

  
  ## Main
  tabPanel('Main',
           sidebarLayout(
             sidebarPanel(
               # Input Panel
               
               tags$h3('Input Panel'),
               
               ## Description
               
               ## load Image
               fileInput('image1', 'Upload Images', multiple = TRUE),
               bsTooltip('image1',
                         'Select and upload the chamber images.',
                         'right', options = list(container = "body")),
               numericInput('roi_num', 'Number of Chamber', 7, 1, 50, 1),
               bsTooltip('roi_num',
                         'Select the number of Sample on a chip. Default is 7',
                         'right', options = list(container = "body")),
               tags$hr(),
               #tags$p(' Crop images'),
               ## Selection Parameters
               #numericInput('x1', 'X', 5, min = 1, max = 10),
               #numericInput('x2', 'X 2', 1.3, min = 1, max = 10),
               #numericInput('y1', 'y', 1.9, min = 1, max = 10),
               #numericInput('y2', 'y 2', 5, min = 1, max = 10),
               #tags$hr(),
               tags$p('Adjust parameters'),
               sliderInput('threshold', 'Threshold', 1, 99, 29, 1),
               bsTooltip('threshold',
                         'Choose a threshold for excluding the background.',
                         'right', options = list(container = "body")),
               sliderInput('shrink', 'Shrink', 1, 10, 1, 1),
               bsTooltip('shrink',
                         'Shrink the selected area by eroding the bounderies around it.',
                         'right', options = list(container = "body")),
               sliderInput('grow', 'Grow', 1, 10, 1, 1),
               bsTooltip('grow',
                         'Grow the selected area by dilating the bounderies around it.',
                         'right', options = list(container = "body")),
               sliderInput('fill', 'Fill', 1, 10, 1, 1),
               bsTooltip('fill',
                         'Remove holes in the selected area.',
                         'right', options = list(container = "body")),
               sliderInput('clean', 'Clean', 1, 10, 1, 1),
               bsTooltip('clean',
                         'Remove small isolated parts in the selected area.',
                         'right', options = list(container = "body")),
               sliderInput('tolerance', 'Tolerance', 0, .99, .1, .1),
               bsTooltip('tolerance',
                         'Set value to determine which two neighboring pixels are in same selected area.',
                         'right', options = list(container = "body"))
               
               ),
             mainPanel(
               # Output Views
               # This part of the app contains the different views of the output
               # It is divided into description and four tabs. The four tabs are:
               #    1. Preview
               #    2. Graph View
               #    3. Data table
               
               fluidRow(
                 
                 ## Description
                 tags$h2('Colorimetric Real time LAMP'),
                 tags$br(),
                 tags$li('Step 1: Upload the images'),
                 tags$li('Step 2: Adjust parameters in the input panel. The auto-selected chambers show in Preview tab with red borders'),
                 tags$li('Step 3: The graphical results show in ReaLAMP analysis tab, the raw data show in Data table tab'),
                 tags$br(),
                 tags$br(),
                 
                 # Tabs
                 tabsetPanel(
                   
                   ## Select ROI
                   tabPanel('Preview',
                            plotOutput("image_plot", height = "6000px")
                           
                   ),
                   
                   ## RGB ratio
                   
                   tabPanel('ReaLAMP analysis',
                  #sliderInput('thresh', 'Threshold value', 0,1.5, 0.95, 0.05),
                            plotOutput('graph_plot', height = "1000")),
                  
            
                   
                  
                   
                   # Graph View
                   tabPanel('Data tables',
                            tableOutput('table'))
                 )
                 )
                 )
             )),
  
  # GitHub
  tabPanel('GitHub',
           "Comments, issues and contributions are welcomed")
  )


# Define server
server <- function(input, output) {
  # intiate interactive values
  values <- reactiveValues()
  
  # load images
  img1 <- reactive({
    image_load(input$image1$datapath)
             
  })
  
  ## calculate the pixset
  px <- reactive({
    roi_select(img1(),
               threshold = input$threshold,
               shrink = input$shrink,
               grow = input$grow,
               fill = input$fill,
               clean = input$clean,
               n = input$roi_num)
  })
  
  ## calculate correlations
  corr <- reactive({
    roi_test(px(), type = 'both')
  })
  b1 <- reactive({
    roi_check(px())  
  })
  # Output Views
  ## Select ROI
  
  # plots
  output$image_plot <- renderPlot({
    req(input$image1)
    
    n <- length(input$image1$name)
    
    par(mfrow=c(n,4), mar = rep(1, 4))
    roi_show(px())
  })
  
  
  ## Pixel Intensities
  output$graph_plot <- renderPlot({
    req(input$image1)
    par(mfrow=c(4,1))
    b1()
  })
  ## Pixel Intensities
  output$table <- renderTable({
    req(input$image1)
    b1()
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
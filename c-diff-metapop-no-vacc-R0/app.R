#
# App Title: R_0 for C. difficile metapopulation model (without vaccination)
#
# This is a Shiny application. You can run the application by clicking
# the 'Run App' button above.
#

library(shiny)

# Define UI for application
ui <- fluidPage(
  withMathJax(),
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(6,
               sliderInput(
                 inputId = "N1",
                 label = helpText("$$N_1$$"),
                 min = 50000, max = 150000, value = 100000, step = 1000,
                 ticks = FALSE
               )
        ),
        column(6,
               sliderInput(
                 inputId = "N2",
                 label = helpText("$$N_2$$"),
                 min = 500, max = 3500, value = 2000, step = 100,
                 ticks = FALSE
               )
        )
      )
    ),
    mainPanel(
      plotOutput("R0")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
}

# Run the application 
shinyApp(ui = ui, server = server)


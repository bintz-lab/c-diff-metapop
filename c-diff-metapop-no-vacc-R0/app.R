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
  fluidRow(
    column(
      width = 4,
      wellPanel(
        selectInput(
          inputId = "parameter_to_vary",
          label = "Parameter to vary continuously?",
          choices = c("v", "alpha", "beta"),
          selected = "v"
        )
      ),
      wellPanel(
        selectInput(
          inputId = "parameter_to_vary",
          label = "Parameter to vary discretely?",
          choices = c("v", "alpha", "beta"),
          selected = "v"
        )
      ),
      wellPanel(
        sliderInput(
          inputId = "param_val_1",
          label = "",
          min = 0, max = 0.6, value = 0, step = 0.1
        ),
        sliderInput(
          inputId = "param_val_2",
          label = "",
          min = 0, max = 0.6, value = 0.3, step = 0.1
        ),
        sliderInput(
          inputId = "param_val_3",
          label = "",
          min = 0, max = 0.6, value = 0.6, step = 0.1
        )
      )
    ),
    column(
      width = 8,
      h4(helpText("$$\\mathscr{R}_0\\,\\text{for } \\textit{C. difficile } \\text{metapopulation model without vaccination}$$")),
      plotOutput("states"),
      actionButton(
        inputId = "update",
        label = "Update"
      )
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


#
# App Title: C. difficile with constant vaccination
#
# This is a Shiny application. You can run the application by clicking
# the 'Run App' button above.
#

library(shiny)
library(tidyverse)
library(deSolve)

# Define UI for application
ui <- fluidPage(
  fluidRow(
    column(
      width = 4,
      wellPanel(
        p("Parameter Group A"),
        sliderInput(
          inputId = "Av", 
          label = "vaccination rate",
          min = 0, max = 0.6, value = 0, step = 0.1
        ),
        sliderInput(
          inputId = "Aalpha", 
          label = "antibiotic prescription rate",
          min = 0, max = 5, value = 0.5, step = 0.1
        ),
        selectInput(
          inputId = "Abeta", 
          label = "transmission coefficient",
          choices = c(0.000001, 0.00001, 0.0001, 0.001, 0.01), 
          selected = 0.000001
        )
      ),
      wellPanel(
        p("Parameter Group B"),
        sliderInput(
          inputId = "Bv", 
          label = "vaccination rate",
          min = 0, max = 0.6, value = 0.3, step = 0.1
        ),
        sliderInput(
          inputId = "Balpha", 
          label = "antibiotic prescription rate",
          min = 0, max = 5, value = 0.5, step = 0.1
        ),
        selectInput(
          inputId = "Bbeta", 
          label = "transmission coefficient",
          choices = c(0.000001, 0.00001, 0.0001, 0.001, 0.01), 
          selected = 0.000001
        )
      )
    ),
    column(
      width = 8,
      h3(em("C. difficile"), "with constant vaccination"),
      plotOutput("states"),
      actionButton(
        inputId = "update",
        label = "Update"
      )
    )
  )
)

# a few functions to run only once
cdiff_vacc <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),
       {
         deltaN <- kR*R + kS*S + kCn*Cn + kCp*Cp + kD*D + kV*V + kCv*Cv
         lambda <- betaC*(Cn+Cp+Cv) + betaD*D
         dR <- aR*deltaN + theta*S - (alpha+kR)*R
         dS <- aS*deltaN + alpha*R + p*eps*D - (theta+lambda+v+kS)*S
         dCn <- aCn*deltaN + (1-f)*lambda*S - (phi+v+kCn)*Cn
         dCp <- aCp*deltaN + f*lambda*S - (v+kCp)*Cp
         dD <- aD*deltaN + phi*Cn - (p*eps+kD)*D
         dV <- aV*deltaN + v*S - (lambda+kV)*V
         dCv <- aCv*deltaN + lambda*V + v*(Cn+Cp) - kCv*Cv
         list(c(dR, dS, dCn, dCp, dD, dV, dCv))
       })
}

fixed_parameters <- c(
  aR = 0.75, aS = 0.09, aCn = 0.06, aCp = 0.09, aD = 0.01, aV = 0.0, aCv = 0.0,
  kR = 0.33, kS = 0.15, kCn = 0.15, kCp = 0.15, kD = 0.068, kV = 0.15, kCv = 0.15,
  theta = 0.033, phi = 0.06, p = 0.8, eps = 0.10, f = 0.6
)
state <- c(R = 112, S = 14, Cn = 9, Cp = 13, D = 2, V = 0, Cv = 0)
t <- seq(0, 30, length = 1000)

# Define server logic
server <- function(input, output) {
  data <- eventReactive(input$update, {
    varied_parameters_A <- c(v = input$Av, alpha = input$Aalpha, betaC = as.numeric(input$Abeta), betaD = as.numeric(input$Abeta))
    varied_parameters_B <- c(v = input$Bv, alpha = input$Balpha, betaC = as.numeric(input$Bbeta), betaD = as.numeric(input$Bbeta))
    parameters_A <- c(fixed_parameters, varied_parameters_A)
    parameters_B <- c(fixed_parameters, varied_parameters_B)
    
    list(A = parameters_A, B = parameters_B) %>% 
      map(~ode(y = state, time = t, func = cdiff_vacc, parms = .) %>%
            unclass() %>% 
            as_tibble() %>% 
            mutate(
              sumC = Cn + Cp + Cv,
              sumV = V + Cv
            ) %>% 
            rename(days = time) %>% 
            gather(state, num, -days) %>% 
            mutate(state = factor(state, levels = c("R", "S", "Cn", "Cp", "sumC", "D", "V", "Cv", "sumV")))
      ) %>%
      bind_rows(.id = "parameter_group")
  }, ignoreNULL = FALSE)
  
  output$states <- renderPlot({
    data() %>% 
      filter(!state %in% c("Cp","V","Cv","deltaN")) %>% 
      ggplot(aes(days, num, col = parameter_group)) + 
      geom_line() +
      facet_wrap(~state, nrow = 3, scales = "free_y")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)


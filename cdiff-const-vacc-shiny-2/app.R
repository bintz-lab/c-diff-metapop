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
        selectInput(
          inputId = "parameter_to_vary",
          label = "Parameter to vary?",
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
         betaC <- beta
         betaD <- beta
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

get_new_par_val_vec <- function(x) {
  for (i in seq_along(x)) {
    names(x)[i] <- as.character(x[i])
  }
  x
}

get_new_param_vec <- function(x, y, p) {
  y[p] <- as.numeric(x)
  y
}

parameters <- c(
  aR = 0.75, aS = 0.09, aCn = 0.06, aCp = 0.09, aD = 0.01, aV = 0.0, aCv = 0.0,
  kR = 0.33, kS = 0.15, kCn = 0.15, kCp = 0.15, kD = 0.068, kV = 0.15, kCv = 0.15,
  theta = 0.033, phi = 0.06, p = 0.8, eps = 0.10, f = 0.6,
  v = 0.3, alpha = 0.5, beta = 0.000001
)
state <- c(R = 112, S = 14, Cn = 9, Cp = 13, D = 2, V = 0, Cv = 0)
t <- seq(0, 30, length = 1000)

# Define server logic
server <- function(input, output, session) {
  observeEvent(input$parameter_to_vary,
    {
      val <- input$parameter_to_vary
      if (val == "v") {
        parmin <- 0
        parmax <- 0.6
        parval <- c(0, 0.3, 0.6)
        parstep <- 0.1
      } else if (val == "alpha") {
        parmin <- 0
        parmax <- 5
        parval <- c(0.5, 2, 5)
        parstep <- 0.1
      } else {
        parmin <- 0.000001
        parmax <- 0.0001
        parval <- c(0.000001, 0.00001, 0.0001)
        parstep <- 0.000001
      }
      updateSliderInput(session,
        inputId = "param_val_1",
        label = paste(val, "value 1"),
        min = parmin, max = parmax, value = parval[1], step = parstep
      )
      updateSliderInput(session,
        inputId = "param_val_2",
        label = paste(val, "value 2"),
        min = parmin, max = parmax, value = parval[2], step = parstep
      )
      updateSliderInput(session,
        inputId = "param_val_3",
        label = paste(val, "value 3"),
        min = parmin, max = parmax, value = parval[3], step = parstep
      )
    }
  )
  
  state_plot <- eventReactive(input$update, {
    vary_par <- input$parameter_to_vary
    par_vals <- c(input$param_val_1, input$param_val_2, input$param_val_3)
    par_vals <- get_new_par_val_vec(par_vals)
    
    par_vals %>% 
      map(get_new_param_vec, parameters, vary_par) %>%
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
      bind_rows(.id = "par_vals") %>% 
      mutate(par_vals = as.factor(par_vals)) %>% 
      filter(!state %in% c("Cp","V","Cv","deltaN")) %>% 
      ggplot(aes(days, num, col = par_vals)) + 
      geom_line() +
      facet_wrap(~state, nrow = 3, scales = "free_y") +
      scale_color_discrete(name=vary_par)
  }, ignoreNULL = FALSE)
  
  output$states <- renderPlot({
    state_plot()
  })
}

# Run the application 
shinyApp(ui = ui, server = server)


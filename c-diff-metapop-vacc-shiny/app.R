#
# App Title: C. difficile metapopulation model with constant vaccination
#
# This is a Shiny application. You can run the application by clicking
# the 'Run App' button above.
#

library(shiny)
library(tidyverse)
library(deSolve)

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
            fluidRow(
                column(
                    width = 6,
                    wellPanel(
                        plotOutput("p1L")
                    )
                ),
                column(
                    width = 6,
                    wellPanel(
                        plotOutput("p2L")
                    )
                )
            ),
            fluidRow(
                column(
                    width = 6,
                    wellPanel(
                        plotOutput("p1H")
                    )
                ),
                column(
                    width = 6,
                    wellPanel(
                        plotOutput("p2H")
                    )
                )
            ),
            fluidRow(
                column(
                    width = 6,
                    wellPanel(
                        plotOutput("p1V0")
                    )
                ),
                column(
                    width = 6,
                    wellPanel(
                        plotOutput("p2V0")
                    )
                )
            ),
            fluidRow(
                column(
                    width = 6,
                    wellPanel(
                        plotOutput("p1V")
                    )
                ),
                column(
                    width = 6,
                    wellPanel(
                        plotOutput("p2V")
                    )
                )
            )
        )
    )
)

# a few functions to run only once
cdiff_metapop_vacc <- function(t, state, parameters) {
    with(as.list(c(state, parameters)),
         {
             sumC1 <- C1L+C1H+C1V0+C1V
             sumD1 <- D1L+D1H+D1V0
             
             lambda1L <- betaC1L*sumC1 + betaD1L*sumD1
             lambda1H <- betaC1H*sumC1 + betaD1H*sumD1
             lambda1V0 <- betaC1V0*sumC1 + betaD1V0*sumD1
             lambda1V <- betaC1V*sumC1 + betaD1V*sumD1
             
             sumC2 <- C2L+C2H+C2V0+C2V
             sumD2 <- D2L+D2H+D2V0
             
             lambda2L <- betaC2L*sumC2 + betaD2L*sumD2
             lambda2H <- betaC2H*sumC2 + betaD2H*sumD2
             lambda2V0 <- betaC2V0*sumC2 + betaD2V0*sumD2
             lambda2V <- betaC2V*sumC2 + betaD2V*sumD2
             
             dR1L <- b - (d1+alpha1L+(delta+d2)*N2/N1)*R1L + theta1L*S1L + xi1L*C1L + (1-etaRL*qRL)*delta*R2L
             dS1L <-   - (d1+theta1L+lambda1L+(delta+d2)*N2/N1)*S1L + alpha1L*R1L + p1L*eps1L*D1L + (1-etaSL*qSL)*delta*S2L
             dC1L <-   - (d1+xi1L+phi1L+(delta+d2)*N2/N1)*C1L + lambda1L*S1L + (1-etaCL*qCL)*delta*C2L
             dD1L <-   - (d1+p1L*eps1L+(delta+d2)*N2/N1)*D1L + phi1L*C1L + delta*D2L
             
             dR1H <-  - (d1+alpha1H+(delta+d2)*N2/N1)*R1H + theta1H*S1H + xi1H*C1H + (1-etaRH*qRH)*delta*R2H
             dS1H <-  - (d1+theta1H+lambda1H+(delta+d2)*N2/N1)*S1H + alpha1H*R1H + p1H*eps1H*D1H + (1-etaSH*qSH)*delta*S2H
             dC1H <-  - (d1+xi1H+phi1H+(delta+d2)*N2/N1)*C1H + lambda1H*S1H + (1-etaCH*qCH)*delta*C2H
             dD1H <-  - (d1+p1H*eps1H+(delta+d2)*N2/N1)*D1H + phi1H*C1H + delta*D2H
             
             dR1V0 <- - (d1+alpha1V0+(delta+d2)*N2/N1+nuR)*R1V0 + theta1V0*S1V0 + xi1V0*C1V0 + delta*(R2V0+etaRL*qRL*R2L+etaRH*qRH*R2H)
             dS1V0 <- - (d1+theta1V0+lambda1V0+(delta+d2)*N2/N1+nuS)*S1V0 + alpha1V0*R1V0 + p1V0*eps1V0*D1V0 + delta*(S2V0+etaSL*qSL*S2L+etaSH*qSH*S2H)
             dC1V0 <- - (d1+xi1V0+phi1V0+(delta+d2)*N2/N1+nuC)*C1V0 + lambda1V0*S1V0 + delta*(C2V0+etaCL*qCL*C2L+etaCH*qCH*C2H)
             dD1V0 <- - (d1+p1V0*eps1V0+(delta+d2)*N2/N1)*D1V0 + phi1V0*C1V0 + delta*D2V0
             
             dR1V <-  - (d1+alpha1V+(delta+d2)*N2/N1)*R1V + theta1V*S1V + xi1V*C1V + delta*R2V + nuR*R1V0
             dS1V <-  - (d1+theta1V+lambda1V+(delta+d2)*N2/N1)*S1V + alpha1V*R1V + delta*S2V + nuS*S1V0
             dC1V <-  - (d1+xi1V+(delta+d2)*N2/N1)*C1V + lambda1V*S1V + delta*C2V + nuC*C1V0
             
             dR2L <-   - (d2+alpha2L+delta)*R2L + theta2L*S2L + xi2L*C2L + (1-rR)*(delta+d2)*N2/N1*R1L
             dS2L <-   - (d2+theta2L+lambda2L+delta)*S2L + alpha2L*R2L + p2L*eps2L*D2L + (1-rS)*(delta+d2)*N2/N1*S1L
             dC2L <-   - (d2+xi2L+phi2L+delta)*C2L + lambda2L*S2L + (1-rC)*(delta+d2)*N2/N1*C1L
             dD2L <-   - (d2+p2L*eps2L+delta)*D2L + phi2L*C2L + (1-rD)*(delta+d2)*N2/N1*D1L
             
             dR2H <-  - (d2+alpha2H+delta)*R2H + theta2H*S2H + xi2H*C2H + (delta+d2)*N2/N1*(rR*R1L+R1H)
             dS2H <-  - (d2+theta2H+lambda2H+delta)*S2H + alpha2H*R2H + p2H*eps2H*D2H + (delta+d2)*N2/N1*(rS*S1L+S1H)
             dC2H <-  - (d2+xi2H+phi2H+delta)*C2H + lambda2H*S2H + (delta+d2)*N2/N1*(rC*C1L+C1H)
             dD2H <-  - (d2+p2H*eps2H+delta)*D2H + phi2H*C2H + (delta+d2)*N2/N1*(rD*D1L+D1H)
             
             dR2V0 <- - (d2+alpha2V0+delta+nuR)*R2V0 + theta2V0*S2V0 + xi2V0*C2V0 + (delta+d2)*N2/N1*R1V0
             dS2V0 <- - (d2+theta2V0+lambda2V0+delta+nuS)*S2V0 + alpha2V0*R2V0 + p2V0*eps2V0*D2V0 + (delta+d2)*N2/N1*S1V0
             dC2V0 <- - (d2+xi2V0+phi2V0+delta+nuC)*C2V0 + lambda2V0*S2V0 + (delta+d2)*N2/N1*C1V0
             dD2V0 <- - (d2+p2V0*eps2V0+delta)*D2V0 + phi2V0*C2V0 + (delta+d2)*N2/N1*D1V0
             
             dR2V <-  - (d2+alpha2V+delta)*R2V + theta2V*S2V + xi2V*C2V + nuR*R2V0 + (delta+d2)*N2/N1*R1V
             dS2V <-  - (d2+theta2V+lambda2V+delta)*S2V + alpha2V*R2V + nuS*S2V0 + (delta+d2)*N2/N1*S1V
             dC2V <-  - (d2+xi2V+delta)*C2V + lambda2V*S2V + nuC*C2V0 + (delta+d2)*N2/N1*C1V
             
             list(c(dR1L, dS1L, dC1L, dD1L, dR1H, dS1H, dC1H, dD1H, dR1V0, dS1V0, dC1V0, dD1V0, dR1V, dS1V, dC1V, dR2L, dS2L, dC2L, dD2L, dR2H, dS2H, dC2H, dD2H, dR2V0, dS2V0, dC2V0, dD2V0, dR2V, dS2V, dC2V))
         })
}

parameters <- c(N1 = 100000, N2 = 2000, d1 = 1/(78.5*365), d2 = 0.0068,
                
                alpha1L = 0.5/50, theta1L = 0.033, xi1L = 0.0165, phi1L = 0.2, 
                p1L = 0.8, eps1L = 0.1, betaC1L = 0.007/50, betaD1L = 0.007/50,
                
                alpha1H = 0.5/50, theta1H = 0.033, xi1H = 0.0165, phi1H = 0.2, 
                p1H = 0.8, eps1H = 0.1, betaC1H = 0.007/50, betaD1H = 0.007/50,
                
                alpha1V0 = 0.5/50, theta1V0 = 0.033, xi1V0 = 0.0165, phi1V0 = 0.2, 
                p1V0 = 0.8, eps1V0 = 0.1, betaC1V0 = 0.007/50, betaD1V0 = 0.007/50,
                
                alpha1V = 0.5/50, theta1V = 0.033, xi1V = 0.0165, 
                betaC1V = 0.007/50, betaD1V = 0.007/50,
                
                alpha2L = 0.5, theta2L = 0.033, xi2L = 0.0165, phi2L = 0.2, 
                p2L = 0.8, eps2L = 0.1, betaC2L = 0.007, betaD2L = 0.007,
                
                alpha2H = 0.5, theta2H = 0.033, xi2H = 0.0165, phi2H = 0.2, 
                p2H = 0.8, eps2H = 0.1, betaC2H = 0.007, betaD2H = 0.007,
                
                alpha2V0 = 0.5, theta2V0 = 0.033, xi2V0 = 0.0165, phi2V0 = 0.2, 
                p2V0 = 0.8, eps2V0 = 0.1, betaC2V0 = 0.007, betaD2V0 = 0.007,
                
                alpha2V = 0.5, theta2V = 0.033, xi2V = 0.0165, 
                betaC2V = 0.007, betaD2V = 0.007,
                
                etaRL = 1, qRL = 0.25, etaSL = 1, qSL = 0.25, etaCL = 1, qCL = 0.25,
                etaRH = 1, qRH = 0.25, etaSH = 1, qSH = 0.25, etaCH = 1, qCH = 0.25,
                rR = 0.25, rS = 0.25, rC = 0.25, rD = 0.25,
                nuR = 2/365, nuS = 2/365, nuC = 2/365,
                delta = 0.135)

N1 <- as.numeric(parameters["N1"])
N2 <- as.numeric(parameters["N2"])
b <- as.numeric(parameters["d1"]*N1 + parameters["d2"]*N2)
parameters <- c(parameters, b = b)

state <- c(R1L = 0.63*N1, S1L = 0.18*N1, C1L = 0.063*N1, D1L = 0.027*N1,
           R1H = 0.07*N1, S1H = 0.02*N1, C1H = 0.007*N1, D1H = 0.003*N1,
           R1V0 = 0, S1V0 = 0, C1V0 = 0, D1V0 = 0,
           R1V = 0, S1V = 0, C1V = 0,
           
           R2L = 0.63*N2, S2L = 0.18*N2, C2L = 0.063*N2, D2L = 0.027*N2,
           R2H = 0.07*N2, S2H = 0.02*N2, C2H = 0.007*N2, D2H = 0.003*N2,
           R2V0 = 0, S2V0 = 0, C2V0 = 0, D2V0 = 0,
           R2V = 0, S2V = 0, C2V = 0)

t <- seq(0, 1100, length = 1000)

cdiff_df <- ode(
    y = state, time = t, 
    func = cdiff_metapop_vacc, 
    parms = parameters) %>%
    unclass() %>% 
    as_tibble() %>%
    gather(func, val, -time) %>% 
    separate(func, into = c("state", "population", "risk"), sep = c(1,2)) %>% 
    mutate(state = factor(state, levels = c("R", "S", "C", "D")),
           population = factor(population, levels = c("1", "2")),
           population = fct_recode(population, catchment = "1", hospitalized = "2"),
           risk = factor(risk, levels = c("L", "H", "V0", "V")))

patchwise_plot <- function(df) {
    df %>% 
        ggplot(aes(time, val)) + 
        geom_line() +
        facet_wrap(vars(state), scales = "free", nrow = 2) +
        theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.5))
}

plotlist <- cdiff_df %>% 
    group_by(population, risk) %>% 
    nest() %>% 
    mutate(
        patchplot = map(data, patchwise_plot) %>% 
            set_names(paste0(population, risk))
    ) %>% select(patchplot)

# Define server logic
server <- function(input, output, session) {
    new_cdiff_df <- eventReactive(input$update, {
        
    }
        
    )
    #TODO: these need to be redone using map and the titles should somehow be static
    output$p1L <- renderPlot({
        plotlist$patchplot$catchmentL +
            labs(title = "Community: Low Risk")
    })
    
    output$p2L <- renderPlot({
        plotlist$patchplot$hospitalizedL +
            labs(title = "Hospital: Low Risk")
    })
    
    output$p1H <- renderPlot({
        plotlist$patchplot$catchmentH +
            labs(title = "Community: High Risk")
    })
    
    output$p2H <- renderPlot({
        plotlist$patchplot$hospitalizedH +
            labs(title = "Hospital: High Risk")
    })
    
    output$p1V0 <- renderPlot({
        plotlist$patchplot$catchmentV0 +
            labs(title = "Community: Vaccination Initialized")
    })
    
    output$p2V0 <- renderPlot({
        plotlist$patchplot$hospitalizedV0 +
            labs(title = "Hospital: Vaccination Initialized")
    })
    
    output$p1V <- renderPlot({
        plotlist$patchplot$catchmentV +
            labs(title = "Community: Vaccinated")
    })
    
    output$p2V <- renderPlot({
        plotlist$patchplot$hospitalizedV +
            labs(title = "Hospital: Vaccinated")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)


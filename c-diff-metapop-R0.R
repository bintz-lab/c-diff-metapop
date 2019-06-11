R_0_cdiff_metapop <- function(state, parameters) {
  with(as.list(c(state, parameters)),
       {
         BC <- matrix(
           c(betaC1L*S1L, 0, betaC1L*S1L, 0,
             0, betaC2L*S2L, 0, betaC2L*S2L,
             betaC1H*S1H, 0, betaC1H*S1H, 0,
             0, betaC2H*S2H, 0, betaC2H*S2H),
           nrow = 4, byrow = TRUE
         )
         
         BD <- matrix(
           c(betaD1L*S1L, 0, betaD1L*S1L, 0,
             0, betaD2L*S2L, 0, betaD2L*S2L,
             betaD1H*S1H, 0, betaD1H*S1H, 0,
             0, betaD2H*S2H, 0, betaD2H*S2H),
           nrow = 4, byrow = TRUE
         )
         
         TC1L <- d1+xi1L+phi1L+(delta+d2)*N2/N1
         TC2L <- d2+xi2L+phi2L+delta
         TC1H <- d1+xi1H+phi1H+(delta+d2)*N2/N1+rhoC
         TC2H <- d2+xi2H+phi2H+delta
         
         TD1L <- d1+p1L*eps1L+(delta+d2)*N2/N1
         TD2L <- d2+p2L*eps2L+delta
         TD1H <- d1+p1H*eps1H+(delta+d2)*N2/N1+rhoD
         TD2H <- d2+p2H*eps2H+delta
         
         TC <- matrix(
           c(TC1L, -delta, -rhoC, 0,
             -(1-rhoC)*(delta+d2)*N2/N1, TC2L, 0, 0,
             0, 0, TC1H, -delta,
             -rhoC*(delta+d2)*N2/N1, 0, -(delta+d2)*N2/N1, TC2H),
           nrow = 4, byrow = TRUE
         )
         
         TD <- matrix(
           c(TD1L, -delta, -rhoD, 0,
             -(1-rhoD)*(delta+d2)*N2/N1, TD2L, 0, 0,
             0, 0, TD1H, -delta,
             -rhoD*(delta+d2)*N2/N1, 0, -(delta+d2)*N2/N1, TD2H),
           nrow = 4, byrow = TRUE
         )
         
         PHI <- diag(c(phi1L, phi2L, phi1H, phi2H), 4, 4)
         eigen((BC + BD%*%solve(TD)%*%PHI)%*%solve(TC), only.values = TRUE)
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
                rhoR = 0.001, rhoS = 0.001, rhoC = 0.001, rhoD = 0.001,
                delta = 0.135)

state <- c(R1L = 0.63*N1, S1L = 0.18*N1, C1L = 0.063*N1, D1L = 0.027*N1,
           R1H = 0.07*N1, S1H = 0.02*N1, C1H = 0.007*N1, D1H = 0.003*N1,
           R1V0 = 0, S1V0 = 0, C1V0 = 0, D1V0 = 0,
           R1V = 0, S1V = 0, C1V = 0,
           
           R2L = 0.63*N2, S2L = 0.18*N2, C2L = 0.063*N2, D2L = 0.027*N2,
           R2H = 0.07*N2, S2H = 0.02*N2, C2H = 0.007*N2, D2H = 0.003*N2,
           R2V0 = 0, S2V0 = 0, C2V0 = 0, D2V0 = 0,
           R2V = 0, S2V = 0, C2V = 0)
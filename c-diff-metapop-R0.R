R_0_cdiff_metapop <- function(S0, parameters) {
  with(as.list(c(S0, parameters)),
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
         FVinv <- (BC + BD %*% solve(TD) %*% PHI) %*% solve(TC)
         spectral_radius(FVinv)
       })
}

spectral_radius <- function(X) {
  eigs <- eigen(X, only.values = TRUE)$values
  eigs[which.max(abs(eigs))]
}

parameters <- c(N1 = 100000, N2 = 2000, d1 = 1/(78.5*365), d2 = 0.0068,
                
                xi1L = 0.0165, phi1L = 0.2, 
                p1L = 0.8, eps1L = 0.1, betaC1L = 0.007/50, betaD1L = 0.007/50,
                
                xi1H = 0.0165, phi1H = 0.2, 
                p1H = 0.8, eps1H = 0.1, betaC1H = 0.007/50, betaD1H = 0.007/50,
                
                xi2L = 0.0165, phi2L = 0.2, 
                p2L = 0.8, eps2L = 0.1, betaC2L = 0.007, betaD2L = 0.007,
                
                xi2H = 0.0165, phi2H = 0.2, 
                p2H = 0.8, eps2H = 0.1, betaC2H = 0.007, betaD2H = 0.007,
                
                rR = 0.25, rS = 0.25, rC = 0.25, rD = 0.25,
                rhoR = 0.001, rhoS = 0.001, rhoC = 0.001, rhoD = 0.001,
                delta = 0.135)

N1 <- as.numeric(parameters["N1"])
N2 <- as.numeric(parameters["N2"])
S0 <- c(S1L = 0.18*N1, S1H = 0.02*N1, S2L = 0.18*N2, S2H = 0.02*N2)

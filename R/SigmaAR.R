#' @title Simulate Covariance Matrix with an Auto-regressive (AR) Model
#'
#' @description Helper function to simulate Simulate Gaussian Data with an Autoregressive (AR) Model
#'
#' @param p Number of variables (rows x columns).
#' @param phi Auto-regression coefficient.
#'
#' @return A list with the values:
#' \item{\code{Sigma}}{A covariance matrix.}
#'
#' @examples
#' library(FastStepGraph)
#' library(MASS) # mvrnorm
#' phi = 0.4
#' p = 100  # Dimension
#' n = 100 # Sample size
#' Sigma = SigmaAR(p, phi) # Simulate Gaussian Data with Autoregressive (AR) Model
#' Omega = solve(Sigma)
#' Omega[abs(Omega) < 1e-5] = 0
#' X = list() # Generate Data from a Gaussian distribution
#' X = mvrnorm(n, mu=rep(0,p), Sigma)
#' X = scale(X)
#'
#' @author Prof. Juan G. Colonna, PhD. \email{juancolonna@icomp.ufam.edu.br}
#' @author Prof. Marcelo Ruiz, PhD. \email{mruiz@exa.unrc.edu.ar}
#'
#' @export
SigmaAR = function(p, phi){
    Sigma = diag(p)
    for (i in 1:p) {
        for (j in 1:p) {
            Sigma[i,j] = phi^(abs(i-j))
        }
    }
    return(Sigma)
}


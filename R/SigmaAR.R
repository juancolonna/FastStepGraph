#' @title Simulate Covariance Matrix with an Auto-regressive (AR) Model
#'
#' @description Helper function to simulate Simulate Gaussian Data with an Autoregressive (AR) Model
#'
#' @param n_rows Number of samples (rows of X).
#' @param p_columns Number of variables (columns of X).
#' @param phi Auto-regression coefficient.
#'
#' @return A list with the values:
#' \item{\code{Sigma}}{A covariance matrix.}
#' \item{\code{Omega}}{A precision matrix.}
#' \item{\code{X}}{A normalized data matrix with Gaussian distribution.}
#'
#' @importFrom MASS mvrnorm
#' 
#' @author Prof. Juan G. Colonna, PhD. \email{juancolonna@icomp.ufam.edu.br}
#' @author Prof. Marcelo Ruiz, PhD. \email{mruiz@exa.unrc.edu.ar}
#'
#' @export
#' 
#' @examples
#' data <- FastStepGraph::SigmaAR(30, 50, 0.4) # Simulate Gaussian Data
SigmaAR <- function(n_rows, p_columns, phi){
    Sigma <- diag(p_columns)
    for (i in seq_len(p_columns)) {
        for (j in seq_len(p_columns)) {
            Sigma[i,j] <- phi^(abs(i-j))
        }
    }
    
    Omega <- solve(Sigma)
    Omega[abs(Omega) < 1e-5] <- 0
    X <- list() # Generate Data from a Gaussian distribution
    X <- MASS::mvrnorm(n_rows, mu=rep(0, p_columns), Sigma)
    X <- scale(X) # data normalization to zero mean and unit variance 
    return(list(Sigma=Sigma, Omega=Omega, X=X))
}


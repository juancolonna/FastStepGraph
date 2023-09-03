#' @title Simulate Covariance Matrix with an Auto-regressive (AR) Model
#'
#' @description Helper function to simulate synthetic data
#'
#' @param p Number of variables (rows x columns).
#' @param phi Auto-regression coefficient.
#'
#' @return A list with the values:
#' \item{\code{Sigma}}{A covariance matrix.}
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


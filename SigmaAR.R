# Simulate Gaussian Data with Autoregressive (AR) Model
SigmaAR = function(p, phi){
    Sigma = diag(p)
    for (i in 1:p) {
        for (j in 1:p) {
            Sigma[i,j] = phi^(abs(i-j))
        }
    }
    return(Sigma)
}


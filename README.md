# FastStepGraph
FastStepGraph

Example of use in RStudio

```{r}
setwd(getwd())
library(MASS)     # mvrnorm
source('cv.FastStepGraph.R')
source('FastStepGraph.R')
source('SigmaAR.R')

# Simulate Gaussian Data with Autoregressive (AR) Model
set.seed(1234567)
phi = 0.4 
p = 200  # Dimension
n = 100 # Sample size

Sigma = SigmaAR(p, phi)
Omega = solve(Sigma)  
Omega[abs(Omega) < 1e-5] = 0  

# Generate Data from a Gaussian distribution 
X = list()
X = mvrnorm(n, mu=rep(0,p), Sigma)
X = scale(X)

alpha_f = 0.22
alpha_b = 0.10
nei.max = 5

t0 <- Sys.time() # INITIAL TIME
G = FastStepGraph(X, alpha_f = alpha_f, alpha_b = alpha_b, nei.max=nei.max)
difftime(Sys.time(), t0, units = "secs")

# G$Omega
```



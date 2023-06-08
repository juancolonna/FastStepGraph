# FastStepGraph

Example of use in RStudio. 

First, we should generate some sample data as follows:
```{r}
setwd(getwd())
library(MASS)     # mvrnorm
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
```

Afterwards, fit the Omega matrix by calling the Fast Step Graph function, like:

```{r}
source('FastStepGraph.R')

alpha_f = 0.22
alpha_b = 0.10
nei.max = 5

t0 <- Sys.time() # INITIAL TIME
G = FastStepGraph(X, alpha_f = alpha_f, alpha_b = alpha_b, nei.max=nei.max)
difftime(Sys.time(), t0, units = "secs")

# G$Omega
```

To find the optimal $\mathbf{\alpha_f}$ and $\mathbf{\alpha_b}$ parameters for the previously generated **X** data, we can perform a cross-validation on a combination grid as follows:

```{r}
source('cv.FastStepGraph.R')

alpha_f_min = 0.01
alpha_f_max = 0.7
n_alpha = 20 # size of the grid search
alpha_b_min = alpha_f_min
alpha_b_max = alpha_f_max

t0 <- Sys.time() # INITIAL TIME
res = cv.FastStepGraph(X, 
                       n_folds = 5, 
                       alpha_f_min = alpha_f_min, 
                       alpha_f_max = alpha_f_max,
                       n_alpha = n_alpha, 
                       nei.max = 5)

G = FastStepGraph(X, res$alpha_f_opt , res$alpha_b_opt, nei.max=nei.max)
difftime(Sys.time(), t0, units = "secs")
```
However, this is not an exhaustive grid search. This is a heuristic that always sets $\mathbf{\alpha_b}$ = $\frac{\mathbf{\alpha_f}}{2}$.

To increase performace you can execute the same cross-validation in parallel, but the **doParallel** dependency must be installed:

```{r}
install.packages("doParallel")
```

The **doParallel** package will require the following dependencies: **foreach**, **iterators** and **parallel**. Make sure you satisfy them. Then call the method setting the parameter **parallel = TRUE**, as follows:

```{r}
source('cv.FastStepGraph.R')

alpha_f_min = 0.01
alpha_f_max = 0.7
n_alpha = 20 # size of the grid search
alpha_b_min = alpha_f_min
alpha_b_max = alpha_f_max

t0 <- Sys.time() # INITIAL TIME
res = cv.FastStepGraph(X, 
                       n_folds = 5, 
                       alpha_f_min = alpha_f_min, 
                       alpha_f_max = alpha_f_max,
                       n_alpha = n_alpha, 
                       nei.max = 5,
                       parallel = TRUE)

G = FastStepGraph(X, res$alpha_f_opt , res$alpha_b_opt, nei.max=nei.max)
difftime(Sys.time(), t0, units = "secs")
```

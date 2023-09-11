To install this package run `devtools::install_github("juancolonna/FastStepGraph")`


First, we should generate some sample data as follows:

```{r}
setwd(getwd())
library(MASS)     # mvrnorm
library(FastStepGraph)

# In case you clone the github repo
# source('FastStepGraph.R')
# source('SigmaAR.R')

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
X = MASS::mvrnorm(n, mu=rep(0,p), Sigma)
X = scale(X)
```

Afterwards, fit the Omega matrix by calling the Fast Step Graph function, like:

```{r}
source('FastStepGraph.R')

t0 <- Sys.time() # INITIAL TIME
G = FastStepGraph(X, alpha_f = 0.22, alpha_b = 0.10)
difftime(Sys.time(), t0, units = "secs")

print(G$Omega)
```

If the `nei.max` argument is omitted, it will be 5. To find the optimal $\mathbf{\alpha_f}$ and $\mathbf{\alpha_b}$ parameters for the previously generated **X** data, we can perform a cross-validation on a combination grid as follows:

```{r}
t0 <- Sys.time() # INITIAL TIME
res = cv.FastStepGraph(X, data_shuffle = TRUE)
difftime(Sys.time(), t0, units = "secs")

print(res$alpha_f_opt)
print(res$alpha_b_opt)
print(res$Omega)
```

The arguments `n_folds = 5`, `alpha_f_min = 0.1`, `alpha_f_max = 0.9`, `n_alpha = 32` (size of the grid search) and `nei.max = 5`, have defaults values and can be omitted. Note that, `cv.FastStepGraph(X)` is not an exhaustive grid search over $\mathbf{\alpha_f}$ and $\mathbf{\alpha_b}$. This is a heuristic that always sets $\mathbf{\alpha_b}$ = $\frac{\mathbf{\alpha_f}}{2}$.

To increase time performance, you can run `cv.FastStepGraph(X)` in parallel. Next, you'll need to install and register a parallel backend. Since Windows does not support forking, the same backend that works in a Linux or OS X environment will not work in Windows. To run on a Linux system the **doParallel** dependency must be installed `install.packages("doParallel")`. To run on a Windows system, **doSNOW** is used `install.packages("doSNOW")`. These parallel packages will also require the following dependencies: **foreach**, **iterators** and **parallel**. Make sure you satisfy them. Then call the method setting the parameter **parallel = TRUE**, as follows:

```{r}
t0 <- Sys.time() # INITIAL TIME
res = cv.FastStepGraph(X, , data_shuffle = TRUE, parallel = TRUE)
difftime(Sys.time(), t0, units = "secs")

print(res$alpha_f_opt)
print(res$alpha_b_opt)
print(res$Omega)
```

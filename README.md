# Fast Step Graph

Fast Step Graph is an optimized implementation in **R** of the Stepwise approach used for discovering high-dimensional Gaussian graphical models. It aims to accurately estimate the $\mathbf{\Omega}$ precision matrix when dealing with datasets where the number of features is significantly larger than the number of samples ($p >> n$), such as in genomics.

This implementation builds upon the original code available at [link](https://jdssv.org/index.php/jdssv/article/view/11), which accompanies the associated research paper. **Fast Step Graph** enhances the computational efficiency of the original code to handle much larger graph structures than previously reported, while reducing the training time. Several improvements have been made, including the elimination of redundant code, utilization of column-wise data structures (better for R), avoid list creation, manipulation and expansion within loops, and the integration of a faster subroutine for regression. Additionally, this implementation addresses a bug introduced in the original code.

Despite these enhancements, the primary bottleneck of Fast Step Graph lies in the requirement of substantial memory resources ($\frac{p(p-1)}{2}$) for storing the entire graph, particularly when $p$ grows.

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

t0 <- Sys.time() # INITIAL TIME
G = FastStepGraph(X, alpha_f = 0.22, alpha_b = 0.10, nei.max=5)
difftime(Sys.time(), t0, units = "secs")

print(G$Omega)
```

To find the optimal $\mathbf{\alpha_f}$ and $\mathbf{\alpha_b}$ parameters for the previously generated **X** data, we can perform a cross-validation on a combination grid as follows:

```{r}
source('cv.FastStepGraph.R')

t0 <- Sys.time() # INITIAL TIME
res = cv.FastStepGraph(X, 
                       n_folds = 5, 
                       alpha_f_min = 0.01, 
                       alpha_f_max = 0.70,
                       n_alpha = 20, # size of the grid search
                       nei.max = 5)
difftime(Sys.time(), t0, units = "secs")

print(res$alpha_f_opt)
print(res$alpha_b_opt)
```
However, this is not an exhaustive grid search. This is a heuristic that always sets $\mathbf{\alpha_b}$ = $\frac{\mathbf{\alpha_f}}{2}$.

To increase performace you can execute the same cross-validation in parallel, but the **doParallel** dependency must be installed:

```{r}
install.packages("doParallel")
```

The **doParallel** package will require the following dependencies: **foreach**, **iterators** and **parallel**. Make sure you satisfy them. Then call the method setting the parameter **parallel = TRUE**, as follows:

```{r}
t0 <- Sys.time() # INITIAL TIME
res = cv.FastStepGraph(X, 
                       n_folds = 5, 
                       alpha_f_min = 0.01, 
                       alpha_f_max = 0.70,
                       n_alpha = 20, # size of the grid search
                       nei.max = 5,
                       parallel = TRUE)
difftime(Sys.time(), t0, units = "secs")

print(res$alpha_f_opt)
print(res$alpha_b_opt)
```

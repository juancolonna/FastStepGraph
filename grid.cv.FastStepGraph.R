#'
#' @title Searches for the optimal combination of alpha_f and alpha_b parameters using Cross-Validation
#'
#' @description \code{cv.FastStepGraph} implements the cross-valiation for the Fast Step Graph algorithm.
#'
#' @param x Data matrix (of size n x p).
#' @param n_folds Number of folds for the cross-validation procedure (default value 5).
#' @param alpha_f_min Minimum alpha_f threshold value for the cross-validation procedure (default value 0.1).
#' @param alpha_f_max Maximum alpha_f threshold value for the cross-validation procedure (default value 0.9).
#' @param n_alpha Number of elements in the grid for the cross-validation (default value 32).
#' @param nei.max Maximum number of variables in every neighborhood (default value 5).
#' @param data_scale Boolean parameter (TRUE or FALSE), when to scale data to zero mean and unit variance (default FALSE).
#' @param return_model Boolean parameter (TRUE or FALSE), when to return the output of FastStepGraph function (default TRUE).
#' @param parallel Boolean parameter (TRUE or FALSE), when to run Cross-Validation in parallel using a multicore architecture (default FALSE).
#'
#' @return A list with the values: \cr \cr
#' \code{alpha_f_opt}: the optimal alpha_f value, \cr \cr
#' \code{alpha_b_opt}: the optimal alpha_f value, \cr \cr
#' \code{CV.loss}: minimum loss, \cr \cr
#' \code{loss_grid}: a loss grid, \cr \cr
#'
#' @export
#'
#' @author Prof. Juan G. Colonna, PhD. \email{juancolonna@icomp.ufam.edu.br}
#'

grid.cv.FastStepGraph = function(x, n_folds=5, alpha_f_min=0.1, alpha_f_max=0.9, n_alpha=32, nei.max = 5, data_scale = FALSE, return_model = TRUE, parallel = FALSE){
    if (n_folds <= 1) { stop('Number of folds must be equal or larger than 2') }
    if (nei.max >= n) { stop('The maximum number of neighbors (nei.max) must be less than n-1.') }
    if (nei.max == 0) { stop('The minimum number of neighbors (nei.max) must be greater than 0.') }
    if (data_scale) { x = scale(x) }

    n = nrow(x)
    p = ncol(x)
 
    # set.seed(123)
    ntest = floor(n/n_folds)
    ntrain = n - ntest
    ind = sample(n)
    alpha_f = seq(alpha_f_min, alpha_f_max, length=n_alpha)
    alpha_b = alpha_f
    loss_grid = matrix(NA, nrow = n_alpha, ncol = n_alpha)

    old_loss = Inf
    new_loss = NA
    alpha_f_opt = NA
    alpha_b_opt = NA
    
    if (parallel) { 
      if(.Platform$OS.type == "unix") {
        library(doParallel) # needed to parallelize on Linux
        cores = detectCores()
        cl = makeCluster(cores[1]-1, type = "FORK") # not to overload your computer "FORK"
        registerDoParallel(cl)
      } 
      else {
        library(doSNOW) # needed to parallelize on Windows
        cores = detectCores()
        cl = makeCluster(cores[1]-1, type="SOCK")
        registerDoSNOW(cl)
      }
      
      # alpha_f_losses <- foreach(f = alpha_f, .combine = 'c', .inorder=TRUE) %dopar% {
      losses_and_alphas <- foreach(i = 1:n_alpha, .combine = 'c', .inorder=TRUE) %dopar% {
      source('FastStepGraph.R')

      # for (i in 1:n_alpha) {
          for (j in 1:n_alpha) {
              if (alpha_f[i] >= alpha_b[j]) {
                  loss = 0
                  for (k in 1:n_folds) {
                      sel = ((k-1)*ntest+1):(k*ntest)
                      x.train = x[ind[-sel], ]
                      x.test = x[ind[sel], ]
                      beta = FastStepGraph(x.train,
                                           alpha_f = alpha_f[i],
                                           alpha_b = alpha_b[j],
                                           nei.max = nei.max)$beta
                      loss = loss + sum(colSums((x.test - x.test%*%beta)^2))
                  }
                  new_loss = loss/n_folds
                  loss_grid[i, j] = new_loss
                  if (new_loss < old_loss) {
                      old_loss = new_loss
                      alpha_f_opt = alpha_f[i]
                      alpha_b_opt = alpha_b[j]
                  }
              }
          }
        c(old_loss, alpha_f_opt, alpha_b_opt)
      }
      stopCluster(cl) #stop cluster

      best_loss = losses_and_alphas[1]
      alpha_f_opt = losses_and_alphas[2]
      alpha_b_opt = losses_and_alphas[3]
      for(i in seq(4, (3*n_alpha)-2, by = 3)){
        if (losses_and_alphas[i] < best_loss) {
          best_loss = losses_and_alphas[i]
          alpha_f_opt = losses_and_alphas[i+1]
          alpha_b_opt = losses_and_alphas[i+2]
        }
      }
    } # if (parallel)
    
    else {
      for (i in 1:n_alpha) {
        for (j in 1:n_alpha) {
          if (alpha_f[i] >= alpha_b[j]) {
            loss = 0
            for (k in 1:n_folds) {
              sel = ((k-1)*ntest+1):(k*ntest)
              x.train = x[ind[-sel], ]
              x.test = x[ind[sel], ]
              beta = FastStepGraph(x.train,
                                   alpha_f = alpha_f[i],
                                   alpha_b = alpha_b[j],
                                   nei.max = nei.max)$beta
              loss = loss + sum(colSums((x.test - x.test%*%beta)^2))
            }
            new_loss = loss/n_folds
            loss_grid[i, j] = new_loss
            if (new_loss < old_loss) {
              old_loss = new_loss
              alpha_f_opt = alpha_f[i]
              alpha_b_opt = alpha_b[j]
            }
          }
        }
      }
    }

    if (return_model) {
        G = FastStepGraph(x, alpha_f = alpha_f_opt, alpha_b = alpha_b_opt, nei.max = nei.max)
        return(list(alpha_f_opt = alpha_f_opt, alpha_b_opt = alpha_b_opt, CV.loss = old_loss, model=G))
    }
    else{
        return(list(alpha_f_opt = alpha_f_opt, alpha_b_opt = alpha_b_opt, CV.loss = old_loss))
    }
}
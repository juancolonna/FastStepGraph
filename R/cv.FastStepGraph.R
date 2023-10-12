#' @title Searches for the optimal combination of alpha_f and alpha_b parameters using Cross-Validation
#'
#' @description \code{cv.FastStepGraph} implements the cross-validation for the Fast Step Graph algorithm.
#'
#' @param x Data matrix (of size n x p).
#' @param n_folds Number of folds for the cross-validation procedure (default value 5).
#' @param alpha_f_min Minimum threshold value for the cross-validation procedure (default value 0.2).
#' @param alpha_f_max Minimum threshold value for the cross-validation procedure  (default value 0.8).
#' @param b_coef This parameter applies the empirical rule alpha_b=b_coef*alpha_f during the initial search for the optimal alpha_f parameter while alpha_b remains fixed, after finding optimal alpha_f, alpha_b is varied to find its optimal value. The default value of b_coef is 0.5.
#' @param n_alpha Number of elements in the grid for the cross-validation (default value 32).
#' @param nei.max Maximum number of variables in every neighborhood (default value 5).
#' @param data_scale Boolean parameter (TRUE or FALSE), when to scale data to zero mean and unit variance (default FALSE).
#' @param data_shuffle Boolean parameter (default TRUE), when samples (rows of X) must be randomly shuffled.
#' @param max.iterations Maximum number of iterations (integer), the defaults values is set to p*(p-1).
#' @param return_model Default FALSE. If set to TRUE, at the end of cross-validation, FastStepGraph is called with the optimal parameters alpha_f and alpha_b, returning \code{vareps}, \code{beta}, \code{Edges} and \code{Omega}.
#' @param parallel Boolean parameter (TRUE or FALSE), when to run Cross-Validation in parallel using a multicore architecture (default FALSE).
#' @param n_cores An 'int' value specifying the number of cores do you want to use if 'parallel=TRUE'. 
#' If n_cores is not specified, the maximum number of cores on your machine minus one will be set automatically.
#' 
#' @return A list with the values:
#' \item{\code{alpha_f_opt}}{the optimal alpha_f value.}
#' \item{\code{alpha_f_opt}}{the optimal alpha_f value.}
#' \item{\code{CV.loss}}{minimum loss.}
#' If return_model=TRUE, then also returns:
#' \item{\code{vareps}}{Response variables.}
#' \item{\code{beta}}{Regression coefficients.}
#' \item{\code{Edges}}{Estimated set of edges.}
#' \item{\code{Omega}}{Estimated precision matrix.}
#'
#' @author Prof. Juan G. Colonna, PhD. \email{juancolonna@icomp.ufam.edu.br}
#' @author Prof. Marcelo Ruiz, PhD. \email{mruiz@exa.unrc.edu.ar}
#'
#' @export
#' @importFrom foreach %dopar%
#' 
#' @examples
#' data <- FastStepGraph::SigmaAR(30, 50, 0.4) # Simulate Gaussian Data
#' res <- FastStepGraph::cv.FastStepGraph(data$X, data_scale=TRUE)
cv.FastStepGraph <- function(x, n_folds = 5, 
                             alpha_f_min = 0.2, 
                             alpha_f_max = 0.8, 
                             b_coef = 0.5,
                             n_alpha = 32, 
                             nei.max = 5, 
                             data_scale = FALSE, 
                             data_shuffle = TRUE, 
                             max.iterations = NULL,
                             return_model = FALSE,
                             parallel = FALSE,
                             n_cores = NULL){
  f <- NULL
  n <- nrow(x)
  p <- ncol(x)

  if (nei.max == 0) { stop('The minimum number of neighbors (nei.max) must be greater than 0.') }
  if (nei.max >= n && n <= p) { stop('Neiborgs must be less than n-1') }
  if (nei.max >= p && p <= n) { stop('Neiborgs must be less than p-1') }
  if (nei.max >= p && p <= n) { stop('Neiborgs must be less than p-1') }
  if (n_alpha <= 3) { stop('n_alpha must be larger than 3') }
  if (n_folds <= 1) { stop('Number of folds must be equal or larger than 2') }
  if (is.numeric(n_cores)) { if (n_cores <= 0) {stop('n_cores must be greater than 0.')}}
  if ((n/n_folds) < 2 ) { stop('Insufficient number of samples to perform cross-validation.') }
  if (data_scale) { x <- scale(x) }
  if (data_shuffle) { x <- x[sample(seq_len(n)),] }
  if (is.null(max.iterations)){ max.iterations <- p*(p-1) }
  
  ntest <- floor(n/n_folds)
  ntrain <- n - ntest
  alpha_f <- seq(alpha_f_min, alpha_f_max, length=n_alpha)

  old_loss <- Inf
  new_loss <- NA
  alpha_f_opt <- NA
  alpha_b_opt <- NA

  if (parallel) {
    if (is.null(n_cores)) { n_cores <- parallel::detectCores()-1 }  # -1 to not overload your computer 
    if(.Platform$OS.type == "unix") {
      cl <- parallel::makeCluster(n_cores[1], type = "FORK") # in linux "FORK"
      doParallel::registerDoParallel(cl)
    }
    else {
      cl <- parallel::makeCluster(n_cores[1], type="SOCK") # in windows "SOCK"
      doParallel::registerDoParallel(cl)
    }

    alpha_f_losses <- foreach::foreach(f = alpha_f, .combine = 'c', .inorder=TRUE) %dopar% {
      loss <- 0
      for (k in seq_len(n_folds)) {
        sel <- ((k-1)*ntest+1):(k*ntest)
        x.train <- x[-sel, ]
        x.test <- x[sel, ]
        beta <- FastStepGraph(x.train,
                             alpha_f = f,
                             alpha_b = b_coef*f,
                             # alpha_b = 0.5*f,
                             nei.max = nei.max,
                             max.iterations = max.iterations)$beta
        loss <- loss + sum(colSums((x.test - x.test%*%beta)^2))
      }
      loss/n_folds
    }

    indx_min_loss <- min(which.min(alpha_f_losses))
    alpha_f_opt <- alpha_f[indx_min_loss]
    alpha_b_opt <- b_coef*alpha_f_opt
    # alpha_b_opt <- 0.5*alpha_f_opt    
    
    alpha_b <- seq(0.1, 0.9*alpha_f_opt, length=10*alpha_f_opt)
    alpha_b_losses <- foreach::foreach(b = alpha_b, .combine = 'c', .inorder=TRUE) %dopar% {
      loss <- 0
      for (k in seq_len(n_folds)) {
        sel <- ((k-1)*ntest+1):(k*ntest)
        x.train <- x[-sel, ]
        x.test <- x[sel, ]
        beta <- FastStepGraph(x.train,
                             alpha_f = alpha_f_opt,
                             alpha_b = b,
                             nei.max = nei.max,
                             max.iterations = max.iterations)$beta
        loss <- loss + sum(colSums((x.test - x.test%*%beta)^2))
      }
      loss/n_folds
    }

    parallel::stopCluster(cl) #stop cluster
    
    if (min(alpha_b_losses) < min(alpha_f_losses)){
      indx_min_loss <- min(which.min(alpha_b_losses))
      alpha_b_opt <- alpha_b[indx_min_loss]
    }
  }
  
  else {
    for (f in alpha_f) {
      loss <- 0
      for (k in seq_len(n_folds)) {
        sel <- ((k-1)*ntest+1):(k*ntest)
        x.train <- x[-sel, ]
        x.test <- x[sel, ]
        beta <- FastStepGraph(x.train,
                             alpha_f = f,
                             alpha_b = b_coef*f,
                             # alpha_b = 0.5*f,
                             nei.max = nei.max,
                             max.iterations = max.iterations)$beta
        loss <- loss + sum(colSums((x.test - x.test%*%beta)^2))
      }
      if (loss/n_folds < old_loss) {
        old_loss <- loss/n_folds
        alpha_f_opt <- f
        alpha_b_opt <- b_coef*f
        # alpha_b_opt <- 0.5*f
      }
    }

    alpha_b <- seq(0.1, 0.9*alpha_f_opt, length=10*alpha_f_opt)
    for (b in alpha_b) {
      loss <- 0
      for (k in seq_len(n_folds)) {
        sel <- ((k-1)*ntest+1):(k*ntest)
        x.train <- x[-sel, ]
        x.test <- x[sel, ]
        beta <- FastStepGraph(x.train,
                             alpha_f = alpha_f_opt,
                             alpha_b = b,
                             nei.max = nei.max,
                             max.iterations = max.iterations)$beta
        loss <- loss + sum(colSums((x.test - x.test%*%beta)^2))
      }
      if (loss/n_folds < old_loss) {
        old_loss <- loss/n_folds
        alpha_b_opt <- b
      }
    }
  }
  
  if (return_model) {
    G <- FastStepGraph(x, alpha_f = alpha_f_opt, alpha_b = alpha_b_opt, nei.max = nei.max)
   return(list(alpha_f_opt = alpha_f_opt,
              alpha_b_opt = alpha_b_opt,
              CV.loss = old_loss,
              vareps = G$vareps,
              beta = G$beta,
              Edges = G$Edges,
              Omega = G$Omega))
  }
  else {
    return(list(alpha_f_opt = alpha_f_opt,
                alpha_b_opt = alpha_b_opt,
                CV.loss = old_loss))
  }
}

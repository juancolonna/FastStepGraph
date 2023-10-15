#' @title Searches for the optimal combination of alpha_f and alpha_b parameters using Cross-Validation
#'
#' @description \code{cv.FastStepGraph} implements the cross-validation for the Fast Step Graph algorithm.
#'
#' @param x Data matrix (of size n x p).
#' @param alpha_f_min Minimum alpha_f value for the cross-validation procedure (example 0.1).
#' @param alpha_f_max Maximum alpha_f value for the cross-validation procedure  (example 0.9).
#' @param n_folds Number of folds for the cross-validation procedure (default value 10). This parameter also accepts the string 'LOOCV' to perform Leave-One-Out cross-validation.
#' @param b_coef This parameter applies the empirical rule alpha_b=b_coef*alpha_f during the initial search for the optimal alpha_f parameter while alpha_b remains fixed, after finding optimal alpha_f, alpha_b is varied to find its optimal value. The default value of b_coef is 0.5.
#' @param n_alpha Number of elements in the grid for the cross-validation (default value 20).
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
#' res <- FastStepGraph::cv.FastStepGraph(data$X, alpha_f_min=0.1, alpha_f_max = 0.9, data_scale=TRUE)
cv.FastStepGraph <- function(x, 
                             alpha_f_min, 
                             alpha_f_max, 
                             n_folds = 10,
                             b_coef = 0.5,
                             n_alpha = 20, 
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

  if (alpha_f_max >= 1) { stop('Please, decrease alpha_f_max.') }
  if (alpha_f_min <= 0) { stop('Please, increase alpha_f_min.') }
  if (alpha_f_min >= alpha_f_max) { stop('alpha_f_min must be less than alpha_f_max.') }
  if (nei.max == 0) { stop('The minimum number of neighbors (nei.max) must be greater than 0.') }
  if (nei.max >= n && n <= p) { stop('Neiborgs must be less than n-1') }
  if (nei.max >= p && p <= n) { stop('Neiborgs must be less than p-1') }
  if (nei.max >= p && p <= n) { stop('Neiborgs must be less than p-1') }
  if (n_alpha <= 3) { stop('n_alpha must be larger than 3') }
  if (n_folds < 2) { stop('Number of folds must be equal or larger than 2') }
  if (is.numeric(n_cores)) { if (n_cores <= 0) {stop('n_cores must be greater than 0.')}}
  if (data_scale) { x <- scale(x) }
  if (data_shuffle) { x <- x[sample(seq_len(n)),] }
  if (is.null(max.iterations)){ max.iterations <- p*(p-1) }
  if (n_folds == 'LOOCV' ) { n_folds = n }
  
  ntest <- floor(n/n_folds)
  ntrain <- n - ntest
  alpha_f <- seq(alpha_f_min, alpha_f_max, length=n_alpha)

  min_loss <- Inf
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
                             nei.max = nei.max,
                             max.iterations = max.iterations)$beta
        loss <- loss + sum(colSums((x.test - x.test%*%beta)^2))
      }
      loss/n_folds
    }
    
    # avoiding NA
    alpha_f <- alpha_f[!is.na(alpha_f_losses)]
    alpha_f_losses <- alpha_f_losses[!is.na(alpha_f_losses)]

    if (length(alpha_f) == 0) { 
      parallel::stopCluster(cl)
      stop("Convergence failed, no alpha_f value found (possibly due to a problem in the input data).
            To solve this issue, try increasing alpha_f_min, decreasing b_coef, and/or increasing n_folds. 
            You can also add Gaussian noise to the data with zero mean and very small variance.") 
    }
    
    indx_min_loss <- min(which.min(alpha_f_losses))
    min_loss = alpha_f_losses[indx_min_loss]
    alpha_f_opt <- alpha_f[indx_min_loss]
    alpha_b_opt <- b_coef*alpha_f_opt

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
    
    # avoiding NA
    alpha_b <- alpha_b[!is.na(alpha_b_losses)]
    alpha_b_losses <- alpha_b_losses[!is.na(alpha_b_losses)]
    
    if (length(alpha_b) == 0) { }
    else if (min(alpha_b_losses) < min(alpha_f_losses)){
      indx_min_loss <- min(which.min(alpha_b_losses))
      min_loss = alpha_b_losses[indx_min_loss]
      alpha_b_opt <- alpha_b[indx_min_loss]
    }
  }
  
  # Not Parallel loop
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
                             nei.max = nei.max,
                             max.iterations = max.iterations)$beta
        loss <- loss + sum(colSums((x.test - x.test%*%beta)^2))
      }
 
      if (is.na(loss)){ }
      else if (loss/n_folds < min_loss) {
        min_loss <- loss/n_folds
        alpha_f_opt <- f
        alpha_b_opt <- b_coef*f
      }
    }
    
    if (is.infinite(min_loss)) {
      stop("Convergence failed, no alpha_f value found (possibly due to a problem in the input data).
            To solve this issue, try increasing alpha_f_min, decreasing b_coef, and/or increasing n_folds. 
            You can also add Gaussian noise to the data with zero mean and very small variance.")
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
      if (is.na(loss)){ }
      else if (loss/n_folds < min_loss) {
        min_loss <- loss/n_folds
        alpha_b_opt <- b
      }
    }
  }
  
  if (return_model) {
    G <- FastStepGraph(x, alpha_f = alpha_f_opt, alpha_b = alpha_b_opt, nei.max = nei.max)
   return(list(alpha_f_opt = alpha_f_opt,
              alpha_b_opt = alpha_b_opt,
              CV.loss = min_loss,
              vareps = G$vareps,
              beta = G$beta,
              Edges = G$Edges,
              Omega = G$Omega))
  }
  else {
    return(list(alpha_f_opt = alpha_f_opt,
                alpha_b_opt = alpha_b_opt,
                CV.loss = min_loss))
  }
}

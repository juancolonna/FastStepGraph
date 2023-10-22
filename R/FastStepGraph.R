#' @title Fast Stepwise Gaussian Graphical Model
#'
#' @description Improved and faster implementation of the Stepwise Gaussian Graphical Algorithm.
#'
#' @param x Data matrix (of size n_samples x p_variables).
#' @param alpha_f Forward threshold  (no default value).
#' @param alpha_b Backward threshold. If alpha_b=NULL, then the rule alpha_b <- 0.5*alpha_f is applied.
#' @param nei.max Maximum number of variables in every neighborhood (default value 5).
#' @param data_scale Boolean parameter (TRUE or FALSE), when to scale data to zero mean and unit variance (default FALSE).
#' @param max.iterations Maximum number of iterations (integer), the defaults values is set to p*(p-1).
#'
#' @return A list with the values:
#' \item{\code{vareps}}{Response variables.}
#' \item{\code{beta}}{Regression coefficients.}
#' \item{\code{Edges}}{Estimated set of edges.}
#' \item{\code{Omega}}{Estimated precision matrix.}
#'
#' @author Prof. Juan G. Colonna, PhD. \email{juancolonna@icomp.ufam.edu.br}
#' @author Prof. Marcelo Ruiz, PhD. \email{mruiz@exa.unrc.edu.ar}
#'
#' @export
#' 
#' @examples
#' data <- FastStepGraph::SigmaAR(30, 50, 0.4) # Simulate Gaussian Data
#' G <- FastStepGraph::FastStepGraph(data$X, alpha_f = 0.22, alpha_b = 0.14, data_scale=TRUE)
FastStepGraph <- function(x, 
                          alpha_f, 
                          alpha_b=NULL, 
                          nei.max=5L, 
                          data_scale=FALSE, 
                          max.iterations=NULL){
  
  .lm.fit = combn = cor = cov = .colvars = .neighbors_of = .compute_omega_and_beta = NULL
  
  if (data_scale) { x <- scale(x) }
  n <- dim(x)[1] # number of rows
  p <- dim(x)[2] # number of columns

  if (alpha_f < alpha_b){ stop("alpha_b must be lower than alpha_f") }
  if (is.null(alpha_b)){ alpha_b <- 0.5*alpha_f }
  if (nei.max == 0L) { stop('The minimum number of neighbors (nei.max) must be greater than 0.') }
  if (nei.max >= n && n <= p) { stop('Neiborgs must be less than n-1') }
  if (nei.max >= p && p <= n) { stop('Neiborgs must be less than p-1') }

  Edges_I <- t(combn(seq_len(p),2)) # Inactive set of ordered pair (i,j)
  Edges_A <- t(matrix(0L,2,dim(Edges_I)[1])) # Active set of ordered pair (i,j)
  N_neighbors <- integer(p) # number of neighbors of each node

  e <- x # (n x p) matrix of regression residuals
  f_ij <- cor(e)
  if (any(is.na(f_ij))) { stop('Columns with zero variance found in input data') }
  lower_tri <- lower.tri(f_ij)
  f_ij <- abs(f_ij[lower_tri])
  b_ij <- rep(2, length(f_ij))

  k <- 1L
  if (is.null(max.iterations)){ max.iterations <- as.integer(p*(p-1)) } # length(Edges_I[,2]) is the total number of edges

  while ((k <= max.iterations)) {
    # Select (i,j) that max(|f_ij|)
    f_ij_indx <- which.max(f_ij)
    f_ij_max <- f_ij[f_ij_indx]

    i_f <- Edges_I[f_ij_indx, 1]
    j_f <- Edges_I[f_ij_indx, 2]

    if ( f_ij_max < alpha_f ) {break} # if |f_ij| lower than alpha, then stop.

    # If the number of neighbors for a given node exceeds the maximum allowed number
    # of neighbors (nei.max), refrain from adding a new neighbor to that node.
    else if ((N_neighbors[i_f]+1L) > nei.max || (N_neighbors[j_f]+1L) > nei.max) { f_ij[f_ij_indx] <- 0 }

    else{
      ################ Forward-Step ######################################
      N_neighbors[i_f] <- N_neighbors[i_f]+1L
      N_neighbors[j_f] <- N_neighbors[j_f]+1L

      Edges_A[f_ij_indx,] <- Edges_I[f_ij_indx,]
      Edges_I[f_ij_indx,] <- c(0L,0L)

      # Update Prediction Errors for (i_f,j_f)
      n_i_f_pos_2 <- which(Edges_A[,2]==i_f)
      n_i_f_pos_1 <- which(Edges_A[,1]==i_f)
      n_i_f <- c(Edges_A[n_i_f_pos_2, 1], Edges_A[n_i_f_pos_1, 2])

      n_j_f_pos_2 <- which(Edges_A[,2]==j_f)
      n_j_f_pos_1 <- which(Edges_A[,1]==j_f)
      n_j_f <- c(Edges_A[n_j_f_pos_2, 1], Edges_A[n_j_f_pos_1, 2])

      # fast LM documentation https://rpubs.com/maechler/fast_lm
      e[,i_f] <- .lm.fit(cbind(1, x[,n_i_f]), x[,i_f])$residuals
      e[,j_f] <- .lm.fit(cbind(1, x[,n_j_f]), x[,j_f])$residuals

      ############### Residuals update ###################################
      # Compute Prediction Errors r_i and r_j for (i,j) in the active set Edges_A
      # L <- unique(c(n_i_f_pos_2, n_i_f_pos_1, n_j_f_pos_1, n_j_f_pos_2)) # check if unique is necessary
      L <- c(n_i_f_pos_2, n_i_f_pos_1, n_j_f_pos_1, n_j_f_pos_2)
      L <- L[-which(L == f_ij_indx)]
      
      for (l in L) { b_ij[l] <- .residuals_update(l, L, Edges_A, x) }
      
      ############### Backward-Step ######################################
      # Select (i,j) that min(|b_ij|) <= alpha_b
      b_ij[f_ij_indx] <- f_ij_max # don't change this line, otherwise it'll cause an infinite loop
      b_ij_indx <- which.min(b_ij)
      b_ij_min <- b_ij[b_ij_indx]

      # Update active set of edges Edges_A
      if (b_ij_min <= alpha_b){

        i_b <- Edges_A[b_ij_indx, 1]
        j_b <- Edges_A[b_ij_indx, 2]

        N_neighbors[i_b] <- N_neighbors[i_b]-1L
        N_neighbors[j_b] <- N_neighbors[j_b]-1L

        Edges_I[b_ij_indx,] <- Edges_A[b_ij_indx,]
        Edges_A[b_ij_indx,] <- c(0L,0L)

        # Update Prediction Errors for (i_b,j_b)
        n_i_b <- .neighbors_of(i_b, Edges_A)
        n_j_b <- .neighbors_of(j_b, Edges_A)

        if(N_neighbors[i_b] > 0L) { e[,i_b] <- .lm.fit(cbind(1, x[,n_i_b]), x[,i_b])$residuals }
        else { e[,i_b] <- x[,i_b] }
        
        if(N_neighbors[j_b] > 0L) { e[,j_b] <- .lm.fit(cbind(1, x[,n_j_b]), x[,j_b])$residuals }
        else { e[,j_b] <- x[,j_b] }

        b_ij[b_ij_indx] <- 2
      }

      f_ij <- cor(e)
      f_ij <- abs(f_ij[lower_tri])
      f_ij[is.na(f_ij)] <- 0 # I have to check this condition
      H <- which(Edges_A[,1] > 0L)
      f_ij[H] <- 0

      ############# steps of main loop (while) ###########################
      k <- k+1L
      # if (k > max.iterations){ message('Maximum number iteration reached') }
    }
  }

  ################## Release some memory #####################################
  rm(k, f_ij, b_ij, Edges_I, N_neighbors, x)

  ################## Compute the precision matrix ############################
  M <- .compute_omega_and_beta(p, e, Edges_A)
  return(list(vareps=e, beta=M$beta, Edges=Edges_A, Omega=M$Omega))
}

.colvars <- function(x) {
  means <- colMeans(x)
  n <- length(x)/length(means)
  centeredX <- x - rep(means, each = n)
  colSums(centeredX^2) / (n-1)
}

.neighbors_of <- function(node, Edges) {
  pos_2 <- which(Edges[,2]==node)
  pos_1 <- which(Edges[,1]==node)
  c(Edges[pos_2, 1], Edges[pos_1, 2])
}

.residuals_update <- function(l, L, Edges_A, x) {
  i <- Edges_A[l, 1]
  j <- Edges_A[l, 2]
  
  n_i <- .neighbors_of(i, Edges_A)
  n_i <- n_i[-which(n_i==j)]
  
  if ( length(n_i) > 0 ){ r_i <- .lm.fit(cbind(1, x[,n_i]), x[,i])$residuals}
  else { r_i <- x[,i] }
  
  n_j <- .neighbors_of(j, Edges_A)
  n_j <- n_j[-which(n_j==i)]
  
  if( length(n_j) > 0 ){ r_j <- .lm.fit(cbind(1, x[,n_j]), x[,j])$residuals }
  else{r_j <- x[,j]}
  
  corr_b <- cor(r_i, r_j)
  ifelse(is.na(corr_b), 0, abs(corr_b)) # I have to check this condition
}

.compute_omega_and_beta <- function(p, e, Edges) {
  col_var <- .colvars(e)
  cor_matrix <- cov(e)
  cor_matrix[is.na(cor_matrix)] <- 0 # I have to check this condition
  Omega <- matrix(0,p,p)
  diag(Omega) <- col_var^(-1)
  Omega[is.infinite(Omega)] <- 1/.Machine$double.eps # I have to check this condition
  beta <- matrix(0,p,p)
  
  for (edge in which(Edges[,2] > 0L)){
    i <- Edges[edge, 1]
    j <- Edges[edge, 2]
    
    Omega[i,j] <- cor_matrix[i,j]*Omega[i,i]*Omega[j,j]
    Omega[j,i] <- Omega[i,j]
    
    if(col_var[j] <= .Machine$double.eps) {beta[i,j] <- 0} # I have to check this condition
    else {beta[i,j] <- -cor_matrix[i,j]/col_var[j]}
    
    if(col_var[i] <= .Machine$double.eps) {beta[j,i] <- 0} # I have to check this condition
    else {beta[j,i] <- -cor_matrix[i,j]/col_var[i]}
  }
  # Omega[abs(Omega) < 1e-5] = 0 # CHEQUEAR ESTA CONDIÇÃO
  return(list(beta=beta, Omega=Omega))
}
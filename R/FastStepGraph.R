#' @title Fast Stepwise Gaussian Graphical Model
#'
#' @description Improved and faster implementation of the Stepwise Gaussian Graphical Algorithm.
#'
#' @param x Data matrix (of size n_samples x p_variables).
#' @param alpha_f Forward threshold  (no default value).
#' @param alpha_b Backward threshold. If alpha_b=0, then the alpha_b=0.5*alpha_f rule is applied.
#' @param nei.max Maximum number of variables in every neighborhood (default value 5).
#' @param data_scale Boolean parameter (TRUE or FALSE), when to scale data to zero mean and unit variance (default FALSE).
#' @param max.iterations Maximum number of iterations (integer), defaults set to three times the combination C(p,2)= p!/(2! * (p-2)!)
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
#' G <- FastStepGraph::FastStepGraph(data$X, alpha_f = 0.22, alpha_b = 0.14)
FastStepGraph <- function(x, alpha_f, alpha_b=0, nei.max=5, data_scale=FALSE, max.iterations=NULL){
  .lm.fit = combn = cor = cov = NULL
  if (data_scale) { x <- scale(x) }
  n <- dim(x)[1] # number of rows
  p <- dim(x)[2] # number of columns

  if (alpha_f < alpha_b){ stop("alpha_b must be lower than alpha_f") }
  if (alpha_b == 0){ alpha_b <- 0.5*alpha_f }
  if (nei.max >= n) { stop('Neiborgs must be less than n-1') }
  if (nei.max == 0 && n <= p){nei.max <- n-1 }
  if (nei.max == 0 && p < n){nei.max <- p }

  Edges_I <- t(combn(seq_len(p),2)) # Inactive set of ordered pair (i,j)
  Edges_A <- t(matrix(0,2,dim(Edges_I)[1])) # Active set of ordered pair (i,j)
  N_neighbors <- integer(p) # number of neighbors of each node

  e <- x # (n x p) matrix of regression residuals
  f_ij <- cor(e)
  lower_tri <- lower.tri(f_ij)
  f_ij <- abs(f_ij[lower_tri])
  b_ij <- rep(2, length(f_ij))

  k <- 1
  if (is.null(max.iterations)){ K <- 2*length(Edges_I[,2]) } # length(Edges_I[,2]) is the total number of edges

  while ((k <= K)) {
    # Select (i,j) that max(|f_ij|)
    f_ij_indx <- which.max(f_ij)
    f_ij_max <- f_ij[f_ij_indx]

    i_f <- Edges_I[f_ij_indx, 1]
    j_f <- Edges_I[f_ij_indx, 2]

    if ( f_ij_max < alpha_f ) {break} # if |f_ij| lower than alpha, then stop.

    # If the number of neighbors for a given node exceeds the maximum allowed number
    # of neighbors (nei.max), refrain from adding a new neighbor to that node.
    else if ((N_neighbors[i_f]+1) > nei.max || (N_neighbors[j_f]+1) > nei.max) { f_ij[f_ij_indx] <- 0 }

    else{
      ################ Forward-Step ######################################
      N_neighbors[i_f] <- N_neighbors[i_f]+1
      N_neighbors[j_f] <- N_neighbors[j_f]+1

      Edges_A[f_ij_indx,] <- Edges_I[f_ij_indx,]
      Edges_I[f_ij_indx,] <- c(0,0)

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

      ############### Backward-Step ######################################
      # Compute Prediction Errors r_i and r_j for (i,j) in the active set Edges_A
      L <- unique(c(n_i_f_pos_2, n_i_f_pos_1, n_j_f_pos_1, n_j_f_pos_2))
      b_ij[f_ij_indx] <- f_ij_max
      L <- L[-which(L == f_ij_indx)]

      ############### Residuals update ###################################
      for (l in L) {
        i <- Edges_A[l, 1]
        j <- Edges_A[l, 2]

        n_i <- c(Edges_A[which(Edges_A[,2] == i), 1], Edges_A[which(Edges_A[,1] == i), 2])
        n_i <- n_i[-which(n_i==j)]

        if ( length(n_i) > 0 ){ r_i <- .lm.fit(cbind(1, x[,n_i]), x[,i])$residuals}
        else { r_i <- x[,i] }

        n_j <- c(Edges_A[which(Edges_A[,2] == j), 1], Edges_A[which(Edges_A[,1] == j), 2])
        n_j <- n_j[-which(n_j==i)]

        if( length(n_j)>0 ){ r_j <- .lm.fit(cbind(1, x[,n_j]), x[,j])$residuals }
        else{r_j <- x[,j]}

        b_ij[l] <- abs(cor(r_i, r_j))
      }

      ################ Edge elimination ######################################
      # Select (i,j) that min(|b_ij|) <= alpha_b
      b_ij_indx <- which.min(b_ij)
      b_ij_min <- b_ij[b_ij_indx]

      # Update active set of edges Edges_A
      if (b_ij_min <= alpha_b){

        i_b <- Edges_A[b_ij_indx, 1]
        j_b <- Edges_A[b_ij_indx, 2]

        N_neighbors[i_b] <- N_neighbors[i_b]-1
        N_neighbors[j_b] <- N_neighbors[j_b]-1

        Edges_I[b_ij_indx,] <- Edges_A[b_ij_indx,]
        Edges_A[b_ij_indx,] <- c(0,0)

        # Update Prediction Errors for (i_b,j_b)
        n_i_b <- c(Edges_A[which(Edges_A[,2] == i_b), 1], Edges_A[which(Edges_A[,1] == i_b), 2])
        n_j_b <- c(Edges_A[which(Edges_A[,2] == j_b), 1], Edges_A[which(Edges_A[,1] == j_b), 2])

        if( N_neighbors[i_b] > 0) { e[,i_b] <- .lm.fit(cbind(1, x[,n_i_b]), x[,i_b])$residuals}
        else{ e[,i_b] <- x[,i_b] }

        if( N_neighbors[j_b] > 0) { e[,j_b] <- .lm.fit(cbind(1, x[,n_j_b]), x[,j_b])$residuals}
        else{ e[,j_b] <- x[,j_b] }

        b_ij[b_ij_indx] <- 2
      }

      f_ij <- cor(e)
      f_ij <- abs(f_ij[lower_tri])
      H <- which(Edges_A[,1]>0)
      f_ij[H] <- 0


      ############# steps of main loop (while) ###########################
      k <- k+1
      if (k > K){ message('Maximum number iteration reached') }
    }
  }

  ################## Release some memory #####################################
  rm(f_ij, b_ij, Edges_I, N_neighbors)

  ################## Compute the precision matrix ############################
  col_var <- colvars(e)
  cor_matrix <- cov(e)
  Omega <- matrix(0,p,p)
  diag(Omega) <- col_var^(-1)
  beta <- matrix(0,p,p)

  for (edge in which(Edges_A[,2] > 0)){
    i <- Edges_A[edge, 1]
    j <- Edges_A[edge, 2]

    Omega[i,j] <- cor_matrix[i,j]*Omega[i,i]*Omega[j,j]
    Omega[j,i] <- Omega[i,j]

    beta[i,j] <- -cor_matrix[i,j]/col_var[j]
    beta[j,i] <- -cor_matrix[i,j]/col_var[i]
  }
  return(list(vareps=e, beta=beta, Edges=Edges_A, Omega=Omega))
}

colvars <- function(x) {
  means <- colMeans(x)
  n <- length(x)/length(means)
  centeredX <- x - rep(means, each = n)
  colSums(centeredX^2) / (n-1)
}

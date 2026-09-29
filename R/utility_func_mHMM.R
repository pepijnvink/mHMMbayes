#' @keywords internal
# Whenever you use C++ code in your package, you need to clean up after yourself
# when your package is unloaded. This function unloads the DLL (H. Wickham(2019). R packages)
.onUnload <- function (libpath) {
  library.dynam.unload("mHMMbayes", libpath)
}

#' @keywords internal
# simple functions used in mHMM
dif_matrix <- function(rows, cols, data = NA){
  return(matrix(data, ncol = cols, nrow = rows))
}

#' @keywords internal
nested_list <- function(n_dep, m){
  return(rep(list(vector("list", n_dep)),m))
}

#' @keywords internal
dif_vector <- function(x){
  return(numeric(x))
}

#' @keywords internal
is.whole <- function(x) {
  return(is.numeric(x) && floor(x) == x)
}

#' @keywords internal
is.mHMM <- function(x) {
  inherits(x, "mHMM")
}

#' @keywords internal
is.mHMM_gamma <- function(x) {
  inherits(x, "mHMM_gamma")
}

#' @keywords internal
is.mHMM_prior_gamma <- function(x) {
  inherits(x, "mHMM_prior_gamma")
}

#' @keywords internal
is.mHMM_prior_emiss <- function(x) {
  inherits(x, "mHMM_prior_emiss")
}

#' @keywords internal
is.cat <- function(x) {
  inherits(x, "cat")
}

#' @keywords internal
is.cont <- function(x) {
  inherits(x, "cont")
}

#' @keywords internal
is.count <- function(x) {
  inherits(x, "count")
}

#' @keywords internal
is.mHMM_pdRW_gamma <- function(x) {
  inherits(x, "mHMM_pdRW_gamma")
}

#' @keywords internal
is.mHMM_pdRW_emiss <- function(x) {
  inherits(x, "mHMM_pdRW_emiss")
}

#' @keywords internal
hms <- function(t){
  paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0"),
        formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"),
        formatC(t %% 60, width = 2, format = "d", flag = "0"),
        sep = ":")
}

#' @keywords internal
depth <- function(x,xdepth=0){
  if(!is.list(x)){
    return(xdepth)
  }else{
    return(max(unlist(lapply(x,depth,xdepth=xdepth+1))))
  }
}

#' @keywords internal
# Calculates the between subject variance from logmu and logvar:
logvar_to_var <- function(logmu, logvar){
  abs(exp(logvar)-1)*exp(2*logmu+logvar)
}

#' @keywords internal
# Use ecr algorithm
ecr <- function(pivot, alloc, m){
  n <- length(pivot) #sequence length
  conf_mat <- table(factor(alloc, levels = 1:m), factor(pivot, levels = 1:m)) # confusion matrix
  cost_mat <- max(conf_mat) - conf_mat # cost matrix to maximize
  # run hungarian algorithm. output: vector length m. element i==j element indicates that if sampled state==i, it should be relabeled to state j
  permutation <- RcppHungarian::HungarianSolver(cost_mat)$pairs[,2]
  is_switched <- !(identical(permutation, 1:m)) # check if relabeling happens
  x_repermute <- permutation[alloc] # old state == i --> take the ith element in permutation
  return(list(switched = is_switched, sequence = x_repermute, permutation = permutation))
}

#' @keywords internal
# Use ecr algorithm
ecr_observed <- function(pivot, alloc, observed, m){
  alloc_observed <- alloc[observed]
  n <- length(pivot)
  conf_mat <- table(factor(alloc_observed, levels = 1:m), factor(pivot, levels = 1:m)) # confusion matrix
  cost_mat <- max(conf_mat) - conf_mat # cost matrix to maximize
  # run hungarian algorithm. output: vector length m. element i==j element indicates that if sampled state==i, it should be relabeled to state j
  permutation <- RcppHungarian::HungarianSolver(cost_mat)$pairs[,2]
  is_switched <- !(identical(permutation, 1:m)) # check if relabeling happens
  x_repermute <- permutation[alloc] # old state == i --> take the ith element in permutation
  return(list(switched = is_switched, sequence = x_repermute, permutation = permutation))
}

#' @keywords internal
# Use PRA algorithm
pra <- function(pivot_emiss, parameters_emiss, parameters_gamma, m, n_dep){
  align_mat <- matrix(0, m, m)

  ## compute dot products (to maximize later). indicates overall similarity of reference state i with sampled parameter j
  for(i in 1:m){
    for(j in 1:m){
      align_mat[i, j] <- sum(pivot_emiss[i,]*parameters_emiss[j,]) # only use emissions for relabeling
    }
  }
  ## transform to cost (needed by hungarian algorithm. also make all elements positive)
  align_mat <- max(align_mat) - align_mat
  ## compute allocations. output is vector of length m. If element i==j, the j'th sampled state parameters will correspond to state i
  permute <- RcppHungarian::HungarianSolver(align_mat)$pairs[,2]
  permute <- round(c(permute), 0) # due to potential rounding issues
  param_emiss_relabel <- parameters_emiss[permute, ]
  param_gamma_relabel <- parameters_gamma[permute, permute]
  is_switched <- !(isTRUE(all.equal(permute, 1:m)))
  return(list(switched = is_switched, emiss_relabeled = param_emiss_relabel, gamma_relabeled = param_gamma_relabel))
}

#' @keywords internal
#' Align pivot sequence to group level when using ECR. Uses the same procedure as PRA
ecr_align_group <- function(
  pivot_emiss, # group_level parameter sequence
  parameters_emiss, # emission parameters
  freq_table, # frequency table with n_vary rows and m columns
  m,
  n_dep
){
    align_mat <- matrix(0, m, m)
  ## compute dot products (to maximize later). indicates overall similarity of reference state i with sampled parameter j
  for(i in 1:m){
    for(j in 1:m){
      align_mat[i, j] <- sum(pivot_emiss[i,]*parameters_emiss[j,]) # only use emissions for relabeling
    }
  }
  ## transform to cost (needed by hungarian algorithm. also make all elements positive)
  align_mat <- max(align_mat) - align_mat
  ## compute allocations. output is vector of length m. If element i==j, the j'th sampled state parameters will correspond to state i. So if permute[1] == 2, second row must become the first row.
  permute <- RcppHungarian::HungarianSolver(align_mat)$pairs[,2]
  permute <- round(c(permute), 0) # due to potential rounding issues
  param_emiss_relabel <- parameters_emiss[permute, ]
  freq_table <- freq_table[, permute]
  is_switched <- !(isTRUE(all.equal(permute, 1:m)))
  return(list(freq_table = freq_table, param_emiss_relabel = param_emiss_relabel))
}

#' @keywords internal
# int_to_prob() without rounding
int_to_prob_noround <- function(int_matrix) {
  if(!is.matrix(int_matrix)){
    stop("int_matrix should be a matrix")
  }
  prob_matrix <- matrix(nrow = nrow(int_matrix), ncol = ncol(int_matrix) + 1)
  for(r in 1:nrow(int_matrix)){
    exp_int_matrix 	<- matrix(exp(c(0, int_matrix[r,])), nrow  = 1)
    prob_matrix[r,] <- exp_int_matrix / as.vector(exp_int_matrix %*% c(rep(1, (dim(exp_int_matrix)[2]))))
  }
  return(prob_matrix)
}



#' @keywords internal
# Reshape the categorical emission probabilities as they are stored in
# PD_subj[[s]]$cat_emiss (one block per dependent variable, within a block the
# q_emiss[q] probabilities of state 1, then of state 2, etc.) into an
# m x sum(q_emiss) matrix, with the states in the rows and the dependent
# variables concatenated over the columns. This is the layout pra_cat() expects,
# and is the categorical counterpart of the m x n_dep matrix of emission means
# that is handed to pra() for continuous observations.
cat_emiss_to_mat <- function(x, m, q_emiss){
  n_dep <- length(q_emiss)
  out   <- matrix(NA_real_, nrow = m, ncol = sum(q_emiss))
  for(q in 1:n_dep){
    in_start  <- sum(c(0, q_emiss)[1:q]) * m
    out_start <- sum(c(0, q_emiss)[1:q])
    out[, (out_start + 1):(out_start + q_emiss[q])] <-
      matrix(x[(in_start + 1):(in_start + m * q_emiss[q])], nrow = m, byrow = TRUE)
  }
  return(out)
}

#' @keywords internal
# Inverse of cat_emiss_to_mat(): flatten an m x sum(q_emiss) matrix back into the
# layout used by PD_subj[[s]]$cat_emiss.
cat_mat_to_emiss <- function(x, m, q_emiss){
  n_dep <- length(q_emiss)
  out   <- numeric(sum(q_emiss) * m)
  for(q in 1:n_dep){
    in_start  <- sum(c(0, q_emiss)[1:q])
    out_start <- sum(c(0, q_emiss)[1:q]) * m
    out[(out_start + 1):(out_start + m * q_emiss[q])] <-
      as.vector(t(x[, (in_start + 1):(in_start + q_emiss[q]), drop = FALSE]))
  }
  return(out)
}

#' @keywords internal
# Use PRA algorithm on categorical emission distributions. Identical to pra(),
# except that pivot_emiss and parameters_emiss are m x sum(q_emiss) matrices of
# emission probabilities instead of m x n_dep matrices of emission means, and
# that the permutation itself is returned. Because a permutation only reorders
# the set of parameter vectors, sum_j ||param_j||^2 is invariant, so maximizing
# the total dot product is equivalent to minimizing the summed squared Euclidean
# distance between the pivot and the permuted emission probability vectors.
pra_cat <- function(pivot_emiss, parameters_emiss, parameters_gamma, m){
  align_mat <- matrix(0, m, m)

  ## compute dot products (to maximize later). indicates overall similarity of reference state i with sampled parameter j
  for(i in 1:m){
    for(j in 1:m){
      align_mat[i, j] <- sum(pivot_emiss[i,]*parameters_emiss[j,]) # only use emissions for relabeling
    }
  }
  ## transform to cost (needed by hungarian algorithm. also make all elements positive)
  align_mat <- max(align_mat) - align_mat
  ## compute allocations. output is vector of length m. If element i==j, the j'th sampled state parameters will correspond to state i
  permute <- RcppHungarian::HungarianSolver(align_mat)$pairs[,2]
  permute <- round(c(permute), 0) # due to potential rounding issues
  param_emiss_relabel <- parameters_emiss[permute, , drop = FALSE]
  param_gamma_relabel <- parameters_gamma[permute, permute]
  is_switched <- !(isTRUE(all.equal(permute, 1:m)))
  return(list(switched = is_switched, emiss_relabeled = param_emiss_relabel,
              gamma_relabeled = param_gamma_relabel, permutation = permute))
}
#' @keywords internal
#' Align the subject level pivot with the group level when using ECR on categorical
#' observations. Identical to ecr_align_group(), except that pivot_emiss and
#' parameters_emiss are m x sum(q_emiss) matrices of emission probabilities instead of
#' m x n_dep matrices of emission means, that the permutation itself is returned, and
#' that the relabeled emission probabilities are returned in the flat layout of
#' PD_subj[[s]]$cat_emiss, so that the running mean of the subject level pivot can be
#' updated with them directly.
ecr_align_group_cat <- function(
  pivot_emiss, # group-level parameter vector, as an m x sum(q_emiss) matrix
  parameters_emiss, # emission parameters, as an m x sum(q_emiss) matrix
  freq_table, # frequency table with n_vary rows and m columns
  m,
  q_emiss
){
  align_mat <- matrix(0, m, m)
  ## compute dot products (to maximize later). indicates overall similarity of reference state i with sampled parameter j
  for(i in 1:m){
    for(j in 1:m){
      align_mat[i, j] <- sum(pivot_emiss[i,]*parameters_emiss[j,]) # only use emissions for relabeling
    }
  }
  ## transform to cost (needed by hungarian algorithm. also make all elements positive)
  align_mat <- max(align_mat) - align_mat
  ## compute allocations. output is vector of length m. If element i==j, the j'th sampled state parameters will correspond to state i. So if permute[1] == 2, second row must become the first row.
  permute <- RcppHungarian::HungarianSolver(align_mat)$pairs[,2]
  permute <- round(c(permute), 0) # due to potential rounding issues
  param_emiss_relabel <- parameters_emiss[permute, , drop = FALSE]
  # the j'th column of the frequency table counts the visits to sampled state j, so
  # the columns follow the same permutation as the rows of the emission parameters
  freq_table <- freq_table[, permute, drop = FALSE]
  is_switched <- !(isTRUE(all.equal(permute, 1:m)))
  return(list(switched = is_switched, freq_table = freq_table,
              param_emiss_relabel = cat_mat_to_emiss(param_emiss_relabel, m = m, q_emiss = q_emiss),
              permutation = permute))
}

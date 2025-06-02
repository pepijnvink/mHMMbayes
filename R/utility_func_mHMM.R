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
  n <- length(pivot)
  conf_mat <- table(factor(alloc, levels = 1:m), factor(pivot, levels = 1:m)) # confusion matrix
  cost_mat <- conf_mat%*%(1-diag(m))
  permutation <- RcppHungarian::HungarianSolver(cost_mat)$pairs[,2]
  is_switched <- !(identical(permutation, 1:m))
  x_repermute <- permutation[alloc]
  return(list(switched = is_switched, sequence = x_repermute))
}

#' @keywords internal
# Use ecr algorithm
ecr2 <- function(pivot, alloc, m){
  n <- length(pivot)
  conf_mat <- table(factor(alloc, levels = 1:m), factor(pivot, levels = 1:m))
  cost_mat <- conf_mat%*%(1-diag(m))
  permutation <- RcppHungarian::HungarianSolver(cost_mat)$pairs[,2]
  return(permutation)
}

#' @keywords internal
# Use ecr algorithm
ecr_observed <- function(pivot, alloc, observed, m){
  alloc_observed <- alloc[observed]
  n <- length(pivot)
  conf_mat <- table(factor(alloc_observed, levels = 1:m), factor(pivot, levels = 1:m)) # confusion matrix
  cost_mat <- conf_mat%*%(1-diag(m))
  permutation <- RcppHungarian::HungarianSolver(cost_mat)$pairs[,2]
  is_switched <- !(identical(permutation, 1:m))
  x_repermute <- permutation[alloc]
  return(list(switched = is_switched, sequence = x_repermute))
}

#' @keywords internal
# Use PRA algorithm
pra <- function(pivot_emiss, pivot_gamma, parameters_emiss, parameters_gamma, m, n_dep){
  align_mat <- matrix(0, m, m)
  for(i in 1:m){
    for(j in 1:m){
      align_mat[i, j] <- sum(pivot_emiss[i,]*parameters_emiss[j,]) + (pivot_gamma[i, i]*parameters_gamma[j,j]) # only use self-transitions for relabeling
    }
  }
  solution <- lpSolve::lp.assign(align_mat, direction = "max")$solution
  permute <- (1:m)%*%t(solution)
  permute <- round(c(permute), 0)
  param_emiss_relabel <- parameters_emiss[permute, ]
  param_gamma_relabel <- parameters_gamma[permute, permute]
  is_switched <- !(identical(permute, 1:m))
  return(list(switched = is_switched, emiss_relabeled = param_emiss_relabel, gamma_relabeled = param_gamma_relabel))
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

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
is.mHMM_cont <- function(x) {
  inherits(x, "mHMM_cont")
}

#' @keywords internal
is.mHMM_vary <- function(x) {
  inherits(x, "mHMM_vary")
}

#' @keywords internal
is.mHMM_gamma <- function(x) {
  inherits(x, "mHMM_gamma")
}

#' @keywords internal
hms <- function(t){
  paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0"),
        formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"),
        formatC(t %% 60, width = 2, format = "d", flag = "0"),
        sep = ":")
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
  return(list(switched = is_switched, sequence = x_repermute))
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
  return(list(switched = is_switched, sequence = x_repermute))
}

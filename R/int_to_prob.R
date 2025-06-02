#' Transforming a set of Multinomial logit regression intercepts to probabilities
#'
#' \code{int_to_prob} transforms a set of Multinomial logit regression
#' intercepts to the corresponding state transition or categorical emission
#' observation probabilities. Note that the first state or category is assumed
#' to be the reference category, hence no intercept is to specified for the
#' first state or category.
#'
#' Designed to ease the specification of informative hyper-prior values for the
#' mean intercepts of the transition probability matrix gamma and categorical
#' emission distribution(s) of the multilevel hidden Markov model through the
#' functions \code{\link{prior_gamma}} and \code{\link{prior_emiss_cat}}. No
#' check is performed on correct specifications of the dimensions.
#'
#' @param int_matrix A matrix with (number of states OR categories - 1) columns
#'   and number of rows to be determined by the user. For obtaining the set of
#'   probabilities of the complete transition probability matrix gamma or
#'   categorical emission distribution matrix, the number of rows equals the
#'   number of states \code{m}. The first state / category is assumed to be the
#'   reference category, no intercept is to be specified for this first
#'   category.
#'
#' @return \code{int_to_prob} returns a matrix containing probabilities with
#'   each row summing to one, with the number of columns equal to the number of
#'   states / categories and the number of rows equal to the number rows
#'   specified in the input matrix.
#'
#' @seealso \code{\link{prob_to_int}} for transforming a set of probabilities to
#'   a set of Multinomial logit regression intercepts, \code{\link{prior_gamma}}
#'   and \code{\link{prior_emiss_cat}} for specifying informative hyper-priors
#'   for the the multilevel hidden Markov model and \code{\link{mHMM}} to fit a
#'   multilevel hidden Markov model.
#'
#' @examples
#'
#' # example for transition probability matrix gamma with 3 states
#' m <- 3
#' gamma_int <- matrix(c(-1, -1,
#'                        3,  0,
#'                        0,  2), ncol = m-1, nrow = m, byrow = TRUE)
#' gamma_prob <- int_to_prob(gamma_int)
#' gamma_prob
#'
#' @export
#' @keywords internal
# computes probabilities from intercepts
int_to_prob <- function(int1) {
  if(is.matrix(int1)){
    prob1 <- matrix(nrow = nrow(int1), ncol = ncol(int1) + 1)
    for(r in 1:nrow(int1)){
      exp_int1 	<- matrix(exp(c(0, int1[r,])), nrow  = 1)
      prob1[r,] <- exp_int1 / as.vector(exp_int1 %*% c(rep(1, (dim(exp_int1)[2]))))
    }
  } else {
    exp_int1 	<- matrix(exp(c(0, int1)), nrow  = 1)
    prob1 		<- exp_int1 / as.vector(exp_int1 %*% c(rep(1, (dim(exp_int1)[2]))))
  }
  return(round(prob1,4))
}

# computes intercepts from probabilities, per row of input matrix
# first catagory is reference catagory
prob_to_int <- function(prob1){
  prob1 <- prob1 + 0.00001
  b0 <- matrix(NA, nrow(prob1), ncol(prob1)-1)
  sum_exp <- numeric(nrow(prob1))
  for(r in 1:nrow(prob1)){
    sum_exp[r] <- (1/prob1[r,1]) - 1
    for(cr in 2:ncol(prob1)){
      #for every b0 except the first collumn (e.g. b012 <- log(y12/y11-y12))
      b0[r,(cr-1)] <- log(prob1[r,cr]*(1+sum_exp[r]))
    }
  }
  return(round(b0,4))
}

obtain_emiss_new <- function (object, level = "group", burn_in = NULL) 
{
  if (!mHMMbayes:::is.mHMM(object) & !mHMMbayes:::is.mHMM_cont(object) & !mHMMbayes:::is.mHMM_vary(object)) {
    stop("The input object used should either be from the class mHMM, mHMM_cont or mHMM_vary obtained by using their respective functions.")
  }
  if (level != "group" & level != "subject") {
    stop("The specification at the input variable -level- should be set to either group or subject")
  }
  input <- object$input
  dep_labels <- input$dep_labels
  n_subj <- input$n_subj
  if (is.null(burn_in)) {
    burn_in <- input$burn_in
  }
  J <- input$J
  if (burn_in >= (J - 1)) {
    stop(paste("The specified burn in period should be at least 2 points smaller\n               compared to the number of iterations J, J =", 
               J))
  }
  m <- input$m
  if (mHMMbayes:::is.mHMM(object)) {
    q_emiss <- input$q_emiss
  }
  if (mHMMbayes:::is.mHMM_vary(object)) {
    data_distr <- input$data_distr
    n_cat <- sum(data_distr %in% "categorical")
    n_cont <- sum(data_distr %in% "continuous")
    which_cat <- which(data_distr %in% "categorical")
    which_cont <- which(data_distr %in% "continuous")
    if (n_cat > 0) {
      q_emiss <- input$q_emiss
    }
  }
  n_dep <- input$n_dep
  if (level == "group") {
    if (mHMMbayes:::is.mHMM(object)) {
      est <- lapply(q_emiss, mHMMbayes:::dif_matrix, rows = m)
      for (j in 1:n_dep) {
        colnames(est[[j]]) <- paste("Category", 1:q_emiss[j])
        rownames(est[[j]]) <- paste("State", 1:m)
      }
      names(est) <- dep_labels
      for (j in 1:n_dep) {
        est[[j]][] <- mHMMbayes:::int_to_prob(matrix(apply(object$emiss_int_bar[[j]][((burn_in + 
                                                                                         1):J), ], 2, median), byrow = TRUE, 
                                                     ncol = q_emiss[j]-1, nrow = m))
      }
    }
    else if (mHMMbayes:::is.mHMM_vary(object)) {
      est <- lapply(q_emiss, mHMMbayes:::dif_matrix, rows = m)
      for (j in which_cat) {
        colnames(est[[j]]) <- paste("Category", 1:q_emiss[j])
        rownames(est[[j]]) <- paste("State", 1:m)
      }
      names(est) <- dep_labels
      if (n_cat > 0) {
        for (q in 1:n_cat) {
          ind <- which_cat[q]
          est[[ind]][] <- mHMMbayes:::int_to_prob(matrix(apply(object$emiss_int_bar[[q]][((burn_in + 
                                                                            1):J), ], 2, median), byrow = TRUE, 
                                 ncol = q_emiss[ind]-1, nrow = m))
        }
      }
      if (n_cont > 0) {
        for (q in 1:n_cont) {
          ind <- which_cont[q]
          est[[ind]] <- matrix(round(c(apply(object$emiss_mu_bar[[q]][((burn_in + 
                                                                          1):J), ], 2, median), apply(object$emiss_var_bar[[q]][((burn_in + 
                                                                                                                                    1):J), ], 2, median)), 3), ncol = 2, nrow = m, 
                               dimnames = list(paste("State", 1:m), c("Mean", 
                                                                      "Variance")))
        }
      }
    }
    else if (mHMMbayes:::is.mHMM_cont(object)) {
      est <- rep(list(matrix(, nrow = m, ncol = 2, dimnames = list(paste("State", 
                                                                         1:m), c("Mean", "Variance")))), n_dep)
      names(est) <- dep_labels
      for (j in 1:n_dep) {
        est[[j]][] <- matrix(round(c(apply(object$emiss_mu_bar[[j]][((burn_in + 
                                                                        1):J), ], 2, median), apply(object$emiss_var_bar[[j]][((burn_in + 
                                                                                                                                  1):J), ], 2, median)), 3), ncol = 2, nrow = m)
      }
    }
    est_emiss <- est
  }
  if (level == "subject") {
    est_emiss <- vector("list", n_dep)
    names(est_emiss) <- dep_labels
    if (mHMMbayes:::is.mHMM(object)) {
      for (j in 1:n_dep) {
        est <- matrix(, ncol = q_emiss[j], nrow = m, 
                      dimnames = list(paste("State", 1:m), paste("Category", 
                                                                 1:q_emiss[j])))
        est_emiss[[j]] <- rep(list(est), n_subj)
        names(est_emiss[[j]]) <- paste("Subject", 1:n_subj)
      }
      start <- c(0, q_emiss * m)
      for (i in 1:n_subj) {
        for (j in 1:n_dep) {
          est_emiss[[j]][[i]][] <- mHMMbayes:::int_to_prob(matrix(apply(object$emiss_int_subj[[i]][burn_in:J, 
                                                                          (sum(start[1:j]) + 1):sum(start[1:(j + 1)])], 
                                                      2, median), byrow = TRUE, ncol = q_emiss[j]-1, 
                                          nrow = m))
        }
      }
    }
    else if (mHMMbayes:::is.mHMM_vary(object)) {
      if (n_cat > 0) {
        start <- c(0, q_emiss * m)
        for (q in 1:n_cat) {
          ind <- which_cat[q]
          est <- matrix(, ncol = q_emiss[ind], nrow = m, 
                        dimnames = list(paste("State", 1:m), paste("Category", 
                                                                   1:q_emiss[ind])))
          est_emiss[[ind]] <- rep(list(est), n_subj)
          names(est_emiss[[ind]]) <- paste("Subject", 
                                           1:n_subj)
          for (i in 1:n_subj) {
            est_emiss[[ind]][[i]][] <- mHMMbayes:::int_to_prob(matrix(apply(object$emiss_int_subj[[i]][[q]][burn_in:J,], 2, median), byrow = TRUE, 
                                              ncol = q_emiss[ind]-1, nrow = m))
          }
        }
      }
      if (n_cont > 0) {
        for (q in 1:n_cont) {
          ind <- which_cont[q]
          est_emiss[[ind]] <- rep(list(matrix(, nrow = m, 
                                              ncol = 2, dimnames = list(paste("State", 
                                                                              1:m), c("Mean", "Variance")))), n_subj)
          names(est_emiss[[ind]]) <- paste("Subject", 
                                           1:n_subj)
          for (i in 1:n_subj) {
            est_emiss[[ind]][[i]][] <- matrix(round(c(apply(object$PD_subj[[i]]$cont_emiss[((burn_in + 
                                                                                               1):J), ((q - 1) * m + 1):((q - 1) * m + 
                                                                                                                           m)], 2, median), apply(object$PD_subj[[i]]$cont_emiss[((burn_in + 
                                                                                                                                                                                     1):J), (n_cont * m + (q - 1) * m + 1):(n_cont * 
                                                                                                                                                                                                                              m + (q - 1) * m + m)], 2, median)), 3), 
                                              ncol = 2, nrow = m)
          }
        }
      }
    }
    else if (mHMMbayes:::is.mHMM_cont(object)) {
      est_emiss <- rep(list(rep(list(matrix(, nrow = m, 
                                            ncol = 2, dimnames = list(paste("State", 1:m), 
                                                                      c("Mean", "Variance")))), n_subj)), n_dep)
      names(est_emiss) <- dep_labels
      for (j in 1:n_dep) {
        names(est_emiss[[j]]) <- paste("Subject", 1:n_subj)
        for (i in 1:n_subj) {
          est_emiss[[j]][[i]][] <- matrix(round(c(apply(object$PD_subj[[i]][((burn_in + 
                                                                                1):J), ((j - 1) * m + 1):((j - 1) * m + 
                                                                                                            m)], 2, median), apply(object$PD_subj[[i]][((burn_in + 
                                                                                                                                                           1):J), (n_dep * m + (j - 1) * m + 1):(n_dep * 
                                                                                                                                                                                                   m + (j - 1) * m + m)], 2, median)), 3), 
                                          ncol = 2, nrow = m)
        }
      }
    }
  }
  return(est_emiss)
}

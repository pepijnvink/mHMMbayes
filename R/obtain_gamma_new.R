obtain_gamma_new <- function(object, level = "group", burn_in = NULL){
    if (!inherits(object, "mHMM") & !inherits(object, "mHMM_vary")){
      stop("The input object used should be from the class mHMM, obtained by using the function mHMM.")
    }
    if (level != "group" & level != "subject"){
      stop("The specification at the input variable -level- should be set to either group or subject")
    }
    input   <- object$input
    n_subj  <- input$n_subj
    if (is.null(burn_in)){
      burn_in <- input$burn_in
    }
    J       <- input$J
    if (burn_in >= (J-1)){
      stop(paste("The specified burn in period should be at least 2 points smaller
               compared to the number of iterations J, J =", J))
    }
    m       <- input$m
    n_dep   <- input$n_dep
    est <- matrix(, ncol = m, nrow = m)
    colnames(est) <- paste("To state", 1:m)
    rownames(est) <- paste("From state", 1:m)
    if (level == "group"){
      gamma_int <- matrix(apply(object$gamma_int_bar[((burn_in + 1): J),], 2, median), byrow = TRUE, ncol = m-1, nrow = m)
      est[]     <- round(mHMMbayes:::int_to_prob(gamma_int),3)
      est_gamma <- est
    }
    if (level == "subject"){
      est_gamma <- rep(list(est), n_subj)
      names(est_gamma) <- paste("Subject", 1:n_subj)
      for(i in 1:n_subj){
        gamma_int <- matrix(apply(object$gamma_int_subj[[i]][((burn_in + 1): J),], 2, median), byrow = TRUE, ncol = m-1, nrow = m)
        est_gamma[[i]][] <- round(mHMMbayes:::int_to_prob(gamma_int),3)
      }
    }
    class(est_gamma) <- append(class(est_gamma), "mHMM_gamma")
    return(est_gamma)
}

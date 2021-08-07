#' Prepare Interval Data
#' 
#' This function can be used to prepare data for fitting
#' 
#' @param dat a data.frame with columes named er, el sl,sr
#' 

prep_interval_data <- function(dat){

    out <- list(
        N = nrow(dat),
        E_L = dat[["el"]],
        E_R = dat[["er"]],
        S_L = dat[["sl"]]
        S_R = dat[["sr"]]
    )

    return(out)

}
#' Prepare Interval Data
#'
#' This function can be used to prepare data for fitting
#'
#' @param dat a data.frame with columes named er, el sl,sr
#' @export
#' @returns a named list for use with `fit_interval_data`

prep_interval_data <- function(dat){

    out <- list(
        N = nrow(dat),
        E_L = dat[["EL"]],
        E_R = dat[["ER"]],
        S_L = dat[["SL"]],
        S_R = dat[["SR"]]
    )

    return(out)

}

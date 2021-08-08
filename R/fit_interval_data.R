#' Fitting Data for Interval Calculation
#' 
#' This function takes data from the `prep_interval_data` function
#' and fits it to a series of speficied distributions
#' 
#' @param interval_data a named list containing the interval data
#' @param distribution a string vector indicating the distributions with which
#'     to fit the data
#' @param return_samples a logical, should  the samples be returned
#' @export 
#' 

fit_interval_data <- function(interval_data, distribution = c("gamma"), return_samples = FALSE){
    staninside::copy_models("intervalcalc")

    fit_to_use <- match.arg(distribution, c("gamma"))

    local_location <- rappdirs::user_cache_dir(appname = "intervalcalc")

    mod <- cmdstanr::cmdstan_model(file.path(local_location, paste0(fit_to_use, ".stan")))

    mod$sample(data = interval_data)

    return(mod)

}
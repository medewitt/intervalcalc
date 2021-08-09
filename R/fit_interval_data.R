#' Fitting Data for Interval Calculation
#'
#' This function takes data from the `prep_interval_data` function
#' and fits it to a series of speficied distributions
#'
#' @param interval_data a named list containing the interval data
#' @param distribution a string vector indicating the distributions with which
#'     to fit the data
#' @param stan_opts a named list with additional options passed to CmdStanR
#' @param return_raw a logical, should  the samples be returned
#' @export
#'

fit_interval_data <- function(interval_data,
                              distribution = c("gamma", "gamma-trunc"),
                              stan_opts = list(),
                              return_raw = FALSE) {

    local_location <- rappdirs::user_cache_dir(appname = this_pkg())

    if (length(list.files(local_location, pattern = ".stan")) > 1) {
        cli::cli_alert_info("Using cached Stan models")
        cli::cli_alert_info(
        "Use `intervalcalc::clear_cache` if you need to refresh")
    } else {
        cli::cli_alert_info("Copying Stan models to cache")
        staninside::copy_models(this_pkg())
        cli::cli_alert_success("Models copied!")
    }


    fit_to_use <- match.arg(distribution,
                            distribution_options(),
                            several.ok = TRUE)

    out <- list()

    for (i in fit_to_use) {

        model_file_path <- file.path(local_location, paste0(i, ".stan"))
        mod <- cmdstanr::cmdstan_model(model_file_path)
        fit <- mod$sample(data = interval_data,
                        parallel_chains = stan_opts[["iter_sampling"]] %||% 2,
                        iter_sampling = stan_opts[["iter_sampling"]] %||% 1000,
                        iter_warmup = stan_opts[["iter_warmup"]] %||% 1000,
                        refresh = stan_opts[["refresh"]] %||% 250,
                        init = stan_opts[["refresh"]] %||% NULL,
                        seed  = stan_opts[["seed"]] %||% 336
                        )

        if (return_raw) {
            out_contents <- fit
        } else {

        out_contents <- list(sumz = fit$summary(c("mean_SI", "sd_SI")),
                             loo = fit$loo(r_eff = TRUE, cores = TRUE))
        }

        out[[i]] <- out_contents
    }

    return(out)
}

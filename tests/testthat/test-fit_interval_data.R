
fit_test <- function(){
	dat <- prep_interval_data(nishiura)
	fits <- fit_interval_data(interval_data = dat, distribution = "gamma",
														stan_opts = list(iter_warmup = 250,
																						 iter_sampling = 250))
	fits
}

result <- suppressWarnings(fit_test())

test_that("Returns a list", {
  expect_equal(class(result), "list")
})

test_that("A loo object is returned", {
	expect_true("loo" %in% class(result$gamma$loo))
})


test_that("Desired variables appear", {
	expect_true(all(c("mean_SI", "sd_SI") %in% result$gamma$sumz$variable))
})

test_that("Stops on invalid distribution", {
	expect_error(fit_interval_data(interval_data = mtcars, distribution = "triangle"))
})

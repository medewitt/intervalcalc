`%||%` <- function(a,b){
    if(is.null(a)) b else a
}

distribution_options <- function() {
    c("gamma", "gamma-trunc",
    "lognormal", "lognormal-trunc", 
     "weibull", "weilbull-trunc")
}

this_pkg <- function(){
    "intervalcalc"
}
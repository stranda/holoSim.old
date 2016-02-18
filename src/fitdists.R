##
##uses the fitdistr function in R to estimate parameters of distributions
## 

fit.weibull <- function(x)
    {
        require(MASS)
        suppressWarnings(fitdistr(x,"weibull",start=list(shape=1,scale=1))$estimate)
    }



fit.gamma <- function(x)
    {
        require(MASS)
        suppressWarnings(fitdistr(x,"gamma",start=list(shape=1,scale=1))$estimate)
    }

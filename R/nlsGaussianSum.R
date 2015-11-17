MAX.COMPONENTS <- 10

g <- function(x, a, s, mu) {
    a * exp(-(x - mu)^2 / (2 * s^2))
}

# Programmatically create the individual sum of gaussian functions.
# Within each iteration we build a string representation of how the function
# should look like.
for (i in seq(MAX.COMPONENTS)) {
    # Create lists of argument names, e.g. c('a2', 's2', 'mu2')
    getArgNames <- function(j) {
        sprintf(c('a%d', 's%d', 'mu%d'), j)
    }
    # Create a string of the sort 'g(x, a2, s2, mu2)'
    jComponentCall <- function(j) {
        arg.names <- getArgNames(j)
        g.call <- sprintf('g(x, %s)', paste(arg.names, collapse=', '))
        g.call
    }
    # Construct the function's arg list, e.g. c('a1', 's2', ...., 's6', 'mu6')
    all.args.str <- paste(c(sapply(seq(i), getArgNames)), collapse=', ')
    # Construct the function's body, e.g 'g1 + g2 + g2', where gj are the
    # individual calls.
    component.calls.str <- sapply(seq(i), jComponentCall)
    full.call.str <- paste(component.calls.str, collapse=' + ')
    # Parse the string and assign it to the variable sumgX
    sumgX <- eval(parse(text=sprintf('
        function(x, %s) { %s }', all.args.str, full.call.str)))
    assign(sprintf('sumg%d', i), sumgX)
}

#' Create a named vector of starting parameters for an nlsGaussiamSum
#' model with `n.components` gaussians.
#' For n.components == 2, the value of this function could look like this:
#' c(a1=40000, s1=4, mu1=30, a2=40000, s2=4, mu2=60)
#' The exact values depend on the maximum intensity `max.y` that was observed
#' in a trace and the maximum value on the horizontal axis.
#' The parameters result in `n.components` gaussian curves that are spread
#' evenly across the range of x-values.
#'
#' @param n.components The number of gaussian curves to use.
#' @param max.x The maximum value in the x-range.
#' @param max.y The maximum value in the y-range.
#' @return A named vector of parameters that can be plugged into the fitting
#' function used internally.
#'
#' @export
makeStartParams <- function(n.components, max.x, max.y) {
    nested.args <- lapply(seq(n.components), function(i) {
        a.start <- max.y
        s.start <- 4
        mu.start <- i * max.x / (n.components + 1)
        args <- c(a.start, s.start, mu.start)
        names(args) <- sprintf(c('a%d', 's%d', 'mu%d'), i)
        args
    })
    # Flatten list of named vectors into a large named vector
    do.call(c, nested.args)
}


#' Fit a nls model using multiple Gaussian functions.
#'
#' @param y A numeric vector of intensity/signal values.
#' @param x An optional vector of x values.
#' @param n The number of components to fit.
#' @param spar A smoothness parameter that is passed to smooth.spline. Used to
#'  presmooth the data.
#' @param with.plot If diagnostic plots should be made.
#' @return The fitted nls model.
#'
#' @export
nlsGaussianSum <- function(y, x=seq_along(y), n=1, spar=0.5, with.plot=FALSE) {
    # Get the right fitting function from the environment dynamically
    sumg <- get(sprintf('sumg%d', n))
    if (with.plot)
        plot(x, y, type='p', lwd=3, col='blue')
    # Fit smoothing spline and get fitted values
    ss <- smooth.spline(x, y, spar=spar)
    ss.pred <- predict(ss, seq(min(x), max(x), by=0.01))
    x <- ss.pred$x
    y <- ss.pred$y
    # How the SS fit looks like
    if (with.plot)
        lines(x, y, type='l', col='green', lwd=3)
    # Create the starting parameters s.t. the initial gaussians
    # are spread over the whole range of possible values
    start.params <- makeStartParams(n, max(x), max(y))
    # Create a string of the form "sumg2(x, a1, s1, mu1, a2, s2, mu2)"
    fit.call.string <- sprintf('sumg%d(x, %s)', n,
                               paste(names(start.params), collapse=', '))
    # Convert this string to a callable expression
    fit.call.expr <- parse(text=fit.call.string)
    print(start.params)
    # Create a y ~ sumg2(...) formula from the string that can later be passed
    # to the nls function
    fit.call.formula <- as.formula(paste0('y ~ ', fit.call.string))
    print(fit.call.formula)
    # Evaluate the fitting expression in the environment of the starting
    # parameters.
    y.fit.initial <- eval(fit.call.expr, envir=as.list(start.params))
    if (with.plot)
        lines(x, y.fit.initial, col='gray')

    # Lower and uppwer bounds for the parameter (requires algorithm == 'port')
    param.limits.lower <- rep(c(0, 3, 0), n)
    param.limits.upper <- rep(c(1.5 * max(y), 30, 80), n)
    cat('Parameter limits:\nlower: ', param.limits.lower,
        '\nupper: ',
        param.limits.upper)

    # Fit the function
    fit <- nls(fit.call.formula,
               start=start.params,
               control=nls.control(maxiter=1000, tol=1E-5, minFactor=1/1240,
                                   warnOnly=T),
               lower=param.limits.lower,
               upper=param.limits.upper,
               # trace=TRUE,
               algorithm='port')

    if (with.plot)
        lines(x, fit$m$fitted(), col='red', lwd=2)
    pars <- fit$m$getAllPars()
    cat('Final params: ', pars, '\n')
    if (with.plot) {
        for (idx in seq(3, length(pars), 3)) {
            print(idx)
            abline(v=pars[idx], lwd=3)
        }
    }

    # Return the fitted model
    fit
}

#' Choose the best number of gaussian components automatically using K-fold
#' cross validation. The best performing model is refitted and returned.
#'
#' @param y The intensity values on which to fit.
#' @param x The associated x values. Optional.
#' @param K The number of cross validation folds to perform. Larger numbers
#' yield better error estimates.
#' @param n.range A numeric vector of n values that should be evaluated.
#' @return The fitted nls model that performed best.
#'
#' @export
nlsGaussianSumCV <- function(y, x=seq_along(y), K=3, n.range=1:3,
                             ...) {
    res <- cvTuning(nlsGaussianSum, x=x, y=y,
                    tuning=list(n=n.range),
                    args=list(...),
                    K=K)
    # Refit the model using the number of components that led to the
    # smallest error as estimated by CV.
    best.n <- res$best
    nlsGaussianSum(y=y, x=x, n=best.n, ...)
}

# Example:
# res = nlsGaussianSumCV(y, K=3, n.range=1:3, with.plot=T)

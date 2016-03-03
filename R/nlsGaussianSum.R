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
nlsGaussianSum <- function(y, x=seq_along(y), n=1, spar=0.5,
                           with.plot=FALSE, quiet=TRUE, with.legend=with.plot) {
    # Get the right fitting function from the environment dynamically
    x.orig <- x
    y.orig <- y
    sumg <- get(sprintf('sumg%d', n))
    if (with.plot)
        plot(x, y, type='o', lwd=1, col='black',
             main=sprintf('Gaussian sum fit (n = %d)', n))
    # Fit smoothing spline and get fitted values
    ss <- smooth.spline(x, y, spar=spar)
    ss.pred <- predict(ss, seq(min(x), max(x), by=0.01))
    x <- ss.pred$x
    y <- ss.pred$y
    # How the SS fit looks like
    if (with.plot)
        lines(x, y, type='l', col='blue', lwd=2, lty=2)
    # Create the starting parameters s.t. the initial gaussians
    # are spread over the whole range of possible values
    start.params <- makeStartParams(n, max(x), max(y))
    # Create a string of the form "sumg2(x, a1, s1, mu1, a2, s2, mu2)"
    fit.call.string <- sprintf('sumg%d(x, %s)', n,
                               paste(names(start.params), collapse=', '))
    # Convert this string to a callable expression
    fit.call.expr <- parse(text=fit.call.string)
    # Create a y ~ sumg2(...) formula from the string that can later be passed
    # to the nls function
    fit.call.formula <- as.formula(paste0('y ~ ', fit.call.string))
    print(fit.call.formula)
    # Evaluate the fitting expression in the environment of the starting
    # parameters.
    y.fit.initial <- eval(fit.call.expr, envir=as.list(start.params))
    if (with.plot)
        lines(x, y.fit.initial, lty='dotted', col='gray')

    # Lower and uppwer bounds for the parameter (requires algorithm == 'port')
    min.intensity.factor <- 0.2 * max(y)
    min.sigma <- 3
    param.limits.lower <- rep(c(min.intensity.factor, min.sigma, 0), n)
    param.limits.upper <- rep(c(1.5 * max(y), 30, 80), n)
    if (!quiet)
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
    if (!quiet)
        cat('Final params: ', pars, '\n')
    if (with.plot) {
        fitted.apexes <- numeric()
        for (idx in seq(3, length(pars), 3)) {
            print(idx)
            fitted.apexes <- c(fitted.apexes, pars[idx])
            abline(v=pars[idx], lwd=2, col='red', lty=3)
        }
    }


    if (with.legend) {
        legend('topleft', c('smoothed protein signal',
                            'gaussian sum w/ starting parameters',
                            'gaussian sum w/ optimized parameters',
                            'gaussian component means (feature apex)',
                            'measured protein signal'),
               lty=c(2, 3, 1, 3), lwd=2,
               col=c('blue', 'gray', 'red', 'red'),
               bg='white')
        legend('topright', c('measured protein signal'),
               pch=1, lty=1, lwd=2,
               col='black',
               bg='white')
    }


    all.values <- list(
        cont=rbind(
            data.frame(x=x.orig, y=y.orig, type='measured signal', n=n),
            data.frame(x=x, y=y, type='smoothed signal', n=n),
            data.frame(x=x, y=fit$m$fitted(), type='fitted values', n=n),
            data.frame(x=x, y=y.fit.initial, type='starting values', n=n)),
        fitted.apexes=data.frame(x=fitted.apexes, n=n)
    )

    gaussian.sum.fit <-
        list(nls=fit,
             sumg=sumg,
             params=pars,
             all.values=all.values)
    class(gaussian.sum.fit) <- 'nlsGaussianSum'

    # Return the fitted model
    gaussian.sum.fit

}

#' Predict using a nlsGaussianSum model.
#' @export
predict.nlsGaussianSum <- function(fit, x) {
    do.call(fit$sumg, c(list(x), as.list(fit$params)))
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
                             quiet=TRUE, ...) {
    res <- cvTuning(nlsGaussianSum, x=x, y=y,
                    tuning=list(n=n.range),
                    cost=function(obs, pred){
                        err = sqrt(mean((obs - pred)^2))
                        if (!quiet)
                            cat('ERROR:\t', err, '\n')
                        err
                    },
                    args=c(list(quiet=quiet), list(...)),
                    K=K)
    # Refit the model using the number of components that led to the
    # smallest error as estimated by CV.
    best.n <- res$best
    nlsGaussianSum(y=y, x=x, n=best.n, quiet=quiet, ...)
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
nlsGaussianSumCV2 <- function(y, x=seq_along(y), K=3, n.range=1:3, ...) {
    folds <- cvFolds(length(y), K=K)
    errors <- matrix(0, nrow=length(n.range), ncol=K)

    for (i in 1:length(n.range)) {
        n <- n.range[i]

        for (k in 1:K) {
            is.training <- folds$which != k
            is.test <- folds$which == k

            y.train <- y[is.training]
            x.train <- x[is.training]
            y.test <- y[is.test]
            x.test <- x[is.test]

            fit <- nlsGaussianSum(y.train, n=n, with.plot=T)
            y.pred <- predict(fit, x.test)
            err <- sqrt(mean((y.pred - y.test)^2))
            errors[i, k] <- err
        }
    }

    errs <- rowMeans(errors)
    best.n <- n.range[which.min(errs)]
    nlsGaussianSum(y, x, n=best.n, with.plot=T)
}


#' Deconvolute a wide format data.table of protein traces using Gaussian
#' deconvolution.
#' @param protein.traces A wide format data.table with the columns 
#'        \itemize{
#'         \item protein_id
#'         \item 1  SEC fraction 1
#'         \item ...
#'         \item n SEC fraction n
#'        }
#'        Note that there is no complex id! These are only protein traces.
#' @return A data.table holding the Gaussian components fitted to each protein.
#' @examples
#' res <- deconvProteinTraces(protein.traces.2peps, n.cores=32)
#' @export
deconvProteinTraces <- function(protein.traces, n.cores=2, K=3, n.range=1:4) {
    require(cvTools)
    registerDoMC(cores=n.cores)
    pfeatures = foreach(i=1:nrow(protein.traces)) %dopar% {
        cat(sprintf('Current iteration: %d / %d', i, nrow(protein.traces)), '\n')
        # Convert row of DT to vector
        ptrace = as.matrix(
            subset(protein.traces[i, ], select=-protein_id)
        )[1, ]
        res <- try({
            nlsGaussianSumCV(ptrace, K=K, n.range=n.range, with.plot=T)
        })
        res
    }

    params.deconv <- lapply(pfeatures, function(res) res['params'])
    names(params.deconv) <- protein.traces$protein_id

    features.deconv <- do.call(rbind, lapply(seq_along(params.deconv), function(i) {
        if (!is.na(params.deconv[[i]])) {
            params <- params[[i]]$params
            n.feats <- length(params) / 3
            params.grouped <- split(params, rep(1:n.feats, each=3))
            params.df <- as.data.frame(do.call(rbind, params.grouped))
            colnames(params.df) <- c('a', 's', 'mu')
            params.df$protein_id <- names(params.deconv)[i]
            params.df$n_features <- n.feats
            params.df
        } else {
            NULL
        }
    }))

    features.deconv
}

# Example:
# res = nlsGaussianSumCV(y, K=3, n.range=1:5, with.plot=T)
# res = nlsGaussianSum(y, n=10, with.plot=T)

# peptraces <- widePepTracesToLong(e4.peptide.traces.wide.filtered)
# prottraces.l <- produceProteinTraces(peptraces)
# prottraces.wc <- annotateProteinTraces(prottraces.l, corum.complex.protein.assoc)
# prottraces <- longProtTracesToWide(prottraces.wc)
# prottraces

# y = as.matrix(subset(longProtTracesToWide(prottraces[protein_id ==
#                                           "P49720"]),
#                      select=-protein_id))[1, ]
# res = nlsGaussianSumCV(y, K=2, n.range=1:3, with.plot=T)
# res = nlsGaussianSum(y, n=3, with.plot=T)

# removeStretchesFilter(e4.peptide.traces.fil[protein_id == "O00267"], 4)

# prot <- 'P06576'
# setkey(prottraces)
# trace <- as.matrix(subset(unique(prottraces[protein_id == prot]),
#                           select=-c(1, 2, 3)))[1, ]


# plot(trace, type='o')

# res = nlsGaussianSumCV(trace, with.plot=T, K=3, n.range=1:4)

# protein.traces <- fread('~/Dev/cprophet/data/e4_peptides_mscore_lt_1percent_no_requant_no_decoy_wide_CC_ONLY_2_PEPS.tsv')
# setnames(protein.traces, 'protein_id', 'complex_id')
# setnames(protein.traces, 'peptide_id', 'protein_id')
# protein.traces.long <- wideProtTracesToLong(protein.traces)
# plotTraces(protein.traces.long[complex_id==192], 'protein_id')
# features.annot[complex_id == 192]


# Generate decoys

# target.complexes <- subset(corum.complex.protein.assoc, select=-complex_name)
# decoy.complexes <- data.table::copy(target.complexes)
# decoy.complexes[, complex_id := sample(complex_id)]
# decoy.complexes[, complex_id := paste0('DECOY_', complex_id)]
# decoy.complexes[, is_decoy := TRUE]
# target.complexes[, is_decoy := FALSE]
# all.complexes <- rbind(target.complexes, decoy.complexes)

#' WIP
saveFeaturesForDC <- function(features.deconv) {
    protein.features <- as.data.table(features.deconv)
    protein.features[, RT := mu]
    protein.features[, leftWidth := mu - 2 * s]
    protein.features[, rightWidth := mu + 2 * s]
    protein.features[, area := sqrt(2) * a * abs(s) * sqrt(pi)]
    protein.features[, relative_area := area / sum(area), by=protein_id]

    protein.concentrations <- fread('~/Dev/cprophet/data/e4_input_proteins_all_imputed.tsv')
    protein.features.w.conc <- 
        merge(protein.features, protein.concentrations, by='protein_id', all.x=T)
    protein.features.w.conc[, relative_concentration := relative_area *
                            protein_concentration]

    protein.features.w.conc[, feature_id := paste('feature', protein_id, 1:.N, sep='_'), by=protein_id]
    protein.features.w.conc[, total_area := sum(area), by=protein_id]

    write.table(protein.features.w.conc,
                '~/Dev/MACode/protein_features_gauss_DC_input.tsv', sep='\t',
                row.names=F)
}



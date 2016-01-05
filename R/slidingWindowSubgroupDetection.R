#' Detect subgroups within a window. This method should be called by
#' `slidingWindow`.
#'
#' @return A vector of the same length as the number of proteins. Each element
#'         corresponds to a integer-based cluster label.
#' @export
detectGroupsWithinWindow <- function(tracemat, corr.cutoff, with.plot=F) {
    if (nrow(tracemat) == 2) {
        corr <- proxy::simil(tracemat, method='correlation')[1]
        if (corr > corr.cutoff) {
            group.assignment <- c(1, 1)
        } else {
            group.assignment <- c(1, 2)
        }
        return(group.assignment)
    }
    corrmat <- cor(t(tracemat))
    # Compute distance between elements as measured by the pearson correlation,
    # i.e., dist = 2 - abs(pearson corr)
    distance <- proxy::dist(corrmat, method='correlation')
    # Cluster correlation vectors hierarchically s.t. proteins that correlate
    # well with a similar group of other proteins cluster together.
    cl <- hclust(distance)
    if (with.plot) {
        plot(cl)
        abline(h=1 - corr.cutoff, col='red')
    }
    # Cut the dendrogram at distance 0.3, i.e. pearson corr == 0.7,
    # this will give a vector of group labels.
    group.assignments <- cutree(cl, h=1 - corr.cutoff)
    group.assignments
}

#' Detect subgroups of proteins within a matrix of protein intensity traces by
#' sliding a window across the SEC dimension. Within each window, the traces
#' that correlate well are clustered together.
#'
#' @param trace.mat A numeric matrix where rows correspond to the different
#' traces.
#' @param start.window.idx Column index where the window starts.
#' @param window.size Size of the window. Numeric.
#' @param score The type of correlation to use. A string with value 'pearson'
#'        or 'diff'.
#' @return The mean correlation between all traces within this window.
#' @export
detectSubgroupsSW <- function(trace.mat,
                              protein.names,
                              corr.cutoff=0.75,
                              window.size=15,
                              with.plot=F) {
    if (nrow(trace.mat) < 2) {
        stop('Trace matrix needs at least 2 proteins')
    }

    end.idx <- ncol(trace.mat) - window.size

    groups.by.window <- sapply(seq(1, ncol(trace.mat)), function(i) {
        start.window.idx <- min(end.idx, i)
        end.window.idx <- start.window.idx + window.size
        window.trace.mat <- trace.mat[, start.window.idx:end.window.idx]
        ### Impute matrix
        # Rows with only 0 will lead to uncomputable standard deviations when
        # computing the correlation. Those rows are imputed with Norm(0, 1) noise.
        is.all.zero.row <- rowSums(window.trace.mat) == 0
        if (any(is.all.zero.row)) {
            effective.window.size <- end.window.idx - start.window.idx + 1
            n.zero.entries <- effective.window.size * sum(is.all.zero.row)
            window.trace.mat[is.all.zero.row, ] <- rnorm(n.zero.entries)
        }

        groups.within.window <-
            detectGroupsWithinWindow(window.trace.mat, corr.cutoff=corr.cutoff,
                                     with.plot=with.plot)

        groups.within.window
    })

    # Check for each SEC position if there is at least one subgroup above the
    # cutoff
    any.subgroups.in.window <- apply(groups.by.window, 2, function(col) any(col != 1:length(col)))

    groups.by.window <- apply(groups.by.window, 2, function(col) {
        lapply(unique(col), function(g) {
            protein.names[g == col]
        })
    })

    groups.dt.list <- lapply(1:length(groups.by.window), function(i) {
        # List of string vectors
        subgroups <- groups.by.window[[i]]
        subgroup.sizes <- sapply(subgroups, function(grp) length(grp))
        # A list of data.tables, each DT describing a subgroup in long list
        subgroups.dt.list <- lapply(subgroups, function(grp) {
            if (length(grp) > 1) {
                rt.dt <- data.table(sec=i,
                                    protein_id=grp,
                                    n_subunits=length(grp),
                                    subgroup=paste(grp, collapse=';'))
                # We want to report a feature RT for each position _within_ the
                # window, where the correlation was high enough. So for example
                # if the window at RT == 20 found some subgroup, then the
                # subgroup should be reported for the interval
                # [20, 20 + window.size].
                # To achieve this, we replicate the data.table and each time
                # change the RT value.
                do.call(rbind, lapply(seq(i, min(i + window.size, ncol(trace.mat))),
                   function(t) {
                       new.rt.dt <- rt.dt
                       new.rt.dt$sec <- t
                       new.rt.dt
                   }))
            } else {
                data.table(sec=integer(length=0),
                           protein_id=character(0),
                           n_subunits=integer(length=0),
                           subgroup=character(0))
            }
        })
        do.call(rbind, subgroups.dt.list)
    })
    groups.dt <- do.call(rbind, groups.dt.list)

    result <- list(subgroups.dt=groups.dt,
                   window.size=window.size,
                   corr.cutoff=corr.cutoff)
    class(result) <- 'swResult'
    result
}

#' Plot the result of the slidingWindow algorithm.
#' @export
plot.swResult <- function(sw.result, protein.traces.long, n.largest=NULL,
                          log=FALSE) {
    subgroups.dt <- sw.result$subgroups.dt
    complex.id <- unique(protein.traces.long$complex_id)
    trace.plot <- plotTraces(protein.traces.long,
                             'protein_id', 'complex_id',
                             paste('Protein traces of complex', complex.id),
                             plot=F,
                             log=log)
    if (!is.null(n.largest)) {
        n.subunits <- unique(subgroups.dt$n_subunits)
        n.subunits.ord <- n.subunits[order(n.subunits, decreasing=TRUE)]
        largest.n <- head(n.subunits.ord, n.largest)
        subgroups.dt <- subgroups.dt[n_subunits %in% largest.n]
    }

    merged.data <- merge(subgroups.dt, protein.traces.long,
                         by=c('protein_id', 'sec'))
    subgroup.plot <- ggplot(merged.data) +
        geom_point(aes(x=sec, y=protein_id, color=protein_id),
                  size=3, shape=15) +
        # geom_area(aes(x=sec, y=2, fill=protein_id), size=2, position='stack') +
        facet_wrap(~ subgroup, nrow=length(unique(subgroups.dt$subgroup))) +
        ggtitle(sprintf('Detected subgroups (window size = %d, correlation cutoff = %f)',
                         sw.result$window.size, sw.result$corr.cutoff)) +
        theme(axis.line.y=element_blank()
              # axis.text.y=element_blank(),
              # axis.ticks.y=element_blank(),
              # axis.title.y=element_blank()
              # panel.background=element_blank(),
              # panel.border=element_blank(),
              # panel.grid.major=element_blank(),
              # panel.grid.minor=element_blank(),
              # plot.background=element_blank()
              )
    p <- multiplot(trace.plot, subgroup.plot)
    print(p)
    p
}

#' Run the sliding window algorithm over a wide format data.table of protein
#' traces.
#' @export
targetedSW <- function(protein.traces.sw, corr.cutoff=0.99, window.size=15,
                       n.cores=detectCores()) {

    setkey(protein.traces.sw, complex_id)

    # Generate decoys
    protein.traces.sw.decoy <- data.table::copy(protein.traces.sw)
    protein.traces.sw.decoy[, complex_id := sample(complex_id)]
    protein.traces.sw.decoy[, complex_id := paste0('DECOY_', complex_id)]
    protein.traces.sw.decoy[, is_decoy := TRUE]
    protein.traces.sw[, is_decoy := FALSE]

    ## Concatenate target and decoy complexes
    protein.traces.sw.all <- rbind(protein.traces.sw.decoy, protein.traces.sw)
    # Throw out all complexes with less than 2 traces
    protein.traces.sw.all[, n_subunits_observed := length(protein_id), by='complex_id']
    protein.traces.sw.all <- protein.traces.sw.all[n_subunits_observed >= 2]
    protein.traces.sw.all[, n_subunits_observed := NULL]

    ## All complexes used for input
    input.complexes.sw <- unique(protein.traces.sw.all$complex_id)

    registerDoMC(n.cores)

    is.complex.validated <- logical(length(input.complexes.sw))

    sw.results <- foreach(i=seq_along(input.complexes.sw)) %dopar% {
        complex.id <- input.complexes.sw[i]
        cat('********************************************************************************\n')
        cat(sprintf('CHECKING RUN:  %d / %d', i, length(input.complexes.sw)), '\n')
        cat('********************************************************************************\n')

        # Extract the protein traces belonging to the current complex
        traces.subs <- protein.traces.sw.all[complex_id == complex.id]
        is.decoy.complex <- traces.subs$is_decoy[1]
        protein.names <- traces.subs$protein_id
        # Convert traces to a matrix
        traces.mat <- as.matrix(subset(traces.subs, select=-c(complex_id, protein_id, is_decoy)))

        # Run the algorithm
        res.sw <- detectSubgroupsSW(traces.mat, protein.names, 0.9, window.size=15)

        # Get the size of the complex as annotated in CORUM
        if (!is.decoy.complex) {
            annotated.complex.size <-
                corum.complex.protein.assoc[complex_id == complex.id, n_subunits_annotated][1]
        } else {
            annotated.complex.size <- nrow(traces.mat)
        }

        if (nrow(res.sw$subgroups.dt) > 0) {
            if (any(res.sw$subgroups.dt$n_subunits >= 0.5 * annotated.complex.size)) {
                is.complex.validated[i] <- TRUE
            } else {
                is.complex.validated[i] <- FALSE
            }
        } else {
            is.complex.validated[i] <- FALSE
        }

        cat('Complex', complex.id, 'was found to be', is.complex.validated[i], '\n')
        cat('Is manually annotated? ==>', complex.id %in% manual.annotations.final$complex_id, '\n')
        if (is.decoy.complex) {
            cat('==> DECOY\n')
        }
    }

    res <- list(sw_results=sw.results,
         input_complexes=input.complexes.sw,
         is.complex.validated=is.complex.validated,
         corr.cutoff=corr.cutoff,
         window.size=window.size)

    class(res) <- 'targetedSWResult'

    res

#     plotTraces(wideProtTracesToLong(subset(protein.traces.sw.decoy, select=-is_decoy))[complex_id == 'DECOY_10'], 'protein_id', 'complex_id')


#     true.complexes <- unique(manual.annotations.final$complex_id)
#     detected.complexes <- input.complexes.sw[is.complex.validated]

#     TP <- sum(detected.complexes %in% true.complexes)
#     FN <- length(true.complexes) - TP

#     FP <- sum(!(detected.complexes %in% true.complexes))
#     negative.complexes <- # == false targets and decoys
#         input.complexes.sw[!(input.complexes %in% true.complexes)]
#     TN <- sum(!(negative.complexes %in% detected.complexes))

#     TPR <- TP / (TP + FN)
#     FPR <- FP / (TN + FP)
}


#' Perform the targeted sliding window detection over a grid of parameters
#' @export
targetedSWGridSearch <- function(protein.traces.sw,
                                 corr.cutoffs=c(0.1, 0.5, 0.9, 0.99),
                                 window.sizes=c(5, 10, 15, 30),
                                 n.cores=detectCores()) {
    grid.size <- length(corr.cutoffs) * length(window.sizes)
    lapply(seq_along(corr.cutoffs), function(corr.cutoff) {
        corr.cutoff <- corr.cutoffs[i]
        lapply(seq_along(window.sizes), function(j) {
            window.size <- window.sizes[j]
            cat(paste(rep('+', 80), collapse=''), '\n\n')
            cat(sprintf('GRID SEARCH ITERATION: %d / %d\n',
                        (i-1) * length(corr.cutoffs) + j,
                        grid.size))
            cat(paste0(rep('+', 80), collapse=''), '\n\n')
            targetedSW(protein.traces.sw, corr.cutoff=corr.cutoff,
                       window.size=window.size,
                       n.cores=n.cores)
        })
    })
}

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

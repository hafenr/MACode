require(devtools)
require('proxy')
devtools::load_all('~/Dev/MACode')


peptraces <- widePepTracesToLong(e4.peptide.traces.wide.filtered)
prottraces <- produceProteinTraces(peptraces)
prottraces.wc <- annotateProteinTraces(prottraces, corum.complex.protein.assoc)


calcNCorrelatingTraces <- function(trace.mat, rt, window.size=14, corr.cutoff=0.7) {
    # Extract the submatrix from the full trace matrix
    window.size.half <- floor(window.size / 2)
    left.border <- max(1, rt - window.size.half)
    right.border <- min(ncol(trace.mat), rt + window.size.half)
    trace.mat.win <- trace.mat[, left.border:right.border]

    n.subunits.tot <- nrow(trace.mat)  # the number of theoretical subunits

    if (n.subunits.tot < 2) {
        return(NA)
    } else if (n.subunits.tot == 2) {
        correl <- cor(t(trace.mat.win))[1, 2]
        if (correl > corr.cutoff) {
            return(1)
        } else {
            return(NA)
        }
    }

    for (i in 1:nrow(trace.mat.win)) {
        all.zero <- all(trace.mat.win[i, ] == 0)
        if (all.zero) {
            trace.mat.win[i, ] <- rnorm(ncol(trace.mat.win))
        }
    }

    corrmat <- cor(t(trace.mat.win))
    # Compute distance between elements as measured by the pearson correlation,
    # i.e., dist = 2 - abs(pearson corr)
    distance <- proxy::dist(corrmat, method='correlation')

    # Cluster correlation vectors hierarchically s.t. proteins that correlate
    # well with a similar group of other proteins cluster together.
    cl <- hclust(distance)
    plot(cl)

    abline(h=1 - corr.cutoff, col='red')
    # Cut the dendrogram at distance 0.3, i.e. pearson corr == 0.7,
    # this will give a vector of group labels.
    group.assignments <- cutree(cl, h=1 - corr.cutoff)
    # Count how often each group label occurs.
    group.sizes <- table(group.assignments)
    max.group.size <- max(group.sizes)
    has.large.enough.subgroup <- any(group.sizes >= 0.5 * n.subunits.tot)

    return(max.group.size)
}


#' Estimate the completness of a peak group that was annotated as being
#' partially complete.
#' @param manual.annotations.partial A manual annotation data.table
#' @param protein.traces Wide format protein traces with columns:
#'   \itemize{
#'     \item complex_id
#'     \item protein_id
#'     \item complex_name
#'     \item additional columns for each sec rt
#'   }
#' @param window.size The window in SEC that is considered around each manually
#'        annotated RT.
#' @param corr.cutoff Above which value a group of traces is considered to be well
#'        correlated.
#' @examples
#' peptraces <- widePepTracesToLong(e4.peptide.traces.wide.filtered)
#' prottraces <- produceProteinTraces(peptraces)
#' prottraces.wc <- annotateProteinTraces(prottraces, corum.complex.protein.assoc)
#' imputePartialPeakgroupCompletness(manual.annotations.partial, prottraces.wc)
#' @export
imputePartialPeakgroupCompletness <- function(manual.annotations.partial,
                                              protein.traces,
                                              window.size=14,
                                              corr.cutoff=0.7) {
    manual.annotations.partial <- manual.annotations.partial  # copy
    manual.annotations.partial$complex_id <-
        as.character(manual.annotations.partial$complex_id)
    complex.ids <- unique(manual.annotations.partial[, complex_id])

    for (complex.id in complex.ids) {
        cat('Complex id:\t', complex.id, '\n')
        rts <- manual.annotations.partial[complex_id == complex.id, rt]
        traces <- protein.traces[complex_id == complex.id, ]
        traces.wide <- longProtTracesToWide(traces)
        is.sec.col <- grepl('\\d+', colnames(traces.wide))
        trace.mat.full <- as.matrix(subset(traces.wide, select=is.sec.col))
        manual.annotations.partial[
            complex_id == complex.id,
            n_subunits_annotated := nrow(trace.mat.full)
        ]

        for (rt in rts) {
            cat('RT:\t', rt, '\n')
            max.group.size <-
                calcNCorrelatingTraces(trace.mat.full, rt,
                                       window.size=window.size,
                                       corr.cutoff=corr.cutoff)
            manual.annotations.partial[
                complex_id == complex.id & rt == rt,
                max_group_size := max.group.size
            ]
        }
    }
    manual.annotations.partial
}

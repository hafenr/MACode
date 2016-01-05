diffScore <- function(trace.mat) {
    n.traces <- nrow(trace.mat)
    combination.idx.matrix <- combn(n.traces, 2)
    correls <- apply(combination.idx.matrix, 2, function(idxs) {
        trace1 <- trace.mat[idxs[1], ]
        trace2 <- trace.mat[idxs[2], ]
        sd(trace1 - trace2)
    })
    mean(correls)
}

pearsonCorrScore <- function(trace.mat) {
    correls <- cor(t(trace.mat), t(trace.mat))
    correls <- correls[upper.tri(correls)]
    mean(correls)
}

detectGroupsWithinWindow <- function(tracemat, protein.names, corr.cutoff) {
    corrmat <- cor(t(tracemat))
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
    group.assignments
}

#' Compute the correlation between rows in a matrix considering only a window.
#' This is a helper method.
#'
#' @param trace.mat A numeric matrix where rows correspond to the different
#' traces.
#' @param start.window.idx Column index where the window starts.
#' @param window.size Size of the window. Numeric.
#' @param score The type of correlation to use. A string with value 'pearson'
#'        or 'diff'.
#' @return The mean correlation between all traces within this window.
#' @export
computeWindowedCorrelation <- function(trace.mat, start.window.idx,
                                       window.size,
                                       score='pearson') {
    start.window.idx <- min(ncol(trace.mat) - window.size, start.window.idx)
    end.window.idx <- start.window.idx + window.size
    window.trace.mat <- trace.mat[, start.window.idx:end.window.idx]

    # Rows with only 0 will lead to uncomputable standard deviations when
    # computing the correlation. Those rows are imputed with Norm(0, 1) noise.
    is.all.zero.row <- rowSums(window.trace.mat) == 0
    if (any(is.all.zero.row)) {
        # window.trace.mat <- window.trace.mat[!is.all.zero.row, ]
        effective.window.size <- end.window.idx - start.window.idx + 1
        n.zero.entries <- effective.window.size * sum(is.all.zero.row)
        window.trace.mat[is.all.zero.row, ] <- rnorm(n.zero.entries)
    }

    # Compute the mean correlation between all the traces within this window.
    # The mean is computed only on the upper triangular matrix.
    correls <- cor(t(window.trace.mat), t(window.trace.mat))
    correls <- correls[upper.tri(correls)]
    if (score == 'pearson') {
        r <- pearsonCorrScore(window.trace.mat)
    } else if (score == 'diff') {
        r <- diffScore(window.trace.mat)
    } else {
        error(paste('Unknown score:', score))
    }
    # The last values that are within one window are set to the last computed
    # value.
    r
}

#' Compute the correlation between rows in a matrix using a sliding window.
#'
#' @param trace.mat A numeric matrix where rows correspond to the different
#' traces.
#' @param window.size Size of the window. Numeric.
#' @return A vector of length `ncol(trace.mat) where each element corresponds to
#' a mean correlation within the window at that position.
#' @export
slidingWindowCorrelation <- function(trace.mat, window.size, score='pearson') {
    end.index <- ncol(trace.mat) - window.size
    corr <- sapply(seq(1, end.index), function(i) {
        computeWindowedCorrelation(trace.mat, i, window.size, score=score)
    })
    all.correls <- c(corr, rep(corr[length(corr)], window.size))
    if (score == 'diff') {
        all.correls <- all.correls / max(all.correls)
    }
    all.correls
}

# protein.traces <-
#     fread('~/Dev/cprophet/data/e4_peptides_mscore_lt_1percent_no_requant_no_decoy_wide_CC.tsv')
# setnames(protein.traces, 'protein_id', 'complex_id')
# setnames(protein.traces, 'peptide_id', 'protein_id')

# protein.traces.subs <- protein.traces[complex_id == 192]
# protein.traces.subs.long <- wideProtTracesToLong(protein.traces.subs)
# plotTraces(protein.traces.subs.long,
#            'protein_id',
#            'complex_id', '')

# samplemat <- as.matrix(subset(protein.traces.subs, select=-c(protein_id, complex_id)))
# detectGroupsWithinWindow(samplemat,
#                          unique(protein.traces.subs.long$protein_id),
#                          0.7)

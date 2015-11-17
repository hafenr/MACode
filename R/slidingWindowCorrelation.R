#' Compute the correlation between rows in a matrix considering only a window.
#' This is a helper method.
#'
#' @param trace.mat A numeric matrix where rows correspond to the different
#' traces.
#' @param start.window.idx Column index where the window starts.
#' @param window.size Size of the window. Numeric.
#' @return The mean correlation between all traces within this window.
#' @export
computeWindowedCorrelation <- function(trace.mat, start.window.idx, window.size) {
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
    correls.upper <- correls[upper.tri(correls)]
    r <- mean(correls.upper)
    # The last values that are within one window are set to the last computed
    # value.
    r
}

#' Compute the correlation between rows in a matrix using a sliding window.
#'
#' @param trace.mat A numeric matrix where rows correspond to the different
#' traces.
#' @param window.size Size of the window. Numeric.
#' @return A vector of length nrow(trace.mat) where each element corresponds to
#' a mean correlation within the window at that position.
#' @export
slidingWindowCorrelation <- function(trace.mat, window.size) {
    corr <- sapply(seq(1, ncol(trace.mat)), function(i) {
        computeWindowedCorrelation(tace.mat, i, window.size)
    })
    corr
}

# protein.traces.with.complex <- produceComplexAnnotatedProteinTraces(
#     '~/Dev/MAScripts/data/corum_complex_protein_assoc.tsv',
#     '~/Dev/MAScripts/data/4_Fixed_weights_wMS1/Subsetted_OSW_Output/Peptide_Matrices/Long_Lists/4_osw_output_mscore_lt_1percent_no_requant_no_decoy_FILTERED.tsv'
# )

# protein.traces.subs <- protein.traces.with.complex[complex_id == complex.id, ]
# protein.traces.subs[, n_proteins_in_complete_complex := NULL]

# protein.traces.subs.wide <-
#     dcast(protein.traces.subs, complex_id + protein_id + complex_name ~ sec,
#           value.var='intensity')
# mat <- as.matrix(subset(protein.traces.subs.wide, select=-c(1, 2, 3)))

# window.size <- 10


# corr.long <- data.table(corr=corr, sec=seq(min(protein.traces.subs$sec),
#                                            max(protein.traces.subs$sec)),
#                         complex_id=complex.id)

# max.intensity <- max(protein.traces.subs$intensity)
# p <- ggplot(protein.traces.subs) +
#             geom_point(aes(x=sec, y=intensity, color=protein_id)) +
#             geom_line(aes(x=sec, y=intensity, color=protein_id)) +
#             geom_line(aes(x=sec, y=corr * max.intensity/2), color='black',
#                       data=corr.long) +
#             geom_abline(intercept=max.intensity/2, linetype=2)

# print(p)





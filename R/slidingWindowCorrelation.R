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
        computeWindowedCorrelation(trace.mat, i, window.size)
    })
    corr
}

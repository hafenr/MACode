removeStretchesFilter <- function(intensity.matrix, cutoff=4) {
    imat <- if (!is.matrix(intensity.matrix))
                as.matrix(intensity.matrix)
            else intensity.matrix
    .C('remove_stretches_filter',
       imat=imat, nrow(m), ncol(m), as.integer(cutoff),
       PACKAGE='MACode')$imat
}

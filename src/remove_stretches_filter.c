void remove_stretches_filter(double * imat,
                      unsigned int * p_nrows,
                      unsigned int * p_ncols,
                      unsigned int * p_strlen_cutoff) {
    unsigned int nrows = *p_nrows;
    unsigned int ncols = *p_ncols;
    unsigned int strlen_cutoff = *p_strlen_cutoff;

    unsigned int strlen;
    unsigned int strend;
    unsigned int idx, i, j, k;
    unsigned int remove_stretch;  // boolean flag

    // for each row
    for (i = 0; i < nrows; i++) {
        strlen = 0;
        remove_stretch = 0;

        for (j = 0; j < ncols; j++) {
            // R saves matrix column by column
            idx = i + j * nrows;

            // Check if there is no measurement at this position
            if (imat[idx] == 0) {
                // If there are measurements at previous adjacent positions
                // and their amount is below the cutoff, remove the stretch.
                if (strlen != 0 && strlen <= strlen_cutoff) {
                    remove_stretch = 1;
                // Otherwise, the stretch was long enough and the counter
                // should be reset.
                } else {
                    strlen = 0;
                }
            // There is a measurement at this position.
            } else {
                // Increase the stretch length
                strlen += 1;
                // Check if we're at the last iteration: if the there is a
                // short stretch that extends all the way to the end, remove it.
                if (j == (ncols - 1) && strlen <= strlen_cutoff) {
                    remove_stretch = 1;
                }
            }

            if (remove_stretch) {
                // In the case that the stretch ends at the border,
                // this will ensure that it is removed nonetheless.
                if (j == (ncols - 1)) {
                    strend = j;
                } else {
                    strend = j - 1;
                }
                k = j - strlen;
                while (k <= j) {
                    imat[i + k * nrows] = 0;
                    k += 1;
                }
                remove_stretch = 0;
                strlen = 0;
            }
        }
    }
}

#include <Rcpp.h>
using namespace Rcpp;

//' Remove consecutive stretches from a matrix of peptide traces.
//'
//' @param x A single integer.
// [[Rcpp::export]]
NumericMatrix removeStretchesFilter(NumericMatrix mat,
                                    int cutoff) {

    mat = clone(mat);

    unsigned int nrows = mat.nrow();
    unsigned int ncols = mat.ncol();
    /* unsigned int cutoff = *p_cutoff; */

    unsigned int strlen;
    unsigned int strend;
    unsigned int idx, i, j, k;
    unsigned int removestretch;  // boolean flag

    // for each row
    for (i = 0; i < nrows; i++) {
        Rcpp::checkUserInterrupt();
        strlen = 0;
        removestretch = 0;

        for (j = 0; j < ncols; j++) {
            // R saves matrix column by column
            idx = i + j * nrows;

            // Check if there is no measurement at this position
            if (mat(i, j) == 0) {
                // If there are measurements at previous adjacent positions
                // and their amount is below the cutoff, remove the stretch.
                if (strlen != 0 && strlen <= cutoff) {
                    removestretch = 1;
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
                if (j == (ncols - 1) && strlen <= cutoff) {
                    removestretch = 1;
                }
            }

            if (removestretch) {
                // In the case that the stretch ends at the border,
                // this will ensure that it is removed nonetheless.
                if (j == (ncols - 1)) {
                    strend = j;
                } else {
                    strend = j - 1;
                }
                k = j - strlen;
                while (k <= j) {
                    mat(i, k) = 0;
                    k += 1;
                }
                removestretch = 0;
                strlen = 0;
            }
        }
    }

    return mat;
}

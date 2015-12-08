#' @param traces.dt A long list style data.table holding the traces.
#' @param subunit.id.col A string giving the name of the subunit column.  
#' @param subunit.id.col A string giving the name of the parent column.  
#' @param subunit.id.col A string giving the name of the parent name column.  
#' @return The ggplot plot object
#' @examples
#' peptraces <- widePepTracesToLong(e4.peptide.traces.wide.filtered)
#' prottraces <- produceProteinTraces(peptraces)
#' prottraces.wc <- annotateProteinTraces(prottraces, corum.complex.protein.assoc)
#' plotTraces(prottraces.wc[complex_id == 635], 'protein_id', 'complex_id',
#'            'Complex 635')
plotTraces <- function(traces.dt, subunit.id.col, parent.id.col, title='') {
    p <- ggplot(traces.dt) +
                geom_point(aes_string(x='sec', y='intensity',
                                      color=subunit.id.col)) +
                geom_line(aes_string(x='sec', y='intensity',
                                     color=subunit.id.col)) +
                ggtitle(title)
    print(p)
    p
}


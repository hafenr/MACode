#' Sum up a list of peptide intensities into a list of protein intensities.
#' @param peptide.traces A long format peptide trace data.table.
#' @return A data.table of protein intensity observations with. 
#'         The DT will have the columns: 'protein_id', 'sec', 'intensity'.
#' @examples
#' peptraces <- widePepTracesToLong(e4.peptide.traces.wide.filtered)
#' prottraces <- produceProteinTraces(peptraces)
#'        
#' @export
produceProteinTraces <- function(peptide.traces) {
    # Sum peptide traces together to produce the protein traces
    protein.traces <- setnames(peptide.traces[, sum(intensity),
                               by=list(protein_id, sec)],
                      'V1', 'intensity')
    protein.traces
}

#' Annotate protein traces with the complex id and complex name to which they
#' might belong.
#'
#' @param protein.traces A long list data.table of protein traces that has the
#'        columns \itemize{
#'          \item protein_id : character
#'          \item sec        : numeric
#'          \item intensity  : numeric
#'        }
#' @param complex.assoc A datatable holding the protein <-> complex
#'        associations. It has the format:
#'        columns \itemize{
#'          \item complex_id   : character
#'          \item complex_name : character
#'          \item protein_id   : character
#'        }
#' @examples
#' prottraces <- produceProteinTraces(peptraces)
#' prottraces.wc <- annotateProteinTraces(prottraces, corum.complex.protein.assoc)
annotateProteinTraces <- function(protein.traces, complex.assoc) {
    merge(protein.traces, complex.assoc, by='protein_id', allow.cartesian=TRUE)
}


#' Convert a wide format peptide trace data.table to long format.
#' This functionc an be used to convert the included data.table
#' 'e4.peptide.traces' to long format.
#' @param dt The wide format data.table.
#' @return The same data.table in long format.
#' @export
widePepTracesToLong <- function(dt) {
    dt.long <-
        melt(dt, id.var=c('protein_id', 'peptide_id'),
             variable.name='sec', value.name='intensity',
             variable.factor=FALSE)
    dt.long$sec <- as.numeric(dt.long$sec)
    dt.long
}

#' Convert a long list peptide trace to a wide format.
#' export
longPepTracesToWide <- function(dt) {
    dcast(dt,
          protein_id + peptide_id ~ sec,
          value.var='intensity')
}

#' Convert a long list protein trace to a wide format.
#' export
longProtTracesToWide <- function(dt) {
    id.cols <- 'protein_id'
    for (maybe.id.col in c('complex_id', 'complex_name')) {
        if (maybe.id.col %in% colnames(dt))
            id.cols <- c(id.cols, maybe.id.col)
    }
    formula.str <- paste(paste(id.cols, collapse='+'), '~', 'sec')
    dcast(dt,
          as.formula(formula.str),
          value.var='intensity')
}

#' Convert a wide format protein trace data.table to long format.
#' @export
wideProtTracesToLong <- function(dt) {
    id.vars <- intersect(colnames(dt), c('protein_id', 'complex_id', 'complex_name'))
    print(id.vars)
    dt.long <- melt(dt,
         id.var=id.vars,
         variable.name='sec',
         value.name='intensity',
         variable.factor=FALSE)
    dt.long$sec <- as.numeric(dt.long$sec)
    dt.long
}

#' Convert a peptide trace matrix as is input to cprophet to the same matrix
#' that associates peptides to complexes, rather than to proteins.
#' The result of this transformation can then be used as the input for a CC
#' workflow that tries to create complex features directly from peptides.
#' @param pc.input The data.table that is normally input to cprophet. It has
#' columns: \itemize{
#'  \item peptide_id
#'  \item protein_id
#'  \item The remaining columns run across the SEC dimension.
#' }
#' @export
convertToPeptideCCData <- function(pc.input) {
    pep.prot.assoc <- pc.input[, list(peptide_id, protein_id)]
    complex.prot.assoc <- corum.complex.protein.assoc[, list(complex_id, protein_id)]
    pep.compl.assoc <-
        subset(merge(pep.prot.assoc, complex.prot.assoc, by='protein_id',
               allow.cartesian=T),
               select=-protein_id)

    pep.traces <- subset(merge(pc.input, complex.prot.assoc, by='protein_id',
                               allow.cartesian=T),
                         select=-protein_id)
    cols <- colnames(pep.traces)
    setcolorder(pep.traces, c(cols[1], cols[length(cols)],
                              cols[2:(ncol(pep.traces)-1)]))
    setnames(pep.traces, 'complex_id', 'protein_id')
    write.table(pep.traces,
                '~/Dev/cprophet/data/e4_peptides_mscore_lt_1percent_no_requant_no_decoy_wide_CC_PEPTIDES.tsv',
                sep='\t',
                row.names=F)
}


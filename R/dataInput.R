#' Combine a list of complex <-> protein associations and observations of
#' peptide intensities into a list of protein intensities that are annotated
#' with the complex they MIGHT belong to.
#' @param peptide.traces A wide format peptide trace data.table.
#' @param filename.corum.complex.assoc The file to the association TSV file.
#'        This table must have the columns: 'complex_id', 'protein_id'.
#' @return A data.table of protein intensity observations with. 
#'         The DT will have the columns: 'protein_id', 'sec', 'intensity',
#'         'complex_id', 'complex_name'.
#' @examples
#' produceComplexAnnotatedProteinTraces(e4.peptide.traces.wide.filtered,
#'                                      'corum_complex_protein_assoc.tsv')
#' @export
produceComplexAnnotatedProteinTraces <- function(peptide.traces,
                                                 filename.corum.complex.assoc) {
    peptide.traces.long <- widePepTracesToLong(peptide.traces)

    # Read corum identifiers
    corum.protein.assoc <-
        fread(filename.corum.complex.assoc, sep='\t',
              stringsAsFactors=FALSE, colClasses=c(complex_id='character'))

    # Sum peptide traces together to produce the protein traces
    protein.traces <- setnames(peptide.traces.long[, sum(intensity),
                               by=list(protein_id, sec)],
                      'V1', 'intensity')

    # Add a new column to the corum proteins that tells how many proteins
    # are theoretically in the corum complex.
    # corum.protein.assoc[
    #     ,
    #     n_proteins_in_complete_complex := length(unique(protein_id)),
    #     by=complex_id
    # ]

    # Merge protein traces together with a list of corum complexes
    # so that each intensity observation is annotated with the complex id
    # it __could potentially__ belong to.
    setkey(protein.traces, protein_id)
    setkey(corum.protein.assoc, protein_id)
    protein.traces.with.complex <-
        merge(protein.traces, corum.protein.assoc, allow.cartesian=TRUE)

    protein.traces.with.complex
}


#' Convert a wide format peptide trace data.table to long format.
#' This functionc an be used to convert the included data.table
#' 'e4.peptide.traces' to long format.
#' @param dt The wide format data.table.
#' @return The same data.table in long format.
#' @export
widePepTracesToLong <- function(dt) {
    melt(dt, id.var=c('protein_id', 'peptide_id'),
         variable.name='sec', value.name='intensity')
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
    dcast(dt,
          protein_id + complex_id + complex_name ~ sec,
          value.var='intensity')
}

#' Convert a wide format protein trace data.table to long format.
#' @export
wideProtTracesToLong <- function(dt) {
    melt(dt,
         id.var=c('protein_id', 'complex_id', 'complex_name'),
         variable.name='sec',
         value.name='intensity')
}

#' Combine a list of complex <-> protein associations and observations of
#' peptide intensities into a list of protein intensities that are annotated
#' with the complex they MIGHT belong to.
#' @param filename.corum.complex.assoc The file to the association TSV file.
#'        This table must have the columns: 'complex_id', 'protein_id'.
#' @param filename.peptide.traces.long The file to the peptide intensity
#'        This table must have the columns: 'peptide_id', 'sec', 'peptide_intensity',
#'        and 'protein_id'. 'sec' corresponds to a numeric time value and 'protein_id'
#'        to the protein that was inferred for this peptide.
#'
#' @return A data.table of protein intensity observations.
#' @export
produceComplexAnnotatedProteinTraces <- function(filename.corum.complex.assoc,
                                                 filename.peptide.traces.long) {
    # Read corum identifiers
    corum.protein.assoc <-
        fread(filename.corum.complex.assoc, sep='\t',
              stringsAsFactors=FALSE, colClasses=c(complex_id='character'))

    # Read peptide traces and produce protein trace
    peptide.traces <- fread(filename.peptide.traces.long)
    # Sum peptide traces together to produce the protein traces
    protein.traces <- setnames(peptide.traces[, sum(peptide_intensity),
                               by=list(protein_id, sec)],
                      'V1', 'intensity')

    # Add a new column to the corum proteins that tells how many proteins
    # are theoretically in the corum complex.
    corum.protein.assoc[
        ,
        n_proteins_in_complete_complex := length(unique(protein_id)),
        by=complex_id
    ]

    # Merge protein traces together with a list of corum complexes
    # so that each intensity observation is annotated with the complex id
    # it __could potentially__ belong to.
    setkey(protein.traces, protein_id)
    setkey(corum.protein.assoc, protein_id)
    protein.traces.with.complex <-
        merge(protein.traces, corum.protein.assoc, allow.cartesian=TRUE)

    protein.traces.with.complex
}

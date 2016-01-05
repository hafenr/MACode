#' Convert the output CC cprophet to a format that looks like it was output by
#' the PC workflow.
#' @export
convertToPCComplexFeatureFormat <- function(complex.features,
                                            corum.complex.protein.assoc) {

    complex.features[, complex_id := gsub('pep_', '', transition_group_id)]

    n_subunits <- sapply(strsplit(complex.features$aggr_Fragment_Annotation, ';'), length)
    complex.features[
         ,
         n_subunits :=
               length(strsplit(aggr_Fragment_Annotation, ';')),
         by=complex_id
    ]
    subunit.counts <- corum.complex.protein.assoc[, length(protein_id), by=complex_id]
    setnames(subunit.counts, 'V1', 'n_subunits_annotated')
    complex.features <- merge(complex.features, subunit.counts, by='complex_id', all.x=TRUE)

    data.table(complex_id=complex.features$complex_id,
               center_rt=complex.features$RT,
               left_boundary_rt=complex.features$leftWidth,
               right_boundary_rt=complex.features$rightWidth,
               subunit_ids=gsub(';', ',', complex.features$aggr_Fragment_Annotation),
               n_subunits=complex.features$n_subunits,
               n_subunits_annotated=complex.features$n_subunits_annotated)
}

#' Create input data wor the CC workflow.
#' All this function essentially does is to take a long list peptide trace
#' data.table, produce the protein traces, annotate them with their complex
#' membership, rename all the columns such that they are compatible with
#' cprophet.
#' @export
createCCInputData <- function(peptide.traces.long,
                              out.file.name,
                              complex.protein.assoc) {
    protein.traces.long <- produceProteinTraces(peptide.traces.long)
    protein.traces.with.complex <-
        merge(protein.traces.long, corum.complex.protein.assoc, by='protein_id',
              allow.cartesian=T)

    protein.traces.wide <- longProtTracesToWide(protein.traces.with.complex)
    protein.traces.wide <- protein.traces.wide[, complex_name := NULL]

    setnames(protein.traces.wide, 'protein_id', 'peptide_id')
    setnames(protein.traces.wide, 'complex_id', 'protein_id')

    write.table(protein.traces.wide, out.file.name, sep='\t', row.names=FALSE)

}


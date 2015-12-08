convertToPCComplexFeatureFormat <- function(complex.features,
                                            corum.complex.protein.assoc) {

    complex.features[, complex_id := gsub('pep_', '', transition_group_id)]

    n_subunits <- sapply(strsplit(complex.features$aggr_Fragment_Annotation, ';'), length)

    subunit.counts <- corum.complex.protein.assoc[, length(protein_id), by=complex_id]
    setnames(subunit.counts, 'V1', 'n_subunits_annotated')

    complex.features <- merge(complex.features, subunit.counts, by='complex_id', all.x=TRUE)

    data.table(complex_id=complex.features$complex_id,
               center_rt=complex.features$RT,
               left_boundary_rt=complex.features$leftWidth,
               right_boundary_rt=complex.features$rightWidth,
               subunit_ids=gsub(';', ',',
                                complex.features$aggr_Fragment_Annotation),
               n_subunits=n_subunits,
               n_subunits_annotated=complex.features$n_subunits_annotated)
}

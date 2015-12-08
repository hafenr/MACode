#' Filter a list of peptide intensities s.t. only those are kept that belong to
#' peptides whose protein has at least `min.peptides` peptides.
#'
#' @param peps A long format data.table of peptide intensities.
#'        The table should have the following columns:
#'        \itemize {
#'          \item protein_id (character)
#'          \item peptide_id (character)
#'          \item sec (integer)
#'          \item intensity (numeric)
#'        }
requireAtLeastXPeptides <- function(peps, min.peptides) {
    peps[, n_peptides := length(unique(peptide_id)), by='protein_id']
    peps.filtered <- peps[n_peptides >= min.peptides, ]
    subset(peps.filtered, select=-n_peptides)
}

#' Filter a list of peptide intensities s.t. only those peptides are kept that
#' belong to the `max.peptides` peptides for each protein that have the highest
#' total intensity.
#'
#' @param peps A long format data.table of peptide intensities.
#'        The table should have the following columns:
#'        \itemize {
#'          \item protein_id (character)
#'          \item peptide_id (character)
#'          \item sec (integer)
#'          \item intensity (numeric)
#'        }
retainXMostIntensePeptides <- function(peps, max.peptides) {
    # Get the peptide ids for the `max.peptides` peptides that have the highest
    # total intensity.
    takeMax <- function(pep.group) {
        order.idxs <- order(pep.group$total_peptide_intensity, decreasing=TRUE)
        take.idxs <- head(order.idxs, max.peptides)
        pep.group[take.idxs, peptide_id]
    }
    total.intensities <- peps[, sum(intensity), by=list(peptide_id, protein_id)]
    setnames(total.intensities, 'V1', 'total_peptide_intensity')
    # Apply function takeMax to the individual groups via the special variable
    # .SD (a data.table without the grouping variable protein_id).
    peptide.ids <- total.intensities[, takeMax(.SD), by='protein_id']
    peps[peptide_id %in% peptide.ids$V1]
}

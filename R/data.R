#' A wide format data.table of peptide trace intensities.
#'
#' The settings applied to the raw data are:
#' \itemize {
#'  \item m_score < 0.01
#'  \item no requant values
#'  \item no decoy values
#' }
#'
#' @format A data.table with 83266 rows and the following columns:
#' \itemize {
#'  \item protein_id: the id of the protein to which this peptide might belong
#'  \item peptide_id: the id of the peptide itself
#'  \item 3: intensity measurements for sec fraction 3
#'  \item 4: intensity measurements for sec fraction 4
#'  \item ...: further measurements
#'  \item 83: intensity measurements for sec fraction 84
#' }
'e4.peptide.traces.wide'


#' A wide format data.table of peptide trace intensities.
#'
#' The settings applied to the raw data are:
#' \itemize {
#'  \item m_score < 0.01
#'  \item no requant values
#'  \item no decoy values
#' }
#'
#' Further filters were applied:
#' \itemize {
#'  \item require_at_least_x_peptides(df, 2)
#'  \item remove_low_correlating_peptides(df)
#'  \item retain_x_most_intense_peptides(df, 10)
#' }
#'
#' @format A data.table with 28925 rows and the following columns:
#' \itemize {
#'  \item protein_id: the id of the protein to which this peptide might belong
#'  \item peptide_id: the id of the peptide itself
#'  \item 3: intensity measurements for sec fraction 3
#'  \item 4: intensity measurements for sec fraction 4
#'  \item ... further measurements
#'  \item 83: intensity measurements for sec fraction 83
#' }
'e4.peptide.traces.wide'

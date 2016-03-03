#' A wide format data.table of peptide trace intensities.
#'
#' The settings applied to the raw data are:
#' \itemize{
#'  \item m_score < 0.01
#'  \item no requant values
#'  \item no decoy values
#' }
#'
#' @format A data.table with 83266 rows and the following columns:
#' \itemize{
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
#' \itemize{
#'  \item m_score < 0.01
#'  \item no requant values
#'  \item no decoy values
#' }
#'
#' @format A data.table with 83266 rows and the following columns:
#' \itemize{
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
#' \itemize{
#'  \item m_score < 0.01
#'  \item no requant values
#'  \item no decoy values
#' }
#'
#' Further filters were applied:
#' \itemize{
#'  \item require_at_least_x_peptides(df, 2)
#'  \item remove_low_correlating_peptides(df)
#'  \item retain_x_most_intense_peptides(df, 10)
#' }
#'
#' @format A data.table with 28925 rows and the following columns:
#' \itemize{
#'  \item protein_id: the id of the protein to which this peptide might belong
#'  \item peptide_id: the id of the peptide itself
#'  \item 3: intensity measurements for sec fraction 3
#'  \item 4: intensity measurements for sec fraction 4
#'  \item ... further measurements
#'  \item 83: intensity measurements for sec fraction 83
#' }
'e4.peptide.traces.wide.filtered'

#' A data.table of associations between corum complexes and their
#' protein subunits. Complexes that have multiple possible subunit
#' compositions are split into separate complexes that are numbered
#' with `-1`, `-2`, etc. suffixes.
'corum.complex.protein.assoc'

#' The merged manual CORUM annotations with only those features that consist of
#' a peak group that was fully observed.
'manual.annotations.full'

#' The merged manual CORUM annotations with only those features that consist of
#' a peak group that was fully observed.
#' In addition this dataset holds only those fully observed peak groups that
#' are also complete in the sense that all protein subunits as annotated in
#' CORUM are actually part of the peak group.
'manual.annotations.full.complete'

#' The merged manual CORUM annotations with those features that were partially
#' observed i.e. some of the proteins that should belong to the CORUM complex
#' did not take part in the peak group.
'manual.annotations.partial'


#' The raw manual complex annotations by rhafen.
#' These annotations are based on the traces
#' summed up protein traces from 
#' e4.peptide.traces.wide.filtered.
#' A already processed and merged  version exists in the data frame
#' \itemize{
#'  \item manual.annotations.partial
#'  \item manual.annotations.full
#'  \item manual.annotations.full.complete
#' }
'manual.annotations.raw.rhafen'

#' The raw manual complex annotations by mheusel
#' These annotations are based on the traces
#' summed up protein traces from 
#' e4.peptide.traces.wide.filtered.
#' A already processed and merged  version exists in the data frame
#' \itemize{
#'  \item manual.annotations.partial
#'  \item manual.annotations.full
#'  \item manual.annotations.full.complete
#' }
'manual.annotations.raw.mheusel'

#' Manual annotations of partial peak groups that whose completeness
#' was estimated by hierchical clustering of the trace correlation profiles.  
'manual.annotations.partial.imputed.win14.corr0.7'

#' A list of all proteins that are found in the traces table
#' e4.peptide.traces.wide. Each protein has a protein_mw that is taken from
#' UniProt and a protein_concentration that is either imputed or computed by
#' aLFQ.
'e4.input.proteins'

#' Summed up protein traces using peptides from experiment 4 with the settings:
#' \itemize{
#'  \item applied corr filter on peptide level
#'  \item m-score less than 1 percent  
#'  \item no requant
#'  \item no decoy
#'  \item only 2 most intense peptides summed up
#' }
'protein.traces.2peps'


#' A data.table produced by deconvProteinTraces using the protein traces
#' protein.traces.2peps as a base.
#' Deconv settings were: K=3, n.range=1:4.
'features.deconv'

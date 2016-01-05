#' Impute concentrations for proteins based on their total intensity across the
#' SEC dimension and a list of proteins for which both SEC intensity
#' information and concentration estimates exist.
#' 
#' Furthermore, missing protein molecular weights are also added.
#'
#' @param peptide.traces A long list format data.table of peptide traces.
#' @param protein.conc.estimates A data.table holding the existing protein
#'        concentrations. It has the format:
#'        \itemize{
#'          \item protein_id
#'          \item protein_mw
#'          \item protein_concentration
#'        }
#' @param protein.ids.file A file path where the protein ids are stored for
#'        which a manual uniprot lookup has to be made.
#' @examples
#' peptide.traces <- widePepTracesToLong(e4.peptide.traces.wide)
#' protein.ids.file <- '~/Desktop/proteinmw.txt'
#' imputeProteinConcentrations(peptide.traces, HEK293.protein.conc.estimates,
#'                             protein.ids.file)
#' @return A data.table holding all proteins with their molecular weight and
#'         concentration.
imputeProteinConcentrations <- function(peptide.traces, protein.conc.estimates,
                                        protein.ids.file) {
    prot.intensity <- peptide.traces[, sum(intensity), by='protein_id']
    setnames(prot.intensity, 'V1', 'total_intensity')
    # Merge the DT holding the estimated protein concentrations and the DT holding
    # the summed up intensities across the SEC dimension.
    # The new DT shouldn't contain protein ids that aren't present in either of the DTs
    # that are to be merged (nomatch=0).
    setkey(prot.intensity, protein_id)
    setkey(protein.conc.estimates, protein_id)
    merged <- prot.intensity[protein.conc.estimates, nomatch=0]

    # Fit a line to the logged values
    fit <- lm(log(protein_concentration) ~ log(total_intensity), merged)
    prots.imputed <- subset(prot.intensity,
                            !(prot.intensity$protein_id %in% protein.conc.estimates$protein_id))
    # Inpute the concentration values for these proteins (log concentrations are converted
    # back to their original unit).
    prots.imputed[, protein_concentration := exp(predict(fit, prots.imputed))]

    # Plot how this stuff looks like
    with(merged, plot(log(total_intensity), log(protein_concentration)))

    p <- smooth.spline(log(merged$total_intensity),
                       log(merged$protein_concentration))
    fitted.vals <- predict(p, log(merged$total_intensity))
    # points(fitted.vals$x, fitted.vals$y, col='blue')
    # with(prots.imputed, points(log(total_intensity), log(protein_concentration),
    #                            col='green', pch=4))

    abline(fit, col='red', lwd=3)
    print(summary(fit))

    # For each of the proteins with an estimated conc. there is also an entry indicating
    # the proteins molecular weight. So for the proteins with imputed concentrations, we also
    # need to get those molecular weight values.
    # List from: http://www.uniprot.org/uploadlists/ by pasting output from
    proteins.without.mws <-  paste0(setdiff(prots.imputed$protein_id, 'iRT_protein'),
                                    collapse=' ')
    write(proteins.without.mws, protein.ids.file)

    cat('Now go to http://www.uniprot.org/uploadlists/ and upload the file ',
        protein.ids.file,
        '. If this doesn\'t work then try to copy paste it into the text area
        and click GO. Then customize the columns and choose "Entry" and "Mass",
        then download the file. Choose "uncompresesd" and "tab-separated file".')

    mw.file <- readline(prompt='Give the path to the file you just downloaded:  ')
    molweights <- fread(mw.file)[, list(Entry, Mass)]
    molweights <- setnames(molweights, c('protein_id', 'protein_mw'))
    # Since uniprot saves their numbers with commas in them they are read as strings and
    # we need to convert them back.
    molweights[, protein_mw := as.numeric(gsub(',', '', protein_mw))]
    # Add the molecular weight to the other proteins with imputed values
    setkey(molweights, protein_id)
    setkey(prots.imputed, protein_id)
    prots.imputed <- prots.imputed[molweights]
    # Delete unused column
    prots.imputed <- prots.imputed[, total_intensity := NULL]

    # Combine the proteins that had a concentration estimates with those that now have an
    # imputed value and write to file.
    all.proteins <- rbind(protein.conc.estimates, prots.imputed)
    all.proteins
}

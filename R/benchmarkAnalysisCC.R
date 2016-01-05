#' Calculate stats for the CC workflow.
#' @export
calculateDetectionStatsCC <- function(gs.run.directory,
                                      manual.annotations,
                                      sec.range=3:83,
                                      input.complexes=NULL,
                                      d.score.cutoffs=c(-2, -1, 0, 1, 2, 3),
                                      window.range=c(2, 4, 8, 12, 18)) {
    files <- list.files(gs.run.directory)
    iter.dirs <- files[grepl('iteration-', files)]
    iter.nums <- rep(as.integer(str_match(iter.dirs, 'iteration-(\\d+)')[, 2]),
                     each=length(d.score.cutoffs))
    iter.params <- data.table()

    site.localization.rate <-
        matrix(0, ncol=length(window.range), nrow=length(iter.dirs))
    FPR <- numeric()
    TPR <- numeric()

    for (i in seq_along(iter.dirs)) {
        cat('********************************************************************************\n')
        cat(sprintf('CHECKING RUN:  %d / %d', i, length(iter.dirs)), '\n')
        cat('********************************************************************************\n')
        d <- iter.dirs[i]
        # Save which parameters were used during this run.
        iter.num <- iter.nums[i]
        params <- fread(file.path(gs.run.directory, d, 'params.csv'))

        input.complexes <- unique(input.complexes)

        feature.file.name <-
            file.path(gs.run.directory, d, 'sec_proteins_with_dscore.csv')

        if (!file.exists(feature.file.name)) {
            TPR <- c(TPR, rep(NA, length(d.score.cutoffs)))
            FPR <- c(FPR, rep(NA, length(d.score.cutoffs)))
            next
        }

        detected.complex.features.cc <- fread(feature.file.name)
        detected.complex.features <- 
            convertToPCComplexFeatureFormat(detected.complex.features.cc,
                                            corum.complex.protein.assoc)
        detected.complex.features$d_score <- detected.complex.features.cc$d_score

        # Only complexes with >= 2 subunits should be considered 'true'
        detected.complex.features <- detected.complex.features[n_subunits >= 2, ]
        # Only complexes with a 50% completeness should pass
        # (should be all of them in the case of the PC workflow, since the pipeline
        # filters this already)
        detected.complex.features <- 
            detected.complex.features[(n_subunits / n_subunits_annotated) >= 0.5, ]
        
        true.complexes <- unique(manual.annotations$complex_id)

        for (d.score.cutoff in d.score.cutoffs) {
            params$d_score_cutoff <- d.score.cutoff
            iter.params <- rbind(iter.params, params)

            detected.complex.features.cut <-
                detected.complex.features[d_score >= d.score.cutoff]
            detected.complexes <- unique(detected.complex.features.cut$complex_id)

            TP <- sum(detected.complexes %in% true.complexes)
            FN <- length(true.complexes) - TP

            FP <- sum(!(detected.complexes %in% true.complexes))
            negative.complexes <- # == false targets and decoys
                input.complexes[!(input.complexes %in% true.complexes)]
            TN <- sum(!(negative.complexes %in% detected.complexes))

            TPR <- c(TPR, TP / (TP + FN))
            FPR <- c(FPR, FP / (TN + FP))
        }
        # Compute the site localization rate
        site.localization.rate[i, ] <- sapply(window.range, function(w) {
            computeSiteLocalizationRate(detected.complex.features,
                                        manual.annotations.full.complete,
                                        w)
        })
    }
    roc.stats <- data.table(iternum=iter.nums, d_score_cutoff=d.score.cutoffs,
                            FPR=FPR, TPR=TPR)

    list(roc.stats=roc.stats,
         site.localization.rate=site.localization.rate,
         params=iter.params)
}


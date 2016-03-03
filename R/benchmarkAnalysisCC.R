#' Calculate stats for the CC workflow.
#' @export
calculateDetectionStatsCC <- function(gs.run.directory,
                                      manual.annotations,
                                      sec.range=3:83,
                                      input.complexes.target,
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
    TPs <- numeric()
    FPs <- numeric()
    TNs <- numeric()
    FNs <- numeric()

    for (i in seq_along(iter.dirs)) {
        cat('********************************************************************************\n')
        cat(sprintf('CHECKING RUN:  %d / %d', i, length(iter.dirs)), '\n')
        cat('********************************************************************************\n')
        d <- iter.dirs[i]
        # Save which parameters were used during this run.
        iter.num <- iter.nums[i]
        params <- fread(file.path(gs.run.directory, d, 'params.csv'))

        feature.file.name <-
            file.path(gs.run.directory, d, 'sec_proteins_with_dscore.csv')
        decoy.complex.assoc.file.name <-
            file.path(gs.run.directory, d, 'decoy_proteins.tsv')
        decoy.complex.protein.assoc <- 
            fread(decoy.complex.assoc.file.name)

        input.complexes.target <- unique(input.complexes.target)
        input.complexes.decoy <- unique(decoy.complex.protein.assoc$protein_id)
        input.complexes <- c(input.complexes.target, input.complexes.decoy)

        # FILTER
        input.complexes <- FINAL.TARGETS

        if (!file.exists(feature.file.name)) {
            TPR <- c(TPR, rep(NA, length(d.score.cutoffs)))
            FPR <- c(FPR, rep(NA, length(d.score.cutoffs)))
            TPs <- c(TPs, rep(NA, length(d.score.cutoffs)))
            FPs <- c(FPs, rep(NA, length(d.score.cutoffs)))
            TNs <- c(TNs, rep(NA, length(d.score.cutoffs)))
            FNs <- c(FNs, rep(NA, length(d.score.cutoffs)))
            next
        }

        detected.complex.features.cc <- fread(feature.file.name)
        detected.complex.features <- 
            convertToPCComplexFeatureFormat(detected.complex.features.cc,
                                            decoy.complex.protein.assoc,
                                            corum.complex.protein.assoc)
        # Only complexes with >= 2 subunits should be considered 'true'
        detected.complex.features <- detected.complex.features[n_subunits >= 2, ]
        # Only complexes with a 50% completeness should pass
        detected.complex.features <- 
            detected.complex.features[(n_subunits / n_subunits_annotated) >= 0.5, ]
        
        true.complexes <- unique(manual.annotations$complex_id)

        for (d.score.cutoff in d.score.cutoffs) {
            params$d_score_cutoff <- d.score.cutoff
            iter.params <- rbind(iter.params, params)

            detected.complex.features.cut <-
                detected.complex.features[d_score >= d.score.cutoff]
            detected.complexes <- unique(detected.complex.features.cut$complex_id)

            # FILTER
            detected.complexes <-
                detected.complexes[detected.complexes %in% FINAL.TARGETS]

            TP <- sum(detected.complexes %in% true.complexes)
            FN <- length(true.complexes) - TP

            FP <- sum(!(detected.complexes %in% true.complexes))
            negative.complexes <- # == false targets and decoys
                input.complexes[!(input.complexes %in% true.complexes)]
            TN <- sum(!(negative.complexes %in% detected.complexes))

            TPR <- c(TPR, TP / (TP + FN))
            FPR <- c(FPR, FP / (TN + FP))
            TNs <- c(TNs, TN)
            FPs <- c(FPs, FP)
            TPs <- c(TPs, TP)
            FNs <- c(FNs, FN)
        }
        # Compute the site localization rate
        # site.localization.rate[i, ] <- sapply(window.range, function(w) {
        #     computeSiteLocalizationRate(detected.complex.features,
        #                                 manual.annotations.full.complete,
        #                                 w)
        # })
    }
    roc.stats <- data.table(iternum=iter.nums, d_score_cutoff=d.score.cutoffs,
                            FPR=FPR, TPR=TPR,
                            FN=FNs, TP=TPs, FP=FPs, TN=TNs)

    list(roc.stats=roc.stats,
         site.localization.rate=site.localization.rate,
         params=iter.params)
}


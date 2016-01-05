require(stringr)
require(hash)

#' Calculate stats for the PC workflow.
#' @export
calculateDetectionStatsPC <- function(gs.run.directory, manual.annotations,
                                      sec.range=3:83,
                                      # with.SEC=FALSE,
                                      window.range=c(2, 4, 8, 12, 18)) {
    files <- list.files(gs.run.directory)
    iter.dirs <- files[grepl('iteration-', files)]
    iter.params.map <- vector('list', length(iter.dirs))
    iter.nums <- as.integer(str_match(iter.dirs, 'iteration-(\\d+)')[, 2])

    site.localization.rate <-
        matrix(0, ncol=length(window.range), nrow=length(iter.dirs))
    FPR <- numeric(length=length(iter.dirs))
    TPR <- numeric(length=length(iter.dirs))

    for (i in seq_along(iter.dirs)) {
        cat('********************************************************************************\n')
        cat(sprintf('CHECKING RUN:  %d / %d', i, length(iter.dirs)), '\n')
        cat('********************************************************************************\n')
        d <- iter.dirs[i]
        # Save which parameters were used during this run.
        iter.num <- iter.nums[i]
        iter.params.map[[as.character(iter.num)]] <- fread(file.path(gs.run.directory, d, 'params.csv'))

        input.complexes.and.proteins <-
            fread(file.path(gs.run.directory, d, 'input_complexes.tsv'))
        input.complexes <- unique(input.complexes.and.proteins$complex_id)
        detected.complex.features <-
            fread(file.path(gs.run.directory, d, 'sec_complexes.tsv'))

        # Only complexes with >= 2 subunits should be considered 'true'
        detected.complex.features <- detected.complex.features[n_subunits >= 2, ]
        # Only complexes with a 50% completeness should pass
        # (should be all of them in the case of the PC workflow, since the pipeline
        # filters this already)
        detected.complex.features <- 
            detected.complex.features[(n_subunits / n_subunits_annotated) >= 0.5, ]
        
        # if (with.SEC) {
        #     annotation.result <-
        #         annotateComplexFeatures(detected.complex.features, manual.annotations)
        #     detected.complex.features <- annotation.result$annotated.features

        #     false.negatives <- annotation.result$false.negatives

        #     # The complexes that wen't into the analysis
        #     theoretically.possible.complexes <- unique(input.complexes$complex_id)

        #     TP <- sum(detected.complex.features$is_true_positive)
        #     FP <- nrow(detected.complex.features) - TP
        #     FN <- nrow(false.negatives)
        #     n.all.possible.features <-
        #         length(sec.range) * length(theoretically.possible.complexes)
        #     TN <- n.all.possible.features - TP - FP - FN

        #     TPR[i] <- TP / (TP + FN)
        #     FPR[i] <- FP / (TN + FP)
        # } else {
        true.complexes <- unique(manual.annotations$complex_id)
        detected.complexes <- unique(detected.complex.features$complex_id)

        TP <- sum(detected.complexes %in% true.complexes)
        FN <- length(true.complexes) - TP

        FP <- sum(!(detected.complexes %in% true.complexes))
        negative.complexes <- # == false targets and decoys
            input.complexes[!(input.complexes %in% true.complexes)]
        TN <- sum(!(negative.complexes %in% detected.complexes))

        TPR[i] <- TP / (TP + FN)
        FPR[i] <- FP / (TN + FP)
        # }

        # Compute the site localization rate
        site.localization.rate[i, ] <- sapply(window.range, function(w) {
            computeSiteLocalizationRate(detected.complex.features,
                                        manual.annotations,
                                        w)
        })
    }
    roc.stats <- data.table(iternum=iter.nums, FPR=FPR, TPR=TPR)

    list(roc.stats=roc.stats,
         site.localization.rate=site.localization.rate,
         params=iter.params.map)
}


#' Compute for a data table of detected features how many of those that belong
#' to a true complex also have a annotated feature that is in the vicinity of
#' them.
computeSiteLocalizationRate <- function(detected.features,
                                        manual.annotations,
                                        feature.vicinity.tol=5) {
    # Subset the features so that only features that belong to complexes
    # that have a manually annotated feature are kept.
    detected.complexes <- unique(detected.features$complex_id)
    true.complexes <- unique(manual.annotations$complex_id)
    true.pos.complexes <-
        detected.complexes[detected.complexes %in% true.complexes]
    detected.features <- detected.features[complex_id %in% true.pos.complexes, ]

    # Check for each of those complexes how many of its features actually have
    # a manual annotation nearby.
    n.close.enough.rts <- sapply(true.pos.complexes, function(cid) {
        cat('Computing site localization rate for complex ', cid, '\n')

        rt.true <- manual.annotations[complex_id == cid, ]$rt
        rt.exp <- detected.features[complex_id == cid, ]$center_rt

        has.close.enough.manual.rt <- sapply(rt.exp, function(rt) {
            any(abs((rt - rt.true)) <= feature.vicinity.tol)
        })

        sum(has.close.enough.manual.rt)
    })
    n.close.enough.rts.sum <- if (length(n.close.enough.rts) != 0) {
        sum(n.close.enough.rts)
    } else {
        0
    }
    
    rate <- n.close.enough.rts.sum / nrow(detected.features)
    rate
}


#' Plot the pseudo ROC curve for the output of `calculateDetectionStats`.
#' @export
plotROC <- function(stats, complex.centric=FALSE, title=NULL) {
    roc.stats <- stats$roc.stats
    distances.to.topleft <- apply(subset(roc.stats, select=c(FPR, TPR)), 1, function(fpr.tpr) {
        sqrt(sum((c(0, 1) - fpr.tpr)^2))
    })
    best.iter <- which.min(distances.to.topleft)

    if (is.null(title)) {
        title <- paste0(sprintf('Grid search results - %s-centric workflow',
                               if (complex.centric) 'complex' else 'protein'))
    }
    p <- ggplot(roc.stats) +
        geom_abline(slope=1) +
        geom_point(aes(x=FPR, y=TPR), alpha=0.5) +
        geom_point(aes(x=FPR, y=TPR), data=roc.stats[best.iter, ], color='red') +
        # geom_segment(aes(xend=FPR, yend=TPR), x=0, y=1, alpha=0.01) +
        geom_segment(aes(xend=FPR, yend=TPR), x=0, y=1, color='red',
                     linetype=2,
                     data=roc.stats[best.iter, ]) +
        xlim(0, 1) +
        ylim(0, 1) +
        # geom_vline(xintercept=0.1, linetype=3) +
        labs(title=title,
             x='false positive rate',
             y='true positive rate')
    print(p)
    p
}

#' Plot the site localization rate for the output of `calculateDetectionStats`.
#' @export
plotSiteLocalizationRate <- function(stats, index) {
    p <- ggplot(data.frame(localization.rate=stats$site.localization.rate[index, ],
                      window.size=c(2, 4, 8, 12, 18))) +
        geom_point(aes(x=window.size, y=localization.rate)) +
        geom_line(aes(x=window.size, y=localization.rate)) +
        labs(x='window size',
             y='ratio of features within window ') +
        ylim(c(0, 1))
    print(p)
    p
}


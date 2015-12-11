# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE


#' Check for each feature in `detected.features` if there is another feature in
#' `true.positive.features` that is close to it. The endresult is a DF of feature
#' retention times that are flagged with either 'FP', 'TP', 'TN', or 'FN'.
#' 
#' @param true.positive.features data.table of true, manually annotated features.
#'        The DF must have the columns: 'complex_id', 'rt'.
#' @param detected.features data.table of detected features.
#'        The DF must have the columns: 'complex_id', 'center_rt'.
#' @param feature.vicinity.tol A number that indicates how close an
#'        experimentally determined feature has to be to a manually annotated one, to
#'        still count as a true positive.
#'        experimentally determined feature has to be to a manually annotated one, to
#'        still count as a true positive.
#' @param all.complexes A list of unique complex ids that should be checked for
#'        features. Every SEC position for each of those
#'        complexes will be checked for a detected feature.
#' @return A data.table with the columns: 'complex_id', 'rt', 'type', where
#'         type is of type character an is either 'FP', 'FN', 'TP', 'TN'.
#' 
#' @examples 
#' manual.annotations.raw <- readManualAnnotationFile(annotations.1.raw)
#' manual.annotations <-
#'     createManualComplexAnnotations(
#'      manual.annotations.1.raw, 'apexes_partially_observed') 
#' detected.features <- fread('cprophet_output/sec_complexes.tsv')
#' assessed.feats <-
#'     assessComplexFeaturesWithSEC(manual.annotations,
#'                                  detected.features, all.complexes)
#' 
#' @export
assessComplexFeaturesWithSEC <- function(true.positive.features,
                                  detected.features,
                                  all.complexes,
                                  min.rt, max.rt,
                                  feature.vicinity.tol=5) {
    complex.ids <- unique(detected.features$complex_id)
    # Helper function to create factors for true positives etc. 

    fac <- function(t) {
        factor(t, levels=c('TN', 'FN', 'TP', 'FP'))
    }
    # All possible RTs where a feature could elute
    possible.rt <- seq(min.rt, max.rt)

    do.call(rbind, lapply(theoretically.possible.complexes,
                          function(cid) {
        # Declare integer vectors of true positives/false positives/false
        # negatives/true negatives/ that are build up in the following loops.
        TPs <- integer(0)
        FPs <- integer(0)
        FNs <- integer(0)
        TNs <- integer(0)

        rt.true <- true.positive.features[complex_id == cid, ]$rt
        rt.exp <- round(detected.features[complex_id == cid, ]$center_rt, 0)
        is.decoy.complex <- grepl('^DECOY', cid)

        features.present <- length(rt.true) != 0
        features.detected <- length(rt.exp) != 0

        # If there are no manual annotations for this complex and no features
        # were detected, all features were correctly identified as being
        # not there, i.e. negative.
        if (!features.present && !features.detected) {
            # TNs <- possible.rt
        }
        # If there are no features present, but the algorithm detected some
        # nonetheless, those are false positives.
        else if (!features.present && features.detected) {
            FPs <- rt.exp
            # TNs <- setdiff(possible.rt, rt.exp)
        }
        # If there are true features, but the algorithm didn't detect any,
        # all of those are automatically false negatives.
        else if (features.present && !features.detected) {
            FNs <- rt.true
        }
        # There are features and some may have been detected,
        # check which ones are close enough to be a TP.
        else {
            # Check for each feature RT...
            for (t.exp in rt.exp) {
                most.proximate.true.rt <- integer(0)
                smallest.delta.encountered <- Inf
                # is the experimental rt close to some of the annotated ones?
                for(t.true in rt.true) {
                    t.delta <- abs(t.true - t.exp)
                    if (t.delta <= feature.vicinity.tol
                            && t.delta < smallest.delta.encountered) {
                        # This is the closest annotated RT to the experimental one.
                        # Save it.
                        most.proximate.true.rt <- t.true
                        smallest.delta.encountered <- t.delta
                    }
                }
                # Check if there was any annotated value that could be assigned to this
                # experimental rt.
                no.corresponding.true.rt.found <- length(most.proximate.true.rt) == 0
                if (no.corresponding.true.rt.found) {
                    # No, this must be a false positive.
                    FPs <- c(FPs, t.exp)
                } else {
                    # Yes, this must be a true positive.
                    TPs <- c(TPs, t.exp)
                    # Remove the annotated vaue from the list so that it won't get
                    # assigned to another feature rt of the same complex.
                    # rt.true <- setdiff(rt.true, most.proximate.true.rt)
                }
            }
            # All those annotated RT values that were not assigned (i.e. were also not
            # removed from the original array) by the setdiff call above are by
            # definition false negatives.
            FNs <- rt.true
        }

        # The true negatives are all other RTs that aret a FP, FN or TP.
        TNs <- setdiff(possible.rt, c(TPs, FPs, FNs))

        # Build a dataframe of the feature rts and add a character indicator
        # flag to what type they belong. 
        # The if expression takes care of the situation where one of the rt
        # vectors has length 0. In that case the whole DF should have 0 rows.
        classifed.rts <- rbind(
            data.frame(rt=FNs, type=(if (length(FNs) > 0) fac('FN') else character(0))),
            data.frame(rt=TPs, type=(if (length(TPs) > 0) fac('TP') else character(0))),
            data.frame(rt=FPs, type=(if (length(FPs) > 0) fac('FP') else character(0))),
            data.frame(rt=TNs, type=(if (length(TNs) > 0) fac('TN') else character(0)))
        )

        classifed.rts$complex_id <- cid

        as.data.table(classifed.rts)
    }))
}


#' Annotate complex features using the manually annotated CORUM
#' complexes.
#' Each row will get an additional column `is_true_positive`.
#' @param detected.features A data.table holding the features as output
#' by cprophet PC.
#' @param manual.annotations A data.table with the columns 'complex_id', 'rt'.
#'        The built in data sets like `manual.annotation.full` can be used
#'        here.
#' @param feature.vicinity.tol A numeric value that specifies how near an
#'                             experimental feature has to be to a manually
#'                             annotated one to count as a true positive
#' @return The same data.table with a new column 'is_true_positive'.
annotateComplexFeatures <- function(detected.features,
                                    manual.annotations,
                                    with.sec=TRUE,
                                    feature.vicinity.tol=5) {
    n.features <- nrow(detected.features)
    feature.annotations <- detected.features
    # Remove features that have a center_rt that is within 1 unit.
    # If those features are not removed, only one of them can be a TP, all
    # others will be FP. Since each SEC number that is not TRUE is FALSE, we
    # have to evaluate RTs on an integer basis. Otherwise the calculation
    # ALL_SEC_OF_ALL_COMPLEXES - TP - FP - FN will not equal TN, since FP might
    # be arbitrarily large.
    feature.annotations[, center_rt := round(center_rt, 0)]
    # TODO: Now to something like, but this doesn't work
    setkey(feature.annotations, c('complex_id', 'center_rt'))
    feature.annotations <- unique(feature.annotations)
    false.negatives <- list()
    for (idx in seq(n.features)) {
        cat(sprintf('Comparing complex feature to manual annotations (%d / %d)', idx, n.features), '\n')
        cid <- detected.features[idx, complex_id]
        rt.true <- manual.annotations[complex_id == cid, ]$rt
        # rt.exp <- round(detected.features[complex_id == cid, ]$center_rt, 0)
        rt.exp <- detected.features[complex_id == cid, ]$center_rt

        features.present <- length(rt.true) != 0
        features.detected <- length(rt.exp) != 0

        # If there are no features present, but the algorithm detected some
        # nonetheless, those are false positives regardless of their RT.
        if (!features.present && features.detected) {
            feature.annotations[complex_id == cid, is_true_positive := FALSE]
        } else {
            # Check for each feature of this complex...
            feats <- detected.features[complex_id == cid, ]
            for (k in seq(nrow(feats))) {
                t.exp <- rt.exp[k]
                most.proximate.true.rt <- integer(0)
                smallest.delta.encountered <- Inf
                # is the experimental rt close to some of the annotated ones?
                for(t.true in rt.true) {
                    t.delta <- abs(t.true - t.exp)
                    if (t.delta <= feature.vicinity.tol
                            && t.delta < smallest.delta.encountered) {
                        # This is the closest annotated RT to the experimental one.
                        # Save it.
                        most.proximate.true.rt <- t.true
                        smallest.delta.encountered <- t.delta
                    }
                }
                # Check if there was any annotated value that could be assigned to this
                # experimental rt.
                no.corresponding.true.rt.found <- length(most.proximate.true.rt) == 0
                if (no.corresponding.true.rt.found) {
                    # No, this must be a false positive.
                    feature.annotations[complex_id == cid & center_rt == t.exp,
                                        is_true_positive := FALSE]
                } else {
                    # Yes, this must be a true positive.
                    feature.annotations[complex_id == cid & center_rt == t.exp,
                                        is_true_positive := TRUE]
                    # Remove the annotated vaue from the list so that it won't get
                    # assigned to another feature rt of the same complex.
                    rt.true <- setdiff(rt.true, most.proximate.true.rt)
                }
            }
            # If there are true RTs left that weren't assigned, then save them
            # since they are false negatives.
            if (length(rt.true) > 0) {
                false.negatives[[cid]] <- data.table(complex_id=cid, rt=rt.true)
            }
        }
        if (any(is.na(feature.annotations[complex_id == cid, is_true_positive]))) {
            browser()
        }
    }
    # Convert the list of false negatives to a DT
    if (length(false.negatives) == 0) {
        false.negatives <- data.table(complex_id=character(0),
                                      rt=numeric(length=0))
    } else {
        false.negatives <- do.call(rbind, false.negatives)
    }
    list(annotated.features=feature.annotations,
         false.negatives=false.negatives)
}

assessComplexFeatures <- function(true.positive.features,
                                  detected.features,
                                  n.all.complexes) {
    detected.features[, is_decoy := grepl('^DECOY', complex_id)]

    complex.ids.true <- unique(true.positive.features$complex_id)
    complex.ids.detected <- unique(detected.features$complex_id)
    complex.ids.detected.target <-
        unique(detected.features[is_decoy == F, complex_id])
    complex.ids.detected.decoy <-
        unique(detected.features[is_decoy == T, complex_id])

    detected.false.targets <-
        complex.ids.detected.target[
            !(complex.ids.detected.target %in% complex.ids.true)
        ]

    # TP are all manually annotated complexes that were identified
    TP <- sum(complex.ids.true %in% complex.ids.detected)
    # FN are all manually annotated complexes that were not identified
    FN <- sum(!(complex.ids.true %in% complex.ids.detected))
    # FP are all decoy complexes that were falsely detected as having a feature
    # and all target complexes that were identified buy weren't manually
    # annotated.
    FP <- sum(!(complex.ids.detected %in% complex.ids.true))
    # TN are complexes that weren't detected as P and are not FN.
    TN <- n.all.complexes - (TP + FP + FN)

    c(TP=TP, FN=FN, FP=FP, TN=TN)
}



#' Read a TSV file where each row corresponds to an annotated complex. 
#' The table must have the columns (in that order):
#' - 'complex_id',
#' - 'apexes_fully_observed'
#' - 'apexes_partially_observed'.
#' - 'n_proteins_in_complex'
#' - 'n_proteins_in_complete_complex'
#' The numbers within the apex columns have to be
#' comma-separated and should not be surrounded by whitespace.
#' 
#' @param fname The filename of the TSV file.
#' @return A list of two data.tables:
#'         The component 'annotations' holds the individual annotations.
#'         Whereas 'complexes' holds information about the complex'
#'         completeness.
#'
#' @export
readManualAnnotationFile <- function(fname) {
    annot <- fread(fname, sep='\t', sep2=',',
                   colClasses=c(rep('character', 3), rep('numeric', 2)))
    setnames(annot, c('complex_id',
                      'apexes_fully_observed',
                      'apexes_partially_observed',
                      'n_proteins_in_complex',
                      'n_proteins_in_complete_complex'))
    list(annotations=annot[, list(complex_id, apexes_fully_observed,
                                  apexes_partially_observed)],
         complexes=annot[, list(complex_id, n_proteins_in_complex,
                                n_proteins_in_complete_complex)])
}


#' Split a list of apexes like '1,2,3,4' and create a data.frame
#' data.frame(rt=c(1, 2, 3, 4), complex.id=X, apex.type=Y).
#' If no split was possible, NULL is returned.
#'
#' @param complex.id A string identifying the complex.
#' @param sep.apexes A string of comma-separated numbers.
#' @param apex.type A string identifying the type of apex.
#' @return Either NULL or a data.frame.
apexStringToDF <- function(complex.id, sep.apexes, apex.type) {
    rt <- as.numeric(strsplit(sep.apexes, ',')[[1]])
    if (length(rt) > 0) {  # if there was something to split on
        data.frame(complex_id=complex.id, rt=rt)
    }
}


#' Given a DF with columns 'complex_id', and another column holding 
#' comma-separated strings of numbers, create a long list style DF
#' where each row corresponds to a number.
#' @param annotations A data.frame with columns 'complex_id' and a second
#'                    column with comma-separated retention times of the
#'                    complex features.
#' @param apex.col.name The name of the column that holds the retention
#'                      times.
#' @return A data.table with the following columns: 'complex_id', 'rt'.
#' @examples
#' manual.annotations.raw <- readManualAnnotationFile('somefile.tsv')
#' manual.annotations <-
#'     createManualComplexAnnotations(manual.annotations.raw$annotations, 
#'                                    'apexes_partially_observed') 
#' 
#' @export
createManualComplexAnnotations <- function(annotations, apex.col.name) {
    dframes.list <- mapply(apexStringToDF, annotations$complex_id,
                           annotations[[apex.col.name]])
    # combine a list of dataframes into one large dataframe.
    # rbind will ignore entries that are NULL.
    apex.df <- as.data.table(do.call(rbind, dframes.list))
    rownames(apex.df) <- NULL
    apex.df
}


#' Merge to vectors of numbers in such a way that the output won't
#' contain numbers of the second vector that are within a interval
#' [i - window, i + window] for each number i in the first vector.
#' As an example, mergeRTs(c(1, 5), c(2, 3)) will result in c(1, 3, 5).
mergeRTs <- function(rts1, rts2, window=1) {
    ref.rts.with.spacings <- c(
        rts1,
        sapply(seq(window), function(i) {
            c(rts1 - i, rts1 + i)
        })
    )
    other.rts <- setdiff(rts2, ref.rts.with.spacings)
    merged.rts <- c(rts1, other.rts)
    merged.rts
}
stopifnot(setequal(mergeRTs(c(1, 5), c(3, 2)), c(1, 5, 3)))
stopifnot(setequal(mergeRTs(c(1, 5), c(3, 2), window=2), c(1, 5)))
stopifnot(setequal(mergeRTs(integer(0), c(3, 2)), c(3, 2)))
stopifnot(setequal(mergeRTs(c(3, 2), integer(0)), c(3, 2)))



#' Merge the RTs for apexes of `apex.type` for two DTs.
#' dt1 is treated as the reference DT.
#' Merge two data.tables where each row corresponds to manually annotated
#' complex feature.
#' The table must have the columns: 'complex_id' and 'rt'.
#' 
#' @param dt1 The first DT.
#' @param dt2 The second DT.
#' 
#' @examples
#' manual.annotations.1.raw <- readManualAnnotationFile(annotations.1.raw)
#' manual.annotations.2.raw <- readManualAnnotationFile(annotations.2.raw)
#' manual.annotations.1 <-
#'     createManualComplexAnnotations(manual.annotations.1.raw$annotations,
#'                                    'apexes_partially_observed') 
#' manual.annotations.2 <-
#'     createManualComplexAnnotations(manual.annotations.2.raw$annotations, 
#'                                    'apexes_partially_observed') 
#' manual.annotations <- mergeManualComplexAnnotations(manual.annotations.1,
#'                                                     manual.annotations.2)
#'
#' @export
mergeManualComplexAnnotations <- function(dt1, dt2) {
    complex.ids <- unique(dt1$complex_id)
    do.call(rbind, lapply(complex.ids, function(cid) {
        ref.rts <- dt1[complex_id == cid, ]$rt
        other.rts <- dt2[complex_id == cid, ]$rt
        merged <- mergeRTs(ref.rts, other.rts)
        if (length(merged) > 0) {
            data.table(complex_id=cid, rt=merged)
        } else {
            NULL
        }
    }))
}


#' Calculate the true positive rate (== recall)
calcTPR <- function(dt) {
    fn <- dt[type == 'FN', .N]
    tp <- dt[type == 'TP', .N]
    tp / (tp + fn)  # TP / P 
}


#' Calculate the false positive rate (== 1 - specificity)
calcFPR <- function(dt) {
    fp <- dt[type == 'FP', .N]
    tn <- dt[type == 'TN', .N]
    fp / (tn + fp)  # FP / N
}


#' @param detected.features A dataframe with the columns: 'complex_id', 'rt'
#' @param true.positive.features A dataframe with the columns: 'complex_id', 'rt'
#' @param cutoffs A numeric vector of apex_mw_fit values that should be used to
#'                define if a feature passes the molecular weight test or not
#' @param feature.vicinity.tol A numeric value that specifies how near an
#'                             experimental feature has to be to a manually
#'                             annotated one to count as a true positive
#' @export
makeROCWithSEC <- function(detected.features,
                           true.positive.features,
                           cutoffs,
                           feature.vicinity.tol=5) {
    fpr <- numeric(length=length(cutoffs))
    tpr <- numeric(length=length(cutoffs))
    for (i in seq_along(cutoffs)) {
        cat(sprintf('Calculating TPR/FPR. Iteration: %d\n', i))
        cval <- cutoffs[i]
        detected.features.filtered <-
            detected.features[apex_apmw_fit < cval, list(complex_id, center_rt)]
        if (nrow(detected.features.filtered) != 0) {
            all.features <- assessComplexFeaturesWithSEC(
                true.positive.features=true.positive.features,
                detected.features=detected.features.filtered,
                all.complexes,
                min.rt=3,
                max.rt=84,
                feature.vicinity.tol=5
            )
            fpr[i] <- calcFPR(all.features)
            tpr[i] <- calcTPR(all.features)
        } else {
            fpr[i] <- 0
            tpr[i] <- 0
        }
    }
    fpr.order <- order(fpr)

    data.frame(FPR=fpr[fpr.order], TPR=tpr[fpr.order], cutoff=cutoffs[fpr.order])
}


#' @param detected.features A dataframe with the columns: 'complex_id', 'rt'
#' @param true.positive.features A dataframe with the columns: 'complex_id', 'rt'
#' @param params A numeric vector of apex_mw_fit values that should be used to
#'                define if a feature passes the molecular weight test or not
#' @param n.all.complexes Number of complexes that were inputto cprophet's 
#'        complex detection (this includes decoys).
#' @param subsetFeaturesFunc A function that should receive a feature DT and a
#'        parameter value. It should return the same DT with some rows removed.
#' @export
makeROC <- function(detected.features,
                    true.positive.features,
                    params,
                    n.all.complexes,
                    subsetFeaturesFunc) {
    fpr <- numeric(length=length(params))
    tpr <- numeric(length=length(params))
    for (i in seq_along(params)) {
        cat(sprintf('Calculating TPR/FPR. Iteration: %d\n', i))
        p <- params[i]
        detected.features.filtered <- subsetFeaturesFunc(detected.features, p)
        detected.features.filtered <-
            detected.features.filtered[, list(complex_id, center_rt)]
        if (nrow(detected.features.filtered) != 0) {
            counts <- assessComplexFeatures(
                true.positive.features=true.positive.features,
                detected.features=detected.features.filtered,
                n.all.complexes
            )
            fpr[i] <- counts['FP'] / (counts['TN'] + counts['FP'])
            tpr[i] <- counts['TP'] / (counts['TP'] + counts['FN'])
        } else {
            fpr[i] <- 0
            tpr[i] <- 0
        }
    }
    fpr.order <- order(fpr)

    data.frame(FPR=fpr[fpr.order], TPR=tpr[fpr.order], param=params[fpr.order])
}


#' Plot a ROC curve.
#' @param df A dataframe as returned by makeROCWithSEC
#' @param dynamic.axis If the axis should be adjusted s.t. only visible values
#'                     are shown.
#' @export
plotROC <- function(df, dynamic.axis=TRUE) {
    p <- ggplot(df) +
        geom_point(aes(x=FPR, y=TPR, size=param), alpha=0.5) +
        geom_line(aes(x=FPR, y=TPR)) +
        geom_abline(slop=1, linetype='dashed') +
        xlab('FPR | (1 - specificity)') +
        ylab('TPR | sensitivity')
    if (!dynamic.axis) {
        p <- p + scale_x_continuous(limits=c(0, 1)) +
                 scale_y_continuous(limits=c(0, 1))
    }
    print(p)
    p
}

#' Annotate complex features by simply checking if there is a manual annotation
#' nearby. This function is faster than `annotateComplexFeatures`, which loops
#' through all the annotations and ensures that only one detected feature can
#' be assigned to one manually annotated feature. 
#' @export
annotateComplexFeaturesAlt <- function(detected.features,
                                       manual.annotations,
                                       feature.vicinity.tol=5) {
    # Check for each of those complexes how many of its features actually have
    # a manual annotation nearby.
    detected.features <- detected.features
    n.close.enough.rts <- sapply(unique(detected.features$complex_id), function(cid) {
        cat('Comparing features to manual annotations. Current complex: ', cid, '\n')

        rt.true <- manual.annotations[complex_id == cid, ]$rt
        rt.exp <- detected.features[complex_id == cid, ]$center_rt

        for (rt in rt.exp) {
            has.close.enough.manual.rt <- any(abs((rt - rt.true)) <= feature.vicinity.tol)
            detected.features[complex_id == cid & center_rt == rt,
                              is_true_positive := has.close.enough.manual.rt]
        }
    })
    detected.features
}


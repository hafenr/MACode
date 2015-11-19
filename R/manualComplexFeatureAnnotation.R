# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE


#' Check for each feature in `detected.features` if there is another feature in
#' `true.positive.features` that is close to it. The endresult is a DF of feature
#' retention times that are flagged with either 'FP', 'TP', 'TN', or 'FN'.
#' 
#' @param true.positive.features data.table of true, manually annotated features.
#' The DF must have the columns: 'complex_id', 'rt'.
#' @param detected.features data.table of detected features.
#' The DF must have the columns: 'complex_id', 'center_rt'.
#' @param feature.vicinity.tol A number that indicates how close an
#' experimentally determined feature has to be to a manually annotated one, to
#' still count as a true positive.
#' @return A data.table with the columns: 'complex_id', 'rt', 'type', where
#' type is of type character an is either 'FP', 'FN', 'TP', 'TN'.
#' 
#' @examples 
#' manual.annotations.raw <- readManualAnnotationFile(annotations.1.raw)
#' manual.annotations <-
#'     createManualComplexAnnotations(manual.annotations.1.raw, 'apexes_partially_observed') 
#' detected.features <- fread('cprophet_output/sec_complexes.tsv')
#' assessed.feats <- assessComplexFeatures(manual.annotations, detected.features)
#' 
#' @export
assessComplexFeatures <- function(true.positive.features,
                                  detected.features,
                                  min.rt, max.rt,
                                  feature.vicinity.tol=5) {
    complex.ids <- unique(detected.features$complex_id)
    # Helper function to create factors for true positives etc. 
    fac <- function(t) {
        factor(t, levels=c('TN', 'FN', 'TP', 'FP'))
    }
    # All possible RTs where a feature could elute
    possible.rt <- seq(min.rt, max.rt)

    do.call(rbind, lapply(complex.ids, function(cid) {
        # Declare integer vectors of true positives/false positives/false
        # negatives/true negatives/ that are build up in the following loops.
        TPs <- integer(0)
        FPs <- integer(0)
        FNs <- integer(0)
        TNs <- integer(0)

        rt.true <- true.positive.features[complex_id == cid, ]$rt
        rt.exp <- round(detected.features[complex_id == cid, ]$center_rt, 0)
        is.decoy.complex <- grepl('^DECOY', cid)

        # If this complex is a decoy complex, all of the detected features are
        # automatically false positives.
        if (is.decoy.complex) {
            FPs <- c(FPs, rt.exp)
        # Not a decoy complex, need to check each detected feature separately.
        } else {
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
            TNs <- setdiff(possible.rt, c(TPs, FPs, FNs))
            # The true negatives would be all theoretical complexes that aren't in TP.
        }
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


#' Read a TSV file where each row corresponds to an annotated complex. 
#' The table must have the columns: 'complex_id', 'apexes_fully_observed',
#' 'apexes_partially_observed'. The numbers within the apex columns have to be
#' comma-separated and should not be surrounded by whitespace.
#' 
#' @param fname The filename of the TSV file.
#' @return A data.table.
#'
#' @export
readManualAnnotationFile <- function(fname) {
    annot <- fread(fname, sep='\t', sep2=',', colClasses=rep('character', 3))
    setnames(annot, c('complex_id', 'apexes_fully_observed',
                      'apexes_partially_observed'))
    annot
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
#' @returns A data.table with the following columns: 'complex_id', 'rt'.
#' @examples
#' manual.annotations.raw <- readManualAnnotationFile('somefile.tsv)
#' manual.annotations <-
#'     createManualComplexAnnotations(manual.annotations, 'apexes_partially_observed') 
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
#'     createManualComplexAnnotations(manual.annotations.1.raw, 'apexes_partially_observed') 
#' manual.annotations.2 <-
#'     createManualComplexAnnotations(manual.annotations.2.raw, 'apexes_partially_observed') 
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
makeROC <- function(detected.features,
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
            all.features <- assessComplexFeatures(
                true.positive.features=true.positive.features,
                detected.features=detected.features.filtered,
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

#' Plot a ROC curve.
#' @param df A dataframe as returned by makeROC
#' @param dynamic.axis If the axis should be adjusted s.t. only visible values
#'                     are shown.
#' @export
plotROC <- function(df, dynamic.axis=TRUE) {
    p <- ggplot(res) +
        geom_point(aes(x=FPR, y=TPR, size=cutoff), alpha=0.5) +
        geom_line(aes(x=FPR, y=TPR)) +
        geom_abline(slop=1, linetype='dashed') +
        xlab('FPR | (1 - specificity)') +
        ylab('TPR | sensitivity')
    if (!dynamic.axis) {
        p <- p +
            scale_x_continuous(limits=c(0, 1)) +
            scale_y_continuous(limits=c(0, 1))
    }
    print(p)
    p
}

# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Check for each feature in `detected.features` if there is another feature in
#' `true.features` that is close to it. The endresult is a DF of feature
#' retention times that are flagged with either 'FP', 'TP', or 'FN'.
#' 
#' @param true.features data.table of true, manually annotated features.
#' The DF must have the columns: 'complex_id', 'rt'.
#' @param detected.features data.table of detected features.
#' The DF must have the columns: 'complex_id', 'rt'.
#' @param feature.vicinity.tol A number that indicates how close an
#' experimentally determined feature has to be to a manually annotated one, to
#' still count as a true positive.
#' @return A data.table with the columns: 'complex_id', 'rt', 'type', where
#' type is of type character an is either 'FP', 'FN', 'TP'.
#' 
#' @export
assessComplexFeatures <- function(true.features, detected.features,
                                  feature.vicinity.tol=5) {
    complex.ids <- unique(detected.features$complex_id)
    do.call(rbind, lapply(complex.ids, function(cid) {
        rt.true <- true.features[complex_id == cid, ]$rt
        rt.exp <- detected.features[complex_id == cid, ]$center_rt
        # Declare integer vectors of true positives/false positives/false
        # negatives that are build up in the following loops.
        TPs <- integer(0)
        FPs <- integer(0)
        FNs <- integer(0)
        # Check for each experimentally found feature, i.e. RT...
        for (t.exp in rt.exp) {
            most.proximate.true.rt <- integer(0)
            smallest.delta.encountered <- Inf
            # Is the experimental rt close to some of the annotated ones?
            for(t.true in rt.true) {
                t.delta <- abs(t.true - t.exp)
                if (t.delta <= feature.vicinity.tol
                        && t.delta < smallest.delta.encountered) {
                    # This is the closest annotated RT to the experimental one.
                    # Save it.
                    most.proximate.true.rt <- t.true
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
                rt.true <- setdiff(rt.true, most.proximate.true.rt)
            }
        }
        # All those annotated RT values that were not assigned (i.e. were also not
        # removed from the original array) by the setdiff call above are by
        # definition false negatives.
        FNs <- rt.true
        # The true negatives would be all theoretical complexes that aren't in TP.

        # Build a dataframe of the feature rts and add a character indicator
        # flag to what type they belong. 
        classifed.rts <- rbind(
            data.frame(rt=FNs, type=(if (length(FNs) > 0) 'FN' else character(0))),
            data.frame(rt=TPs, type=(if (length(TPs) > 0) 'TP' else character(0))),
            data.frame(rt=FPs, type=(if (length(FPs) > 0) 'FP' else character(0)))
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
        data.frame(complex_id=complex.id, rt=rt,
                   apex_type=apex.type,
                   stringsAsFactors=F)
    }
}

#' Given a DF with columns 'complex_id', and another column holding 
#' comma-separated strings of numbers, create a long list style DF
#' where each row corresponds to a number.
createApexDF <- function(annotations, apex.col.name) {
    dframes.list <- mapply(apexStringToDF, annotations$complex_id,
                           annotations[[apex.col.name]], apex.col.name)
    # combine a list of dataframes into one large dataframe.
    # rbind will ignore entries that are NULL.
    apex.df <- do.call(rbind, dframes.list)
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

#' Merge the RTs for apexes of `apex.type` for two DTs.
#' dt1 is treated as the reference DT.
#' 
#' @export
createMergedList <- function(dt1, dt2, apex.type) {
    complex.ids <- unique(dt1$complex_id)
    do.call(rbind, lapply(complex.ids, function(cid) {
        ref.rts <- dt1[complex_id == cid & apex_type == apex.type, ]$rt
        other.rts <- dt2[complex_id == cid & apex_type == apex.type, ]$rt
        merged <- mergeRTs(ref.rts, other.rts)
        if (length(merged) > 0) {
            data.table(complex_id=cid, rt=merged, apex_type=apex.type)
        } else {
            NULL
        }
    }))
}

#' @export
mergeManualComplexAnnotations <- function(annotations.fname.1,
                                          apex.col.name.1,
                                          annotations.fname.2,
                                          apex.col.name.2) {
    annotations.1 <- readManualAnnotationFile(annotations.fname.1)
    annotations.2 <- readManualAnnotationFile(annotations.fname.2)

    apex.df.1 <- createApexDF(annotations.1, apex.col.name.1)
    apex.df.2 <- createApexDF(annotations.2, apex.col.name.2)
    apex.dt.1 <- as.data.table(apex.df.1, key='complex_id')
    apex.dt.2 <- as.data.table(apex.df.2, key='complex_id')

    apex.dt.merged <- rbind(
        createMergedList(apex.dt.1, apex.dt.2, 'apexes_fully_observed'),
        createMergedList(apex.dt.1, apex.dt.2, 'apexes_partially_observed')
    )

    apex.dt.merged
}

stopifnot(setequal(mergeRTs(c(1, 5), c(3, 2)), c(1, 5, 3)))
stopifnot(setequal(mergeRTs(c(1, 5), c(3, 2), window=2), c(1, 5)))
stopifnot(setequal(mergeRTs(integer(0), c(3, 2)), c(3, 2)))
stopifnot(setequal(mergeRTs(c(3, 2), integer(0)), c(3, 2)))

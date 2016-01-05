#' Plot traces for a single complex or protein.
#'
#' @param traces.dt A long list style data.table holding the traces.
#' @param subunit.id.col A string giving the name of the subunit column.  
#' @param subunit.id.col A string giving the name of the parent column.  
#' @param subunit.id.col A string giving the name of the parent name column.  
#' @return The ggplot plot object
#' @examples
#' peptraces <- widePepTracesToLong(e4.peptide.traces.wide.filtered)
#' prottraces <- produceProteinTraces(peptraces)
#' prottraces.wc <- annotateProteinTraces(prottraces, corum.complex.protein.assoc)
#' plotTraces(prottraces.wc[complex_id == 635], 'protein_id', 'complex_id',
#'            'Complex 635')
#' @export
plotTraces <- function(traces.dt, subunit.id.col, parent.id.col, title='',
                       plot=T, log=F) {
    p <- ggplot(traces.dt) +
                geom_point(aes_string(x='sec', y='intensity',
                                      color=subunit.id.col)) +
                geom_line(aes_string(x='sec', y='intensity',
                                     color=subunit.id.col)) +
                ggtitle(title) +
                xlab('SEC fraction') +
                ylab('intensity')
    if (log) {
        p <- p + scale_y_log10('log(intensity)')
    }
    if (plot) print(p)
    p
}

#' Plot traces and features for CC wf.
#'
#' @param detected.features.cc A datatable of the output from openms.
#' @param complex.assoc The complex <-> protein associations.
#' @param traces.cc The data.file that was input into the cprophet CC wf.
#' @param pdf.file.loc The location where to store the pdf.
#' @export
plotTracesWithFeaturesCC <- function(detected.features.cc,
                                     complex.assoc=corum.complex.protein.assoc,
                                     traces.cc,
                                     pdf.file.loc,
                                     nplots=NULL) {
    ### Create trace plots with features
    complex.features.cc <- convertToPCComplexFeatureFormat(detected.features.cc,
                                                           complex.assoc)
    #  Copy all score columns back into the complex features data.frame
    score.column.names <- colnames(detected.features.cc)[
        grep('var_', colnames(detected.features.cc))
    ]
    for (col in c(score.column.names, 'd_score')) {
        complex.features.cc[, col := detected.features.cc[[col]], with=F]
    }
    complex.features.cc <- 
        complex.features.cc[order(complex.features.cc$d_score, decreasing=T), ]
    traces.cc <- widePepTracesToLong(traces.cc)
    setnames(traces.cc, 'protein_id', 'complex_id')
    setnames(traces.cc, 'peptide_id', 'protein_id')
    traces.cc <- subset(traces.cc, select=-complex_id)
    setkey(traces.cc)
    traces.cc <- unique(traces.cc)

    traces <- traces.cc
    # Need the information which decoy belongs to which protein in the CC
    # workflow!!!

    # prots <- fread(file.path(gs.run.directory, 'iteration-0000/sec_proteins.csv'))
    # strsplit(prots$aggr_Fragment_Annotation, ';')[[1]]
    # lappl

    # Decoys dropped for now
    complex.features.cc <- complex.features.cc[!grepl('DECOY', complex_id)]

    nplots <- if (is.null(nplots)) nrow(complex.features.cc) else nplots
    pdf(pdf.file.loc)
    for (i in 1:nplots) {
        feat <- complex.features.cc[i, ]
        cat('DECOY: ', feat$is_decoy == T, '\n')
        cat('D_SCORE: ', feat$d_score, '\n')
        subunits <- strsplit(feat$subunit_ids, ',')[[1]]
        trace.data <- traces[protein_id %in% subunits, ]
        complex.id <- feat$complex_id
        feat$max_intensity <- max(trace.data$intensity)

        p <- ggplot(aes(x=sec, y=intensity, color=protein_id),
               data=trace.data) +
            geom_line() +
            geom_point() +
            geom_vline(xintercept=feat$left_boundary_rt) +
            geom_vline(xintercept=feat$center_rt, linetype=3) +
            geom_vline(xintercept=feat$right_boundary_rt) +
            ggtitle(complex.id)

        max.sec <- max(trace.data$sec)
        max.intensity <- max(trace.data$intensity)
        p <- p + annotate('text', x=0.88*max.sec, y=0.96*max.intensity,
                          label=round(feat$d_score, 2), size=10)
        # p <- p + annotate('text', x=0.88*max.sec, y=0.86*max.intensity,
        #                   label='TRUE' if feat$is_true_positive else 'FALSE',
        #                   size=8)
        print(p)
        # readline(prompt = "Pause. Press <Enter> to continue...")
    }
    dev.off()
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

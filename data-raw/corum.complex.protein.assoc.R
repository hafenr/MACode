corum.complex.protein.assoc <-
    fread('data-raw/corum_complex_protein_assoc.tsv',
          stringsAsFactors=F, colClasses=c(complex_id='character'))
devtools::use_data(corum.complex.protein.assoc, overwrite=TRUE)

HEK293.protein.conc.estimates <-
    fread('data-raw/HEK293_protein_conc_estimates.csv',
          stringsAsFactors=F)
devtools::use_data(HEK293.protein.conc.estimates, overwrite=TRUE)

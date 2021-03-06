% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/dataConversion.R
\name{annotateProteinTraces}
\alias{annotateProteinTraces}
\title{Annotate protein traces with the complex id and complex name to which they
might belong.}
\usage{
annotateProteinTraces(protein.traces, complex.assoc)
}
\arguments{
\item{protein.traces}{A long list data.table of protein traces that has the
columns \itemize{
  \item protein_id : character
  \item sec        : numeric
  \item intensity  : numeric
}}

\item{complex.assoc}{A datatable holding the protein <-> complex
associations. It has the format:
columns \itemize{
  \item complex_id   : character
  \item complex_name : character
  \item protein_id   : character
}}
}
\description{
Annotate protein traces with the complex id and complex name to which they
might belong.
}
\examples{
prottraces <- produceProteinTraces(peptraces)
prottraces.wc <- annotateProteinTraces(prottraces, corum.complex.protein.assoc)
}


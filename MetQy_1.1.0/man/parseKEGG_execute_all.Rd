% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parseKEGG_execute_all.R
\name{parseKEGG_execute_all}
\alias{parseKEGG_execute_all}
\title{Execute all parseKEGG parent functions to format KEGG databases into data frames}
\usage{
parseKEGG_execute_all(KEGG_path, ...)
}
\arguments{
\item{KEGG_path}{- string pointing to the location of the KEGG database parent folder.}

\item{...}{- further arguments, such as \code{outDir}, for \link{parseKEGG_file},
\link{parseKEGG_file.list} and database-specific functions (below).}
}
\description{
Execute all parseKEGG parent functions to format specific KEGG databases into data frames.
}
\examples{
KEGG_path <- "~/KEGG" # MODIFY!
parseKEGG_parseKEGG_execute_all(KEGG_path)
# multiple reference_table objects in workspace and .txt files written to
#    output/ (relative to current working directory)

}
\seealso{
\link{parseKEGG_compound}, \link{parseKEGG_enzyme},
         \link{parseKEGG_genome}, \link{parseKEGG_module},

\link{parseKEGG_ko}, \link{parseKEGG_reaction},
         \link{parseKEGG_ko_enzyme}, \link{parseKEGG_ko_reaction}
}

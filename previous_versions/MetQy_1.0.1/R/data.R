#' Example of a vector containing KEGG Ortholog (KO) identifiers
#'
#' A dataset containing an example of a vector containing KEGG Orthologs (K numbers or KOs)
#' that can be used with the function \code{misc_geneVector_module()}
#' or by the function \code{query_genomes_to_modules()} after formatting it into a data frame (see function description).
#'
#' @format A character vector with 2896 entries
#' \describe{
#'   \item{KOs}{K00005, K00009, K00012, K00013, K00014, ...}
#' }
#' @source \url{http://www.kegg.jp/kegg-bin/get_htext?eco00001} and expand the subsections to see the KOs.
#' @examples
#' # Load data
#' data("data_example_KOnumbers_vector")
#' head(data_example_KOnumbers_vector)
#' # [1] "K00005" "K00009" "K00012" "K00013" "K00014" "K00024"
#'
#' @seealso \link{data_example_ECnumbers_vector}

"data_example_KOnumbers_vector"

#' Example of a vector containing Enzyme Commission (EC) numbers.
#'
#' A dataset containing an example of a vector containing Enzyme Commission (EC) numbers
#' that can be used with the function \code{misc_geneVector_module()}
#' or by the function \code{query_genomes_to_modules()} after formatting it into a data frame (see function description).
#'
#' When generating an EC number vector, make sure that all entries have the 4 nomenclature positions
#' (".-" denotes an unspecified field, e.g. "1.1.1.-").
#'
#' @format A character vector with 568 entries
#' \describe{
#'   \item{ECs}{1.1.1.10, 1.1.1.102, 1.1.1.105, 1.1.1.12, 1.1.1.14, 1.1.1.153, ...}
#' }
#' @source \url{http://www.kegg.jp/kegg-bin/get_htext?eco00001} and expand the subsections to see the ECs.
#' @examples
#' # Load data
#' data("data_example_ECnumbers_vector")
#' head(data_example_ECnumbers_vector)
#' # [1] "1.1.1.10"  "1.1.1.102" "1.1.1.105" "1.1.1.12"  "1.1.1.14"  "1.1.1.153"
#'
#' @seealso \link{data_example_KOnumbers_vector}


"data_example_ECnumbers_vector"

#' Example of a data frame with multiple datasets.
#'
#' A data frame containing an example of multiple datasets that could be analysed with the function \code{misc_geneVector_module()}
#' or by the function \code{query_genomes_to_modules()} after formatting it into a data frame (see function description).
#' The dataset entries are to illustrate what the KEGG genome data looks like.
#'
#'@format
#' A data frame with the following columns:\preformatted{
#'  ID       - KEGG genome ID, T number (e.g. "T00001")
#'  ORG_ID   - KEGG organism ID,3-4 letter code (e.g. "eco")
#'  ORGANISM - Organism name (e.g. "Escherichia coli K-12 MG1655")
#'  KOs       - Concatenated string with the K numbers (e.g. "K00013;K00014;K00018;...").
#'  ECs       - Concatenated string with the EC numbers 
#'              (e.g. "1.1.1.1;1.1.1.100;1.1.1.130;...")
#' }
#' @details
#' The columns 'KOs' or 'ECs' need to be pointed at to do the module mapping and is specified by using the argument 'mapBy' in \code{query_genomes_to_modules()} and
#' The default is to use the column named 'KOs', but it can be named differently and specified in the argument field.
#'
#' Specify the character to be used to split the string with the K/EC numbers using the argument 'split_vector_by' in \code{query_genomes_to_modules()}
#' @seealso \link{query_genomes_to_modules}
#' @examples
#' # Load data
#' data("data_example_multi_ECs_KOs")
#' head(data_example_multi_ECs_KOs)

"data_example_multi_ECs_KOs"

#' Data frame used to populate the short name in module_reference_table
#'
#' This data frame contains the following columns and is used to provide a manually abbreviated name to the 
#' modules to ease plotting.
#' @format
#' A data frame with the following columns:\preformatted{
#'  ID         - KEGG module ID, M number (e.g. "M00001")
#'  NAME       - KEGG module name
#'  NAME_SHORT - Manually abbreviated KEGG module name (unique entries)
#' }
#' @seealso \link{parseKEGG_module}
#' @examples
#' # Load data
#' data("module_shortName_mapping")
#' head(module_shortName_mapping)
"data_module_shortName_mapping"

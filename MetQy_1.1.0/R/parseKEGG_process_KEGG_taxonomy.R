#' Retrieve the taxonomy listed as part of the KEGG genome database.
#'
#' @details
#' After processing the KEGG genome database, this function can extract the taxonomic information.
#' This is obtained from splitting the \code{TAXONOMY} column after the "LINEAGE" tag into "words".
#' This information should specify KINGDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS (6), but there are ocasionally more or less words than expected.
#' Therefore NOTE that it is incomplete and inconsistent, most likely because it is derived from multiple sources.
#' Use with caution.
#'
#' @param genome_reference_table - table containing the different genome entries (rows) and data (columns). See \link{parseKEGG_genome}.
#'
#' @param taxonomy_header        - optional. Character with header name or number indicating column index for taxonomic information. Default ("TAXONOMY").
#'
#' @param org_header             - optional. Character with header name or number indicating column index for organism name or ID. Use to name the output rownames. Default (1; first column).
#'
#' @param ...  - further arguments (currently unsupported)
#'
#' @return Data frame with 6 columns containing the taxonomic information (columns: KINGDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS).
#' The rownames contain the organism information (name or ID) as specified with \code{org_header}.
#' \code{UNKNOWN} label is added where no information is available.
#'
#' @examples
#' # Generate KEGG's genome database table
#' genome_reference_table <- parseKEGG_genome("~/KEGG")
#'
#' # Extract taxonomic information from genome_referene_table
#' TAXONOMY_table <- parseKEGG_process_KEGG_taxonomy(genome_reference_table)
#' TAXONOMY_table[1,]
#' #            KINGDOM    PHYLUM             CLASS                  ORDER             FAMILY             GENUS
#' # T00001     Bacteria    Proteobacteria    Gammaproteobacteria    Pasteurellales    Pasteurellaceae    Haemophilus
#'
#' @seealso \link{parseKEGG_genome}
#'
#' @export

############################################################################################################################################

parseKEGG_process_KEGG_taxonomy <- function(genome_reference_table,taxonomy_header = "TAXONOMY",org_header = 1,...){

  ## MANAGE INPUT ----
  if(!is.character(org_header)||!is.numeric(org_header)||length(org_header)>1) org_header <- 1
  ## FIND ORGANISM COLUMNS ----
  organism_col <- ifelse(is.numeric(org_header),org_header,which(names(genome_reference_table)==org_header))
  ORGANISM     <- genome_reference_table[,organism_col]

  #
  TAXONOMY    <- genome_reference_table$TAXONOMY
  TAXONOMY    <- gsub(".*LINEAGE ","",TAXONOMY)

  # INCONSISTENT: KINGDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES
  TAXONOMY_table <- data.frame(matrix("",nrow = nrow(genome_reference_table),ncol = 6),stringsAsFactors = F)
  names(TAXONOMY_table) <- c('KINGDOM', 'PHYLUM', 'CLASS', 'ORDER', 'FAMILY', 'GENUS')
  for(S in 1:length(TAXONOMY)){
    temp <- strsplit(TAXONOMY[S],split = "; ")[[1]]
    TAXONOMY_table[S,1:ifelse(length(temp)>ncol(TAXONOMY_table),ncol(TAXONOMY_table),length(temp))] <- temp[1:ifelse(length(temp)>ncol(TAXONOMY_table),ncol(TAXONOMY_table),length(temp))]
  }
  # head(TAXONOMY_table)
  for(C in 1:ncol(TAXONOMY_table)) TAXONOMY_table[which(TAXONOMY_table[,C]==""),C] <- "UNKNOWN"

  # binary <- matrix(as.numeric(as.logical(TAXONOMY_table!="UNKNOWN")),nrow = nrow(TAXONOMY_table),ncol = ncol(TAXONOMY_table))
  # head(binary)

  # ADD ORGANISM NAME AS ROWNAME ----
  if(is.null(organism_col)){
    warning("No matching 'org_header' column found in 'genome_reference_table'. Output rownames not determined.")
  } else {
    duplicated_rows <- unique(ORGANISM[which(duplicated(ORGANISM))])
    if(length(duplicated_rows)>0){
      warning("Duplicated dataset names provided. Appending \\\\s**[A-Z]\\\\d to distinguish them ")

      for(D in duplicated_rows){
        index <- which(ORGANISM==D)
        ORGANISM[index] <- paste(ORGANISM[index]," **",misc_create_labels(length(index)),sep="")
      }
    }
    rownames(TAXONOMY_table) <- ORGANISM
  }

  return(TAXONOMY_table)

}

#' Find the genome(s) that contain a set of genes.
#'
#' The genes can be either an enzyme, given by its Enzyme Classification (EC) number (e.g. "1.1.1.1"),
#' or a KEGG ortholog identifier (K number, e.g. "K00001").
#'
#' @param genes - character vector of KO or EC numbers.
#'
#' @param use_genome_reference_table - optional. Provide a data frame with updated KEGG genome database.
#'                                     Default (\code{NULL}; inbuilt data used). See Details.
#'
#' @param use_ko_reference_table     - optional. Provide a data frame with updated KEGG ortholog database.
#'                                     Default (\code{NULL}; inbuilt data used). See Details.
#'
#' @param use_enzyme_reference_table - optional. Provide a data frame with updated KEGG enzyme database.
#'                                     Default (\code{NULL}; inbuilt data used). See Details.
#'
#' @return Data frame containing the \code{genes} (rows, specified in the rownames) and the T number of the KEGG genomes (columns) with a binary indicator for presence of gene in given genome.
#'
#' @details
#' The \code{use_} set of arguments allow users with KEGG FTP access to provide the updated data from the 
#' KEGG databases in the form of reference tables AND/OR for advanced users to provide custom-made modules (see below).
#' These reference tables can be generated with the \code{parseKEGG} family of functions and need to have a 
#' specific format (see function descriptions for details on format).
#'
#' If providing \code{use_genome_reference_table}, make sure that the \code{parseKEGG_genome} function is run with arguments \code{addECs = T} and/or \code{addKOs = T}
#' to include genes that make up the genomes.
#'
#' \code{ko_reference_table} and \code{enzyme_reference_table} are used to verify that the \code{genes} provided exist within the KEGG database.
#'
#' @examples
#' genomes <- query_genes_to_genomes("K00844")
#'
#' genomes <- query_genes_to_genomes("1.1.1.1")
#'
#' genomes <- query_genes_to_genomes(genes = paste("K0000",1:3,sep=""))
#'
#' @seealso \link{parseKEGG_genome}, \link{ko_reference_table}, \link{enzyme_reference_table}
#' @export

query_genes_to_genomes <- function(genes,use_genome_reference_table = NULL,use_ko_reference_table = NULL,use_enzyme_reference_table = NULL){

  # MANAGE INPUT ----
  stopifnot(is.character(genes))

  # CHECK FOR GENOME REFERENCE TABLE
  if(!is.null(use_genome_reference_table)&&is.data.frame(use_genome_reference_table)){
    GENOME_REFERENCE_TABLE <- use_genome_reference_table
  } else if(!is.null(use_genome_reference_table)&&!is.data.frame(use_genome_reference_table)){
    warning("'use_genome_reference_table' should be a data frame. Using inbuilt data.")
    GENOME_REFERENCE_TABLE <- genome_reference_table
  }else{
    GENOME_REFERENCE_TABLE <- genome_reference_table
  }

  nKOs <- length(grep("^K\\d{5}$",genes))
  nECs <- length(grep("^\\d{1}[.][[:punct:] [:digit:]]+$",genes))

  # STOP if not KOs or ECs
  if(nKOs==0&&nECs==0) stop("'genes' should contain an EC numbers (e.g. '1.-.-.-') or a KEGG ortholog ID (e.g. 'K00001')")

  ## KO ENTRIES ----
  if(nKOs>nECs){
    GENOME_REF <- GENOME_REFERENCE_TABLE$KOs

    # CHECK FOR KO REFERENCE TABLE ----
    if(!is.null(use_ko_reference_table)&&is.data.frame(use_ko_reference_table)){
      KO_REFERENCE_TABLE <- use_ko_reference_table
    } else if(!is.null(use_ko_reference_table)&&!is.data.frame(use_ko_reference_table)){
      warning("'use_ko_reference_table' should be a data frame. Using inbuilt data.")
      KO_REFERENCE_TABLE <- ko_reference_table
    }else{
      KO_REFERENCE_TABLE <- ko_reference_table
    }

    # CHECK ENTRIES ARE IN DATA BASE ----
    remove <- NULL
     for(i in 1:length(genes)){
       if(length(grep(genes[i],KO_REFERENCE_TABLE$ID))==0) remove <- c(remove,i)
     }

     if(length(remove)>0){
       warning(paste(paste(genes[remove],collapse = ","),"- not valid KEGG ortholog entry"))
       # genes <- genes[-remove]
     }
    rm(remove)

    ## EC ENTRIES ----
  }else{
    GENOME_REF <- GENOME_REFERENCE_TABLE$ECs

    # CHECK FOR ENZYME REFERENCE TABLE ----
    if(!is.null(use_enzyme_reference_table)&&is.data.frame(use_enzyme_reference_table)){
      ENZYME_REFERENCE_TABLE <- use_enzyme_reference_table
    } else if(!is.null(use_enzyme_reference_table)&&!is.data.frame(use_enzyme_reference_table)){
      warning("'use_enzyme_reference_table' should be a data frame. Using inbuilt data.")
      ENZYME_REFERENCE_TABLE <- enzyme_reference_table
    }else{
      ENZYME_REFERENCE_TABLE <- enzyme_reference_table
    }

    # CHECK ENTRIES ARE IN DATA BASE ----
    remove <- NULL
    for(i in 1:length(genes)){
      if(is.na(match(genes[i],ENZYME_REFERENCE_TABLE$ID))==0) remove <- c(remove,i)
    }

    if(length(remove)>0){
      warning(paste(paste(genes[remove],collapse = ","),"- not valid KEGG enzyme entry"))
      # genes <- genes[-remove]
    }
    rm(remove)

  }

  # FIND GENOMES ----
  TABLE <- NULL
  for(i in 1:length(genes)){
    genomes_index     <- grep(pattern = genes[i],GENOME_REF,fixed = T)
    if(length(genomes_index)>0){
      TABLE           <- rbind(TABLE,cbind(rep(genes[i],length(genomes_index)),GENOME_REFERENCE_TABLE$ID[genomes_index]))
    } else {
      if(is.null(dim(TABLE))){
        TABLE           <- c(genes[i],"no_match") # ADD AS PADDING TO REGISTER NO GENES
      } else{
        TABLE           <- rbind(TABLE,cbind(genes[i],"no_match"))
      }
    }
  }

  TABLE <- as.data.frame(TABLE,stringsAsFactors = F)
  TABLE$V3 <- 1

  # CAST DATA.FRAME TO MATRIX ----
  genomes           <- reshape2::dcast(TABLE,formula = V1~V2,value.var = "V3",drop = F,fill = 0)
  rownames(genomes) <- genomes[,1]
  genomes           <- genomes[,-1]

  # REMOVE PADDING - no genes found
  ROWS <- rownames(genomes)
  if(length(which(colnames(genomes)=="no_match"))>0) genomes <- genomes[,-which(colnames(genomes)=="no_match")]

  return(genomes)

}

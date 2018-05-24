#' Find the genome(s) that contain a set of genes.
#'
#' The genes can be either an enzyme, given by its Enzyme Classification (EC) number (e.g. "1.1.1.1"),
#' or a KEGG ortholog identifier (K number, e.g. "K00001").
#'
#' @param genes - character vector of KO or EC numbers. See Details for more info on the use of ECs.
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
#' @param ...  - further arguments (currently unsupported)
#'
#' @return Data frame containing the \code{genes} (rows, specified in the rownames) and
#' the T number of the KEGG genomes (columns) with a binary indicator for presence of gene in given genome.
#'
#' @details
#' When providing EC numbers, the user can provide a full EC number as described above or the first three values/positions (e.g. "1.1.1" and "1.1.1.-" are both allowed).
#' In the latter case, the function evaluates all enzymes (EC numbers) that match the three values/positions provided
#' (e.g., using “1.1.1” would cause the function to evaluate all enzymes starting with that EC number combination).
#' In other words, using a three number entry would act as a wildcard functioning for the 4th position of the EC notation.
#'
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

query_genes_to_genomes <- function(genes,use_genome_reference_table = NULL,use_ko_reference_table = NULL,use_enzyme_reference_table = NULL,...){

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
    REMOVE <- NULL
     for(i in 1:length(genes)){
       if(length(grep(genes[i],KO_REFERENCE_TABLE$ID))==0) REMOVE <- c(REMOVE,i)
     }

     if(length(REMOVE)>0){
       warning(paste(paste(genes[REMOVE],collapse = ","),"- not valid KEGG ortholog entry"))
       # genes <- genes[-REMOVE]
     }
    rm(REMOVE)

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

    # CHECK IF IT's A COMPLETE ENZYME ID ----
    add_genes    <- NULL
    remove_index <- NULL
    for(i in 1:length(genes)){
      enzyme_parts <- strsplit(genes[i],"[.]")[[1]]
      if(length(enzyme_parts)==4){

        if(enzyme_parts[4]=="-"){
          enzyme_parts <- enzyme_parts[1:3]
        }else{
          next
        }
      }

      # CHECK FOR SPACE
      if(enzyme_parts[length(enzyme_parts)]==" ") enzyme_parts <- enzyme_parts[1:(length(enzyme_parts)-1)]
      # ADD UNDEFINED ENTRIES
      if(length(enzyme_parts)==2||enzyme_parts[3]=="-"){
        warning("Only two levels defined. Please specify at least three level (e.g. '1.1.1).'")
        remove_index <- c(remove_index,i)
      }
      if(length(enzyme_parts)==3){
        enzyme <- paste(enzyme_parts,collapse = "[.]")
        specific_enzymes <- ENZYME_REFERENCE_TABLE$ID[grep(paste("^",enzyme,sep=""),ENZYME_REFERENCE_TABLE$ID)]
        if(length(specific_enzymes)>0){
          add_genes <- c(add_genes,specific_enzymes)
          warning(paste("Enzyme",genes[i],"is either undefinied or incomplete. Adding",length(specific_enzymes),"enzyme IDs to search."))
        }
        remove_index <- c(remove_index,i)

        # NO UNDEFINED ENZYMES IN KEGG DATABASE (ONLY IN MODULE DEFINITIONS)

      }
    } # for genes

    # CHECK THERE ARE VALIDLY DEFINED ENZYMES
    if(!is.null(add_genes)){
      add_genes <- unique(add_genes)
      genes     <- sort(unique(c(genes,add_genes)))
    }
    if(!is.null(remove_index)) genes <- genes[-remove_index]
    if(length(genes)==0) stop("Please provide enzyme with at least three levels defined (e.g. '1.1.1')")

    # CHECK ENTRIES ARE IN DATA BASE ----
    REMOVE <- NULL
    for(i in 1:length(genes)){
      if(is.na(match(genes[i],ENZYME_REFERENCE_TABLE$ID))) REMOVE <- c(REMOVE,i)
    }

    if(length(REMOVE)>0){
      genes <- genes[-REMOVE]
      if(length(genes)==0) stop("No valid KEGG enzyme provided")
      warning(paste(paste(genes[REMOVE],collapse = ","),"- not valid KEGG enzyme entry"))
    }
    rm(REMOVE)

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

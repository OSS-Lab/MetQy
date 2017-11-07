#' Identify missing genes from a KEGG module
#'
#' Identify missing genes from a KEGG module given a set of genes or a genome ID to obtain a complete module.
#'
#' @param GENOME          - character vector containing a single genome identifier or set of genes or enzymes that define a [meta]genome. See Details.
#'
#' @param MODULE_ID       - KEGG module ID to be analysed (e.g. \code{"M00001"}).
#'
#' @param PRINT_TO_SCREEN - logical. Should a print-friendly result be displayed? Default(\code{TRUE}).
#'
#' @param use_module_reference_table - optional. Provide a data frame with updated KEGG module database OR with costum-made modules.
#'                                     Default (\code{NULL}; inbulit data used). See Details.
#'
#' @param use_genome_reference_table - optional. Provide a data frame with updated KEGG genome database.
#'                                     Default (\code{NULL}; inbulit data used). See Details.
#'
#' @details
#' \code{GENOME} can be a genome identifier (T0 number or a 3 or 4 letter code; e.g. "T00001" or "eco", respectively) OR
#' a character vector contianing a set of genes (i.e. either EC or K numbers; e.g. "1.1.1.1" or "K00001", respectively).
#' Examples of the latter can be found in the following objects:
#' \link{data_example_KOnumbers_vector} or \link{data_example_ECnumbers_vector}.
#'
#' Organism name NOT supported. Use the KEGG database website (\url{http://www.genome.jp/kegg/catalog/org_list.html}) to determine the genome identifier.
#'
#' Note that the pipe ('|') in the DEFINITION indicates an OR operation
#' This means that there are several possibile genes that carry out the same reaction or function and only one is required.
#'
#' The \code{use_} set of arguments allow users with KEGG FTP access to provide the updated data from the KEGG databases in the form of reference tables AND/OR for advanced users to provide contum-made modules (see below).
#' These reference tables can be generated with the \code{parseKEGG} family of functions and need to have a specific format (see function descriptions for details on format).
#'
#' The module definition (contained in \code{module_reference_table}) describes the relationship between genes and modules
#' and is used to identify the modules in which \code{gene} is involved.
#' The user can provide costum-made module definitions that use the logical expression format (however, the table format must be conserved!).
#'
#' @return Data frame containing the following columns:
#' \preformatted{
#'	 BLOCK_DEF     - the KEGG module DEFINITION of each block, with missing genes flagged by '*';
#'	 PRESENT       - binary indicator of the automatic evaluation;
#'	 MISSING_GENES - list of missing genes.
#' }
#' @examples
#' # Load data
#' data("data_example_KOnumbers_vector")
#' OUT <- query_missingGenes_from_module(data_example_KOnumbers_vector,"M00001")
#'
#' @export

query_missingGenes_from_module  <- function(GENOME,MODULE_ID, PRINT_TO_SCREEN = TRUE,use_genome_reference_table = NULL,use_module_reference_table = NULL){

  # MANAGE INPUT ----
  stopifnot(is.character(GENOME))
  stopifnot(length(MODULE_ID)==1, MODULE_ID!="")

  # CHECK FOR GENOME REFERENCE TABLE
  if(!is.null(use_genome_reference_table)&&is.data.frame(use_genome_reference_table)){
    GENOME_REFERENCE_TABLE <- use_genome_reference_table
  } else if(!is.null(use_genome_reference_table)&&!is.data.frame(use_genome_reference_table)){
    warning("'use_genome_reference_table' should be a data frame. Using inbuilt data.")
    GENOME_REFERENCE_TABLE <- genome_reference_table
  }else{
    GENOME_REFERENCE_TABLE <- genome_reference_table
  }

  # CHECK FOR MODULE REFERENCE TABLE
  if(!is.null(use_module_reference_table)&&is.data.frame(use_module_reference_table)){
    MODULE_REFERENCE_TABLE <- use_module_reference_table
  } else if(!is.null(use_module_reference_table)&&!is.data.frame(use_module_reference_table)){
    warning("'use_module_reference_table' should be a data frame. Using inbuilt data.")
    MODULE_REFERENCE_TABLE <- module_reference_table
  }else{
    MODULE_REFERENCE_TABLE <- module_reference_table
  }
  index <- grep(MODULE_ID,MODULE_REFERENCE_TABLE$ID)
  if(length(index)==0) stop(paste("No Module ID matches:", MODULE_ID))

  ### DEFINE GENOME GENES ----
  T_ID      <- gsubfn::strapplyc(pattern = "^T\\d{5}$",X = GENOME,simplify = T)
  letter_ID <- gsubfn::strapplyc(pattern = "^[a-z]{3,4}$",X = GENOME,simplify = T)
  KOs       <- gsubfn::strapplyc(pattern = "^K\\d{5}$",X = GENOME,simplify = T)
  ECs       <- gsubfn::strapplyc(pattern = "^\\d{1}[.][[:punct:] [:digit:]]+$",X = GENOME,simplify = T)

  # CHECK THAT ONLY ONE TYPE OF ENTRY WAS GIVEN AND GET GENES
  if(!is.list(T_ID)){
    if(length(T_ID)>1) warning("Multiple T identidiers provided. Only first will be used")
    ref_index   <- grep(paste("^",T_ID[1],"$",sep=""),GENOME_REFERENCE_TABLE$ID)
    if(length(ref_index)==0){
      stop(paste(T_ID,"(T number identifier) does not match any in KEGG genome database"))
    }else{
      gene_vector <- strsplit(GENOME_REFERENCE_TABLE$KOs[ref_index],split = "[;]")[[1]]
      KO          <- T
    }
  }
  if(!is.list(letter_ID)){
    if(length(letter_ID)>1) warning("Multiple 3 or 4 letter code identidiers provided. Only first will be used")
    ref_index   <- grep(paste("^",letter_ID[1],"$",sep=""),GENOME_REFERENCE_TABLE$ORG_ID)
    if(length(ref_index)==0){
      stop(paste(letter_ID,"(3 or 4 letter code identifier) does not match any in KEGG organism database"))
    }else{
      gene_vector <- strsplit(GENOME_REFERENCE_TABLE$KOs[ref_index],split = "[;]")[[1]]
      KO          <- T
    }
  }
  if(!is.list(KOs)){
    gene_vector <- GENOME
    KO          <- T
  }
  if(!is.list(ECs)){
    gene_vector <- GENOME
    KO          <- F
  }
  if(!exists('gene_vector')) stop("No GENOME_INFO matched.")

  # DETERMINE KOs or ECs ----
   if(KO){
    DEFINITION <- MODULE_REFERENCE_TABLE$DEFINITION_KOs
    PATTERN    <- "K\\d{5}"
  } else {
    DEFINITION <- MODULE_REFERENCE_TABLE$DEFINITION_ECs
    PATTERN    <- "\\d{1}[.][[:punct:] [:digit:]]+"
  }

  ### SEARCH BLOCKS ----
  block_defs  <- strsplit(DEFINITION[index], split = " ")[[1]]
  nBlocks     <- length(block_defs)

  TABLE <- data.frame("block_No"      = 1:nBlocks,
                      "PRESENT"       = numeric(nBlocks),
                      "BLOCK_DEF"     = character(nBlocks),
                      "MISSING_GENES" = character(nBlocks),
                      stringsAsFactors = F)

  for(B in 1:nBlocks){
    thisBlock         <- block_defs[B]

    TABLE$PRESENT[B]  <- misc_evaluate_block(gene_vector,thisBlock)

    # FIND MISSING GENES
    genes_needed  <- gsubfn::strapplyc(thisBlock,PATTERN,simplify = T)
    gene_index    <- na.omit(match(genes_needed,gene_vector))

    genes_present <- gene_vector[gene_index]

    # OR  <- gsubfn::strapplyc(thisBlock,"\\|",simplify = T)
    # FIND OUT STRUCTURE
    # AND <- gsubfn::strapplyc(thisBlock,"\\+",simplify = T)

    REPLACE <- setdiff(genes_needed,genes_present)
    for(R in REPLACE) thisBlock <- gsub(R,paste("*",R,"*",sep=""),thisBlock)

    # IF THERE ARE NO "AND" OPERATORS
    # if(!is.list(AND))
    if(length(REPLACE)>0) TABLE$MISSING_GENES[B] <- paste(REPLACE,collapse = ";")

    block_defs[B] <- thisBlock
  }

  TABLE$BLOCK_DEF <- block_defs
  present_index <- which(TABLE$PRESENT==1)
  if(length(present_index)>0) TABLE$MISSING_GENES[present_index] <- ""

  # PRINT TO SCREEN
  if(PRINT_TO_SCREEN){
    cat(MODULE_ID, fill = T)
    for(D in TABLE$BLOCK_DEF) cat("\t",D, fill = T)
    cat('\n\t\tFor more info visit:',  paste('http://www.genome.jp/kegg-bin/show_module?',MODULE_ID,sep=""),fill=T)
  }

 return(TABLE)

}

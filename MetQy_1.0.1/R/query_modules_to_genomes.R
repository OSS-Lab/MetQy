#' Find a KEGG genome that has a given KEGG module complete
#'
#' Find a KEGG genome that has a given KEGG module complete (to a \code{threshold}).
#' The module ‘completeness’ is based on the fraction of complete module blocks (given the KEGG module logical definition).
#' Thus, a threshold of 0.5 would mean that the function would return all genomes that contain at least half of the blocks of the given module.
#' See \link{parseKEGG_module} for further details on the KEGG module database.
#'
#' @param MODULE_ID - KEGG module identifier (M number, e.g. "M00001").
#'
#' @param threshold - optional. completness fraction desired (greater or equal). Default (\code{1}). Used as \code{fraction>=threshold}.
#'
#' @param use_matrix_dataframe - optional. Provide a matrix or data frame with the output from \code{query_genomes_to_modules} using all KEGG genome entries with updated data.
#'                                     Default (\code{NULL}; inbulit data used). See Details.
#'
#' @param use_module_reference_table - optional. Provide a data frame with updated KEGG module database OR with costum-made modules.
#'                                     Default (\code{NULL}; inbulit data used). See Details.
#'
#' @return Vector containing the KEGG genome identifiers (T number) that express the module to the set \code{threshold}.
#'
#' @details
#' The \code{use_} set of argument allows users with KEGG FTP access to provide the updated data from the KEGG databases in the form of reference tables AND/OR for advanced users to provide contum-made modules (see below).
#' These reference tables can be generated with the \code{parseKEGG} family of functions and need to have a specific format (see function descriptions for details on format).
#'
#' The user can then run the function \code{query_genomes_to_modules} with the updated KEGG databases information AND/OR with
#' costum-made module definitions and use the matrix that is outputed to run this function.
#'
#' The module definition (contained in \code{module_reference_table}) describes the relationship between genes and modules
#' and is used to identify the modules in which \code{gene} is involved.
#' The user can provide costum-made module definitions that use the logical expression format (however, the table format must be conserved!).
#'
#' FTP users - this function does not currently support using updated data.
#' This funciton can be replicated manually by first running \code{query_genomes_to_modules} on all genomes and modules and then applying the desired threshold (see function code).
#'
#' @format \code{use_matrix_dataframe}  rows - genome identifiers, columns - module IDs.
#'
#' @seealso \link{parseKEGG_module}
#'
#' @examples
#' genomes <- query_modules_to_genomes("M00001")
#'
#' genomes <- query_modules_to_genomes(MODULE_ID = c("M00001","M00002"),,threshold = 0.75)
#'
#' @export
#'

query_modules_to_genomes <- function(MODULE_ID,threshold = 1, use_matrix_dataframe = NULL,use_module_reference_table = NULL){

  #### MANAGE INPUT ----
  stopifnot(is.character(MODULE_ID),length(MODULE_ID)>0,threshold>=0,threshold<=1)

  ## REFERENCE TABLE
  # all_genomes_all_modules_matrix - contains all KEGG genomes (rows) and modules (cols)

  # CHECK FOR REFERENCE TABLE ----
  if(is.matrix(use_matrix_dataframe)||is.data.frame(use_matrix_dataframe)){
    GENOME_MODULE_MATRIX <- use_matrix_dataframe
  } else if(!is.null(use_matrix_dataframe)&&!(is.matrix(use_matrix_dataframe)&&is.data.frame(use_matrix_dataframe))){
    warning("'use_matrix_dataframe' should be a matrix or data frame. Using inbuilt data.")
  }else{
    GENOME_MODULE_MATRIX <- all_genomes_all_modules_matrix
  }

  # CHECK FOR MODULE REFERENCE TABLE ----
  if(!is.null(use_module_reference_table)&&is.data.frame(use_module_reference_table)){
    MODULE_REFERENCE_TABLE <- use_module_reference_table
  } else if(!is.null(use_module_reference_table)&&!is.data.frame(use_module_reference_table)){
    warning("'use_module_reference_table' should be a data frame. Using inbuilt data.")
    MODULE_REFERENCE_TABLE <- module_reference_table
  }else{
    MODULE_REFERENCE_TABLE <- module_reference_table
  }

  # CHECK ENTRIES ARE IN DATA BASE ----
  remove <- NULL
  for(i in 1:length(MODULE_ID)){
    if(length(grep(MODULE_ID[i],MODULE_REFERENCE_TABLE$ID))==0) remove <- c(remove,i)
  }

  if(length(remove)>0){
    warning(paste(paste(MODULE_ID[remove],collapse = ","),"- not matched with KEGG module. Format: M#####, e.g. 'M00001'"))
    MODULE_ID <- MODULE_ID[-remove]
  }
  rm(remove)

  if(length(MODULE_ID)==0) stop("No valid entries provided. Check MODULE_ID input")

  # FIND MODULES ----
  MATRIX <- matrix(as.numeric(GENOME_MODULE_MATRIX>=threshold),dim(GENOME_MODULE_MATRIX))
  colnames(MATRIX) <- colnames(GENOME_MODULE_MATRIX)
  rownames(MATRIX) <- rownames(GENOME_MODULE_MATRIX)

  # match_index          <- which(names(GENOME_MODULE_MATRIX)==MODULE_ID[i])
  match_index  <- match(MODULE_ID,colnames(GENOME_MODULE_MATRIX))

  genome_match <- MATRIX[,match_index]
  # REMOVE EMPTY MODULE COLUMNS AND ZERO ROWS (GENOMES)
  if(is.matrix(genome_match)){
    if(length(which(apply(genome_match,1,sum)==0))>0) genome_match <- genome_match[-which(apply(genome_match,1,sum)==0),]
    if(length(which(apply(genome_match,2,sum)==0))>0) genome_match <- genome_match[,-which(apply(genome_match,2,sum)==0)]
  }else{
    if(length(which(genome_match==0))>0) genome_match <- genome_match[-which(genome_match==0)]
  }

  # RETURN
  return(genome_match)

}

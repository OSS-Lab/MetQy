#' Given a set of genes, find the modules they are involved in.
#'
#' The genes can be either enzymes, given by its Enzyme Classification (EC) number (e.g. "1.1.1.1"),
#' or KEGG ortholog identifiers (K number, e.g. "K00001").
#' Using EC or K numbers might give different results, as an EC number might map to multiple K numbers or none.
#'
#' @param genes - character vector (length 1). K or EC number.
#'
#' @param use_module_reference_table - optional. Provide a data frame with updated KEGG module database OR with costum-made modules.
#'                                     Default (\code{NULL}; inbulit data used). See Details.
#'
#' @param use_ko_reference_table     - optional. Provide a data frame with updated KEGG ortholog database.
#'                                     Default (\code{NULL}; inbulit data used). See Details.
#'
#' @param use_enzyme_reference_table - optional. Provide a data frame with updated KEGG enzyme database.
#'                                     Default (\code{NULL}; inbulit data used). See Details.
#'
#' @return Data frame containing a binary indicator for \code{genes} (rows, specified by rownames) involved in modules (columns, Module IDs specified in names).
#' A column called '\code{no_match}' is return if one or more of the genes is not involved in any modules.
#'
#' \code{NULL} is returned when there are no valid entries to evaluate. Note that enzyme entries are check for 4 position completeness.
#'
#' @details
#' The \code{use_} set of arguments allow users with KEGG FTP access to provide the updated data from the KEGG databases in the form of reference tables AND/OR for advanced users to provide contum-made modules (see below).
#' These reference tables can be generated with the \code{parseKEGG} family of functions and need to have a specific format (see function descriptions for details on format).
#'
#' The module definition (contained in \code{module_reference_table}) describes the relationship between genes and modules
#' and is used to identify the modules in which \code{genes} is involved.
#' The user can provide costum-made module definitions that use the logical expression format (however, the table format must be conserved!).
#'
#' If enzyme identifiers are provided note that as there are lingering unspecified enzymes (e.g. '1.1.1.-') in the KEGG module EC-based definition, this function also returns a row entry for
#' the unspecified-versions of enzymes provided involved in one or more modules.
#'
#' @examples
#' modules <- query_genes_to_modules("K00844")
#'
#' modules <- query_genes_to_modules("1.1.1.1")
#'
#' @seealso \link{parseKEGG_module}
#' @export

query_genes_to_modules <- function(genes,use_module_reference_table = NULL,use_ko_reference_table = NULL,use_enzyme_reference_table = NULL){

  # MANAGE INPUT ----
  stopifnot(is.character(genes))
  # if(length(genes)>1){
  #   genes <- genes[1]
  #   warning("Multiple genes provided. Only first used.")
  # }

  # CHECK FOR MODULE REFERENCE TABLE
  if(!is.null(use_module_reference_table)&&is.data.frame(use_module_reference_table)){
    MODULE_REFERENCE_TABLE <- use_module_reference_table
  } else if(!is.null(use_module_reference_table)&&!is.data.frame(use_module_reference_table)){
      warning("'use_module_reference_table' should be a data frame. Using inbuilt data.")
      MODULE_REFERENCE_TABLE <- module_reference_table
  }else{
    MODULE_REFERENCE_TABLE <- module_reference_table
  }

  nKOs <- length(grep("^K\\d{5}$",genes))
  nECs <- length(grep("^\\d{1}[.][[:punct:] [:digit:]]+$",genes))

  if(nKOs==0&&nECs==0){
    stop("'genes' should contain valid EC numbers (e.g. '1.-.-.-') or KEGG ortholog IDs (e.g. 'K00001')")
  }

  if(nKOs>nECs){
    DEFINITIONS <- MODULE_REFERENCE_TABLE$DEFINITION_KOs
    EC <- FALSE

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

  }else{
    DEFINITIONS <- MODULE_REFERENCE_TABLE$DEFINITION_ECs
    EC <- TRUE

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
    REMOVE <- NULL
    for(i in 1:length(genes)){
      if(is.na(match(genes[i],ENZYME_REFERENCE_TABLE$ID))) REMOVE <- c(REMOVE,i)
    }

    if(length(REMOVE)>0){
      warning(paste(paste(genes[REMOVE],collapse = ","),"- not valid KEGG enzyme entry"))

      check_enzymes <- genes[REMOVE]
      uncomplete_enzyme <- NULL
      for(U in length(check_enzymes):1){
        if(length(strsplit(check_enzymes[U], split = "[.]")[[1]])==3){
          uncomplete_enzyme <- c(uncomplete_enzyme,REMOVE[U])
          genes[REMOVE[U]] <- paste(genes[REMOVE[U]],".-",sep="")
        }
      }
      # REMOVE_REMOVE <- setdiff(REMOVE,uncomplete_enzyme)
      # if(length(REMOVE_REMOVE)) genes <- genes[-REMOVE_REMOVE]

      if(length(genes)==length(REMOVE)&&is.null(uncomplete_enzyme)) return(NULL)
      rm(REMOVE,uncomplete_enzyme,check_enzymes,U)
    }

  }

  # ADD UNSPECIFIC ENZYMES TO REPERTOIR
  if(EC) {
    enzyme            <- t(as.data.frame(strsplit(genes, split = "[.]")))
    enzyme[,4]        <- "-"
    unspecific        <- unique(apply(enzyme,1,paste,collapse = "."))
    genes             <- c(genes,unspecific)
    genes             <- unique(genes)
    # unspecific_index  <- (length(genes)-length(unspecific)+1):length(genes)
  }

  # FIND MODULES ----
  TABLE <- NULL
  for(i in 1:length(genes)){
    match_index     <- grep(pattern = genes[i],DEFINITIONS,fixed = T)
    if(length(match_index)>0){
      TABLE           <- rbind(TABLE,cbind(rep(genes[i],length(match_index)),MODULE_REFERENCE_TABLE$ID[match_index]))
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
  modules           <- reshape2::dcast(TABLE,formula = V1~V2,value.var = "V3",drop = F,fill = 0)
  rownames(modules) <- modules[,1]
  modules           <- modules[,-1]

  # REMOVE UNSPECIFIC GENES WITHOUT MATCHES ----
  if(EC){
    no_match_col <- which(names(modules)=="no_match")
    if(length(no_match_col)>0){
      unspecific_index <- grep("-",rownames(modules),fixed = T)
      if(length(unspecific_index)>0){
        REMOVE_index <- which(modules[unspecific_index,no_match_col]==1)
        if(length(REMOVE_index)>0) modules <- modules[-unspecific_index[REMOVE_index],]
      }
    }
  }
  no_match_col <- which(names(modules)=="no_match")
  if(length(no_match_col)>0&&sum(modules[,no_match_col])==0){
    COLS <- names(modules); ROWS <- rownames(modules)
    if(length(COLS)<2){
      modules <- NULL
    }else if(length(COLS)<3){
      modules         <- data.frame(modules[,-no_match_col],row.names = ROWS)
      names(modules)  <- COLS[-no_match_col]
    }else{
      modules <-  modules[,-no_match_col]
    }
  }

  return(modules)
}

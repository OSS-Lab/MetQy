#' Evaluates the KEGG modules presence given genome information.
#'
#' This function returns the ‘completeness’ of KEGG modules in the provided genome(s)
#' (either as a genome identifier or as a lists of Enzyme Comission (EC) number or KEGG ortholog identifier, K number).
#' The user can define which modules the function should check by providing a single or set of modules under ‘\code{MODULE_ID}’.
#' If this is left blank, the function returns completeness of all KEGG modules (excluding modules defined in terms of other modules).
#' The function can also be restricted to a subset of all KEGG modules by using the \code{SEARCH_NAME} and \code{SEARCH_CLASS} arguments. See Details.
#'
#' @param GENOME_INFO     - character vector containing genome identifier(s) or organism name(s) OR data frame containing genome identifier/names(s) and gene list. See Details.
#'
#' @param GENOME_ID_COL   - optional. Column NAME or NUMBER containing genome NAME or IDENTIFIER. Default (\code{1}; <first column>). See Details.
#'
#' @param GENES_COL       - optional. Column NAME or NUMBER pointing to the GENES. Default (\code{2}; <second column>). See Details.
#'
#' @param splitBy         - string indicating the split pattern for the data contained in the column indicated by GENES_COL.
#'                          Default (" [;] ": uses ' \code{;} ' , the ' [ ] ' indicates that it is NOT a regular expression).
#'
#' @param MODULE_ID        - optional. Character vector listing specific KEGG module IDs (e.g. "M00001"). Default ("").
#'
#' @param SEARCH_NAME     - optional. Character vector listing terms to search in KEGG module \code{NAME} field (case-insensitive). Default ("").
#'
#' @param SEARCH_CLASS_I  - optional. Character vector listing terms to search in KEGG module \code{CLASS_I} field (case-insensitive). Default ("").
#'
#' @param SEARCH_CLASS_II - optional. Character vector listing terms to search in KEGG module \code{CLASS_II} field (case-insensitive). Default ("").
#'
#' @param SEARCH_CLASS_III - optional. Character vector listing terms to search in KEGG module \code{CLASS_III} field (case-insensitive). Default ("").
#'
#' @param EXCLUDE_NAME     - optional. Character vector listing terms that if matched in KEGG module \code{NAME} field will be excluded (case-insensitive).
#'                           Default ("").
#'
#' @param OUT_MODULE_NAME  - optional (logical). Should the column names of \code{MATRIX} be the module IDs (M numbers) or module names?
#'                           Default (\code{FALSE}; return matrix with M numbers)
#'
#' @param META_OUT         - optional (logical). Should the KEGG module metadata be outputed? Default (\code{FALSE}).
#'
#' @param ADD_OUT          - optional (logical). Should additional information be outputed? Default (\code{FALSE}).
#'
#' @param use_genome_reference_table - optional. Provide a data frame with updated KEGG module database.
#'                                     Default (\code{NULL}; inbulit data used). See Details.
#'
#' @param ...              - further arguments passed to \link{misc_geneVector_module}, such as \code{KO_in_DEF_EC} or \code{use_module_reference_table}. See Details.
#'
#' @details
#' This function processes the \code{GENOME_INFO} by pasing each in turn to \code{misc_geneVector_module()} and collating all the output.
#' \code{GENOME_ID_COL} and \code{GENES_COL} are only used if \code{GENOME_INFO} is a data frame.
#' Post-analysis of this output can be carried out with \code{analysis_genomes_module_output()}.
#'
#' Appropriate \code{GENOME_INFO} input can be either a character vector or a data frame:
#' \preformatted{
#' *character vector*
#'    options:
#'	     (1) KEGG taxonomy identifier (T0 number; e.g. "T00001") AND/OR KEGG organism identifier (3 or 4 letter code; e.g. "eco").
#'	     (2) Organism's scientific name (e.g. "Escherichia coli"; case-insensitive). Multiple strains might be matched in the KEGG organisms database and all will be processed.
#'	         Strain information or full organism name can be added to reduce the search results; use the 'Definition' entry in KEGG \url{http://www.genome.jp/kegg/catalog/org_list.html}.
#'
#'	         NOTE: do NOT combine NAMES with IDENTIFIERS
#'	         WARNING issued when there is no matching identifier in the KEGG genome/organism databases.
#' *data frame*
#'    > column I  - genome name/identifier (pointed to by \code{GENOME_ID_COL}).
#'    > column II - concatenated string of either EC or K numbers using \code{splitBy} as delimiter (pointed to by \code{GENES_COL}).
#' }
#'
#' \code{KO_in_DEF_EC} - If the genes given are EC numbers, should K numbers present in the KEGG module definition be assumed to be present or not? Default (\code{FALSE}).
#'
#' The \code{use_} set of argument allows users with KEGG FTP access to provide the updated data from the KEGG databases in the form of reference tables.
#' These reference tables can be generated with the \code{parseKEGG} family of functions and need to have a specific format (see function descriptions for details on format).
#'
#' @return Returns a list containing the following objects:
#'
#' \code{$MATRIX} - matrix of the datasets (rows) and the modules searched (columns) containing the fraction completeness.
#'
#' \code{$QUERIES} - data frame listing the SEARCH_TERMS and ARGUMENTS used.
#'
#' \code{$METADATA} - data frame containing the metadata from the modules analysed (if \code{META_OUT == TRUE}).\preformatted{
#'  Columns:
#'   (1) MODULE_ID     (4) CLASS_I     (7) DEFINITION
#'   (2) MODULE_NAME   (5) CLASS_II    (8) OPTIONAL
#'   (3) NAME_SHORT    (6) CLASS_III
#'  }
#'  An OPTIONAL entry of 'NA' indicates that there are no optional K numbers in that module.
#'
#' \code{$ADD_INFO} - data frame containing additional information of the analysis (if \code{ADD_OUT == TRUE}).\preformatted{
#'  Columns:
#'   (1) GENOME_ID     (3) NAME_SHORT   (5) nBLOCKS    (7) OPTIONAL_PRESENT
#'   (2) MODULE_ID     (4) FRACTION     (6) COVERAGE
#'  }
#'   where 'COVERAGE' refers to the genes provided that are involved in the given module and genome.
#'
#'   See \link{misc_geneVector_module} for additional details on the output.
#'
#' @seealso \link{misc_geneVector_module}, \link{analysis_genomes_module_output}, \link{plot_heatmap}, \link{data_example_multi_ECs_KOs}
#'
#' @examples
#' ## USE T numbers
#' T_NUMEBERS <- paste("T0000",1:5,sep="")
#' OUT <- query_genomes_to_modules(T_NUMEBERS,MODULE_ID = paste("M0000",1:5,sep=""),META_OUT = T, ADD_OUT = T)
#'
#' ## USE SPECIES NAMES
#' names <- c("escherichia coli","heliobacter")
#' OUT <- query_genomes_to_modules(names,MODULE_ID = paste("M0000",1:5,sep=""),META_OUT = T, ADD_OUT = T)
#'
#' ## USE USER-SPECIFIED GENE SETS
#' data(data_example_multi_ECs_KOs) # load example data set
#' names(data_example_multi_ECs_KOs)
#'    # "ID"       "ORG_ID"   "ORGANISM" "KOs"      "ECs"
#' OUT <- query_genomes_to_modules(data_example_multi_ECs_KOs,GENOME_ID_COL = "ID",GENES_COL = "KOs",MODULE_ID = paste("M0000",1:5,sep=""),META_OUT = T,ADD_OUT = T)
#'
#' # Using EC numbers (less accurate)
#' OUT <- query_genomes_to_modules(data_example_multi_ECs_KOs,GENOME_ID_COL = "ID",GENES_COL = "ECs",MODULE_ID = paste("M0000",1:5,sep=""),META_OUT = T,ADD_OUT = T)
#'
#' @export

############################################################################################################################################

query_genomes_to_modules <- function(GENOME_INFO,splitBy = "[;]",GENOME_ID_COL = 1,GENES_COL = 2,
                                     MODULE_ID="",SEARCH_NAME="",SEARCH_CLASS_I="",SEARCH_CLASS_II="",SEARCH_CLASS_III="",
                                     EXCLUDE_NAME ="",OUT_MODULE_NAME = FALSE,META_OUT = FALSE,ADD_OUT = FALSE,use_genome_reference_table = NULL,...){

  ##### MANAGE INPUT ----
  stopifnot(is.character(GENOME_INFO)||is.data.frame(GENOME_INFO))
  stopifnot(is.character(splitBy),length(splitBy)==1)

  # CHECK FOR MODULE REFERENCE TABLE
  if(!is.null(use_genome_reference_table)&&is.data.frame(use_genome_reference_table)){
    GENOME_REFERENCE_TABLE <- use_genome_reference_table
  } else if(!is.null(use_genome_reference_table)&&!is.data.frame(use_genome_reference_table)){
    warning("'use_genome_reference_table' should be a data frame. Using inbuilt data.")
    GENOME_REFERENCE_TABLE <- genome_reference_table
  }else{
    GENOME_REFERENCE_TABLE <- genome_reference_table
  }

  ## GENOME INFO GIVEN ----
  if(is.character(GENOME_INFO)){
    # USE GENOME_REFERENCE_TABLE
    if(mean(nchar(GENOME_INFO))<7){
      Ts <- grep("T\\d{5}",GENOME_INFO)              # T number identifier
      Ls <- grep("[[a-z]{3}|[a-z]43}]",GENOME_INFO)  # 3 or 4 letter code identifier
      if(length(Ts)==0&&length(Ls)==0) stop("GENOME_INFO - unacceptable genome identifier. Must be a T number (eg 'T00001') or a 3 or 4 letter code (eg 'eco').")
      # GET GENES
      GENES     <- NULL
      ref_index <- NULL
      if(length(Ts)>0){
        T_index <- na.omit(match(GENOME_INFO[Ts],GENOME_REFERENCE_TABLE$ID))
        T_missing <- attr(T_index,which = "na.action")
        if(length(T_missing)>0) warning(paste("'",paste(GENOME_INFO[T_missing],collapse = ","),"' NOT matched to KEGG genome database. Exact match only possible.",sep=""))
        ref_index <- c(ref_index,T_index)
      }
      if(length(Ls)>0){
        L_index <- na.omit(match(GENOME_INFO[Ls],GENOME_REFERENCE_TABLE$ORG_ID))
        L_missing <- attr(L_index,which = "na.action")
        if(length(L_missing)>0) warning(paste(paste(GENOME_INFO[L_missing],collapse = ","),"NOT matched to KEGG genome database. Exact match only possible."))
        ref_index <- c(ref_index,L_index)
      }

      # ID COLUMN
      id_col    <- ifelse(length(Ts)>length(Ls),1,2)

      # FIND

    }else{
      # ASSUME NAMES
      ref_index <- NULL
      for(N in GENOME_INFO){
        INDEX <- grep(N,GENOME_REFERENCE_TABLE$ORGANISM,ignore.case = T)
        if(length(INDEX)>0){
          ref_index <- c(ref_index,INDEX)
        }else{
          warning(paste('No matches (partial or otherwise) of name:',N))
        }
      }
      id_col <- 1
    }

    if(length(ref_index)>0){
      ref_index        <- unique(ref_index)
      GENOME_INFO_DATA <- cbind(GENOME_REFERENCE_TABLE$ID[ref_index],
                                GENOME_REFERENCE_TABLE$ORG_ID[ref_index],
                                GENOME_REFERENCE_TABLE$ORGANISM[ref_index],
                                GENOME_REFERENCE_TABLE$KOs[ref_index])
      colnames(GENOME_INFO_DATA) <- c("ID","ORG_ID","ORGANISM","KOs")

            genes_col <- 4
      splitBy   <- "[;]"

      # CHECK DATA AVAILABILITY
      no_info_index <- which(GENOME_INFO_DATA[,genes_col]=="")
      if(length(no_info_index)>0) GENOME_INFO_DATA <- GENOME_INFO_DATA[-no_info_index,]
    }else{
      stop('No matching identifiers or names found.')
    }
  }

  # USE USER DATA
  if(is.data.frame(GENOME_INFO)){

    # FIND GENOME ID COLUMNS
    id_col <- ifelse(is.numeric(GENOME_ID_COL),GENOME_ID_COL,which(names(GENOME_INFO)==GENOME_ID_COL))
    if(length(id_col)==0) stop("No matching 'GENOME_ID_COL' column found in 'GENOME_INFO'. 'GENOME_ID_COL' must be the column number or name containing the genome/organism ID.")
    if(length(id_col)>1)  stop("Multiple matching 'GENOME_ID_COL' column found in 'GENOME_INFO'. 'GENOME_ID_COL' must refer to a single column number or name containing the genome/organism ID.")

    # FIND GENES COLUMN IF NEEDED
    genes_col <- ifelse(is.numeric(GENES_COL),GENES_COL,which(names(GENOME_INFO)==GENES_COL))
    if(length(genes_col)==0) stop("No matching 'GENES_COL' column found in 'GENOME_INFO'. 'GENES_COL' must be the column number or name containing the KO/EC numbers.")
    if(length(genes_col)>1)  stop("Multiple matching 'GENES_COL' column found in 'GENOME_INFO'. 'GENES_COL' must refer to a single column number or name containing the KO/EC numbers.")

    GENOME_INFO_DATA <- GENOME_INFO[,c(id_col,genes_col)]
    id_col    <- 1
    genes_col <- 2

    # CHECK DATA AVAILABILITY
    no_info_index <- which(GENOME_INFO_DATA[,genes_col]=="")
    if(length(no_info_index)>0){
      GENOME_INFO_DATA <- GENOME_INFO_DATA[-no_info_index,]
      warning(paste("Exlcuding ROWS",paste(no_info_index,collapse = ","),"in GENOME_INFO: no gene information"))
    }
  }

  #### PROCESS INPUT DATA ----
  # GET [CONCATENATED] LISTS OF GENES
  # geneString_vector <- GENOME_INFO_DATA[,genes_col]

  # CHECK IF DUPLICATED DATASET NAMES PROVIDED
  DATASETS <- misc_check_duplicate_names(GENOME_INFO_DATA[,id_col])

  # CHECK IF '/' PRESENT IN DATASET NAMES AND SUBSTITUTE
  DATASETS <- gsub("/","--",DATASETS)


  ### DETERMINE QUERY ####
  # MODULE_ID <- SEARCH_NAME <- SEARCH_CLASS_I <- SEARCH_CLASS_II <- SEARCH_CLASS_III <- EXCLUDE_NAME <- ""

  query_terms <- c("MODULE_ID", "SEARCH_NAME", "SEARCH_CLASS_I", "SEARCH_CLASS_II", "SEARCH_CLASS_III", "EXCLUDE_NAME")
  QUERIES     <- data.frame("SEARCH_TERMS" = query_terms,"ARGUMENTS" = "",stringsAsFactors = F)

  QUERIES$SEARCH_TERMS[1] <- "MODULE_ID";        QUERIES$ARGUMENTS[1] <- paste(MODULE_ID,collapse = ",");
  QUERIES$SEARCH_TERMS[2] <- "SEARCH_NAME";      QUERIES$ARGUMENTS[2] <- paste(SEARCH_NAME,collapse = ",");
  QUERIES$SEARCH_TERMS[3] <- "SEARCH_CLASS_I";   QUERIES$ARGUMENTS[3] <- paste(SEARCH_CLASS_I,collapse = ",");
  QUERIES$SEARCH_TERMS[4] <- "SEARCH_CLASS_II";  QUERIES$ARGUMENTS[4] <- paste(SEARCH_CLASS_II,collapse = ",");
  QUERIES$SEARCH_TERMS[5] <- "SEARCH_CLASS_III"; QUERIES$ARGUMENTS[5] <- paste(SEARCH_CLASS_III,collapse = ",");
  QUERIES$SEARCH_TERMS[6] <- "EXCLUDE_NAME";     QUERIES$ARGUMENTS[6] <- paste(EXCLUDE_NAME,collapse = ",");

  if(sum(QUERIES$ARGUMENTS[1:5]=="")==5) QUERIES<- rbind(QUERIES,c("* DEFAULT SEARCH","ALL MODULES *"))
  if(QUERIES$ARGUMENTS[6]=="")           QUERIES<- rbind(QUERIES,c("* NO MODULES","EXCLUDED *"))

  QUERIES$ARGUMENTS[which(QUERIES$ARGUMENTS=="")] <- "< >"

  #### ANALYSIS ----
  MATRIX        <- NULL
  METADATA      <- NULL
  ADD_INFO      <- NULL
  remove_index  <- NULL

  for(i in 1:nrow(GENOME_INFO_DATA)){
    # Retrieve vector of genes
    gene_vector <- strsplit(GENOME_INFO_DATA[i,genes_col],split = splitBy)[[1]]

    # Analyse vector of genes
    try(OUT <- misc_geneVector_module(gene_vector,MODULE_ID,SEARCH_NAME,SEARCH_CLASS_I,SEARCH_CLASS_II,SEARCH_CLASS_III,EXCLUDE_NAME,...))
    # try(OUT <- misc_geneVector_module(gene_vector,MODULE_ID,SEARCH_NAME,SEARCH_CLASS_I,SEARCH_CLASS_II,SEARCH_CLASS_III,EXCLUDE_NAME))

    if(!exists("OUT")){
      remove_index <- c(remove_index,i)
      next
    }

    ## COLLECT ALL INFORMATION ----
    # FRACTION MATCHED
    MATRIX   <- rbind(MATRIX,OUT$FRACTION)
    # ADDITIONAL INFO
    if(ADD_OUT){
      ADD_INFO   <- rbind(ADD_INFO,cbind(rep(GENOME_INFO_DATA[i,id_col],nrow(OUT$ADD_INFO)),OUT$ADD_INFO))
    }
    # FORMAT/RETRIEVE
    if(i==1){
      # MATRIX
      if(OUT_MODULE_NAME){
        colnames(MATRIX)  <- OUT$METADATA$NAME_SHORT
      }else{
        colnames(MATRIX)  <- colnames(OUT$FRACTION)
      }

      # METADATA
      if(META_OUT) METADATA <- OUT$METADATA
    }
    rm(OUT)
  }
  ### PROCESS OUTPUTs ----
  # ASIGN COL AND ROW NAMES
  if(length(remove_index)>0){
    rownames(MATRIX)        <- DATASETS[-remove_index]
  }else{
    rownames(MATRIX)        <- DATASETS
  }
  if(ADD_OUT)
    colnames(ADD_INFO)[1] <- "GENOME_ID"

  ## Format data as data frame
  if(META_OUT){
    METADATA           <- as.data.frame(METADATA, stringAsFactor = F)
  }

  if(ADD_OUT){
    ADD_INFO           <- as.data.frame(ADD_INFO, stringAsFactor = F)
    ADD_INFO$FRACTION  <- as.numeric(ADD_INFO$FRACTION)
  }

  ### OUTPUT ----
  output <- list("MATRIX"   = MATRIX,
                 "QUERIES"  = QUERIES)

  if(META_OUT) output$METADATA = METADATA
  if(ADD_OUT)  output$ADD_INFO = ADD_INFO

  # PROVIDE RELATIONSHIP OR INPUT DATA USED
  if(is.character(GENOME_INFO)){
    GENOME_INFO_DATA <- as.data.frame(GENOME_INFO_DATA,stringsAsFactors = F)
    GENOME_INFO_DATA <- dplyr::select(GENOME_INFO_DATA,-KOs)
  }
  output$GENOME_INFO_DATA <- GENOME_INFO_DATA

  return(output)
}

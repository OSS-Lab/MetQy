#' Map a list of EC or K numbers to KEGG modules.
#'
#' @description
#' This function maps the list of Enzyme Commission (EC) numbers or KEGG orthologs (specified as K numbers) given by \code{gene_vector} to the
#' KEGG modules using an in-build reference table (\code{module_reference_table}) of all the KEGG modules in the database at the time of release.
#'
#' @details
#' The modules to be analysed are determined by search queries.
#' These can search across the module NAME, CLASS (I-III) or by specifying the module ID (M number). See the Argument section.
#' When none of the optional arguments are specified, the default is to search all modules.
#' An additional agument, \code{EXCLUDE_NAME}, can be used to exclude modules with a certain term in the \code{NAME} field.
#' For example, specifying \code{EXCLUDE_NAME} = "biosynthesis" would return the search across all modules,
#' except those that contain "biosynthesis" in the name.
#' For instance, module "M00005" named "PRPP biosynthesis, ribose 5P => PRPP" and 168 others would be excluded form the search.
#'
#' @param gene_vector      - character vector listing EC or K numbers. If more K numbers are given, the K number-based definition will be used. See Value.
#'
#' @param MODULE_ID        - character vector listing specific module IDs (e.g. M00001).
#'
#' @param SEARCH_NAME      - character vector listing terms to search in \code{NAME} field (case-insensitive).
#'
#' @param SEARCH_CLASS_I   - character vector listing terms to search in \code{CLASS_I} field (case-insensitive).
#'
#' @param SEARCH_CLASS_II  - character vector listing terms to search in \code{CLASS_II} field (case-insensitive).
#'
#' @param SEARCH_CLASS_III - character vector listing terms to search in \code{CLASS_III} field (case-insensitive).
#'
#' @param EXCLUDE_NAME     - character vector listing terms that if matched in \code{NAME} field will be excluded (case-insensitive).
#'
#' @param use_module_reference_table - optional. Provide a data frame with updated KEGG module database OR with costum-made modules.
#'                                     Default (\code{NULL}; inbulit data used). See Details.
#'
#' @param ...              - further arguments such as \code{KO_in_DEF_EC}. See Details.
#'
#' @details If none of the search terms are specified, the default is to analyse all modules.
#'
#' The \code{use_} argument allow users with KEGG FTP access to provide the updated data from the KEGG databases in the form of reference tables AND/OR for advanced users to provide contum-made modules (see below).
#' These reference table can be generated with \code{parseKEGG_module} and need to have a specific format (see function descriptions for details on format).
#'
#' The module definition (contained in \code{module_reference_table}) describes the relationship between genes and modules
#' and is used to identify the modules in which \code{gene} is involved.
#' The user can provide costum-made module definitions that use the logical expression format (however, the table format must be conserved!).
#'
#' @return List:
#' \code{$FRACTION} matrix containing the completeness fraction for the modules matched in the query along the columns.
#'
#' \code{$METADATA} data frame with the metadata for the modules matched in the query.
#'  Columns:\preformatted{
#' 	     (1)     MODULE_ID         - M number (e.g. M00001);
#' 	     (2)     MODULE_NAME       - the KEGG given name for the module;
#' 	     (3)     MODULE_NAME_SHORT - manually shortened module name (for plotting purposes);
#' 	 (4 - 6)     CLASS_I - III     - hierarchical module classes;
#' 	     (7)     DEFINITION        - KEGG module definition in terms of K or EC numbers (without optional K or EC numbers);
#' 	     (8)     OPTIONAL          - the optional K or EC numbers that are part of the KEGG module definition (NA otherwise);
#' 	     (9)     KOs_IN_DEF        - logical flag indicating whether there are K numbers involved in that module definition (if EC numbers are given);
#' }
#'
#' \code{$ADD_INFO}  data frame with additional output information.
#'  Columns:\preformatted{
#' 	 (1) MODULE_ID           - M number (e.g. M00001);
#' 	 (2) MODULE_NAME_SHORT   - manually shortened module name (for plotting purposes);
#' 	 (3) FRACTION            - the number of complete blocks divided by the total number of blocks;
#' 	 (4) nBLOCKS             - the number of blocks that make up the KEGG module;
#' 	 (6) COVERAGE            - the K numbers or ECs that are present and involved in the KEGG module definition;
#' 	 (7) OPTIONAL_PRESENT    - the K numbers or ECs that are present and are listed as optional;
#' }
#'
#' See \link{parseKEGG_module} or visit \href{http://www.genome.jp/kegg/module.html}{http://www.genome.jp/kegg/module.html} for more information.
#'
#' @examples
#' output_module_ECs <- misc_geneVector_module(data_exampleECnumbers_vector,MODULE_ID=paste("M0000",1:5,sep=""))
#'
#' output_module_KOs <- misc_geneVector_module(data_exampleKOnumbers_vector,MODULE_ID=paste("M0000",1:5,sep=""))
#' @seealso
#' \link{parseKEGG_module}
#'
#' @export

# Previously called module_mapping and query_ECs_KOs_to_module
############################################################################################################################################

misc_geneVector_module  <- function(gene_vector,MODULE_ID,SEARCH_NAME,SEARCH_CLASS_I,SEARCH_CLASS_II,SEARCH_CLASS_III,EXCLUDE_NAME,use_module_reference_table = NULL,...){

#### MANAGE INPUT ----
  stopifnot(is.character(gene_vector))
  if(length(gene_vector)==1) warning("Single value in gene_vector! Check the 'split_vector_by' argument or use 'query_gene_to_module.")

  if(!exists("MODULE_ID"))        MODULE_ID        <- NULL
  if(!exists("SEARCH_NAME"))      SEARCH_NAME      <- NULL
  if(!exists("SEARCH_CLASS_I"))   SEARCH_CLASS_I   <- NULL
  if(!exists("SEARCH_CLASS_II"))  SEARCH_CLASS_II  <- NULL
  if(!exists("SEARCH_CLASS_III")) SEARCH_CLASS_III <- NULL
  if(!exists("EXCLUDE_NAME"))     EXCLUDE_NAME     <- NULL

  # CHECK FOR MODULE REFERENCE TABLE
  if(!is.null(use_module_reference_table)&&is.data.frame(use_module_reference_table)){
    MODULE_REFERENCE_TABLE <- use_module_reference_table
  } else if(!is.null(use_module_reference_table)&&!is.data.frame(use_module_reference_table)){
    warning("'use_module_reference_table' should be a data frame. Using inbuilt data.")
    MODULE_REFERENCE_TABLE <- module_reference_table
  }else{
    MODULE_REFERENCE_TABLE <- module_reference_table
  }

  ### DEAL WITH M numbers in DEFINITION ----
  # POSSIBLE FUTURE IMPLEMEMTATION TO BE ABLE TO EVALUATE MODULES DEFINED IN TERMS OF OTHER MODULES
  exclude_M_def_modules <-  T

  # REMOVE MODULES THAT ARE DEFINED IN TERMS OF OTHER MODULES FROM REFERENCE TABLES
  #     USING DEFINITION_KOs AS IT IS THE KEGG DEFINED DEFINITION
  if(exclude_M_def_modules){
    module_def_index <- grep("M\\d{5}",MODULE_REFERENCE_TABLE$DEFINITION_KOs)
    if(length(module_def_index)>0) MODULE_REFERENCE_TABLE <- MODULE_REFERENCE_TABLE[-module_def_index,]
  }

#### GET MATCHING INDICES ----

  # IF NO SEARCH TERMS HAVE BEEN SPECIFY, DO THE ANALYIS FOR ALL MODULES
  if((MODULE_ID==""||is.null(MODULE_ID))&&(SEARCH_NAME==""||is.null(SEARCH_NAME))&&(SEARCH_CLASS_I==""||is.null(SEARCH_CLASS_I))&&
     (SEARCH_CLASS_II==""||is.null(SEARCH_CLASS_II))&&(SEARCH_CLASS_III==""||is.null(SEARCH_CLASS_III))){
    match_index <- 1:nrow(MODULE_REFERENCE_TABLE)
  } else {
    match_index <- NULL

    if(length(MODULE_ID)>1||MODULE_ID!=""){
      # GET THE INDEX OF MODULE ENTRIES THAT MATCH THE MODULE ID
      for(MOD in MODULE_ID){
        if(length(grep(MOD,MODULE_REFERENCE_TABLE$ID))==0){
          warning(paste("No Module ID matches:", MOD))
          next
        }
        match_index <- c(match_index,grep(MOD,MODULE_REFERENCE_TABLE$ID))
      }
    }

    # GET THE INDEX OF MODULE NAME ENTRIES THAT MATCH SEARCH_NAME
    if(length(SEARCH_NAME)>1||SEARCH_NAME!=""){
      for(S in SEARCH_NAME){
        if(length(grep(S,MODULE_REFERENCE_TABLE$NAME,ignore.case = T))==0){
          warning(paste("No Module NAME matches:", S))
          next
        }
        match_index <- c(match_index,grep(S,MODULE_REFERENCE_TABLE$NAME,ignore.case = T))
      }
    }

    # GET THE INDEX OF MODULE CLASS_II ENTRIES THAT MATCH SEARCH_CLASS
    if(length(SEARCH_CLASS_I)>1||SEARCH_CLASS_I!=""){
      for(S in SEARCH_CLASS_I){
        if(length(grep(S,MODULE_REFERENCE_TABLE$CLASS_I,ignore.case = T))==0){
          warning(paste("No Module CLASS_I matches:", S))
          next
        }
        match_index <- c(match_index,grep(S,MODULE_REFERENCE_TABLE$CLASS_I,ignore.case = T))
      }
    }
    if(length(SEARCH_CLASS_II)>1||SEARCH_CLASS_II!=""){
      for(S in SEARCH_CLASS_II){
        if(length(grep(S,MODULE_REFERENCE_TABLE$CLASS_II,ignore.case = T))==0){
          warning(paste("No Module CLASS_II matches:", S))
          next
        }
        match_index <- c(match_index,grep(S,MODULE_REFERENCE_TABLE$CLASS_II,ignore.case = T))
      }
    }
    if(length(SEARCH_CLASS_III)>1||SEARCH_CLASS_III!=""){
      for(S in SEARCH_CLASS_III){
        if(length(grep(S,MODULE_REFERENCE_TABLE$CLASS_III,ignore.case = T))==0){
          warning(paste("No Module CLASS_III matches:", S))
          next
        }
        match_index <- c(match_index,grep(S,MODULE_REFERENCE_TABLE$CLASS_III,ignore.case = T))
      }
    }
  }
  # EXCLUDE DUPLICATED MODULE IDs
  match_index <- unique(sort(match_index))

  # EXCLUDE_NAME
  if(length(EXCLUDE_NAME)>1||EXCLUDE_NAME!=""){
    exclude <- NULL
    for(E in EXCLUDE_NAME){

      exclude <- c(exclude,grep(E,MODULE_REFERENCE_TABLE$NAME,ignore.case = T))
    }
    match_index <- setdiff(match_index,exclude)
  }

  if(length(match_index)==0) stop("No module matches")

#### DETERMINE KO or EC AND GET DEFINITION AND OPTIONAL ----

  ## DETERMINE IF INPUT IS KO OR EC NUMBERS ----
  nKOs <- length(grep("^K\\d{5}$",gene_vector))
  nECs <- length(grep("^\\d{1}[.][[:punct:] [:digit:]]+$",gene_vector))


  if(nKOs==0&&nECs==0){
    stop("gene_vector should contain a vector with EC numbers (e.g. 1.-.-.-) or with KEGG ortholog IDs (e.g. K00001)")
  }

  #GIVE A WARNING IF THERE ARE ANY INVALID GENE IDS
  if(nKOs + nECs < length(gene_vector)) {
    error_gene <- gene_vector[!grepl("(^K\\d{5}$|^\\d{1}[.][[:punct:] [:digit:]]+$)",gene_vector)][1] #get the first invalid gene ID
    warning(paste("'",error_gene,"'","is not a valid KO or EC number. Genes should be EC numbers (e.g. 1.-.-.-) or KEGG ortholog IDs (e.g. K00001) "), sep=" ")
  }

  if(nKOs>=nECs){
    # DETERMINE LOOK UP SETS
    DEFINITION <- MODULE_REFERENCE_TABLE$DEFINITION_KOs
    OPTIONAL   <- MODULE_REFERENCE_TABLE$OPTIONAL_KOs
    EC         <- F
  } else {
    # REMOVE MODULES WHOSE EC DEFINITION HAS NO ECs PRESENT
    noECs_index <- setdiff(1:nrow(MODULE_REFERENCE_TABLE),grep("EC\\d",MODULE_REFERENCE_TABLE$DEFINITION_ECs))
    if(length(noECs_index)>0) MODULE_REFERENCE_TABLE <- MODULE_REFERENCE_TABLE[-noECs_index,]

    # DETERMINE LOOK UP SETS
    DEFINITION <- MODULE_REFERENCE_TABLE$DEFINITION_ECs
    OPTIONAL   <- MODULE_REFERENCE_TABLE$OPTIONAL_ECs
    EC         <- T

    # ADD UNSPECIFIED ENZYMES TO LIST
    EC_unspecified  <- t(as.data.frame(strsplit(gene_vector,split = "[.]")))
    EC_unspecified  <- unique(as.character(apply(EC_unspecified[,1:3],1,paste,collapse = ".")))
    gene_vector    <- c(gene_vector,EC_unspecified)         # add unspecified enzymes as having only the first 3 positions
    if(length(grep("-",gene_vector)>0)) gene_vector    <- gene_vector[-grep("-",gene_vector)]  # remove enzymes wih ".-"
  }

  ## FORMAT DEFINITION ----
  # TRIM END SPACES AND REMOVE OTHER CHARACTERS
  DEFINITION <- gsub("\\s+$","",DEFINITION)
  DEFINITION <- gsub("[;]","",DEFINITION)

### STORAGE MODULE METADATA AND ADDITIONAL INFO ----
  metadata  <- data.frame(MODULE_ID        = MODULE_REFERENCE_TABLE$ID[match_index],
                            MODULE_NAME        = MODULE_REFERENCE_TABLE$NAME[match_index],
                            NAME_SHORT         = MODULE_REFERENCE_TABLE$NAME_SHORT[match_index],
                            CLASS_I            = MODULE_REFERENCE_TABLE$CLASS_I[match_index],
                            CLASS_II           = MODULE_REFERENCE_TABLE$CLASS_II[match_index],
                            CLASS_III          = MODULE_REFERENCE_TABLE$CLASS_III[match_index],
                            DEFINITION         = DEFINITION[match_index],
                            OPTIONAL           = OPTIONAL[match_index],
                            stringsAsFactors   = F)
  if(EC){
    metadata$KOs_IN_DEF <- 0
    metadata$KOs_IN_DEF[grep("K\\d{5}",DEFINITION[match_index])] <- 1
  }

  additional_info  <- data.frame(MODULE_ID           = MODULE_REFERENCE_TABLE$ID[match_index],
                                  NAME_SHORT         = MODULE_REFERENCE_TABLE$NAME_SHORT[match_index],
                                  FRACTION           = numeric(length(match_index)),
                                  nBLOCKS            = character(length(match_index)),
                                  COVERAGE           = character(length(match_index)),
                                  OPTIONAL_PRESENT   = character(length(match_index)),
                                  stringsAsFactors   = F)

#### PREALLOCATE QUERY OUTPUT ----
  FRACTION           <- matrix(nrow = 1,ncol = length(match_index))
  colnames(FRACTION) <- MODULE_REFERENCE_TABLE$ID[match_index]


#### MODULE ANALYSIS ----
  for(M in 1:length(match_index)){

    if(is.na(DEFINITION[match_index[M]])) next

  #### GET PERCENTAGE MATCHING ----
    blocks  <- strsplit(gsub(" +"," ",DEFINITION[match_index[M]]),split = " ")[[1]] # remove multiple spaces
    present <- numeric(length(blocks))

    for(B in 1:length(blocks)) present[B] <- misc_evaluate_block(gene_vector,blocks[B],...)
    # for(B in 1:length(blocks)) present[B] <- misc_evaluate_block(gene_vector,blocks[B])

  ### GET FRACTION ----
    FRACTION[1,M] <- sum(present)/length(present)

  ### RECORD BLOCK INFO ----
    additional_info$nBLOCKS[M] <- length(blocks)

  ### FIND OPTIONALS PRESENT ----
    if(!is.na(OPTIONAL[match_index[M]])){
      thisOptions       <- strsplit(OPTIONAL[match_index[M]],split = ";")[[1]]
      OPTIONAL_PRESENT  <- paste(thisOptions[is.na(match(thisOptions,gene_vector))==F],collapse = ";")
      additional_info$OPTIONAL_PRESENT[M] <- OPTIONAL_PRESENT
    }

  ### GET COVERAGE ----
    if(EC){
      EC_data_names               <- unique(gsubfn::strapplyc(X = DEFINITION[match_index[M]],pattern = "(\\d{1}[.][[:punct:]||\\d]+)",simplify = T))
      additional_info$COVERAGE[M] <- paste(sort(gene_vector[na.omit(match(EC_data_names,gene_vector))]),collapse = ";")
    }else{
      KO_data_names               <- unique(gsubfn::strapplyc(X = DEFINITION[match_index[M]],pattern = "(K\\d{5})",simplify = T))
      additional_info$COVERAGE[M] <- paste(sort(gene_vector[na.omit(match(KO_data_names,gene_vector))]),collapse = ";")
    }
  } # end module matching

  # RECORD FRACTION
  additional_info$FRACTION <- FRACTION[1,]

#### OUTPUT ----
  OUTPUT <- list(FRACTION = FRACTION)
  OUTPUT$METADATA <- metadata
  OUTPUT$ADD_INFO <- additional_info
  # OUTPUT$QUERY <- QUERY
  return(OUTPUT)
}

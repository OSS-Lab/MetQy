#' Evaluate a KEGG module block definition.
#'
#' Evaluate a KEGG module block definition given a vector of genes (KEGG Orthologs -KOs-, identified with K number).
#' Some KOs can be mapped to Enzyme Commission (EC) numbers.
#'
#' @param gene_vector  - vector containing either K or EC numbers. See Details.
#'
#' @param BLOCK        - KEGG module DEFINITION BLOCK. See \link{parseKEGG_module}.
#'
#' @param KO_in_DEF_EC - logical. If enzyme IDs are given, should lingering K numbers in the module DEFINITION
#'                       be assumed to be present? Default (\code{FALSE}).
#'
#' @param ...  - further arguments (currently unsupported)
#'
#' @details \code{gene_vector} must contain either K or EC numbers.
#'
#' @export

misc_evaluate_block <- function(gene_vector, BLOCK, KO_in_DEF_EC = FALSE,...){

  EC_KO_data_names <- BLOCK

  # DETERMINE IF GENES ARE KOs or ECs
  nKOs <- length(grep("K\\d{5}",gene_vector))
  nECs <- length(grep("\\d{1}[.][[:punct:] [:digit:]]+",gene_vector))

  if(nKOs==0&&nECs==0)
    stop("gene_vector should contain a vector with EC numbers (e.g. 1.-.-.-) or with KEGG ortholog IDs (e.g. K00001)")

  if(nKOs>=nECs){
    EC <- F
  }else{
    EC <- T
  }

  # FORMAT STRING TO BE ABLE TO EXTRACT EC numbers
  if(EC){
    stopifnot(is.logical(KO_in_DEF_EC)&&length(KO_in_DEF_EC)==1)
    # EXTRACT DEFINITION ELEMENTS
    EC_KO_data_names <- gsub("&"," ",gsub("|"," ",EC_KO_data_names,fixed = T),fixed = T)
    EC_KO_data_names <- gsub("("," ",gsub(")"," ",EC_KO_data_names,fixed = T),fixed = T)

    # GET THE UNIQUE LIST OF ENTRIES
    EC_data_names <- unique(gsubfn::strapplyc(X = EC_KO_data_names,pattern = "(\\d{1}[.][[:punct:]||\\d]+)",simplify = T))
    KO_data_names <- unique(gsubfn::strapplyc(X = EC_KO_data_names,pattern = "(K\\d{5})",simplify = T))

    # MANAGE ABSENT KOs or ECs IN DEFINITION
    if(is.list(EC_data_names)){
      EC_data_names <- NULL
    }else {
      EC_data_names <- paste("EC",EC_data_names,sep="")
    }
    if(is.list(KO_data_names)){
      KO_data_names <- NULL
      KO_present    <- FALSE
    }else{
      KO_present    <- TRUE
    }

    # MAKE DATA FRAME WITH UNIQUE LIST OF ENTRIES AS VARIABLE NAME FOR LOGICAL TEST
    if(is.null(EC_data_names)){
      EC_data <-NULL
    }else{
      EC_data           <- matrix(nrow = 1,ncol = length(EC_data_names))
      colnames(EC_data) <- EC_data_names
      # DEAL WITH UNSPECIFIED ENZYMES IN DEFINITION
      colnames(EC_data) <- gsub(".-","",colnames(EC_data))
    }

    # POPULATE DATA FRAME WITH LOGIC VALUES CORESSPONDING TO THE PRESCENCE OF THE KO OR EC NUMBER (ENTRIES)
    EC_data[1,] <- is.na(match(gsub("EC","",colnames(EC_data)),gene_vector))==F

    # ADD KEGG ORTHOLOGS THAT DO NOT HAVE AN ASSOCIATED EC number to the data frame
    if(is.null(KO_data_names)){
      KO_data <- NULL
    } else {
      KO_data <- matrix(data = KO_in_DEF_EC,nrow = 1,ncol = length(KO_data_names))
      colnames(KO_data) <- KO_data_names
    }

    EC_KO_data <- cbind(EC_data,KO_data)
    EC_KO_data <- as.data.frame(EC_KO_data,stringsAsFactors = F)

    # EVALUATE LOGICAL EXPRESSION TO DETERMINE WHETHER THE MODULE IS PRESENT OR NOT
    # ASSUME KEGG ORTHOLOGS THAT ARE PART OF A DEFINITION ARE PRESENT (IF KOs are present, otherwised ignored)
    present <- eval(expr = parse(text = gsub(".-","",BLOCK)),EC_KO_data)

    # # ASSUME KEGG ORTHOLOGS THAT ARE PART OF A DEFINITION ARE NOT PRESENT
    # if(KO_present){
    #   EC_KO_data[1,grep("K\\d{5}",colnames(EC_KO_data))] <- FALSE
    #   present_F[B] <- eval(expr = parse(text = gsub(".-","",BLOCK)),EC_KO_data)
    # } else {
    #   # EVALUATE LOGICAL EXPRESSION TO DETERMINE WHETHER THE MODULE IS PRESENT OR NOT
    #   # ASSUME KEGG ORTHOLOGS THAT ARE PART OF A DEFINITION ARE NOT PRESENT
    #   present_F[B] <- present
    # }

    # # DETERMINE WHICH ECs ARE PRESENT
    # ECs_needed  <- gsubfn::strapplyc(EC_KO_data_names,pattern = "(\\d{1}[.][[:punct:]||\\d]+)",simplify = T)
    #
    # if(is.list(ECs_needed)==F&&length(ECs_needed)>0){
    #   ECs_needed  <- unique(gsub(".-","",ECs_needed,fixed = T))
    #   ECs_matched <- ECs_needed[as.logical(is.na(match(ECs_needed,gene_vector))==F)]
    #   if(length(ECs_matched)>0) ECs_present <- c(ECs_present,ECs_matched)
    # }

    # FINISH if(EC)
  } else {
    ### K numbers - based definition
    # EC_KO_data_names  <- as.character(unique(gsubfn::strapplyc(X = EC_KO_data_names,pattern = "(K\\d{5})",simplify = T)))
    EC_KO_data_names  <- unique(gsubfn::strapplyc(X = EC_KO_data_names,pattern = "(K\\d{5})")[[1]])
    EC_KO_data        <- data.frame(matrix(nrow = 1,ncol = length(EC_KO_data_names)),stringsAsFactors = F)
    names(EC_KO_data) <- EC_KO_data_names
    EC_KO_data[1,]    <- is.na(match(EC_KO_data_names,gene_vector))==F

    # EVALUATE LOGICAL EXPRESSION TO DETERMINE WHETHER THE MODULE IS PRESENT OR NOT
    names(EC_KO_data) <- gsub(".-","",names(EC_KO_data))  # MANAGE EC NUMBERS WITH '.-'
    present           <- eval(expr = parse(text = gsub(".-","",BLOCK)),EC_KO_data)

    # # DETERMINE WHICH KOs ARE PRESENT
    # # KOs_needed  <- gsubfn::strapplyc(EC_KO_data_names,pattern = "(K\\d{5})",simplify = T)
    # KOs_needed  <- gsubfn::strapplyc(EC_KO_data_names,pattern = "(K\\d{5})")[[1]]
    #
    # if(is.list(KOs_needed)==F&&length(KOs_needed)>0){
    #   KOs_needed  <- unique(KOs_needed)
    #   KOs_matched <- KOs_needed[as.logical(is.na(match(KOs_needed,gene_vector))==F)]
    #   if(length(KOs_matched)>0) KOs_present <- c(KOs_present,KOs_matched)
    # }
  } # end -->  if KO
  present <- as.numeric(present)
  return(present)
}

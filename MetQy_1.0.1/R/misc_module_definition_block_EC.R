#' KEGG module definition processing: subblocks
#'
#' Format the KEGG module database - Format the KEGG module DEFINITION for ease of analysis by excluding optional KEGG orthologs.
#' Used by \code{parseKEGG_module()}.
#'
#' @param BLOCK  - string containing a block of a KEGG module DEFINITION (logical expression using K numbers).
#'
#' @param ORTHOLOGS  - vector listing the K numbers and related EC numbers.
#'
#' @details The KEGG module definition uses optional KEGG orthologs, indicated by a "-". These are removed and the corresponding EC number stored.
#'
#' @seealso \link{parseKEGG_module}
#' @export

############################################################################################################################################

misc_module_definition_block_EC <- function(BLOCK,ORTHOLOGS){

  BLOCK_EC <-BLOCK
  # REPLACE KOs with EC numbers
  KOs   <- gsubfn::strapplyc(BLOCK_EC,"(K\\d{5})",simplify = T)
  for(K in 1:length(KOs)){
    thisEC   <- unique(gsubfn::strapplyc(ORTHOLOGS[grep(KOs[K],ORTHOLOGS)],"(\\d{1}[.][[:punct:]||\\d]+)",simplify = T))

    # IF THERE ARE CORRESPONDING EC NUMBERS ELSE THE KO ID will remain in definition
    if(is.list(thisEC)==F){
      BLOCK_EC <- gsub(KOs[K],paste(paste("EC",thisEC,sep=""),collapse = "|"),BLOCK_EC,fixed = T)
    }
  }

  # OUTPUT
  return(BLOCK_EC)
}

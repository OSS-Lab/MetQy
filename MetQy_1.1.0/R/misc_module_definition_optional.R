#' KEGG module definition processing: optional KEGG ortholog exclusion
#'
#' Format the KEGG module database - Format the KEGG module definition for ease of analysis by excluding optional KEGG orthologs.
#' Used by \code{parseKEGG_module()}.
#'
#' The KEGG module definition uses optional KEGG orthologs, indicated by a "-". These are removed and the corresponding EC number stored.
#'
#' @param BLOCK  - string containing a block of a KEGG module DEFINITION (logical expression using K numbers).
#'
#' @param ORTHOLOGS  - vector listing the K numbers and related EC numbers.
#'
#' @param ...  - further arguments (currently unsupported)
#'
#' @return A list with the formatted block (\code{$thisBlock}) and any optional K and EC numbers (\code{$thisOptional_KO} and \code{$thisOptional_EC})
#'
#' @seealso \link{parseKEGG_module}
#' @export

############################################################################################################################################

misc_module_definition_optional <- function(BLOCK,ORTHOLOGS,...){

  # EXCLUDE optional KOs (indicated by -K\\d{5})
    thisOptional_KO <- NULL
    thisOptional_EC <- NULL

    if(length(grep("-K",BLOCK,fixed = T))>0){
      replace_KOs_1 <- gsub("-","",gsubfn::strapplyc(BLOCK,"(-K\\d{5})",simplify = T))
      for(Opt in 1:length(replace_KOs_1)){
        thisOptional_KO <- c(thisOptional_KO,replace_KOs_1[Opt])
        thisOptEC <- unique(gsubfn::strapplyc(ORTHOLOGS[grep(replace_KOs_1[Opt],ORTHOLOGS)],"(\\d{1}[.][[:punct:]||\\d]+)",simplify = T))

        if(is.character(thisOptEC)) thisOptional_EC <- c(thisOptional_EC,thisOptEC)
      }
    }
    ## OPTAIN THIS BLOCK
    thisBlock <- gsub("-K\\d{5}","",BLOCK)

    # DEAL WITH ADDITIONAL OPTIONAL SUBGROUPS
    if(length(grep("[-]",thisBlock))>0){
        # Deal with nested groups
        if(length(grep("[(]",thisBlock))>0){
          optional_groups <- gsubfn::strapplyc(thisBlock,"(-[(].*?[)])",simplify = T)
          nOptSubgroups <- length(optional_groups)

          for(OSG in 1:nOptSubgroups){
            thisOptional_group <- optional_groups[OSG]
            # REMOVE OPTIONAL GROUP FROM BLOCK
            thisBlock           <- gsub(thisOptional_group,"",thisBlock,fixed = T)

            # ADD OPTIONALS TO VECTORS
            replace_KOs_2   <- gsubfn::strapplyc(thisOptional_group,"(K\\d{5})",simplify = T)
            for(R in 1:length(replace_KOs_2)){
              thisOptional_KO <- c(thisOptional_KO,replace_KOs_2[R])
              thisOptEC2   <- unique(gsubfn::strapplyc(ORTHOLOGS[grep(replace_KOs_2[R],ORTHOLOGS)],"(\\d{1}[.][[:punct:]||\\d]+)",simplify = T))

              if(is.character(thisOptEC2)) thisOptional_EC <- c(thisOptional_EC,thisOptEC2)
            } # end second KO replacement
          } # end subgroups
        } # end if grep "("
    } # end if grep "-"

  # Unique set
  thisOptional_EC <- unique(thisOptional_EC)
  thisOptional_KO <- unique(thisOptional_KO)

  output <- list(thisBlock = thisBlock, thisOptional_KOs = thisOptional_KO,thisOptional_ECs = thisOptional_EC)
  return(output)
}

#' Subgroup check - formatting KEGG module definition.
#'
#' This functions helps format the KEGG module DEFINITION by checking all the bracket-delimitted BLOCKS
#' to then remove flanking brackets (unnecessary).
#' @param DEFINITION - KEGG module definition
#' @seealso \link{parseKEGG_module},\link{misc_module_definition_check}
#'
#' @export

############################################################################################################################################

misc_module_subgroup_indexing <- function(DEFINITION){

  all_chars   <- strsplit(DEFINITION,split="+")[[1]] # split all characters (split by regular expression)
  Group_index <- numeric(length(all_chars)+1) # LEAVE ONE POSITION BEFORE
                          # zeros
  for(C in 1:length(all_chars)){
    if(all_chars[C]=="("){
      Group_index[C+1] <- Group_index[C] + 1
    } else if(C>1&&all_chars[C-1]==")"){
      Group_index[C+1] <- Group_index[C] - 1  # If the previous one was a ')', change group, unless its a '('
    } else {
      Group_index[C+1] <- Group_index[C]
    }
  }
  # Remove first zero
  Group_index <- Group_index[-1]

  return(Group_index)
}

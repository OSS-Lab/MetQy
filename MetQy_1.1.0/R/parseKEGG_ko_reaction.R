#' Map KEGG orthologs (KOs) to Reaction IDs
#'
#' Map KEGG orthologs (KOs) to Reaction IDs and format it into a matrix with binary indicator for mapping/relationship.
#' Generates 'ko_reaction_map' (.txt & .rda).
#'
#' @param KEGG_path - string pointing to the location of the KEGG database parent folder. The path to the required file is contained within the function.
#'
#' @param outDir    - string pointing to the output folder. Default ("output/").
#'
#' @param verbose   - logical. Should progress be printed to the screen? Default (\code{TRUE}).
#'
#' @return Data frame establishing the relationship between KO numbers and reactions (R number) (binary).
#'\preformatted{
#' > ko_reaction_map[1:3,1:3]
#'
#'               R00005 R00006 R00008
#'     K00001      0      0      0
#'     K00002      0      0      0
#'     K00003      0      0      0
#' }
#'
#' @examples
#' KEGG_path     <- "~/KEGG" # MODIFY TO KEGG PARENT FOLDER!
#' # The parent folder should contain the following (KEGG FTP structure):
#' #	brite/
#' #	genes/
#' #	ligand/
#' #	medicus/
#' #	module/
#' #	pathway/
#' #	README.kegg
#' #	RELEASE
#' #	xml/
#'
#' ko_reaction_map <- parseKEGG_ko_reaction(KEGG_path)
#'  # A .txt file (tab separated) is written to output/ (relative to current working directory)
#'
#' @seealso \link{parseKEGG_file.list}
#'
#' @export

############################################################################################################################################


parseKEGG_ko_reaction <- function(KEGG_path, outDir = "output", verbose = T){

  ####  MANAGE INPUT ----
  # CHECKS
  stopifnot(is.character(KEGG_path),length(KEGG_path)==1,dir.exists(KEGG_path))
  stopifnot((is.character(outDir)&&length(outDir)==1)||is.null(outDir))
  stopifnot(is.logical(verbose))

  # PATHS
  KEGG_path <- gsub("/$","",KEGG_path) # end file path WITHOUT ending "/"
  if(!is.null(outDir)) outDir    <- gsub("/$","",outDir)

  # OUTPUT FOLDER
  if(!is.null(outDir)) if(dir.exists(outDir)==F) dir.create(outDir)

  #### PROCESS FILE ----
  if(verbose) cat("\tko_reaction mapping...",fill = T)
  start <- Sys.time()

  ko_reaction_file_path <- paste(KEGG_path,"/genes/ko/ko_reaction.list",sep="")
  ko_reaction_map       <- parseKEGG_file.list(ko_reaction_file_path)

  #### WRITE MATRIX ----
  if(!is.null(outDir)){
    write.table(ko_reaction_map,file = paste(outDir,"/ko_reaction_map.txt",sep=""),sep="\t",quote = F)
    save(ko_reaction_map,file = paste(outDir,"/ko_reaction_map.rda",sep=""))
  }
  if(verbose) cat("\t Completed in",format(difftime(Sys.time(),start,units = "mins"),digits = 4),"\n",fill = T)

  return(ko_reaction_map)
}

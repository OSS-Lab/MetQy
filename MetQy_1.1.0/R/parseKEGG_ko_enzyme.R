#' Map KEGG orthology (KO) to Enzyme Commission (EC) numbers
#'
#' Map KEGG orthologs (K numbers) to EC numbers and format it into a matrix with binary indicator for mapping/relationship.
#' Generates ko_enzyme_map (.txt & .rda)
#'
#' @param KEGG_path - string pointing to the location of the KEGG database parent folder. The path to the required file is contained within the function.
#'
#' @param outDir    - string pointing to the output folder. Default ("output/").
#'
#' @param verbose   - logical. Should progress be printed to the screen? Default (\code{TRUE})
#'
#' @return Data frame establishing the relationship between K numbers and enzymes (binary).
#'\preformatted{
#' > ko_enzyme_map[1:3,1:3]
#'
#'           1.1.1.1 1.1.1.10 1.1.1.100
#' K00001       1        0         0
#' K00002       0        0         0
#' K00003       0        0         0
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
#' ko_enzyme_map <- parseKEGG_ko_enzyme(KEGG_path)
#'  # A .txt file (tab separated) is written to output/ (relative to current working directory)
#'
#' @seealso \link{parseKEGG_file.list}
#'
#' @export

############################################################################################################################################

parseKEGG_ko_enzyme <- function(KEGG_path, outDir = "output", verbose = T){

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
  if(verbose) cat("\tko_enzyme mapping...",fill = T)
  start <- Sys.time()

  ko_enzyme_file_path <- paste(KEGG_path,"/genes/ko/ko_enzyme.list",sep="")
  ko_enzyme_map       <- parseKEGG_file.list(ko_enzyme_file_path)

  #### WRITE MATRIX ----
  if(!is.null(outDir)){
    write.table(ko_enzyme_map,file = paste(outDir,"/ko_enzyme_map.txt",sep=""),sep="\t",quote = F)
    save(ko_enzyme_map,file = paste(outDir,"/ko_enzyme_map.rda",sep=""))
  }
  if(verbose) cat("\t Completed in",format(difftime(Sys.time(),start,units = "mins"),digits = 4),"\n",fill = T)

  ### RETURN REFERENCE TABLE ---
  return(ko_enzyme_map)
}

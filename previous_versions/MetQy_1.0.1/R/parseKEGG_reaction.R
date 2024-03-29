#' Parse KEGG reaction database
#'
#' Read and format the KEGG reaction database text file into a reference table.
#' 
#' @details 
#' The columns are automatically generated by the \code{parseKEGG_file} function into variables, 
#' which are further formatted specifically for the KEGG reaction database.
#' 
#' The text file used is "\code{KEGG_path}/ligand/reaction/reaction". 
#' 
#' It decompresses "\code{KEGG_path}/ligand/reaction.tar.gz" if needed.
#'
#' @param KEGG_path - string pointing to the location of the KEGG database parent folder
#'
#' @param outDir    - string pointing to the output folder. Default ("output/"). \code{NULL} overwrites existing files.
#'
#' @param verbose   - logical. Should progress be printed to the screen? Default (\code{TRUE})
#'
#' @param ...       - other arguments for \code{parseKEGG_file()}.
#'
#' @return Generates reaction_reference_table (.txt & .rda; saved to \code{'outDir'}) and returns
#' a data frame with as many rows as entries and the following columns (or variables):
#' \preformatted{
#' 	 (1) ID         - R number identifier (e.g. "R00001");
#' 	 (2) NAME       - reaction or enzyme name;
#' 	 (3) DEFINITION - reaction definition using compound's names;
#' 	 (4) EQUATION   - reaction definition using compound's IDs;
#' 	 (5) ENZYME     - Enzyme Commission (EC) number (e.g. "1.1.1.1");
#' 	 (6) COMMENT;  (7) RCLASS;
#' }
#' In all instances, multiple entries in a given column are separated by '[;]'.
#' @export
#'
#' @seealso \link{parseKEGG_file}

############################################################################################################################################

parseKEGG_reaction <- function(KEGG_path, outDir = "output", verbose = T,...){

  ####  MANAGE INPUT ----
  # CHECKS
  stopifnot(is.character(KEGG_path),length(KEGG_path)==1,dir.exists(KEGG_path))
  stopifnot((is.character(outDir)&&length(outDir)==1)||is.null(outDir))
  stopifnot(is.logical(verbose))

  # PATHS
  KEGG_path <- gsub("/$","",KEGG_path)
  if(!is.null(outDir)) outDir    <- gsub("/$","",outDir)

  # LOCAL FOLDER
  if(!is.null(outDir)) if(!dir.exists(outDir)) dir.create(outDir)

  #### READ IN FILE ----
  if(verbose) cat("\treaction processing...",fill = T)
  start <- Sys.time()

  #### UNTAR FILE ----
  if(!file.exists(paste(KEGG_path,"/ligand/reaction/reaction",sep=""))) {
    if(verbose) cat("\n\t\tDecompressing reaction file...",fill = T)
    untar(paste(KEGG_path,"/ligand/reaction.tar.gz",sep=""),files = "reaction/reaction",exdir = paste(KEGG_path,"/ligand/",sep="")) # extract only the desired file
  }

  ####  PROCESS FILE ----
  reaction_file           <- paste(KEGG_path,"/ligand/reaction/reaction",sep="")
  reaction_reference_table <- parseKEGG_file(reaction_file,verbose = verbose,...)

  ####  FORMAT TABLE                    ----
  reaction_reference_table$ENTRY      <- gsubfn::strapplyc(pattern = "R\\d{5}",reaction_reference_table$ENTRY,simplify = T)
  names(reaction_reference_table)[1]  <- "ID"

  reaction_reference_table$ENZYME   <- gsub(" ",";",reaction_reference_table$ENZYME)

  ####  WRITE TABLE  (.txt, .rda) ----
  reaction_reference_table_file <- paste(outDir,"/reaction_reference_table.txt",sep="")
  if(!is.null(outDir)){
    write.table(reaction_reference_table,file = reaction_reference_table_file,sep = "\t",row.names = F,quote = F)
    save(reaction_reference_table,file = gsub(".txt",".rda",reaction_reference_table_file))
  }
  if(verbose) cat("\t Completed in",format(difftime(Sys.time(),start,units = "mins"),digits = 4),"\n",fill = T)

  ### RETURN ---
  return(reaction_reference_table)
}

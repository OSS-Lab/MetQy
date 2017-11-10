#' Parse any '.list' KEGG file
#'
#' Reads the KEGG  database text files with '.list' extension (e.g. 'ko_enzyme.list', 'ko_reaction.list')
#' and formats it into a matrix with a binary indicator or relationships or mappings.
#'
#' @param FILE_PATH - string pointing to the location of '.list' file.
#'
#' @return MATRIX
#'
#' @export
#'
#' @seealso \link{parseKEGG_ko_enzyme}, \link{parseKEGG_ko_reaction}

################################ STEPS ########################################################
# 1.    MANAGE INPUT  - format paths and check/make outDir
# 2.    PROCESS FILE  - read table and convenrt into relationship matrix
# 3.    RETURN MATRIX
##############################################################################################

parseKEGG_file.list <- function(FILE_PATH){

  #### MANAGE INPUT ----
  # CHECK FILE EXISTS
  file_name               <- strsplit(FILE_PATH,split = "/")[[1]]
  file_name               <- file_name[length(file_name)]
  file_parent_path        <- gsub(paste("/",file_name,"$",sep=""),"",FILE_PATH)
  file_parent             <- strsplit(file_parent_path,split = "/")[[1]]
  file_parent             <- file_parent[length(file_parent)]
  file_parent_parent_path <- gsub(paste(file_parent,"$",sep=""),"",file_parent_path)

  # UNTAR FILE
  if(!file.exists(FILE_PATH)){
      untar(paste(file_parent_parent_path,file_parent,".tar.gz",sep=""),files = paste(file_parent,"/",file_name,sep=""),exdir = file_parent_parent_path)
  }

  #### PROCESS FILE ----
  TABLE    <- read.delim(FILE_PATH,header = F,stringsAsFactors = F)
  TABLE$V1 <- gsub(".*[:]","",TABLE$V1)
  TABLE$V2 <- gsub(".*[:]","",TABLE$V2)
  TABLE$V3 <- 1

  # CAST DATA.FRAME TO MATRIX ----
  MAP           <- reshape2::dcast(TABLE,formula = V1~V2,value.var = "V3",drop = F,fill = 0)
  rownames(MAP) <- MAP[,1]
  MAP           <- MAP[,2:ncol(MAP)]

  return(MAP)
  # View(MAP)
}

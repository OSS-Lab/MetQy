#' Parse any KEGG file without extension
#'
#' Read the KEGG database text file without extension (e.g. 'module', 'enzyme', 'genome') and format it into a
#' reference table.
#' Generates DATABASE_reference_table (data frame).
#'
#' @details
#' File must be decompressed before being processed.
#' The functions that use this function (listed below in \code{see also}) perform this step if necessary.
#'
#' @param FILE          - string pointing to the location of the file WITHOUT EXTENSION (uncompressed).
#'
#' @param split_pattern - string to use to identify start of new section/entry. Default ("ENTRY").
#'
#' @param pathway_trim  - logical. Should the file 'KEGG_path/pathway/pathway' be trimmed to only include 'map' entries? (i.e. exclude ko, ec and organism-specific pathways).
#'                          Default (\code{TRUE}).
#'
#' @param verbose       - logical. Should progress be printed to the screen? Default (\code{TRUE}).
#'
#' @param ...  - further arguments (currently unsupported)
#'
#' @return A data frame with the formatted data.See the file-specific functions.
#'
#' @seealso \link{parseKEGG_compound},\link{parseKEGG_enzyme},\link{parseKEGG_genome}, \link{parseKEGG_module},
#'
#' @seealso \link{parseKEGG_reaction}, \link{parseKEGG_execute_all}
#'
#' @examples
#' compound_file_path <- "~/KEGG/ligand/compound/compound" # MODIFY!
#' reference_table <- parseKEGG_file(compound_file_path)
#' # WITHOUT DATABASE SPECIFIC FORMATTTING!
#'
#' @export
#'

################################ STEPS ########################################################
# 1.    PROCESS FILE  - read table and convenrt into relationship matrix
# 1.1     Process an entry at a time
# 1.1.1     Skip obsolete enzymes or past module entries
# 1.1.2     Identify new section headings
# 1.1.3     Populate reference table
# 2.    RETURN REFERENCE TABLE
##############################################################################################


parseKEGG_file <- function(FILE_PATH,split_pattern = "ENTRY",pathway_trim=T,verbose = T,...){

  ####  MANAGE INPUT ----
  # CHECKS
  stopifnot(is.character(FILE_PATH),length(FILE_PATH)==1,file.exists(FILE_PATH))
  stopifnot(is.character(split_pattern),length(split_pattern)==1)
  stopifnot(is.logical(verbose))

  #### READ FILE ----
  FILE    <- scan(FILE_PATH,what = "character", sep="\n", quiet = T) # read in each line
  SPLIT   <- grep(split_pattern,FILE)
  LENGTH  <- length(SPLIT)

  if(pathway_trim) if(length(grep("/pathway/pathway$",FILE_PATH))>0){
    index <- grep("map",FILE[SPLIT])
    SPLIT[index[length(index)]+1]
    FILE <- FILE[1:(SPLIT[index[length(index)]+1]-2)]
    SPLIT   <- grep(split_pattern,FILE)
    LENGTH  <- length(SPLIT)
  }

  #### PROCESS FILE ----
  for(S in 1:LENGTH){
    if(verbose) if(S==1||S%%500==0||S==LENGTH) cat("\t",S,"\tof",LENGTH,"chunks", fill = T)
    if(S<LENGTH){
      CHUNK <- FILE[SPLIT[S]:(SPLIT[S+1]-2)]
    }else {
      CHUNK <- FILE[SPLIT[S]:length(FILE)]
    }

    ### SKIP ----
    ## IF READING OBSOLETE ENZYME
    if(length(grep("Obsolete",CHUNK[1]))>0) next

    # SKIP
    if(length(grep("[a-z]+[_]{1}[[:alpha:]]+\\d{5}",CHUNK[1]))>0) next

    ### START REFERENCE TABLE ----
    if(S==1) {
      headers         <- unique(gsubfn::strapplyc(CHUNK[grep("^[A-Z]+",CHUNK)],"(^[A-Z]+[[:punct:]]*[A-Z]*)",simplify = T))
      reference_table <- matrix("",LENGTH,length(headers))
      reference_table <- matrix("",LENGTH,length(headers))
      colnames(reference_table) <- headers
    }

    ### IDENTIFY SECTIONS ----
    section_index         <- grep("^[A-Z]",CHUNK)
    sections_present      <- gsubfn::strapplyc(CHUNK[section_index],"(^[A-Z]+[[:punct:]]*[A-Z]*)",simplify = T)

    ##    ID NEW SECTIONS ----
    more_sections_index   <- which(is.na(match(sections_present,headers)))
    # ADD MISSING SECTIONS
    if(length(more_sections_index)>0){
      missing_sections    <- unique(sections_present[more_sections_index])
      add_cols            <- matrix("",nrow(reference_table),length(missing_sections))
      colnames(add_cols)  <- missing_sections
      reference_table     <- cbind(reference_table,add_cols)
      headers             <- c(headers,missing_sections)
      rm(add_cols)
    }

    ### POPULATE REFERENCE TABLE ----
    for(H in 1:length(headers)){
      where_index   <- which(sections_present==headers[H])
      # HEADER SECTION NOT PRESENT IN CHUNK
      # if(is.na(match(headers[H],sections_present))){
      if(length(where_index)==0){
        # cat(S,H,"no header match", sep="\t",fill = T)
        next
      }

      # RETRIEVE INFO
      if(length(where_index)>1){
        section_text <- NULL
        for(ii in where_index){
          if(ii<length(section_index)){
            to_index      <- section_index[ii+1]-1
          } else{
            to_index      <- length(CHUNK)
          }
          this_section_text <- CHUNK[section_index[ii]:to_index]
          section_text <- c(section_text, paste(gsub("[.]$","",gsub("\\s+"," ",gsub("^\\s*","",gsub(headers[H],"",this_section_text)))),collapse = ","))
        }
        section_text <- paste(section_text,collapse = ";")
      }else{
        if(where_index<length(section_index)){
          to_index      <- section_index[where_index+1]-1
        } else{
          to_index      <- length(CHUNK)
        }
        section_text  <- CHUNK[section_index[where_index]:to_index]
        section_text  <- paste(gsub("\\s+"," ",gsub("^\\s*","",gsub("\\s*$","",gsub(headers[H],"",section_text)))),collapse = ";") # trim white spaces and remove SECTION HEADER
      }
      # WHERE TO STORE INFO
      this_col <- which(colnames(reference_table)==headers[H])
      # S is row?
      reference_table[S,this_col] <- gsub(";;",";",section_text)

      rm(section_text)
    }
  }# end CHUNK loop
  ### ELIMINATE EMPTY ROWS ----
  remove_index    <- which(apply(reference_table,1,paste,collapse = "")=="")
  if(length(remove_index)>0) reference_table <- reference_table[-remove_index,]
  # AS DATA FRAME
  reference_table  <- as.data.frame(reference_table,stringsAsFactors = F)

  return(reference_table)
  # View(reference_table)
}

#' Execute all parseKEGG parent functions to format KEGG databases into data frames
#'
#' Execute all parseKEGG parent functions to format specific KEGG databases into data frames.
#'
#' @param KEGG_path - string pointing to the location of the KEGG database parent folder.
#'
#' @param ...       - further arguments, such as \code{outDir}, for \link{parseKEGG_file},
#'                     \link{parseKEGG_file.list} and database-specific functions (below).
#'
#' @seealso \link{parseKEGG_compound}, \link{parseKEGG_enzyme},
#'          \link{parseKEGG_genome}, \link{parseKEGG_module},
#'
#' @seealso \link{parseKEGG_ko}, \link{parseKEGG_reaction},
#'          \link{parseKEGG_ko_enzyme}, \link{parseKEGG_ko_reaction}
#'
#' @examples
#' KEGG_path <- "~/KEGG" # MODIFY!
#' parseKEGG_parseKEGG_execute_all(KEGG_path)
#' # multiple reference_table objects in workspace and .txt files written to
#' #    output/ (relative to current working directory)
#'
#' @export

############################################################################################################################################

parseKEGG_execute_all <- function(KEGG_path, ...){

  ####  MANAGE INPUT ----
  stopifnot(is.character(KEGG_path),length(KEGG_path)==1,dir.exists(KEGG_path))

  ### EXECUTE FILE-SPECIFIC FUNCTIONS ----
  # COMPOUND
  try(compound_reference_table    <- parseKEGG_compound(KEGG_path, ...))

  # ENZYME
  try(enzyme_reference_table      <- parseKEGG_enzyme(KEGG_path, ...))

  # GENOME
  try(genome_reference_table      <- parseKEGG_genome(KEGG_path, ...))

  # MODULE
  try(module_reference_table      <- parseKEGG_module(KEGG_path, ...))

  # REACTION
  try(reaction_reference_table    <- parseKEGG_reaction(KEGG_path, ...))

  # KO
  try(ko_reference_table          <- parseKEGG_ko(KEGG_path, ...))

  # KO - EC MAPPING
  try(ko_enzyme_map               <- parseKEGG_ko_enzyme(KEGG_path, ...))

  # KO - REACTION MAPPING
  try(ko_reaction_map             <- parseKEGG_ko_reaction(KEGG_path,...))

}

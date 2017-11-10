# test_parse.R

library(MetQy)
context("Parse output")

test_that("correct output is produced", {
  KEGG_path <- "data/test_data"
  # COMPOUND
  expect_true(is.data.frame(parseKEGG_compound(KEGG_path,outDir = NULL,verbose = F)))

  # ENZYME
  expect_true(is.data.frame(parseKEGG_enzyme(KEGG_path,outDir = NULL,verbose = F)))


  # GENOME
  expect_true(is.data.frame(parseKEGG_genome(KEGG_path,outDir = NULL,verbose = F,addECs = F,addKOs = F)))


  # MODULE
  expect_true(is.data.frame(parseKEGG_module(KEGG_path,outDir = NULL,verbose = F)))


  # REACTION
  expect_true(is.data.frame(parseKEGG_reaction(KEGG_path,outDir = NULL,verbose = F)))


  # KO
  expect_true(is.data.frame(parseKEGG_ko(KEGG_path,outDir = NULL,verbose = F)))


  # KO - EC MAPPING
  expect_true(is.data.frame(parseKEGG_ko_enzyme(KEGG_path,outDir = NULL,verbose = F)))


  # KO - REACTION MAPPING
  expect_true(is.data.frame(parseKEGG_ko_reaction(KEGG_path,outDir = NULL,verbose = F)))


})



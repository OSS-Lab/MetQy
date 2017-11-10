# test_queries.R

library(MetQy)
library(testthat)
context("Queries output")

test_that("correct output is produced", {
  expect_true(is.data.frame(query_genes_to_genome("4.1.3.1")))

  expect_true(is.data.frame(query_genes_to_modules("4.1.3.1")))

  expect_true(is.list(query_genomes_to_modules("T00001",MODULE_ID = paste("M0000",1:5,sep=""))))

  expect_true(is.data.frame(query_missingGenes_from_module("T00001","M00001",PRINT_TO_SCREEN = F)))

  expect_true(is.numeric(query_modules_to_genomes("M00001")))
})
test_that("Examples run without error", {
  expect_error(is.data.frame(query_genes_to_genome("4.1.3.1")),NA)
  # query_genes_to_genome ----
  expect_error(query_genes_to_genome("K00844"),NA)
  expect_error(query_genes_to_genome("1.1.1.1"),NA)
  expect_error(query_genes_to_genome(genes = paste("K0000",1:3,sep="")),NA)

  # query_genes_to_modules ----
  expect_error(query_genes_to_modules("K00844"),NA)
  expect_error(query_genes_to_modules("1.1.1.1"),NA)

  # query_genomes_to_modules ----
    expect_error(query_genomes_to_modules(paste("T0000",1:5,sep=""),MODULE_ID = paste("M0000",1:5,sep=""),META_OUT = T, ADD_OUT = T),NA)
    expect_error(query_genomes_to_modules(c("escherichia coli","heliobacter"),MODULE_ID = paste("M0000",1:5,sep=""),META_OUT = T, ADD_OUT = T),NA)
    expect_error(query_genomes_to_modules(data_example_multi_ECs_KOs,GENOME_ID_COL = "ID",GENES_COL = "KOs",MODULE_ID = paste("M0000",1:5,sep=""),META_OUT = T,ADD_OUT = T),NA)
    expect_error(query_genomes_to_modules(data_example_multi_ECs_KOs,GENOME_ID_COL = "ID",GENES_COL = "ECs",MODULE_ID = paste("M0000",1:5,sep=""),META_OUT = T,ADD_OUT = T),NA)

  # query_missingGenes_from_module ----
  expect_error(query_missingGenes_from_module(data_example_KOnumbers_vector,"M00001",PRINT_TO_SCREEN = FALSE),NA)

  # query_modules_to_genomes ----
  expect_error(query_modules_to_genomes("M00001"),NA)
  expect_error(query_modules_to_genomes(MODULE_ID = c("M00001","M00002"),threshold = 0.75),NA)
})

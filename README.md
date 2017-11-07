# MetQy - metabolic query 

MetQy is a R package to ease interfacing with the Kyoto Encyclopedia of Genes and Genomes (KEGG) database to query metabolic functions of genes and genomes. Given a list of genes or an organism, it will return the set of functional modules possessed by an organism and their completeness.

## Getting started

MetQy requires that you have the follwing pacakes installed:

    dplyr, ggplot2, gsubfn, reshape

To install the package, clone the repository, launch R and run:

    install.packages('path_to_directory/MetQy_1.0.1.tar.gz',repos=NULL)
    library(MetQy)
    
## Example of use

The function genomes_to_modules takes as input a genome or a dataframe containing genomes, and returns the completeness (proportion of module 'blocks' contained in the genome(s)) of each KEGG module for each genome. These genomes can be expressed as either:

* KEGG species names
* a list of KEGG ortholog IDs e.g. K00001
* a list of enzyme commison (EC) numbers e.g. 1.1.1.1

Assuming you have a datafile with a list of organisms and KEGG orthologs of the following form:

    organism_ID,genes
    ORG1,K00001;K00002;K00005;K19234 ...
    ORG2,K00002;K00003;K00004;K10023 ...
    ...
    
Run the genomes to modules query as follows:

    genome_table = read.table(filename, stringsAsFactors=FALSE, header=TRUE)
    query_output = query_genomes_to_modules(genome_table, splitBy='[;]')
    modules_table = query_output$MATRIX
    
The output of this is a data frame where the rows are organism IDs, the columns are KEGG modules, and the entries are the completeness (proportion of complete blocks) in that module.  

## Other functionality

MetQy can also perform a variety of other queries:

  * _query_gene_to_modules_:          find metabolic modules a given gene is involved in
  * _query_module_to_genomes_:        find KEGG organisms expressing a given module
  * _query_genes_to_genome_:          find organisms that express a given set of genes
  * _query_missingGenes_from_module_: given a genome and a module, return the set of genes missing from that module, if any.
  
For details on the use of these functions, consult their help functions.

# MetQy - metabolic query 

MetQy is a R package to ease interfacing with the Kyoto Encyclopedia of Genes and Genomes (KEGG) database to query metabolic functions of genes and genomes. Given a list of genes or an organism, it will return the set of functional modules possessed by an organism and their completeness.

## Using MetQy
### CITE US
Please cite Martinez-Vernon et al. MetQy – an R package to query metabolicfunctions of genes and genomes. _In preparation_.

### Commerical use
MetQy is a free software for academic purposes. If interested in commerical use, please read the LICENCE and contact [Warwick Ventures](mailto:ventures@warwick.ac.uk)

### Getting started
#### Installation
There are two ways of installing the package:

1. Download 'MetQy_1.0.1.tar.gz' or 'MetQy_1.0.1.tgz' and run the follwing command withing the R environment:
    ```
    install.packages('path_to_directory/MetQy_1.0.1.tgz',repos=NULL)  # QUICKER
    library(MetQy)
    ```
    _OR_
    ```
    install.packages('path_to_directory/MetQy_1.0.1.tar.gz',repos=NULL)
    library(MetQy)
    ```
2. Run the following command within the R enviornment. Requires package `devtools`
    ```
    devtools::install_github('OSS-Lab/MetQy',subdir = 'MetQy_1.0.1')
    ```

Note that MetQy requires the follwing pacakes, but these should be installed along with tha package (if necessary):

    dplyr, ggplot2, gsubfn, reshape

    
#### Example of use

The function `query_genomes_to_modules` takes as input a genome of set of genes making up a genome and returns the module completeness fraction (_mcf_) (proportion of module 'blocks' contained in the genome(s)) of specified KEGG module for each genome (see https://github.com/OSS-Lab/MetQy/wiki/KEGG-databases-description). These genomes can be expressed as either:

* KEGG genome identifier (T number or 3-4 letter code, e.g. 'T0001' or 'eco')
* Organism name (e.g. 'Escherichia coli')
* a list of KEGG ortholog IDs (e.g. 'K00001')
* a list of Enzyme Commison (EC) numbers (e.g. '1.1.1.1')

The last two options allow the user to provide their own [meta]genome as determine by the genes or enzymes provided.

Assuming you have a data file with a list of organisms and KEGG orthologs of the following form:

    organism_ID,genes
    ORG1,K00001;K00002;K00005;K19234 ...
    ORG2,K00002;K00003;K00004;K10023 ...
    ...
    
    genome_table = read.table(filename, stringsAsFactors=FALSE, header=TRUE)
    
Alternatively, the package has inbuilt data examples to provide users with examples:
    
    # Format compatible with query functions:
    data(data_example_multi_ECs_KOs)
    names(data_example_multi_ECs_KOs)
    # [1] "ID"       "ORG_ID"   "ORGANISM" "KOs"      "ECs"     

Run the genomes to modules query as follows:
    
    query_output = query_genomes_to_modules(genome_table,MODULE_ID = "M0001",splitBy='[;]',GENOME_ID_COL = 1,GENES_COL = 2)
    # OR
    query_output = query_genomes_to_modules(data_example_multi_ECs_KOs,MODULE_ID = "M0001",splitBy='[;]',GENOME_ID_COL = 1,GENES_COL = 4) # Use KOs
    # OR
    query_output = query_genomes_to_modules(data_example_multi_ECs_KOs,MODULE_ID = "M0001",splitBy='[;]',GENOME_ID_COL = 1,GENES_COL = 5) # Use ECs

    modules_table = query_output$MATRIX
    
The output of this is a data frame where the rows are organism IDs, the columns are KEGG modules, and the entries are the module completeness fraction (proportion of complete blocks) for each module. In the examples avobe, there would only be one column as there was only one modules specified.

#### Other functionality

MetQy can also perform a variety of other queries:

 * `query_genes_to_genome` - find organisms that express a given set of genes
 * `query_genes_to_modules` - find metabolic modules a given gene is involved in
 * `query_modules_to_genomes` - find KEGG genomes expressing a given module
 * `query_missingGenes_from_module` - given a genome and a module, return the set of genes missing from that module, if any.
  
For details on the use of these functions, consult their help functions either through R or in the pdf documentation document.

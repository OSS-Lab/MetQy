# MetQy - metabolic query 

MetQy is a R package to ease interfacing with the Kyoto Encyclopedia of Genes and Genomes (KEGG) database to query metabolic functions of genes and genomes. As discussed in detail in the associated paper (see above), MetQy can be used to convert key parts of the data contained in the KEGG database into R data structures, and perform analyses on this data. Below we provide a simple example of how MetQy can be used to analyze a set of given genomes for the presence of a set of given metabolic capabilities (i.e. modules).

## Using MetQy
#### CITE US
Please cite:
```
    Andrea Martinez-Vernon, Fred Farrell, Orkun Soyer. MetQy - an R package to query metabolic functions of genes and genomes. 
        Bioinformatics, Volume 34, Issue 23, 01 December 2018, Pages 4134–4137, https://doi.org/10.1093/bioinformatics/bty447   
```
https://doi.org/10.1093/bioinformatics/bty447  

#### Commercial use
MetQy is a free software for academic purposes. If interested in commercial  use, please read the [LICENCE](https://github.com/OSS-Lab/MetQy/blob/master/LICENCE) and contact [Warwick Ventures](mailto:ventures@warwick.ac.uk)

## Installation
There are three ways of installing the package:

1. Run the following command within the R environment. Requires package `devtools`
    ```
    # > install.packages("devtools") # Uncomment if not previously installed 
    
    # Install MetQy - required dependent packages are installed automatically
    > devtools::install_github('OSS-Lab/MetQy',subdir = 'MetQy_1.1.0',dependencies = TRUE) # 
    ```
    
2. Download 'MetQy_1.1.0.tar.gz' and run the following commands within the R environment:
    ```
    # Manually install dependent packages - remove any that is already installed. Installing will replace local library.
    > install.packages(c("dplyr","ggplot2","gsubfn","reshape2","xtable"))
    
    # Install MetQy
    > install.packages('<path_to_directory>/MetQy_1.1.0.tar.gz',repos=NULL)
    > library(MetQy)
    ```

3. Download 'MetQy_1.1.0.tgz' and run the following commands within the R environment:
    ```
    ## Note that this option might only be compatible with OS X. 
    
    # Manually install dependent packages - remove any that is already installed. Installing will replace local library.
    > install.packages(c("dplyr","ggplot2","gsubfn","reshape2","xtable"))
    
    # Install MetQy
    > install.packages('<path_to_directory>/MetQy_1.1.0.tgz',repos=NULL)  # QUICKER
    > library(MetQy)
    ```

*** 
# More Info

Detailed information on MetQy functions, along with usage examples and a sample analyses can be found on associated [Wiki page](https://github.com/OSS-Lab/MetQy/wiki).

## Wiki menu
* [About MetQy and installation][home]
* [KEGG databases descriptions][kegg]
* [MetQy functions - Query][fun_query]
* [MetQy functions - Parsing][fun_parse]
* [MetQy functions - Visualization][fun_vis]
* [MetQy functions - Analysis][fun_analysis]
* [MetQy example - biological analysis][bio_eg]

[home]: https://github.com/OSS-Lab/MetQy/wiki/About-MetQy-and-installation
[kegg]: https://github.com/OSS-Lab/MetQy/wiki/KEGG-databases-description
[fun_query]: https://github.com/OSS-Lab/MetQy/wiki/MetQy-functions-and-usage-examples-%E2%80%93-Query-functions
[fun_parse]: https://github.com/OSS-Lab/MetQy/wiki/MetQy-functions-and-usage-examples-%E2%80%93-Parsing-functions
[fun_vis]: https://github.com/OSS-Lab/MetQy/wiki/MetQy-functions-and-usage-examples-%E2%80%93-Visualization-functions
[fun_analysis]: https://github.com/OSS-Lab/MetQy/wiki/MetQy-functions-and-usage-examples-%E2%80%93-Analysis-functions
[bio_eg]: https://github.com/OSS-Lab/MetQy/wiki/Example-use-of-MetQy-for%E2%80%93biological-analysis

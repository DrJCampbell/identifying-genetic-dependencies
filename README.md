# identifying-genetic-dependencies

Provides a simple example of running tests of associations between gene mutation or copy number status and genetic dependency on kinome genes.

To get started, look in the 'R-scripts' directory and open the R script named 'identifying_CGDs_RB1_osteosarcoma.R'. Running the code in this file should generate the statistical associations table and image file found in the 'results' directory.

The other script in the 'R-scripts' directory contains a library of functions that abstract the process of combining mutation and siRNA data sets, running association tests on these data and generating an image to summarise the results. Reading through the three functions in this file will provide details on occurs during the entire process.

The 'mutation-data' and 'siRNA-data' directories contain tables with mutation data and Z-scores from an siRNA screen of the kinome. The kinome screen data were extracted from a published data set ([http://www.ncbi.nlm.nih.gov/pubmed/26947069 Campbell, Ryan, et al., 2016]) and provide the siRNA Z-scores for 18 osteosarcoma cell lines screened using a library of siRNAs targeting 714 kinases and kinases-related genes (the kinome). The mutation data sets were generated using the [http://cancer.sanger.ac.uk/cell_lines Cosmic Cell Lines Project (v76)] mutation and copy number data sets processed using [https://github.com/GeneFunctionTeam/cell_line_functional_annotation publically available scripts].




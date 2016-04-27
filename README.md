chromstaR
=========

This package implements functions for combinatorial and differential analysis of ChIP-seq data. It includes uni- and multivariate peak-calling, export to genome browser viewable files, and functions for enrichment analyses.

Installation
------------

### Stable release version from Bioconductor
To install the *current stable* version from Bioconductor, please visit http://bioconductor.org/packages/chromstaR/ and follow the provided instructions.

### Development version from Github
To install the *development* version from Github, follow the steps given below. The installation has only been tested on Ubuntu so far, if you need to install on Windows or Mac additional steps might be necessary (e.g. installation of Rtools from https://cran.r-project.org/bin/windows/Rtools/)

1. Install a recent version of R (3.3.0) from https://www.r-project.org/
2. Optional: For ease of use, install Rstudio from https://www.rstudio.com/
3. Open R and install all dependencies. Please ensure that you have writing permissions to install packages. Execute the following lines one by one:

   install.packages("devtools")  
	 source("http://bioconductor.org/biocLite.R")  
	 biocLite(c("GenomicRanges","GenomicAlignments"))  
	 library(devtools)  
	 install_github("ataudt/chromstaRData")  
	 install_github("ataudt/chromstaR")  
	 # Or alternatively if the above line doesn't work:  
	 install_git("git://github.com/ataudt/chromstaRData.git", branch = "master")  
	 install_git("git://github.com/ataudt/chromstaR.git", branch = "master")


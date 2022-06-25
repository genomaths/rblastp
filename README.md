<!-- README.md is generated from README.Rmd. Please edit that file -->

# rblastp: An R interface to run command line blastp

Robersy Sanchez  
Department of Biology. Eberly College of Science.  
Pennsylvania State University, University Park, PA 16802  
<rus547@psu.edu>  
[ORCID:
orcid.org/0000-0002-5246-1453](https://orcid.org/0000-0002-5246-1453)

# Overview

This an R package to run command line
[blastp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins) from R,
permitting a soft downstream continuation of a pipeline analysis in R.

# Get Started

## Command-line Blast

Blast must be installed in your machine

    sudo apt install ncbi-blast+

# Installation

Be sure that both the R and bioconductor packages are up to date. To get
the latest version of Bioconductor by starting R and entering the
commands:

    if (!requireNamespace("BiocManager")) install.packages("BiocManager")
    BiocManager::install()

Install R dependencies:

    install.packages(seqinr)

    BiocManager::install("BiocParallel", "Biostrings", "S4Vectors"))

You can install ‘rblastp’ from GitHub:

    install.packages("devtools")
    # Master stable version (0.3.2.3)
    devtools::install_git("https://github.com/genomaths/rblastp.git")

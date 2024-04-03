# E.PathDash

This repository contains the R shiny application code, data collection and analysis scripts, and cleaned data files used to create the E.PathDash web application. 

### Introduction

E.PathDash is a Shiny application that facilitates re-analysis of gene expression data from pathogens clinically relevant to Cystic Fibrosis. The application runs pathway activation analysis of KEGG pathways and gene ontology (GO) terms for a set of RNA-seq datasets compiled from the [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/). The datasets cover the pathogens *Pseudomonas aeruginosa*, *Staphylococcus aureus*, *Streptococcus sanguinis*, and *Bacteroides thetaiotaomicron*. Users can rapidly toggle between experimental comparisons and different studies of the same phenomenon, enabling them to judge the extent to which observed responses are reproducible.

In addition to exploring analysis results within the application interface, users can download high quality images that communicate pathway activation, differential gene expression tables, and raw count data. The application is freely available at scangeo.dartmouth.edu/EPathDash/. 

### Related applications

This application was inspired by two other applications developed by our group, [ESKAPE Act PLUS](https://github.com/DartCF/ESKAPE_Act) and [CF-Seq](https://github.com/DartCF/cf-seq). Pathway activation analysis was implemented using logic devleoped for ESKAPE Act PLUS and the original set of RNA-seq datasets was extracted from the database compiled for CF-Seq. 

### Directory Contents

**/App_Files**

This folder contains all application files except for the large R data file that is used to load the pathway activation analysis and RNA-seq data objects into the application. 

- **app.R** implements the UI and server logic for the application

- **www** contains a CSS file for additional styling elements

- **EPathDashUserGuide.pdf** user guide that users can download from the application

- **Supported_Species.csv** table of supported bacterial species users can download from the application

**/GEO_Datasets**

This folder contains the files used to construct the data objects for E.PathDash. The folder has a subdirectory for each species. Within each species subdirectory are *Count_Tables*, *Design_Matrices*, and *Metadata* folders for each of the three table types (all as CSV files). All count, design, and metadata files are named with their GEO accession identifier. 

**/Query_GEO**

This folder contains the bash scripts used to query the NCBI's GEO Datasets database for RNA-seq datasets from the bacterial species and strains of interest. These scripts were used to collect the raw data for E.PathDash. 

- **id_extract.sh** gets the UIDs for datasets returned by a particular query string

- **ftp_extract.sh** gets the FTP locations for supplementary files for a list of UIDs (which define GEO datasets)

- **merge.sh** creates a table of UIDs and their corresponding GEO accession numbers 

- **get_metadata.sh** gets the metadata for given GEO datasets

- **suppl_data_download.sh** downloads supplmentary files stored at given FTP locations

**data_collection.R**

This script was used to pull relevant datasets from CF-Seq, combine these with the new datasets identified using the scripts in the **Query_GEO** folder, and translate all gene identifiers to a consistent encoding schema (UniprotKB) in preparation for pathway activation analysis. 

**pathway_activation_analysis.R**

This script was used to run pathway activation analysis on the componedium of RNA-seq datasets cleaned and compiled by the data_collection script. The final data obejct produced by this script serves as the database that powers E.PathDash.  


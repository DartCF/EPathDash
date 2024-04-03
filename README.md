# E.PathDash

### Introduction

E.PathDash is an Shiny application that facilitates re-analysis of gene expression data from pathogens clinically relevant to cystic fibrosis. The application runs pathway activation analysis of KEGG pathways and gene ontology (GO) terms for a set of RNA-seq datasets compiled from the [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/). The datasets cover the pathogens Pseudomonas aeruginosa, Staphylococcus aureus, Streptococcus sanguinis, and Bacteroides thetaiotaomicron. Users can rapidly toggle between experimental comparisons and different studies of the same phenomenon, enabling them to judge the extent to which observed responses are reproducible.

In addition to exploring analysis results within the application interface, users can download high quality images that communicate pathway activation, differential gene expression tables, and raw count data. The application is freely available at scangeo.dartmouth.edu/EPathDash/. 

### Directory Contents

**/App Files**

This folder contains all application files except for the large R data file that is used to load the pathway activation analysis and RNA-seq data objects into the application. 

- **app.R** implements the UI and server logic for the application

- **www** contains a CSS file for additional styling elements

- **EPathDashUserGuide.pdf** user guide that users can download from the application

- **Supported_Species.csv** table of supported bacterial species users can download from the application


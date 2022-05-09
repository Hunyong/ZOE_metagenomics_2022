# ZOE_metagenomics_2022  
Pathobiont-mediated spatial structuring enhances biofilm virulence in childhood oral disease

# Instruction
`scripts` folder contains the R code for implementing the linear model fitting (C1.R), ...
The `F-` files are the helper functions.

One the data are prepared in the following format, the C11, C21, C22, C31, C32, C41, all ending with `.R`, can be run to obtain the desired outputs and figures. The leading numbers corresponds to the flow of the analyses.  


# Data
Data-processed includes the following collated datasets
* `data.bracken2.full.DRNA.ZOE2.rds`  
* `data.bracken2.full.DRNA.ZOE1.rds`  
* `data.humann3.path.joint.DRNA.ZOE2.rds`   
* `data.humann3.path.marginal.DRNA.ZOE2.rds`  
* `data.geneRPK.full.DRNA.ZOE2.RNASEQ.rds`  

Each file is an R list object that contains `otu`, `meta`, `taxa`, where, in case of `data.bracken2.full.DRNA.ZOE2.rds`,
`otu` is a three-dimensional array of bracken-based species abundance for each species (p = 6412 rows), each subjects (n = 302 columns), and metagenomics/metatranscriptomics (3rd dimension = 2) of ZOE 2.0 main study. 
`meta` is an $n\times q$ data.frame of meta data that include id, phenotypes, age, race, sequencing dates, etc.
`taxa` is a $p\times 1$ data.frame that includes taxonomic information.  

`data.bracken2.full.DRNA.ZOE1.rds` is a same structured data for ZOE 2.0 pilot data with $n=118, p = 5703$.

`data.humann3.path.joint.DRNA.ZOE2.rds` and `data.humann3.path.marginal.DRNA.ZOE2.rds` are same structured data for pathway-species combinations and marginal pathways processed according to HUMMANn 3 for ZOE 2.0 data. Each row of the pathway-species data is in the form of `(gene names) (species names)`.  

`data.geneRPK.full.DRNA.ZOE2.RNASEQ.rds` is a same structured data for RNA sequencing data for the TOP 4 species for ZOE 2.0 data.

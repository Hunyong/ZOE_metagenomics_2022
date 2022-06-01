# ZOE_metagenomics_2022  
This is the analysis code for the paper "Pathobiont-mediated spatial structuring enhances 
biofilm virulence in childhood oral disease" by Cho, Ren, Divaris, et al. (2022). 
This code was prepared in the R computing language (R 4.1.2).  


# Instruction
`scripts` folder contains the R code. See ***scripts*** below.  

Once the data are prepared in the following format using `C01--C03.R` files, the C11, C21, C22, C31, C32, C41, all ending with `.R`, can be run to obtain the desired outputs and figures. The leading numbers corresponds to the flow of the analyses. I.e., the latter numbers may depend on the earlier numbers, but not the other way around.    


# Data
The `Data/` folder contains elements in `CDR_ZOEmicro_20220519.zip` in Carolina Digital Repository (https://cdr.lib.unc.edu/concern/data_sets/5d86p890x).  
* `ZOE2_targeted_RNAseq/`  
* `ZOE2_pilot_Bracken/`  
* `ZOE2_main_Bracken/`
* `ZOE2_h3_RNA_path/`
* `data.outcome.ZOE2_pilot_20200527.xlsx`  
* `data.outcome.ZOE2_main_20200527.xlsx`  

# Data-processed  
The `Data-processed` folder is obtained by implementing the `C01.R` code.  
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

# scripts   
The scripts folder contains the R code files.  
* `C11.R`: demographics (Extended Data Table 5).  
* `C21--C23.R`: discovery and validation of the significant species (Extended Data Table and Figure 1, Figure 2).  
* `C31--C32.R`: associated pathways identification (Extended Table and Figure 2).  
* `C41.R`: associated gene identification (Extended Tables 3 and 4).  
* `C51.R`: correlation analysis of the significant species  (Figure 3)  
* The `F_.R` files are the helper functions.  

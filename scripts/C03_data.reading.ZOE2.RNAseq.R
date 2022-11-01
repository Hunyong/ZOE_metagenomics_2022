### C03.data.reading.ZOE2.RNAseq.R
### Data for the targeted species gene expression analysis
### Created by Hunyong Cho
### input: Data/ZOE2_targeted_RNAseq, Data/data.outcome.ZOE2_xxx.xlsx
### output: Data-processed/data.geneRPK.full.RNA.ZOE2.RNASEQ.rds

### 0. library
    source("scripts/F_file.R")
    source("scripts/F_jointAnalysis.R")
    source("scripts/F_generic.R")
    library(dplyr); library(magrittr); library(tidyr)
    library(readxl)
    if (!dir.exists("Data-processed")) dir.create("Data-processed")

### 1. Data directory
    data.dir <- ("Data/ZOE2_targeted_RNAseq")
    
    rnaseq <- list()
    dir.all <- list.dirs(path = data.dir, recursive=FALSE)  # unnecessary to remove "./" ?
    species.nm <- 
      gsub(data.dir, "", dir.all) %>% 
      gsub("\\/", "", .) %>% 
      gsub("\\_.*", "", .)
    
    for (i.tmp in 1:length(dir.all)) {
      path.tmp <- list.dirs(dir.all[i.tmp], recursive=FALSE)
      rnaseq[[i.tmp]] <- data.frame(path = paste0(path.tmp, "/quant.sf"), stringsAsFactors=FALSE)
      rnaseq[[i.tmp]]$no <- gsub(".*/(.*)\\.SALMON/.*", "\\1", rnaseq[[i.tmp]]$path)  #extracting the current folder names
    }
    names(rnaseq) <- species.nm
    rm(dir.all, data.dir, i.tmp, path.tmp)
    print(rnaseq)

### 2. gene RNA
    
    # S. Mutans
    gene.full.RPK.Mutans <- jAnalyze(direction="long", study = "Mutans", id = "all", type="quant", 
                                     info.file=rnaseq, includeReads = TRUE)
    gene.full.RPK.Mutans[["bacteria"]] = "g__Streptococcus.s__Streptococcus_mutans"
    
    # P. Salivae
    gene.full.RPK.Salivae <- jAnalyze(direction="long", study = "Salivae", id = "all", type="quant", 
                                     info.file=rnaseq, includeReads = TRUE)
    gene.full.RPK.Salivae[["bacteria"]] = "g__Prevotella.s__Prevotella_salivae"
    
    # S. Sputigena
    gene.full.RPK.Sputigena <- jAnalyze(direction="long", study = "Sputigena", id = "all", type="quant", 
                                      info.file=rnaseq, includeReads = TRUE)
    gene.full.RPK.Sputigena[["bacteria"]] = "g__Selenomonas.s__Selenomonas_sputigena"
    
    # L. Wadei
    gene.full.RPK.Wadei <- jAnalyze(direction="long", study = "Wadei", id = "all", type="quant", 
                               info.file=rnaseq, includeReads = TRUE)
    gene.full.RPK.Wadei[["bacteria"]] = "g__Leptotrichia.s__Leptotrichia_wadei"
    
    gene.full.RPK.RNA <- rbind(gene.full.RPK.Mutans, gene.full.RPK.Salivae, gene.full.RPK.Sputigena, gene.full.RPK.Wadei)
    rm(gene.full.RPK.Mutans, gene.full.RPK.Salivae, gene.full.RPK.Sputigena, gene.full.RPK.Wadei); gc()
    # gene.full.RPK.RNA <- gene.full.RPK.RNA[gene.full.RPK.RNA$id != "Blank", ]
    
    gene.full.RPK.RNA$TPM %<>% as.character %>% as.numeric()
    gene.full.RPK.RNA$reads %<>% as.character %>% as.numeric()
    dim(gene.full.RPK.RNA) #  2,703,591 x 7 (gene, bacteria, TPM, reads, eff.length, length, id)
    

### 3. outcome
    outcome <- read_xlsx("Data/data.outcome.ZOE2_main_20200626.xlsx")
    outcome =
      mutate(outcome, id = gsub("IDMG", "subj", outcome$id_MTG)) %>%  # defining a common id across MTG and MTX.
      dplyr::select(id, everything())
    
### 4. Wide form and array of DNA and RNA together.
    
# geneRPK.RNA.wide.fn <- "../Data-processed/data.geneRPK.full.RNA.wide.ZOE2.RNASEQ.rds"
# geneRPK.full.fn     <- "../Data-processed/data.geneRPK.full.DRNA.ZOE2.RNASEQ.rds"
# geneRPK.full.fn.RNA <- "../Data-processed/data.geneRPK.full.RNA.ZOE2.RNASEQ.rds"
# geneRPK.marginal.fn <- "../Data-processed/data.geneRPK.marginal.DRNA.ZOE2.RNASEQ.rds"
# geneRPK.joint.fn    <- "../Data-processed/data.geneRPK.joint.DRNA.ZOE2.RNASEQ.rds"
# geneRPK.RNA.wide.sort.fn <- "../Data-processed/data.geneRPK.full.RNA.wide.sort.ZOE2.RNASEQ.rds"
    
    # 0.1 how many geneRPKs x species x subject? (# long form indexed by id(and ECC, age))
    gene.full.RPK.RNA %>% dim  # 2,703,591 x 7
    gene.full.RPK.RNA = gene.full.RPK.RNA %>% mutate(id = gsub("IDMT", "subj", id))
    
    # Into wide form: geneRPK.bact / geneRPK / bacteria / count.ID1 / count.ID2 / ...
    gene.full.RPK.RNA.wide = 
      gene.full.RPK.RNA %>% 
      transmute(gene.bact = paste0(gene, " ", bacteria), gene, bacteria, TPM, id) %>%
      spread(key = id, value = TPM, fill = 0)
    
    gene.full.RPK.RNA.wide.reads = 
      gene.full.RPK.RNA %>% 
      transmute(gene.bact = paste0(gene, " ", bacteria), gene, bacteria, reads, id) %>%
      spread(key = id, value = reads, fill = 0)
    
    gene.full.RPK.RNA.wide %>% dim  # 9,103    303 
    rm(gene.full.RPK.RNA)
    
    gene.full.RPK.RNA.wide.sort = gene.full.RPK.RNA.wide
    gene.full.RPK.RNA.wide.sort.reads = gene.full.RPK.RNA.wide.reads
    
    
    
    geneRPK.full.DRNA <- array(NA, dim = c(dim(gene.full.RPK.RNA.wide.sort)[1],
                                           dim(gene.full.RPK.RNA.wide.sort)[2] - 3, 2), 
                               dimnames = list(gene.full.RPK.RNA.wide.sort[,1], 
                                               colnames(gene.full.RPK.RNA.wide.sort)[-(1:3)], 
                                               c("DNA", "RNA")))
    geneRPK.full.DRNA.reads <- geneRPK.full.DRNA
    
    geneRPK.full.DRNA[, , 1] <- NA
    geneRPK.full.DRNA.reads[, , 1] <- NA
    # rm(geneRPK.full.DNA.wide.sort); gc()
    # gene.full.RPK.RNA.wide.sort <- readRDS("../Data-processed/tmp.data.gene.full.RPK.RNA.wide.sort.ZOE2.rds")
    geneRPK.full.DRNA[, , 2] <- gene.full.RPK.RNA.wide.sort[, -(1:3)] %>% as.matrix
    geneRPK.full.DRNA.reads[, , 2] <- gene.full.RPK.RNA.wide.sort.reads[, -(1:3)] %>% as.matrix
    row.names(gene.full.RPK.RNA.wide.sort) <- NULL
    row.names(gene.full.RPK.RNA.wide.sort.reads) <- NULL
    id.index = sapply(colnames(gene.full.RPK.RNA.wide.sort)[-(1:3)], function(x) which(outcome$id == x))
    geneRPK.full.DRNA <- list(otu = geneRPK.full.DRNA,
                              taxa = gene.full.RPK.RNA.wide.sort[,1:3],
                              meta = outcome[id.index, ],
                              reads = geneRPK.full.DRNA.reads)
    # This is actually a joint data. We haven't sum this up for each species.
    # Since this only has four species, marginal gene data are not available.
    saveRDS(geneRPK.full.DRNA, "Data-processed/data.geneRPK.full.DRNA.ZOE2.RNASEQ.rds")
    

### C02.data.reading.humann3_2021_path.R
### Created by Hunyong Cho
### input: Data/ZOE2_h3_RNA_path, Data/data.outcome.ZOE2_xxx.xlsx
### output: Data-processed/data.humann3.path.marginal.DRNA.ZOE2.rds, 
###         Data-processed/data.humann3.path.joint.DRNA.ZOE2.rds

### 0. library
    source("scripts/F_file.R")
    source("scripts/F_jointAnalysis.R")
    source("scripts/F_generic.R")
    library(dplyr); library(magrittr); library(tidyr)
    library(readxl)
    if (!dir.exists("Data-processed")) dir.create("Data-processed")

### 1. Data directory
    data.dir <- ("Data/ZOE2_h3_RNA_path")
    humann3 <- list()
    dir.all <- 
      list.dirs(path = data.dir, recursive=FALSE) %>%   # unnecessary to remove "./" ?
      grep("HUMANN3-(D|R)NA", ., value = TRUE)

    study.nm <- gsub(".*\\/", "", dir.all)

    
    for (i.tmp in 1:length(dir.all)) {
      path.tmp <- list.dirs(dir.all[i.tmp], recursive=T)
      path.tmp <- grep("(D|R)NA/[0-9A-Z]*$", path.tmp, value = TRUE) # removing the root directory
      humann3[[i.tmp]] <- data.frame(path = path.tmp, stringsAsFactors=FALSE)
      humann3[[i.tmp]]$no = gsub(".*/", "", humann3[[i.tmp]]$path)  #extracting the current folder names
      humann3[[i.tmp]] = humann3[[i.tmp]][!humann3[[i.tmp]]$no %in% c("Blank", "00000"),]
    }
    names(humann3) <- study.nm

    
### 2. Humann3 pathway data concatenation
    
    ### 2.1 RNA
      ### RNA
      print("path-RNA")
      path.fn.zoe2 <- "Data_propcessed/data.humann3.path.full.RNA.ZOE2.rds"
      
      ### ZOE 2.0
      path.full.RNA <- jAnalyze(direction="long", study = "HUMANN3-RNA", id = "all", type="pathabundance", info.file=humann3)
      print(dim(path.full.RNA)) # 718,881 x 4
      path.full.RNA$id = gsub("IDM.", "subj", path.full.RNA$id)
      
    ### 2.2 DNA
      # path.fn.zoe2 <- gsub("RNA", "DNA", path.fn.zoe2)
      # 
      # ### ZOE 2.0
      # path.full.DNA <- jAnalyze(direction="long", study = "HUMANN3-DNA", id = "all", type="pathabundance", info.file=humann3)
      # path.full.DNA$id %<>% as.integer
      # print(dim(path.full.DNA))
      
      
### 3. outcome files
      outcome.ZOE2.trait <- read_xlsx("Data/data.outcome.ZOE2_main_20200626.xlsx")
      # outcome.ZOE1.trait <- read_xlsx("Data/data.outcome.ZOE2_pilot_20200527.xlsx")
      
      
### 4. Wide form and array of DNA and RNA together.
      outcome = outcome.ZOE2.trait
      outcome =
        mutate(outcome, id = gsub("IDMG", "subj", outcome$id_MTG)) %>%  # defining a common id across MTG and MTX.
        dplyr::select(id, everything())
      
      path.full.fn = "Data-processed/data.humann3.path.full.DRNA.ZOE2.rds"
      path.marginal.fn <- "Data-processed/data.humann3.path.marginal.DRNA.ZOE2.rds"
      path.joint.fn    <- "Data-processed/data.humann3.path.joint.DRNA.ZOE2.rds"
          
      # How many humann3 x species x subject? (# long form indexed by id(and ECC, age))
      path.full.RNA %>% dim  # 718,881
          
      # Into wide form: bracken.bact / bracken / bacteria / count.ID1 / count.ID2 / ...
      path.full.RNA.wide.RPK = 
        path.full.RNA %>% 
        transmute(pathbact = paste0(path, " ", bacteria), 
                  path, bacteria, count, id) %>%
        spread(key = id, value = count, fill = 0)
          
      path.full.RNA.wide.RPK %>% dim  # 10,046   300
          
      tmp <- names(path.full.RNA.wide.RPK)
          
      # ordering the subjects - RNA
      tmp <- matrix(NA, nrow = dim(path.full.RNA.wide.RPK)[1], ncol = length(outcome$id) + 3)
      tmp <- as.data.frame(tmp, row.names = path.full.RNA.wide.RPK[, 1])
      names(tmp) <- c("pathbact", "path", "bacteria", outcome$id)
      tmp$pathbact = path.full.RNA.wide.RPK$pathbact
      tmp$path = path.full.RNA.wide.RPK$path
      tmp$bacteria = path.full.RNA.wide.RPK$bacteria
      for (i in outcome$id) {
        j = i %>% as.character
        if (j %in% names(path.full.RNA.wide.RPK)) {
          tmp[, j] = path.full.RNA.wide.RPK[, j]
        }
      }
      path.full.RNA.wide.RPK <- tmp; rm(tmp)
      
      for (i in 1:3) path.full.RNA.wide.RPK[, i] %<>% as.character
          
    
      # row ordering
      path.full.RNA.wide.sort <- path.full.RNA.wide.RPK %>% arrange(pathbact)
          
      ### putting DNA and RNA #########
      # make an array of (overlapped humann3) x (sample) x (DNA+RNA)
      path.full.DRNA <- array(NA, dim = c(dim(path.full.RNA.wide.sort)[1],
                                             dim(path.full.RNA.wide.sort)[2] - 3, 2), 
                                 dimnames = list(path.full.RNA.wide.sort[, 1], 
                                                 outcome$id, 
                                                 c("DNA", "RNA")))
      path.full.DRNA[, , 2] <- path.full.RNA.wide.sort[, -(1:3)] %>% as.matrix
      
      row.names(path.full.RNA.wide.sort) <- NULL
      path.full.DRNA <- list(otu = path.full.DRNA,
                                taxa = as.data.frame(path.full.RNA.wide.sort[, 1:3]),
                                meta = outcome)
      
      # splitting marginal & joint data from full data
      path.full.DRNA %$% otu[taxa$bacteria == "(TOTAL)",, ] -> path.marginal.DRNA
      path.full.DRNA %$% otu[taxa$bacteria != "(TOTAL)" & taxa$bacteria != "(GAP)",,] -> path.joint.DRNA
      path.marginal.DRNA <- list(otu = path.marginal.DRNA,
                                    taxa = path.full.DRNA$taxa[path.full.DRNA$taxa$bacteria == "(TOTAL)", ],
                                    meta = outcome)
      path.joint.DRNA <- list(otu = path.joint.DRNA,
                                 taxa = path.full.DRNA$taxa[
                                   path.full.DRNA$taxa$bacteria != "(TOTAL)"& 
                                     path.full.DRNA$taxa$bacteria != "(GAP)", ],
                                 meta = outcome)
      rownames(path.marginal.DRNA$taxa) <- rownames(path.joint.DRNA$taxa) <- NULL
      saveRDS(path.marginal.DRNA, path.marginal.fn)
      saveRDS(path.joint.DRNA,    path.joint.fn)
      rm(path.full.RNA.wide.sort, path.full.fn, path.full.DRNA, path.marginal.fn, path.joint.fn,
         path.fn.zoe2, path.tmp)
    
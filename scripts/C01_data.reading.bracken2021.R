### C01.data.reading.bracken2021.R
### Created by Hunyong Cho
### input: Data/ZOE2_xxx_Bracken, Data/data.outcome.ZOE2_xxx.xlsx
### output: Data-processed/data.bracken2.full.DRNA.ZOE2.rds

### 0. library
    source("scripts/F_file.R")
    source("scripts/F_jointAnalysis.R")
    source("scripts/F_generic.R")
    library(dplyr); library(magrittr); library(tidyr)
    library(readxl)
    if (!dir.exists("Data-processed")) dir.create("Data-processed")
    
### 1. Data directory
    data.dir2 <- ("Data/ZOE2_main_Bracken")
    data.dir1 <- ("Data/ZOE2_pilot_Bracken")

    bracken2 <- list()
    dir.2 <- list.dirs(path = data.dir2, recursive=FALSE)  # unnecessary to remove "./" ?
    dir.1 <- list.dirs(path = data.dir1, recursive=FALSE)  # unnecessary to remove "./" ?
    dir.all  <- c(dir.2, dir.1)
    
    study.nm <- gsub(".*\\/", "", dir.all)
    
    for (i.tmp in 1:length(dir.all)) {
      path.tmp <- list.dirs(dir.all[i.tmp])
      file.tmp <- list.files(dir.all[i.tmp])
      bracken2[[i.tmp]] <- 
        data.frame(path = paste0(path.tmp, "/", file.tmp), stringsAsFactors=FALSE) %>% 
        subset(!grepl("Blank", path))
      bracken2[[i.tmp]]$no <- gsub("(.*/)(.*)(\\.bracken.out$)", "\\2", bracken2[[i.tmp]]$path)  #extracting the current folder names
    }
    names(bracken2) <- study.nm
    rm(dir.all, dir.1, dir.2, data.dir1, data.dir2, i.tmp, file.tmp, path.tmp)
    print(bracken2)
    sapply(bracken2, dim)
  
    
### 2. bracken data concatenation
    ### 2.1 Concatenation of files into a matrix
    bracken2_2017D <- jAnalyze(direction="long", study = "2017-DNA", id = "all", type="bracken", info.file=bracken2)
    bracken2_2018D <- jAnalyze(direction="long", study = "2018-DNA", id = "all", type="bracken", info.file=bracken2)
    bracken2_2019D <- jAnalyze(direction="long", study = "2019-DNA", id = "all", type="bracken", info.file=bracken2)
    
    bracken2_2017R <- jAnalyze(direction="long", study = "2017-RNA", id = "all", type="bracken", info.file=bracken2)
    bracken2_2018R <- jAnalyze(direction="long", study = "2018-RNA", id = "all", type="bracken", info.file=bracken2)
    bracken2_2019R <- jAnalyze(direction="long", study = "2019-RNA", id = "all", type="bracken", info.file=bracken2)
    
    ### 2.2 RNA
    # ZOE pilot: bracken2_2017R
    # ZOE 2.0 main: bracken2_2018R, bracken2_2019R
    bracken.full.RNA.wv <- rbind(bracken2_2018R, bracken2_2019R) # with virus
    if (!all(is.numeric(bracken.full.RNA.wv$reads))) {bracken.full.RNA.wv$reads %<>% as.character %>% as.numeric}
    if (!all(is.numeric(bracken2_2017R$reads))) {bracken2_2017R$reads %<>% as.character %>% as.numeric}
    print(dim(bracken.full.RNA.wv)) # 512,968 x 4
    print(dim(bracken2_2017R)) # 112,172 x 4
    
    # Filtering viruses
    virus.list <- grep(virus.keyword, bracken.full.RNA.wv$bacteria, ignore.case = TRUE)
    bracken.full.RNA.ZOE2 <- bracken.full.RNA.wv[-virus.list, ]
    print(dim(bracken.full.RNA.ZOE2)) # 510,011 x 4
    
    virus.list <- grep(virus.keyword, bracken2_2017R$bacteria, ignore.case = TRUE)
    bracken.full.RNA.ZOE1 <- bracken2_2017R[-virus.list, ]
    print(dim(bracken.full.RNA.ZOE1)) # 111,651 x 4
    
    rm(bracken2_2017R, bracken2_2018R, bracken2_2019R, bracken.full.RNA.wv)
    
    ### 2.3 DNA
    bracken.full.DNA.wv <- rbind(bracken2_2018D, bracken2_2019D)
    if (!all(is.numeric(bracken.full.DNA.wv$reads))) {bracken.full.DNA.wv$reads %<>% as.character %>% as.numeric}
    if (!all(is.numeric(bracken2_2017D$reads))) {bracken2_2017D$reads %<>% as.character %>% as.numeric}
    print(dim(bracken.full.DNA.wv)) # 436,760 x 4
    print(dim(bracken2_2017D))   # 116,455 x 4
    
    virus.list <- grep(virus.keyword, bracken.full.DNA.wv$bacteria, ignore.case = TRUE)
    bracken.full.DNA.ZOE2 <- bracken.full.DNA.wv[-virus.list, ]
    virus.list <- grep(virus.keyword, bracken2_2017D$bacteria, ignore.case = TRUE)
    bracken.full.DNA.ZOE1 <- bracken2_2017D[-virus.list, ]
    
    rm(bracken2_2017D, bracken2_2018D, bracken2_2019D, bracken.full.DNA.wv)
    rm(virus.list, virus.keyword)
  
    
    
### 3. outcome files
    outcome.ZOE2.trait <- read_xlsx("Data/data.outcome.ZOE2_main_20200626.xlsx")
    outcome.ZOE1.trait <- read_xlsx("Data/data.outcome.ZOE2_pilot_20200527.xlsx")
    
  
    
### 4. Wide form and array of DNA and RNA together.
    
    for (zoe in 2:1) {
      cat("ZOE", zoe, "\n")
      
      outcome <- switch(zoe, `1` = outcome.ZOE1.trait, `2` = outcome.ZOE2.trait)
      outcome =
        mutate(outcome, id = gsub("IDMG", "subj", outcome$id_MTG)) %>%  # defining a common id across MTG and MTX.
        dplyr::select(id, everything())
      bracken.full.RNA <- switch(zoe, `1` = bracken.full.RNA.ZOE1, `2` = bracken.full.RNA.ZOE2)
      bracken.full.DNA <- switch(zoe, `1` = bracken.full.DNA.ZOE1, `2` = bracken.full.DNA.ZOE2)
      bracken.full.fn  <- sprintf("Data-processed/data.bracken2.full.DRNA.ZOE%s.rds", zoe)
      
        
      # Into wide form: bracken.bact / bracken / bacteria / count.ID1 / count.ID2 / ...
      bracken.full.RNA.wide.reads = 
        bracken.full.RNA %>% 
        transmute(bacteria, reads, id = gsub("IDM.", "subj", id)) %>%
        spread(key = id, value = reads, fill = 0)
      
      bracken.full.DNA.wide.reads = 
        bracken.full.DNA %>% 
        transmute(bacteria, reads, id = gsub("IDM.", "subj", id)) %>%
        spread(key = id, value = reads, fill = 0)
      
      bracken.full.RNA.wide.reads %>% dim  # 6268   300
      bracken.full.DNA.wide.reads %>% dim  # 5972   303
      
      #rm(bracken.full.RNA, bracken.full.DNA)
      
      # filling in missing subjects: Not needed for bracken. Same subjects.
      tmp <- names(bracken.full.DNA.wide.reads)
      tmp <- tmp[which(!tmp %in% names(bracken.full.RNA.wide.reads))]
      bracken.full.RNA.wide.reads[, tmp] <- NA
      
      tmp <- names(bracken.full.RNA.wide.reads)
      tmp <- tmp[which(!tmp %in% names(bracken.full.DNA.wide.reads))]
      bracken.full.DNA.wide.reads[, tmp] <- NA
      
      
      # ordering the subjects - DNA
      tmp <- matrix(NA, nrow = dim(bracken.full.DNA.wide.reads)[1], ncol = length(outcome$id) + 1)
      tmp <- as.data.frame(tmp, row.names = bracken.full.DNA.wide.reads[, 1])
      names(tmp) <- c("bacteria", outcome$id)
      tmp$bacteria = bracken.full.DNA.wide.reads$bacteria
      for (i in outcome$id) {
        if (i %in% names(bracken.full.DNA.wide.reads)) {
          tmp[, as.character(i)] = bracken.full.DNA.wide.reads[, as.character(i)]
        }
      }
      bracken.full.DNA.wide.reads <- tmp; rm(tmp)
    
      # ordering the subjects - RNA
      tmp <- matrix(NA, nrow = dim(bracken.full.RNA.wide.reads)[1], ncol = length(outcome$id) + 1)
      tmp <- as.data.frame(tmp, row.names = bracken.full.RNA.wide.reads[, 1])
      names(tmp) <- c("bacteria", outcome$id)
      tmp$bacteria = bracken.full.RNA.wide.reads$bacteria
      for (i in outcome$id) {
        if (i %in% names(bracken.full.RNA.wide.reads)) {
          tmp[, as.character(i)] = bracken.full.RNA.wide.reads[, as.character(i)]
        }
      }
      bracken.full.RNA.wide.reads <- tmp; rm(tmp)
      all(names(bracken.full.RNA.wide.reads)[-1] == names(bracken.full.DNA.wide.reads)[-1])
      
      bracken.full.RNA.wide.reads[, 1] %<>% as.character
      bracken.full.DNA.wide.reads[, 1] %<>% as.character
      
      # how many overlapped?
      bracken.bact.RNAinDNA = which (bracken.full.RNA.wide.reads[,1] %in% bracken.full.DNA.wide.reads[,1])
      bracken.bact.names = bracken.full.RNA.wide.reads[bracken.bact.RNAinDNA, 1]
      bracken.bact.DNAinRNA = which (bracken.full.DNA.wide.reads[,1] %in% bracken.bact.names)
      bracken.bact.names.D = bracken.full.DNA.wide.reads[bracken.bact.DNAinRNA, 1]
      
      bracken.bact.RNAonly = which (!bracken.full.RNA.wide.reads[,1] %in% bracken.full.DNA.wide.reads[,1])
      bracken.bact.RNAonly.names = bracken.full.RNA.wide.reads[bracken.bact.RNAonly, 1]
      bracken.bact.DNAonly = which (!bracken.full.DNA.wide.reads[,1] %in% bracken.bact.names)
      bracken.bact.DNAonly.names = bracken.full.DNA.wide.reads[bracken.bact.DNAonly, 1]
      
      ###
      length(bracken.bact.names)   # 5829
      length(bracken.bact.names.D) # 5829
      length(bracken.bact.RNAonly.names)   # 439
      length(bracken.bact.DNAonly.names)   # 143
      
      # 1.0 data manipulation
      bracken.full.RNA.wide.sort <- bracken.full.RNA.wide.reads[bracken.bact.RNAinDNA,]
      bracken.full.RNA.wide.sort %<>%  arrange(bacteria)
      bracken.full.DNA.wide.sort <- bracken.full.DNA.wide.reads[bracken.bact.DNAinRNA,]
      bracken.full.DNA.wide.sort %<>%  arrange(bacteria)
      bracken.bact.names.sort <- bracken.full.RNA.wide.sort[,1]
      # make sure the subjects and brackens are in the same order
      all(names(bracken.full.RNA.wide.sort)[-1] == names(bracken.full.DNA.wide.sort)[-1])
      all(as.character(bracken.full.RNA.wide.sort[,1]) == as.character(bracken.full.DNA.wide.sort[,1]))
      all(names(bracken.full.RNA.wide.sort)[-1] == outcome$id)
      
      if (length(bracken.bact.DNAonly)) {
        tmp <- bracken.full.DNA.wide.reads[bracken.bact.DNAonly, ]
        tmp[, -1] <- 0
      } else {
        tmp <- NULL
      }
      bracken.full.RNA.wide.sort <- rbind(bracken.full.RNA.wide.sort,                 # shared
                                          bracken.full.RNA.wide.reads[bracken.bact.RNAonly, ], # RNA only
                                          tmp)                                     # DNA only (zeros)
      
      if (length(bracken.bact.RNAonly)) {
        tmp <- bracken.full.RNA.wide.reads[bracken.bact.RNAonly, ]
        tmp[, -1] <- 0
      } else {
        tmp <- NULL
      }
      bracken.full.DNA.wide.sort <- rbind(bracken.full.DNA.wide.sort,                 # shared
                                          tmp,                                     # RNA only (zeros)
                                          bracken.full.DNA.wide.reads[bracken.bact.DNAonly,])  # DNA only
      
      # sorting once more, after putting them altogether
      index <- order(bracken.full.RNA.wide.sort[,1])
      bracken.full.RNA.wide.sort <- bracken.full.RNA.wide.sort[index, ]
      bracken.full.DNA.wide.sort <- bracken.full.DNA.wide.sort[index, ]
      
      all(bracken.full.DNA.wide.sort[,1] == bracken.full.RNA.wide.sort[,1])
      rm(bracken.full.DNA.wide.reads, bracken.full.RNA.wide.reads, bracken.bact.RNAinDNA, bracken.bact.DNAinRNA,
         tmp, index, bracken.bact.DNAonly, bracken.bact.DNAonly.names,
         bracken.bact.RNAonly, bracken.bact.RNAonly.names,
         bracken.bact.names.sort, bracken.bact.names, bracken.bact.names.D)
      gc()
      
      
      # make an array of (overlapped bracken) x (sample) x (DNA+RNA)
      bracken.full.DRNA <- array(NA, dim = c(dim(bracken.full.DNA.wide.sort)[1],
                                             dim(bracken.full.DNA.wide.sort)[2] - 1, 2), 
                                 dimnames = list(bracken.full.DNA.wide.sort[,1], 
                                                 outcome$id, 
                                                 c("DNA", "RNA")))
      bracken.full.DRNA[, , 1] <- bracken.full.DNA.wide.sort[, -1] %>% as.matrix
      bracken.full.DRNA[, , 2] <- bracken.full.RNA.wide.sort[, -1] %>% as.matrix
      row.names(bracken.full.DNA.wide.sort) <- NULL
      row.names(bracken.full.RNA.wide.sort) <- NULL
      bracken.full.DRNA <- list(otu = bracken.full.DRNA,
                                taxa = data.frame(bacteria = bracken.full.DNA.wide.sort[,1]),
                                meta = outcome)
      
      # splitting marginal & joint data from full data
      saveRDS(bracken.full.DRNA, bracken.full.fn)
      rm(bracken.full.DNA.wide.sort, bracken.full.RNA.wide.sort,
         bracken.full.fn, bracken.full.DRNA)
    }
    
    
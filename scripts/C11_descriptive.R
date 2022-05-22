### C11 descriptive.R
### Author: Hunyong Cho
### output: Demographics corresponds to Ext. Dat. Table 5.
  
  bracken = TRUE; humann = 2
  type = "bact"; DR.no = 1; 
  pheno = "t3c"
  nrm = TRUE #normalization
  brk = ifelse(bracken, "bracken", "humman2")
  trait <- c("t3c", "micro_t3c_b")
  
  library(tidyverse); library(magrittr); library(energy)
  source("scripts/F_generic.R")  # typeDRNA(), sample.exclude, taxa.exclude
  source("scripts/F_eco.R")  # typeDRNA(), sample.exclude, taxa.exclude
  
  for (ZOE in 1:2) {
    cat("############ ", ifelse(ZOE == 2, "ZOE 2.0 main", "ZOE pilot"), "############\n")  
    
    # Subjects filtered / Virus not included / TPM-Normalized (=rescaled) / noninformative taxa excluded
    .initialize(type = type, ZOE = ZOE, DR.no = DR.no,
                test = "none", add.epsilon = F, 
                filter = "both", screen = FALSE, humann = humann, bracken = bracken)
    
    cat("\n\n### 1. phenotypes (mean, sd)\n")  
    dat$meta[, trait] %>% 
      apply(2, function(x) c(mean = mean(x, na.rm = TRUE),
                             sd = sd(x, na.rm = TRUE))) %>% 
      print
    
    cat("\n\n### 2. age (mean, sd)\n")  
    dat$meta$agemo %>% 
      {c(mean = mean(., na.rm = TRUE),
         sd = sd(., na.rm = TRUE))} %>% 
      print
    
    cat("\n\n### 3. race (#, %)\n")  
    dat$meta$race %>%
      summary %>% print %>%
      {./sum(.)} %>% round(2)
    
    cat("\n\n### 3. sex (#, %)\n")  
    dat$meta$sex %>%
      table %>% print %>%
      {./sum(.)} %>% round(2)
    
  }

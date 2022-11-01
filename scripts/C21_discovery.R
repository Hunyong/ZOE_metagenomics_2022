### C21 test_main.R
### Author: Hunyong Cho
### Identify the species with significant association with the ECC.


############################################################################################
### 0.1 library
  library(dplyr)
  source("scripts/F_generic.R")  # typeDRNA(), sample.exclude, taxa.exclude

### 0.2 parameters
  ZOE = 2  # ZOE == 2: ZOE 2.0 main study, ZOE == 1: ZOE 2.0 pilot study
  type =  "bact"
  bracken = TRUE
  nrm = TRUE #normalization
  

### Tests for discovery
  for (DR.no in 1:2) { # DR.no == 1: DNA, 2: RNA
    for (pheno in c("t3c", "microt3cb")) { # t3c = person level quantitative trait,
                                           # microt3cb = localized quantitative trait
      
      ### 1. Test prep
      ## 1.0 data preparation (screening, filtering, etc)
      .initialize(type = type, ZOE = ZOE, DR.no = DR.no, pheno = pheno,
                  test = "LN", screen = TRUE, prev.threshold = 0.1, 
                  avg.detect = 0.00001, detect.rel.abund = TRUE, nrm = nrm,
                  bracken = bracken, humann = humann,
                  screen.among.tested.DNAs = TRUE)
      
      ## 1.1 Model
      cat("ZOE ", ZOE, ", ", DR.no, "", pheno,"\n")
      binaryDNA = ifelse(DR.no == 2, "+ binaryDNA", "")
      model1 = paste0("logUnit ~ batcheffect + phenotype", binaryDNA, " + agemo + race")
      # if there is no batches, drop the term.
      if (na.omit(dat.reg)$batcheffect %>% unique %>% length %>% "<"(2)) {
        model1 <- gsub("batcheffect \\+ ", "", model1) 
      }
      model1 %<>% as.formula
      
      ## 1.2 Skeleton    
      dat.reg$logUnit = log2(as.numeric(dat$otu[1, , DR.no]))
      full <- lm(model1, data = dat.reg)
      nm <- names(full$coefficients)
      result.coef <- result.pval <- 
        matrix(NA, nrow = n.taxa, ncol = length(nm), dimnames = list(dat$taxa[, 1], nm)) %>% 
        as.data.frame()
      # #names of interest
      # nm.int <- grep("phenotype", nm)
          
      nm.common = paste0("bracken_", DRNA, "-", pheno, "-ZOE", ZOE)
      output.nm.final = paste0("output/C2_LM_", nm.common,  ".rds")
          
      ### 2 Testing
      ## 2.1 testing
      for (i in 1:n.taxa) {
        if (!i %% 30) cat(i, " out of ", n.taxa, "  ")
        dat.reg$Unit = as.numeric(dat$otu[i, , DR.no])
        dat.reg$logUnit = log2(as.numeric(dat$otu[i, , DR.no]))
        dat.reg$binaryDNA = ifelse(as.numeric(dat$otu[i, , 1]), 1, 0)
        full <- lm(model1, data = dat.reg)
        full.summary <- full %>% summary %>% coef
        result.coef[i, nm %in% rownames(full.summary)] <- full.summary[, "Estimate"]
        result.pval[i, nm %in% rownames(full.summary)] <- (full %>% summary %>% coef)[, "Pr(>|t|)"]
      }
      result <- list(screen = dat$screen,
                     epsilon = dat$epsilon,
                     coef = result.coef,
                     pval = result.pval)
            
      ## 2.2 p-adjusting
      dim.int = dim(result$pval[, "phenotype", drop = FALSE])
      dimnames.int = dimnames(result$pval[, "phenotype", drop = FALSE])
      result$pval.adj <-
        result$pval %>% 
        apply(2, function(s) p.adjust(s, method = "BH")) %>% 
        as.data.frame
      # BH adjustment (and take the min of both ECC1 and 2)
      result$pval.adj[, "phenotype", drop = FALSE] %>% 
        apply(1, min, na.rm = TRUE) -> result$pval.adj$min   #for CF, apply(1, min) is just identity
      result$sigFeature = dat$taxa[result$pval.adj$min < 0.05, 1]
      
      saveRDS(result, output.nm.final)
      
    }
  }

  
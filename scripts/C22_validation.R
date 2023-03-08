### C22_validation.R
### Author: Hunyong Cho
### Validating with the ZOE pilot data.


### 0. library
  library(dplyr)
  
  source("scripts/F_generic.R")
  source("scripts/F_aesthetics.R")
  res = function(DR = "DNA", pheno = "t3c", ZOE = 2, cutoff = 0.05) {
    result = 
      sprintf("output/C2_LM_bracken_%s-%s-ZOE%d.rds", DR, pheno, ZOE) %>% 
      readRDS()
    sig = result$pval.adj$pheno <= cutoff
    data.frame(species = rownames(result$pval.adj) [sig], 
               coef = result$coef[sig, "phenotype"], 
               pval = result$pval[sig, "phenotype"],
               pval.adj = result$pval.adj[sig, "phenotype"])
  }
  .two_sided_ci = function(stat, p) {
    se = abs(stat)/qnorm(1 - p/2)
    ci = stat + c(-1,1) * qnorm(0.975) * se
    return(ci)
  }
  two_sided_ci = Vectorize(.two_sided_ci)
  

##### 1. Identifying the 23 species all FDR significant, discovered by the main data
  type = "bact"
  bracken = TRUE
  
  key.bact23 = Reduce(intersect, list(res("DNA", "t3c")$species,
                                      res("DNA", "microt3cb")$species,
                                      res("RNA", "t3c")$species,
                                      res("RNA", "microt3cb")$species)) %>% sort

##### 2. validation of the 23 species. Identifying the 16 species significant (at least one of the four) in the replication sample
  ZOE = 1
  counter = 1  # initializing the counter
  coef.plot = list()
  raw_data_file = tibble()
  for (DR.no in 1:2) {
    DRNA = c("DNA", "RNA")[DR.no]
    
    for (pheno in c("microt3cb", "t3c")) {
      plot.nm = paste0(DRNA, "_", pheno)
      cat(plot.nm, "\n")
      
      ## 1. Reading data
      .initialize(type = type, ZOE = ZOE, DR.no = DR.no, pheno = pheno,
                  test = "LN", screen = TRUE, prev.threshold = 0.1,
                  avg.detect = 0.00001, detect.rel.abund = TRUE, nrm = TRUE,
                  humann = 2, bracken = TRUE, screen.among.tested.DNAs = TRUE)
      sig.index = sapply(key.bact23, function(x) {tmp = which(dat$taxa[,1] == x); if (length(tmp)) tmp else NA})
      n.sig = length(sig.index)
      
      
      ## 2. reading the main analysis results
      ## 2.1 setting file names
      output.nm.validation = paste0("output/C2_LM_bracken_", DRNA, "-", pheno, "-validation", ".rds")
      output.nm.combined = gsub("-validation\\.rds", "-combined.rds", output.nm.validation)
      output.nm.combined_csv = gsub("\\.rds", ".csv", output.nm.combined)

      ## 2. running the validation tests
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
        matrix(NA, nrow = n.sig, ncol = length(nm), dimnames = list(dat$taxa[sig.index, 1], nm)) %>% 
        as.data.frame()
      
      for (k in 1:n.sig) {
        i = sig.index[k]
        if (is.na(i)) next
        dat.reg$Unit = as.numeric(dat$otu[i, , DR.no])
        dat.reg$logUnit = log2(as.numeric(dat$otu[i, , DR.no]))
        dat.reg$binaryDNA = ifelse(as.numeric(dat$otu[i, , 1]), 1, 0)
        full <- lm(model1, data = dat.reg)
        full.summary <- full %>% summary %>% coef
        result.coef[k, nm %in% rownames(full.summary)] <- full.summary[, "Estimate"]
        result.pval[k, nm %in% rownames(full.summary)] <- (full %>% summary %>% coef)[, "Pr(>|t|)"]
      }
      result <- list(epsilon = dat$epsilon,
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
      result$sigFeature = key.bact23[result$pval.adj$min < 0.05]
      
      saveRDS(result, output.nm.validation)
      
      
      ZOE2.output = res(DR = DRNA, pheno = pheno, ZOE = 2, cutoff = 1) # all outputs
      
      sig.index.ZOE2 = 
        sapply(key.bact23, 
               function(x) {a = which(ZOE2.output$species == x); ifelse(length(a)==0, NA, a)})
      
      combined = # Combining ZOE 2 main + pilot test information
        tibble(species = dat$taxa[sig.index, 1],
               genus = factor(gsub(" .*", "", dat$taxa[sig.index, 1])),
               zoe2 = ZOE2.output$coef[sig.index.ZOE2],
               zoe2.p = ZOE2.output$pval[sig.index.ZOE2],
               zoe2.q = ZOE2.output$pval.adj[sig.index.ZOE2],
               zoe1 = result$coef[, "phenotype"],
               zoe1.p = result$pval[, "phenotype"],
               zoe1.q = result$pval.adj[, "phenotype"])
      ci2 = two_sided_ci(combined$zoe2, combined$zoe2.p)
      combined$zoe2_lb = ci2[1, ]
      combined$zoe2_ub = ci2[2, ]
      ci1 = two_sided_ci(combined$zoe1, combined$zoe1.p)
      combined$zoe1_lb = ci1[1, ]
      combined$zoe1_ub = ci1[2, ]
      
      saveRDS(combined, output.nm.combined)
      
      # Supplementary tables
      raw_data_file = rbind(raw_data_file, combined %>% mutate(type = DRNA, phenotype = pheno))
      
      counter = counter + 1
    }
  }
  
  library(openxlsx)
  raw_data_file %>% 
    transmute(
      type, phenotype,
      species, genus, 
      zoe2_coef = zoe2, 
      zoe2_lb, zoe2_ub, 
      zoe2_p = zoe2.p,
      zoe2_q = zoe2.q,
      zoe1_coef = zoe1, 
      zoe1_lb, zoe1_ub,
      zoe1_p = zoe1.p,
      zoe1_q = zoe1.q
    ) %>% 
    write.xlsx("figure/raw_data_for_Fig2_and_EFig1.xlsx", sheetName = "Fig2_and_EFig1")
    
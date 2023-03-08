### C32_test_pathMarginal.R
### Author: Hunyong Cho
### Identify the pathways with significant association with the ECC 
###   among the pathways that contain the 16 significant species.
### output: E. Table 2, E. Fig 2

### 0.1 library
  library(dplyr); library(tibble)
  library(ggplot2); library(ggrepel)
  source("scripts/F_aesthetics.R")  # For color pallette
  source("scripts/F_generic.R")  # typeDRNA(), sample.exclude, taxa.exclude

  type =  "path"
  DR.no = 2 ;  ZOE <- 2;
  model.nm = "cov"
  bracken = FALSE
  humann = 3     
  nrm = TRUE #normalization
  
  # The final list of 16 species.
  key.bact16 = read.csv("output/ETable1_C2_coef_table_bracken.csv")$species
  
  
# from C31_path_top.R
  tab2b = readRDS(sprintf("output/C31_pathway_composition_humann3_RNA_sigTaxa16_Marginal.rds"))
  tab2c = read.csv(sprintf("output/C31_pathway_top30_humann3_RNA_sigTaxa16_Marginal.csv"))
  path.with.key.bact = tab2c$pathway

  ### 1. test of marginal pathways.
  for (pheno in c("t3c", "microt3cb")) {
    .initialize(type = "path", ZOE = 2, DR.no = 2, pheno = pheno,
                bracken = FALSE, humann = 3,
                test = "LN", add.epsilon = TRUE, screen = TRUE, 
                nrm = TRUE, nrmScale = 400000, # original scale 383K
                prev.threshold = 0.1,  avg.detect = 1e-8, detect.rel.abund = TRUE)
    
    model1 = logUnit ~ batcheffect + phenotype + agemo + race
    
    dat.reg$logUnit = log2(as.numeric(dat$otu[1, , 2]))
    full <- lm(model1, data = dat.reg)
    nm <- names(full$coefficients)
    result.coef <- result.pval <- result.se <-
      matrix(NA, nrow = length(path.with.key.bact), ncol = length(nm), 
             dimnames = list(path.with.key.bact, nm)) %>% 
      as.data.frame()
    
    nm.common = sprintf("humann3_sigTaxa16_Marginal_path-RNA-%s-ZOE2", pheno)
    output.nm.final = paste0("output/C32_LM_", nm.common,  ".rds")
    
    ## 2.1 testing!
    for (i in path.with.key.bact) {  ### testing only important bacterias
      cat(i,  "  ")
      dat.reg$Unit = as.numeric(dat$otu[i, , 2])
      dat.reg$logUnit = log2(as.numeric(dat$otu[i, , 2]))
      dat.reg$binaryDNA = ifelse(as.numeric(dat$otu[i, , 1]), 1, 0)
      full <- lm(model1, data = dat.reg)
      full.summary <- full %>% summary %>% coef
      result.coef[i, nm %in% rownames(full.summary)] <- full.summary[, "Estimate"]
      result.se[i, nm %in% rownames(full.summary)] <- full.summary[, "Std. Error"]
      result.pval[i, nm %in% rownames(full.summary)] <- (full %>% summary %>% coef)[, "Pr(>|t|)"]
    }
    result <- list(screen = dat$screen,
                   epsilon = dat$epsilon,
                   coef = result.coef,
                   se = result.se,
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
    result$sigFeature = path.with.key.bact[result$pval.adj$min < 0.5]
    
    saveRDS(result, output.nm.final)
    
  }
  
### 2. Summary of test results
  pth.names = data.frame(
    pathway = c("LACTOSECAT-PWY", "PWY-1042", "PWY0-1296", "PWY-7111", "PWY-5097", "PWY-2942", "PWY-6121", "UDPNAGSYN-PWY", "PWY-6386", "PWY-6737", "PWY-6122", "PWY-6277", "PWY-6163", "COMPLETE-ARO-PWY", "PWY-6387", "GLYCOGENSYNTH-PWY", "ARGSYNBSUB-PWY", "ARO-PWY", "PWY-7219", "PWY-6700", "PWY-5100", "PWY-6609", "VALSYN-PWY", "HSERMETANA-PWY", "PWY-7221", "DTDPRHAMSYN-PWY", "BRANCHED-CHAIN-AA-SYN-PWY", "OANTIGEN-PWY", "TRNA-CHARGING-PWY", "PWY-5103"),
    pathname = c("lactose and galactose degradation I", "glycolysis IV"  , "purine ribonucleosides degradation", "pyruvate fermentation to isobutanol (engineered)"  , "L-lysine biosynthesis VI", "L-lysine biosynthesis III"  , "5-aminoimidazole ribonucleotide biosynthesis I", "UDP-N-acetyl-D-glucosamine biosynthesis I"  , "UDP-N-acetylmuramoyl-pentapeptide biosynthesis II (lysine-containing)", "starch degradation V"  , "5-aminoimidazole ribonucleotide biosynthesis II", "superpathway of 5-aminoimidazole ribonucleotide biosynthesis"  , "chorismate biosynthesis from 3-dehydroquinate", "superpathway of aromatic amino acid biosynthesis"  , "UDP-N-acetylmuramoyl-pentapeptide biosynthesis I", "glycogen biosynthesis I (from ADP-D-Glucose)"  , "L-arginine biosynthesis II (acetyl cycle)", "chorismate biosynthesis I"  , "adenosine ribonucleotides de novo biosynthesis", "queuosine biosynthesis I (de novo)"  , "pyruvate fermentation to acetate and lactate II", "adenine and adenosine salvage III"  , "L-valine biosynthesis", "L-methionine biosynthesis III"  , "guanosine ribonucleotides de novo biosynthesis", "dTDP-β-L-rhamnose biosynthesis"  , "superpathway of branched chain amino acid biosynthesis", "O-antigen building blocks biosynthesis (E. coli)"  , "tRNA charging", "L-isoleucine biosynthesis III"))
  
### 3. pathway-marginal test results (involving the key species)
  # Outputs from C05.01.CovLM.pathMarginal.R
  DR.test = "RNA"
  
  pheno.test = "t3c"
  result = readRDS(sprintf("output/C32_LM_humann3_sigTaxa16_Marginal_path-%s-%s-ZOE2.rds", DR.test, pheno.test))
  result_t3c = tibble(pathway =  result$coef %>% rownames, 
                      coef.ind = result$coef$phenotype, 
                      se.ind = result$se$phenotype,
                      p.ind = result$pval$phenotype,
                      q.ind = result$pval.adj$phenotype,
                      star.ind = ifelse(q.ind  <= 0.05, "*", "")) 
  pheno.test = "microt3cb"
  result = readRDS(sprintf("output/C32_LM_humann3_sigTaxa16_Marginal_path-%s-%s-ZOE2.rds", DR.test, pheno.test))
  result_mt3c = tibble(pathway =  result$coef %>% rownames, 
                       coef.loc = result$coef$phenotype, 
                       se.loc = result$se$phenotype,
                       p.loc = result$pval$phenotype,
                       q.loc = result$pval.adj$phenotype,
                       star.loc = ifelse(q.loc  <= 0.05, "*", "")) 
  tab3 =
    left_join(result_t3c, result_mt3c) %>% 
    left_join(tab2c, .) %>% 
    left_join(pth.names) %>% 
    transmute(pathway, description = pathname, p = p.loc, beta = coef.loc,  se = se.loc,
              ci_lb = beta - se * qnorm(0.975),  ci_ub = beta + se * qnorm(0.975),
              `expression rank` = rank.overall, 
              `% sig species` = X..key.species, `sig species` = key.species) %>% 
    arrange(p)
  write.csv(tab3, sprintf("output/ETable2_C32_composition_in_path_humann3_RNA_sigTaxa16_Marginal.csv"))
  
  tab_for_EFig2 =
    tab2b %>% 
    left_join(result_mt3c) %>% 
    left_join(pth.names) %>% 
    mutate(pathway = ifelse(!is.na(pathname), sprintf("%s (%s) ", pathway, pathname), pathway),
           `Expression level` = total,
           rank.overall = as.integer(rank.overall),
           `Log-normal model coefficient` = coef.loc,
           `% significant species` = `% key species`,
           `-log10 FDR-adjust p-values` = -log10(`q.loc`),
           pway = sprintf("%s\n [%s] ", pathway, `key species`)) 
  tab_for_EFig2 %>% 
    ggplot(aes(y = `Log-normal model coefficient`,
               x = `% significant species`,
               size = `-log10 FDR-adjust p-values`,
               alpha = `% significant species`,
               col = `Expression level`)) +
    scale_color_continuous(high = "#132B43", low = "#56B1F7") +
    scale_alpha_continuous(range = c(0.5, 1)) +
    geom_point() +
    guides(alpha = "none") +
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_text_repel(data = . %>% filter(q.loc <= 0.05),
                    aes(label = pway), alpha = 1, hjust = 1, direction = "y",
                    nudge_x = 35, nudge_y = -0.01, segment.size = 0.2, segment.alpha = 1,
                    size = 3) +
    geom_text_repel(data = . %>% filter(q.loc > 0.05),
                    aes(label = pathway), col = "grey", alpha = 1, hjust = 1, direction = "y",
                    nudge_x = -35, nudge_y = -0.01, segment.size = 0.2, segment.alpha = 0.5,
                    size = 2, max.overlaps = getOption("ggrepel.max.overlaps", default = 2)) +
    geom_text(data = . %>% filter(q.loc <= 0.05), 
              aes(x =`% significant species` + 0.5,
                  y = `Log-normal model coefficient` + 0.003,
                  label = sprintf("#%d", rank.overall)),
              size = 3, col = "black")
  ggsave(sprintf("figure/EFig2_C32_composition_in_path_humann3_RNA_sigTaxa16_Marginal.pdf"),
         width = 10, height = 6)
  
  library(openxlsx)
  write.xlsx(tab_for_EFig2, "figure/raw_data_for_EFig2.xlsx", sheetName = "EFig2")

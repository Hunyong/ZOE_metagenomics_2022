### C41 RNAseq_four_species.R
### Author: Hunyong Cho
### Identify the pathobiont genes for the four species.

### 0. library and preprocessing
  library(dplyr)
  source("scripts/F.generic.R")
  dat = readRDS("Data-processed/data.geneRPK.full.DRNA.ZOE2.RNASEQ.rds")
  reads.gene = dat$reads[,, 2] %>% apply(1, mean)
  gene.quality = rownames(dat$reads)[reads.gene > 20] # 542 gene-species
  key = c("mutans", "sputigena", "salivae", "wadei")
  gene.quality.by.species = 
    lapply(key, function(x) grep(x, gene.quality, val = T))
  names(gene.quality.by.species) = key
  
  dat = dat %>% pick.rows(row.index = gene.quality)
  
  ### Mutans data (p = 47)
  rna.mutans = t(dat$otu[gene.quality.by.species$mutans,,2])
  geneSm = colnames(rna.mutans) = gsub("(.*) g__.*mutans$" , "Sm_\\1", colnames(rna.mutans))
  
  ### Mutans data (p = 39)
  rna.sputigena = t(dat$otu[gene.quality.by.species$sputigena,,2])
  geneSs = colnames(rna.sputigena) = gsub("(.*) g__.*sputigena$" , "Ss_\\1", colnames(rna.sputigena))
  
  
  ### Salivae data (p = 8)
  rna.salivae = t(dat$otu[gene.quality.by.species$salivae,,2])
  genePs = colnames(rna.salivae) = gsub("(.*) g__.*salivae$" , "Ps_\\1", colnames(rna.salivae))
  
  ### Wadei data (p = 448)
  rna.wadei = t(dat$otu[gene.quality.by.species$wadei,,2])
  geneLw = colnames(rna.wadei) = gsub("(.*) g__.*wadei$" , "Lw_\\1", colnames(rna.wadei))
  
  dat.tmp = dat$meta
  race <- factor(dat$meta$race2)
  levels(race) <- c("1black", "3other", "2white", "3other", "3other")
  race <- factor(race, levels = c("1black", "2white", "3other"))



### 1. Main effects

  for (k in 1:4) {
    genes = list(geneSm, geneSs, genePs, geneLw)[[k]]
    rna.dat = list(rna.mutans, rna.sputigena, rna.salivae, rna.wadei)[[k]]
    nm  = c("Sm", "Ss", "Ps", "Lw")[k]
    nm.full = c("S. Mutans", "S. Sputigena", "P. Salivae", "L. Wadei")[k]
    
    outMain = data.frame(gene = NA, 
                         coef = rep(NA, length(genes)), 
                         pval = NA, qval = NA)
    
    cat(nm, "\n")
    i = 1
    for (g1 in genes) {
      cat(g1, " \n")
      dat.tmp$rna = rna.dat[, g1]
      out = lm(micro_t3c_b ~ agemo + race2 + rna + batch.RNA, data = dat.tmp) %>% summary %>% coef
      outMain[i, c("gene")] = g1
      outMain[i, c("coef", "pval")] = out["rna", c("Estimate", "Pr(>|t|)")] %>% print
      i = i + 1
    }
    outMain[, "qval"] = p.adjust(outMain[, "pval"], method = "BH")
    outMain %>% filter(qval < 0.2)
    # Extended Data Table 3
    outMain %>% arrange(qval) %>% write.csv(sprintf("output/C41_RNASeq_%s_lm_main.csv", nm))
    sig.genes = outMain %>% filter(qval <= 0.05) %>% "$"("gene")
    
  }



### 2. Interaction effects


  outInteraction = array(NA, dim = c(length(geneSm), length(geneSs), 3),
                         dimnames = list(geneSm, geneSs, c("coef", "pval", "qval")))
  outInteraction2 = data.frame(geneSm = NA,geneSs = NA, 
                               coef = rep(NA, length(geneSm) * length(geneSs)), coef.Sm = NA, coef.Ss = NA,
                               pval = NA, pval.Sm = NA, pval.Ss = NA, 
                               qval = NA, qval.Sm = NA, qval.Ss = NA,
                               sig.interaction = "")
  i = 1
  for (g1 in geneSm) {
    for (g2 in geneSs) {
      cat(g1, " ", g2, " \n")
      dat.tmp$Sm = rna.mutans[, g1]
      dat.tmp$Ss = rna.sputigena[, g2]
      out = lm(micro_t3c_b ~ agemo + race2 + Sm * Ss + batch.RNA, data = dat.tmp) %>% summary %>% coef
      outInteraction2[i, c("geneSm", "geneSs")] = c(g1, g2)
      outInteraction2[i, c("coef", "pval")] = 
        outInteraction[g1, g2, 1:2] = 
        out["Sm:Ss", c("Estimate", "Pr(>|t|)")] %>% print
      outInteraction2[i, c("coef.Sm", "pval.Sm")] = 
        out["Sm", c("Estimate", "Pr(>|t|)")] %>% print
      outInteraction2[i, c("coef.Ss", "pval.Ss")] = 
        out["Ss", c("Estimate", "Pr(>|t|)")] %>% print
      i = i + 1
    }
  }
  outInteraction2[, "qval"] = p.adjust(outInteraction2[, "pval"], method = "BH")
  outInteraction2[, "qval.Sm"] = p.adjust(outInteraction2[, "pval.Sm"], method = "BH")
  outInteraction2[, "qval.Ss"] = p.adjust(outInteraction2[, "pval.Ss"], method = "BH")
  
  outInteraction2[outInteraction2$qval <= 0.05, "sig.interaction"] = "sig"
  outInteraction2 = outInteraction2 %>% arrange(qval) 
  # Extended Data Table 4
  outInteraction2 %>% write.csv("output/C41_RNASeq_four_sepcies_lm_interaction.csv")

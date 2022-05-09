### C23_test_validation.R
### Author: Hunyong Cho
### Validating with the ZOE 2.0 pilot data.


### 0. library
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  
  source("scripts/F.generic.R")
  source("scripts/F.aesthetics.R")
  res = function(DR = "DNA", pheno = "t3c", ZOE = 2, cutoff = 0.05) {
    result = 
      sprintf("output/C2_LM_bracken_cov-bact-%s-%s-ZOE%d.rds", DR, pheno, ZOE) %>% 
      readRDS()
    sig = result$pval.adj$pheno <= cutoff
    data.frame(species = rownames(result$pval.adj) [sig], 
               coef = result$coef[sig, "phenotype"], 
               pval = result$pval[sig, "phenotype"],
               pval.adj = result$pval.adj[sig, "phenotype"])
  }
  res.validation = function(DR = "DNA", pheno = "t3c", cutoff = 0.05, digits.coef = 2) {
    # result = 
    sprintf("output/C2_validation-LM_cov-bact-%s-%s-validation_tab.rds", DR, pheno) %>% 
      readRDS() %>% 
      transmute(species,
                valid = ifelse(zoe1.q <= 0.05, "FDR", 
                               ifelse(zoe1.p <= 0.05, "sig", 
                                      ifelse(zoe1 * zoe2 > 0 , "SD", ""))),
                valid = ifelse(zoe1 * zoe2 < 0, "Opposite", valid),
                coef = sprintf("%s (%s, %s)", round(zoe2, digits.coef), 
                               format(zoe2.p, digits = 2, scientific = T), 
                               format(zoe2.q, digits = 2, scientific = T)),
                coef = gsub("e-0", "x10-", coef))
  }
  
  

##### 1. Key 23 species discovered by the main data and 16 validated sig (at least one of the four) in the replication sample
  type1 = "bact"
  model.nm = "cov"
  bracken = TRUE
  
  int.set = Reduce(intersect, list(res("DNA", "t3c")$species,
                                   res("DNA", "microt3cb")$species,
                                   res("RNA", "t3c")$species,
                                   res("RNA", "microt3cb")$species)) %>% sort
  key.bact23 = int.set %>% gsub("([^ ]*) (.*)", "g__\\1.s__\\1_\\2", .) %>% gsub("sp\\. oral", "sp_oral", .) %>% gsub(" ", "_", .)
  key.bact.m59 = res("DNA", "microt3cb")$species  %>% gsub("([^ ]*) (.*)", "g__\\1.s__\\1_\\2", .) %>% gsub("sp\\. oral", "sp_oral", .) %>% gsub(" ", "_", .)
  setdiff(key.bact.m59, key.bact23)
  ### maculosa (6) is not included any more
  ### longum (15) and periodontii (16) are newly included.

##### 2. validation
  ZOE = 1
  ll = 1
  coef.plot = list()
  for (DR.no in 1:2) {
    for (pheno in c("microt3cb", "t3c")) {
      plot.nm = paste0(c("DNA", "RNA")[DR.no], "_", pheno)
      cat(plot.nm, "\n")
      
      .initialize(test = "LN", screen = TRUE, prev.threshold = 0.1, 
                  avg.detect = 0.00001, detect.rel.abund = TRUE, nrm = TRUE,
                  humann = 2, bracken = TRUE, screen.among.tested.DNAs = TRUE)
      
      ## 1. reading the analysis results
      summary.fn <- paste0("output/C2_summary_bracken_cov-", type, "-", DRNA, "-ZOE2.rds")
      mainAnalysis <- readRDS(summary.fn)
      lm.fn <- paste0("output/C2_LM_bracken_cov-", type, "-", DRNA, "-", pheno, "-ZOE2.rds")
      output.nm.validation = paste0("output/C2_validation-LM_", model.nm, "-", type, "-", DRNA, "-", pheno, "-validation", ".rds")
      fig.coef.nm.validation = paste0("figure/C2_validation_LM_", model.nm, "-", type, "-", DRNA, "-", pheno, ".png")
      comb.table.fn = gsub("\\.rds", "_tab.rds", output.nm.validation)

      sig.taxa  = mainAnalysis$species
      dat.taxa.nm = gsub(".*\\.s__", "", dat$taxa[,1]) %>% gsub("_", " ", .) 
      # sig.index = gsub(".*\\.s__", "", dat$taxa[,1]) %>% gsub("_", " ", .) %>% {which(. %in% sig.taxa)}
      sig.index = sapply(sig.taxa, function(x) {tmp = which(dat.taxa.nm == x); if (length(tmp)) tmp else NA})
      n.sig = length(sig.taxa)
      
      
      
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
      result$sigFeature = sig.taxa[result$pval.adj$min < threshold.pval]
      
      saveRDS(result, output.nm.validation)
      
      
      lm.output <- readRDS(lm.fn)
      
      taxa.nm.lm = 
        gsub(".*\\.s__", "", lm.output$coef %>% rownames) %>% 
        gsub("_", " ", .)
      sig.index.lm = 
        sapply(sig.taxa %>% gsub(".*\\.s__", "", .) %>% gsub("_", " ", .) %>% gsub("sp oral", "sp. oral", .), 
               function(x) {a = which(taxa.nm.lm == x); ifelse(length(a)==0, NA, a)})
      
      rng.zoe2 = lm.output$coef[sig.index.lm, "phenotype"] %>% range(na.rm = TRUE)
      xlim1 <- c(rng.zoe2[1] - (rng.zoe2[2] - rng.zoe2[1]) * 0.2,
                 rng.zoe2[2] + (rng.zoe2[2] - rng.zoe2[1]) * 0.9)
      
      key.bact4 =
        c("g__Streptococcus.s__Streptococcus_mutans", 
          "g__Selenomonas.s__Selenomonas_sputigena",
          "g__Leptotrichia.s__Leptotrichia_wadei",
          "g__Prevotella.s__Prevotella_salivae")
      key.bact16 =   # From C05.09 summary_for_validation.R significant species for at least one level in the replication sample.
        c("g__Prevotella.s__Prevotella_salivae",
          "g__Prevotella.s__Prevotella_oulorum",
          "g__Prevotella.s__Prevotella_melaninogenica", 
          "g__Prevotella.s__Prevotella_sp_oral_taxon_306",
          "g__Prevotella.s__Prevotella_veroralis",
          "g__Stomatobaculum.s__Stomatobaculum_longum",
          "g__Streptococcus.s__Streptococcus_mutans",
          "g__Selenomonas.s__Selenomonas_sputigena",
          "g__Leptotrichia.s__Leptotrichia_wadei",
          "g__Leptotrichia.s__Leptotrichia_sp_oral_taxon_847",
          "g__Leptotrichia.s__Leptotrichia_sp_oral_taxon_223",
          "g__Leptotrichia.s__Leptotrichia_sp_oral_taxon_498",
          "g__Lachnoanaerobaculum.s__Lachnoanaerobaculum_saburreum",
          "g__Lachnospiraceae.s__Lachnospiraceae_bacterium_oral_taxon_082",
          "g__Veillonella.s__Veillonella_atypica",
          "g__Centipeda.s__Centipeda_periodontii" 
        )
      
      combined = 
        tibble(species = dat$taxa[sig.index, 1] %>% gsub("sp\\.", "sp", .),
               genus = factor(gsub(" .*", "", dat$taxa[sig.index, 1])),
               zoe2 = lm.output$coef[sig.index.lm, "phenotype"],
               zoe2.p = lm.output$pval[sig.index.lm, "phenotype"],
               zoe2.q = lm.output$pval.adj[sig.index.lm, "phenotype"],
               zoe1 = result$coef[, "phenotype"],
               zoe1.p = result$pval[, "phenotype"],
               zoe1.q = result$pval.adj[, "phenotype"]) %>%
        mutate(
          # q-value based label
          note = ifelse(zoe1.q <= 0.10, "~'\u2020 '", "~'   '"),
          note = ifelse(zoe1.q <= 0.05, "~'* '", note),
          # p-value based label
          note.p = ifelse(zoe1.p <= 0.10, "~'*  '", "~'    '"),
          note.p = ifelse(zoe1.p <= 0.05, "~'** '", note.p),
          note.p = ifelse(zoe1.p <= 0.01, "~'***'", note.p),
          key4 = species %in% {key.bact4 %>% gsub(".*s__", "", .) %>% gsub("_", " ", .)},
          # p-value based label
          label = sprintf("%s%s%s~%s", 
                          ifelse(key4, "underline(", ""), gsub(" ", "~", species), ifelse(key4, ")", ""),
                          note.p))
      
      saveRDS(combined, comb.table.fn)
      
      coef.plot[[ll]] = 
        combined %>% 
        # arrange(zoe1) %>% 
        # mutate(species = factor(species)) %>% 
        ggplot(aes(zoe2, zoe1, col = genus, size = -log(zoe1.q))) + geom_point() + 
        geom_text_repel(data = combined %>% filter(species %in% {key.bact16 %>% gsub(".*s__", "", .) %>% gsub("_", " ", .)}),
                        aes(label = label), hjust = 1, direction = "y", parse = TRUE,
                        nudge_x = 10, segment.size = 0.2, segment.alpha = 0.5, #segment.colour = "navyblue",
                        size = 3) +
        xlim(xlim1) + guides(col = "none", size = guide_legend(title= "-log (FDR-adjusted p-value) in the replication sample")) +
        geom_abline(slope = 1, intercept = 0, col = "gray50") + 
        scale_size_continuous(range = size.log.p) +
        scale_color_manual(values = col.genera) +
        # scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8), labels = c(-2, 0, 2, 4, "", "")) +
        xlab("Discovery sample log-normal model coefficient") + ylab("Replication sample log-normal model coefficient") +
        theme_bw() + theme(legend.position = "bottom")
      names(coef.plot)[ll] = plot.nm
      
      ggsave(fig.coef.nm.validation, plot = coef.plot[[ll]], width = 6, height = 6)
      ll = ll + 1
    }
  }

### 1B. Supp Table B.
int.set.validation = list(mD = res.validation("DNA", "microt3cb", digits.coef = 2) %>% 
                            rename(mD = valid, mD.c = coef),
                          mR = res.validation("RNA", "microt3cb", digits.coef = 2) %>% 
                            rename(mR = valid, mR.c = coef),
                          tD = res.validation("DNA", "t3c", digits.coef = 3) %>% 
                            rename(tD = valid, tD.c = coef),
                          tR = res.validation("RNA", "t3c", digits.coef = 3) %>% 
                            rename(tR = valid, tR.c = coef))
comb.tab.all4 = 
  Reduce(left_join, int.set.validation) %>% 
  mutate(FDR = (mD == "FDR") + (mR == "FDR") + (tD == "FDR") + (tR == "FDR"),
         sig = (mD %in% c("FDR", "sig")) + (mR  %in% c("FDR", "sig")) + 
           (tD  %in% c("FDR", "sig")) + (tR  %in% c("FDR", "sig")),
         SD = (mD != "Opposite") + (mR != "Opposite") + (tD != "Opposite") + (tR != "Opposite")) %>% 
  left_join(sprintf("output/C2_validation-LM_cov-bact-%s-%s-validation_tab.rds", "DNA", "microt3cb") %>% 
              readRDS() %>% transmute(species, p = zoe2.p), 
            by = "species") %>% 
  arrange(desc(FDR), desc(sig), desc(SD), desc(p))
comb.tab.all4 %>% filter(sig > 0) %>% dplyr::select(species, ends_with(".c"), everything()) %>% View
comb.tab.all4 %>% filter(sig > 0) %>% write.csv("output/_C2_coef_table_bracken_combined-cov_(DnR)_ZOE2.csv")
comb.tab.all4 %>% filter(sig > 0) %>% "$"("species")


### 2. Combined plot: t3c / micro_t3c_b x DNA / RNA   for only 16 key species.     Oct 11, 2021
### Supp Fig B.

coef.plot.legend = get_legend(coef.plot[[1]] + theme(legend.direction = "horizontal", legend.position = "bottom"))

coef.plot.label = 
  gsub("(.*NA)_(microt3cb)", "Localized caries experience, \\1", names(coef.plot)) %>% 
  gsub("(.*NA)_(t3c)", "Person-level caries experience, \\1", .) %>% 
  gsub(", DNA", ", metagenome", .) %>% 
  gsub(", RNA", ", metatranscriptome", .)

coef.plot.main = 
  plot_grid(plotlist = 
              lapply(1:4, function(x) coef.plot[[x]] + theme(legend.position="none", axis.title = element_blank()) +
                       # annotate("text", x = -0.05, y = 1, label = coef.plot.label[x]) + 
                       ggtitle(coef.plot.label[x])), 
            nrow = 2, ncol = 2)+
  draw_label("Discovery sample log-normal model coefficient", x=0.5, y=  0, vjust= 0.5, angle=  0) +
  draw_label("Replication sample log-normal model coefficient", x=0, y=  0.5, vjust= 0, angle= 90) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 

p.final = 
  plot_grid(coef.plot.main, coef.plot.legend, ncol = 1, rel_heights = c(1, 0.05)) 

##### Fig 2.  #####
save_plot(paste0("figure/_C2_validation_plot_2x2_bracken_combined-cov_ZOE2_validation.png"), 
          p.final, base_height = 12, base_width = 12)





### 3. DNA-RNA plot for only 16 key species.     Oct 11, 2021
### Extended Data Fig. 1
  ## 2. lm coefficents D vs R
  pheno = "microt3cb"
  ZOE = 2
  model1 = paste0("~ batcheffect + phenotype", binaryDNA, " + agemo + race")
  
  ll = 1
  DR.plot = list()
  for (ZOE in 2:1) {
    for (pheno in c("microt3cb", "t3c")) {
      plot.nm = paste0("ZOE", ZOE, "_", pheno)
      cat(plot.nm, "\n")
      
      pheno.name <-
        c(t3c = "Person-level caries experience", microt3cb = "Localized caries experience")[pheno]
      ZOE.name <- c("1" = "replication sample", "2" = "discovery sample")[ZOE]
      
      fn2 = paste0("output/C2_LM_bracken_cov-bact-", "DNA", "-", pheno, "-ZOE", ZOE, ".rds")
      res.lm <- readRDS(fn2)
      res.lm2 <- readRDS(fn2 %>% gsub("DNA", "RNA", .))
      tib.D <- 
        tibble(taxa = res.lm$coef %>% rownames(),
               lm.coef    = res.lm$coef$phenotype,
               lm.p    = res.lm$pval$phenotype,
               lm.q    = res.lm$pval.adj$phenotype,
               sig = paste0(ifelse(!is.na(lm.q) & lm.q <= 0.05, "sig", "")),
               direction = paste0(ifelse(!is.na(lm.q) & lm.q <= 0.05, sign2(lm.coef), " "))) %>% 
        mutate(genus = gsub("g\\_\\_(.*)\\.s.*", "\\1", taxa),
               species = gsub(".*\\.s\\_\\_(.*)$", "\\1", taxa),
               species = gsub("\\_", " ", species)) %>% 
        dplyr::filter(taxa %>% gsub("sp\\.", "sp", .) %in% {key.bact16 %>% gsub(".*s__", "", .) %>% gsub("_", " ", .)})
      tib.R <- 
        tibble(taxa = res.lm2$coef %>% rownames(),
               lm.coef    = res.lm2$coef$phenotype,
               lm.p    = res.lm2$pval$phenotype,
               lm.q    = res.lm2$pval.adj$phenotype,
               sig = paste0(ifelse(!is.na(lm.q) & lm.q <= 0.05, "sig", "")),
               direction = paste0(ifelse(!is.na(lm.q) & lm.q <= 0.05, sign2(lm.coef), " "))) %>% 
        mutate(genus = gsub("g\\_\\_(.*)\\.s.*", "\\1", taxa),
               species = gsub(".*\\.s\\_\\_(.*)$", "\\1", taxa),
               species = gsub("\\_", " ", species)) %>% 
        dplyr::filter(taxa %>% gsub("sp\\.", "sp", .) %in% {key.bact16 %>% gsub(".*s__", "", .) %>% gsub("_", " ", .)})
      tib.DR <- 
        full_join(tib.D, tib.R, by = c("taxa", "genus", "species"), suffix = c(".D", ".R")) %>% 
        mutate(sig = paste0(ifelse(is.na(sig.D) | sig.D == "", "", "D"),
                            ifelse(is.na(sig.R) | sig.R == "", "", "R")),
               sig = ifelse(sig == "", " ", sig)) %>% 
        mutate(
          # q-value based label
          note = ifelse(lm.q.D <= 0.10, "~'\u2020 '", "~'   '"),
          note = ifelse(lm.q.D <= 0.05, "~'* '", note),
          # # p-value based label
          # note.p = ifelse(lm.p.D <= 0.05, "~'\u2020 '", "~'   '"),
          # note.p = ifelse(lm.p.D <= 0.05, "~'* '", note.p),
          # p-value based label
          note.p = ifelse(lm.p.D <= 0.10, "~'*  '", "~'   '"),
          note.p = ifelse(lm.p.D <= 0.05, "~'** '", note.p),
          note.p = ifelse(lm.p.D <= 0.01, "~'***'", note.p),
          # note.p = ifelse(ZOE == 2, "~'   '", note.p), # if ZOE==2, no annotation.
          note.pR = ifelse(lm.p.R <= 0.10, "~'*  '", "~'   '"),
          note.pR = ifelse(lm.p.R <= 0.05, "~'** '", note.pR),
          note.pR = ifelse(lm.p.R <= 0.01, "~'***'", note.pR),
          # note.pR = ifelse(ZOE == 2, "~'   '", note.pR), # if ZOE==2, no annotation.
          key4 = species %in% {key.bact4 %>% gsub(".*s__", "", .) %>% gsub("_", " ", .)},
          label = sprintf("%s%s%s~%s/%s",
                          ifelse(key4, "underline(", ""), gsub(" ", "~", species), ifelse(key4, ")", ""),
                          note.p, note.pR)) %>% 
        mutate(genus = factor(gsub(" .*", "", species)))
      y.rng = range(tib.DR$lm.coef.D)
      
      DR.plot[[ll]] =
        tib.DR %>%  
        ggplot(aes(lm.coef.D, lm.coef.R, col = genus, size = -log(lm.q.D)), alpha = 0.5) +
        geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
        geom_abline(slope = 1, intercept = 0, col = "gray50") + 
        scale_x_continuous(limits = c(NA, y.rng[2] + min(y.rng[2], y.rng[2] - y.rng[1]))) +
        scale_size_continuous(range = size.log.p) +
        scale_color_manual(values = col.genera) +
        geom_point() + 
        guides(color = "none", size = guide_legend(title= "-log (FDR-adjusted p-value) in the metagenome sample")) +
        theme_bw() + theme(legend.position = "bottom") +
        # scale_color_manual(values = c("D" = "red", "R" = "royalblue4", "DR" = "slateblue", " " = "gray")) +
        xlab("log-normal model coefficient (metagenome)") +
        ylab("log-normal model coefficient (metatranscriptome)") + 
        # ggtitle(paste0("Estimated log-normal coefficients for species (DNA vs RNA)\n", pheno.name," ", model1)) +
        ggtitle(sprintf("%s in %s", pheno.name, ZOE.name)) +
        geom_text_repel(data = tib.DR, 
                        aes(label = label, col = genus), parse = TRUE,
                        # aes(label = species %>% paste0(" ")), 
                        hjust = 1, direction = "y",
                        nudge_x = 2, segment.size = 0.2, segment.alpha = 0.5, #segment.colour = "navyblue",
                        size = 3)
      names(DR.plot)[ll] = plot.nm
      ll = ll + 1
    }
  }
  DR.plot.legend = get_legend(DR.plot[[1]] + theme(legend.direction = "horizontal", legend.position = "bottom"))
  
  DR.plot.label = 
    gsub("(ZOE[1-2])_(microt3cb)", "Localized caries experience, \\1", names(DR.plot)) %>% 
    gsub("(ZOE[1-2])_(t3c)", "Person-level caries experience, \\1", .) %>% 
    gsub(", ZOE1", " in replication sample", .) %>% 
    gsub(", ZOE2", " in discovery sample", .)
  
  DR.plot.main = 
    plot_grid(plotlist = 
                lapply(1:4, function(x) DR.plot[[x]] + theme(legend.position="none", axis.title = element_blank()) +
                         ggtitle(DR.plot.label[x])), 
              nrow = 2, ncol = 2)+
    draw_label("Log-normal model coefficient in the metagenome sample", x=0.5, y=  0, vjust= 0.5, angle=  0) +
    draw_label("Log-normal model coefficient in the metatranscriptome sample", x=0, y=  0.5, vjust= 0, angle= 90) +
    theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
  
  p.final = 
    plot_grid(DR.plot.main, DR.plot.legend, ncol = 1, rel_heights = c(1, 0.05)) 
  
  ##### Extended Data Fig. 1 ##### 
  save_plot(paste0("figure/_C2_DR_plot_2x2_bracken_", "combined-", model.nm, ".png"), 
            p.final, base_height = 12, base_width = 12)



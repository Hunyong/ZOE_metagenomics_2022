### C23_reporting.R
### Author: Hunyong Cho
### Tables and Figures for the 16 and four species.
### output: Ext. Dat. Table 1, Fig. 2, Ext. Dat. Fig. 1


### 0. library
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  
  source("scripts/F_generic.R")
  source("scripts/F_aesthetics.R")
  res.validation = function(DR = "DNA", pheno = "t3c", cutoff = 0.05, digits.coef = 2) {
    # result = 
    sprintf("output/C2_LM_bracken_%s-%s-combined.rds", DR, pheno) %>% 
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


### Reporting
### 1 Extended Data Table 1
  int.set.validation = list(mD = res.validation("DNA", "microt3cb", digits.coef = 2) %>% 
                              rename(mD.val = valid, mD.c = coef),
                            mR = res.validation("RNA", "microt3cb", digits.coef = 2) %>% 
                              rename(mR.val = valid, mR.c = coef),
                            tD = res.validation("DNA", "t3c", digits.coef = 3) %>% 
                              rename(tD.val = valid, tD.c = coef),
                            tR = res.validation("RNA", "t3c", digits.coef = 3) %>% 
                              rename(tR.val = valid, tR.c = coef))
  comb.tab.all4 = 
    Reduce(left_join, int.set.validation) %>% 
    mutate(FDR = (mD.val == "FDR") + (mR.val == "FDR") + (tD.val == "FDR") + (tR.val == "FDR"),
           sig = (mD.val %in% c("FDR", "sig")) + (mR.val  %in% c("FDR", "sig")) + 
             (tD.val  %in% c("FDR", "sig")) + (tR.val  %in% c("FDR", "sig")),
           SD = (mD.val != "Opposite") + (mR.val != "Opposite") + (tD.val != "Opposite") + (tR.val != "Opposite")) %>% 
    left_join(sprintf("output/C2_LM_bracken_%s-%s-combined.rds", "DNA", "microt3cb") %>% 
                readRDS() %>% transmute(species, p = zoe2.p), 
              by = "species") %>% 
    mutate(rank.FDR = ifelse(FDR > 1, FDR * 10, FDR * 10 + sig)) %>% 
    arrange(desc(rank.FDR), p) %>% 
    dplyr::select(-rank.FDR)
  comb.tab.all4 %>% filter(sig > 0) %>% dplyr::select(species, ends_with(".c"), everything()) %>% View
  comb.tab.all4 %>% filter(sig > 0) %>% write.csv("output/ETable1_C2_coef_table_bracken.csv")
  
  
### 2. Figure 2: Panels of D/RNA & phenotypes
  
  ## 2.0 Some parameters for annotation
  # 16 species being annotated
  key.bact16 = comb.tab.all4 %>% filter(sig > 0) %>% "$"("species")
  # Four species being highlighted
  top.bact4 = c("Streptococcus mutans", "Selenomonas sputigena", 
                "Leptotrichia wadei", "Prevotella salivae")
  
  ## 2.1 Generating each panel over D/RNA & phenotypes
  counter = 1  # initializing the counter
  coef.plot = list()
  for (DR.no in 1:2) {
    DRNA = c("DNA", "RNA")[DR.no]
    for (pheno in c("microt3cb", "t3c")) {
      plot.nm = paste0(DRNA, "_", pheno)
      output.nm.combined = paste0("output/C2_LM_bracken_", DRNA, "-", pheno, "-combined.rds")
      fig.coef.nm.combined = gsub("output/(.*)\\.rds", "figure/\\1.pdf", output.nm.combined)
      
      combined2 = 
        readRDS(output.nm.combined) %>% 
        mutate(
          # q-value based label
          note = ifelse(zoe1.q <= 0.10, "~'\u2020 '", "~'   '"),
          note = ifelse(zoe1.q <= 0.05, "~'* '", note),
          # p-value based label
          note.p = ifelse(zoe1.p <= 0.10, "~'*  '", "~'    '"),
          note.p = ifelse(zoe1.p <= 0.05, "~'** '", note.p),
          note.p = ifelse(zoe1.p <= 0.01, "~'***'", note.p),
          top4 = species %in% top.bact4,
          species.short = gsub("^([a-zA-Z])([a-zA-Z]*) (.*)", "\\1. \\3", species),
          # p-value based label
          label = sprintf("italic(%s%s%s~(~%2.3f))~%s", 
                          ifelse(top4, "underline(", ""), 
                          gsub(" ", "~", species.short), 
                          ifelse(top4, ")", ""),
                          zoe1.p,
                          note.p),
          label = gsub("\\(~0\\.000\\)", "( '<'~0.001)", label),
          label = gsub("(\\(~0\\.\\d\\d)0\\)", "\\10)", label))
      
      # Graphical parameters
      rng.zoe2 = combined2$zoe2 %>% range(na.rm = TRUE)
      xlim1 <- c(0, #rng.zoe2[1] - (rng.zoe2[2] - rng.zoe2[1]) * 0.2,
                 rng.zoe2[2] + (rng.zoe2[2] - rng.zoe2[1]) * 0.9)
      
      coef.plot[[counter]] = 
        combined2 %>% 
        ggplot(aes(zoe2, zoe1, col = genus, size = -log(zoe1.q))) + geom_point() + 
        geom_text_repel(data = combined2 %>% filter(species %in% key.bact16),
                        aes(label = label), hjust = 1, direction = "y", parse = TRUE,
                        nudge_x = 10, segment.size = 0.2, segment.alpha = 0.5, #segment.colour = "navyblue",
                        size = 4) +
        xlim(xlim1) + guides(col = "none", size = guide_legend(title= "-log (FDR-adjusted p-value) in the replication sample")) +
        geom_abline(slope = 1, intercept = 0, col = "gray50") + 
        scale_size_continuous(range = size.log.p) +
        scale_color_manual(values = col.genera) +
        # scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8), labels = c(-2, 0, 2, 4, "", "")) +
        xlab("Discovery sample log-normal model coefficient") + ylab("Replication sample log-normal model coefficient") +
        theme_bw() + theme(legend.position = "bottom")
      names(coef.plot)[counter] = plot.nm
      
      ggsave(fig.coef.nm.combined, plot = coef.plot[[counter]], width = 6, height = 6)
      counter = counter + 1
    }
  }
  
  ## 2.2 Combining across: DNA / RNA and phenotypes
  coef.plot.legend = get_legend(coef.plot[[1]] + theme(legend.direction = "horizontal", legend.position = "bottom"))
  
  coef.plot.label = 
    gsub("(.*NA)_(microt3cb)", "localized caries experience, \\1", names(coef.plot)) %>% 
    gsub("(.*NA)_(t3c)", "person-level caries experience, \\1", .) %>% 
    gsub("(.*), DNA", "MTG, \\1", .) %>% 
    gsub("(.*), RNA", "MTX, \\1", .)
  
  coef.plot.main = 
    plot_grid(plotlist = 
                lapply(1:4, function(x) coef.plot[[x]] + theme(legend.position="none", axis.title = element_blank())),# +
                         # annotate("text", x = -0.05, y = 1, label = coef.plot.label[x]) + 
                         # ggtitle(coef.plot.label[x])), 
              nrow = 2, ncol = 2)+
    draw_label("Discovery sample log-normal model coefficient", x=0.5, y=  0, vjust= 0.5, angle=  0) +
    draw_label("Replication sample log-normal model coefficient", x=0, y=  0.5, vjust= 0, angle= 90) +
    theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
  
  p.final = 
    plot_grid(coef.plot.main, coef.plot.legend, ncol = 1, rel_heights = c(1, 0.05)) 
  
  save_plot(paste0("figure/Fig2_C2_validation_DRNA_phenotype.pdf"), 
            p.final, base_height = 12, base_width = 12)
  
  
  
### 3. Ext. Data Fig 1. DNA-RNA plot for only 16 key species.
  counter = 1
  DR.plot = list()
  for (ZOE in 2:1) {
    for (pheno in c("microt3cb", "t3c")) {
      plot.nm = paste0("ZOE", ZOE, "_", pheno)
      cat(plot.nm, "\n")
      
      pheno.name <-
        c(t3c = "Person-level caries experience", microt3cb = "Localized caries experience")[pheno]
      ZOE.name <- c("1" = "replication sample", "2" = "discovery sample")[ZOE]
      
      fnD = sprintf("output/C2_LM_bracken_DNA-%s-%s.rds", pheno, ifelse(ZOE == 2, "ZOE2", "validation"))
      fnR = gsub("DNA", "RNA", fnD)
      result.D <- readRDS(fnD)
      result.R <- readRDS(fnR)
      tib.D <- 
        tibble(taxa    = result.D$coef %>% rownames(),
               lm.coef = result.D$coef$phenotype,
               lm.p    = result.D$pval$phenotype,
               lm.q    = result.D$pval.adj$phenotype,
               sig = paste0(ifelse(!is.na(lm.q) & lm.q <= 0.05, "sig", "")),
               direction = paste0(ifelse(!is.na(lm.q) & lm.q <= 0.05, sign2(lm.coef), " "))) %>% 
        mutate(genus = gsub("g\\_\\_(.*)\\.s.*", "\\1", taxa),
               species = gsub(".*\\.s\\_\\_(.*)$", "\\1", taxa),
               species = gsub("\\_", " ", species)) %>% 
        dplyr::filter(taxa %in% key.bact16)
      tib.R <- 
        tibble(taxa = result.R$coef %>% rownames(),
               lm.coef    = result.R$coef$phenotype,
               lm.p    = result.R$pval$phenotype,
               lm.q    = result.R$pval.adj$phenotype,
               sig = paste0(ifelse(!is.na(lm.q) & lm.q <= 0.05, "sig", "")),
               direction = paste0(ifelse(!is.na(lm.q) & lm.q <= 0.05, sign2(lm.coef), " "))) %>% 
        mutate(genus = gsub("g\\_\\_(.*)\\.s.*", "\\1", taxa),
               species = gsub(".*\\.s\\_\\_(.*)$", "\\1", taxa),
               species = gsub("\\_", " ", species)) %>% 
        dplyr::filter(taxa %in% key.bact16)
      tib.DR <- 
        full_join(tib.D, tib.R, by = c("taxa", "genus", "species"), suffix = c(".D", ".R")) %>% 
        mutate(sig = paste0(ifelse(is.na(sig.D) | sig.D == "", "", "D"),
                            ifelse(is.na(sig.R) | sig.R == "", "", "R")),
               sig = ifelse(sig == "", " ", sig)) %>% 
        mutate(
          # q-value based label
          note = ifelse(lm.q.D <= 0.10, "~'\u2020 '", "~'   '"),
          note = ifelse(lm.q.D <= 0.05, "~'* '", note),
          
          # p-value based label
          note.p = ifelse(lm.p.D <= 0.10, "~'*  '", "~'   '"),
          note.p = ifelse(lm.p.D <= 0.05, "~'** '", note.p),
          note.p = ifelse(lm.p.D <= 0.01, "~'***'", note.p),
          
          note.pR = ifelse(lm.p.R <= 0.10, "~'*  '", "~'   '"),
          note.pR = ifelse(lm.p.R <= 0.05, "~'** '", note.pR),
          note.pR = ifelse(lm.p.R <= 0.01, "~'***'", note.pR),
          species.short = gsub("^([a-zA-Z])([a-zA-Z]*) (.*)", "\\1. \\3", species),
          
          top4 = species %in% {top.bact4 %>% gsub(".*s__", "", .) %>% gsub("_", " ", .)},
          label = sprintf("italic(%s%s%s)~(~%2.3f~','~%2.3f )~%s/~%s",
                          ifelse(top4, "underline(", ""), 
                          gsub(" ", "~", species.short), 
                          ifelse(top4, ")", ""),
                          lm.p.D, lm.p.R, note.p, note.pR),
          label = gsub("0\\.000", " '<'~0.001", label)
          ) %>% 
        mutate(genus = factor(gsub(" .*", "", species)))
      y.rng = range(tib.DR$lm.coef.D)
      
      DR.plot[[counter]] =
        tib.DR %>%  
        ggplot(aes(lm.coef.D, lm.coef.R, col = genus, size = -log(lm.q.D)), alpha = 0.5) +
        geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
        geom_abline(slope = 1, intercept = 0, col = "gray50") + 
        scale_x_continuous(limits = c(NA, y.rng[2] + min(y.rng[2], y.rng[2] - y.rng[1]))) +
        scale_size_continuous(range = size.log.p) +
        scale_color_manual(values = col.genera) +
        geom_point() + 
        guides(color = "none", size = guide_legend(title= "-log (FDR-adjusted p-value) in MTG")) +
        theme_bw() + theme(legend.position = "bottom") +
        xlab("log-normal model coefficient (metagenome)") +
        ylab("log-normal model coefficient (metatranscriptome)") + 
        ggtitle(sprintf("%s in %s", pheno.name, ZOE.name)) +
        geom_text_repel(data = tib.DR, 
                        aes(label = label, col = genus), parse = TRUE,
                        # aes(label = species %>% paste0(" ")), 
                        hjust = 1, direction = "y",
                        nudge_x = 2, segment.size = 0.5, segment.alpha = 0.5, #segment.colour = "navyblue",
                        size = 4)
      DR.plot[[counter]]
      names(DR.plot)[counter] = plot.nm
      counter = counter + 1
    }
  }
  DR.plot.legend = 
    get_legend(DR.plot[[1]] + 
                 theme(legend.direction = "horizontal", legend.position = "bottom"))
  
  DR.plot.label = 
    gsub("(ZOE[1-2])_(microt3cb)", "Localized caries experience, \\1", names(DR.plot)) %>% 
    gsub("(ZOE[1-2])_(t3c)", "Person-level caries experience, \\1", .) %>% 
    gsub(", ZOE1", " in the replication sample", .) %>% 
    gsub(", ZOE2", " in the discovery sample", .)
  
  DR.plot.main = 
    plot_grid(plotlist = 
                lapply(1:4, function(x) 
                  DR.plot[[x]] + 
                    theme(legend.position="none", axis.title = element_blank()) +
                    ggtitle("")),
                    #ggtitle(DR.plot.label[x])), 
              nrow = 2, ncol = 2)+
    draw_label("Log-normal model coefficient in MTG", 
               x=0.5, y=  0, vjust= 0.5, angle=  0) +
    draw_label("Log-normal model coefficient in MTX", 
               x=0, y=  0.5, vjust= 0, angle= 90) +
    theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 
  
  p.final = 
    plot_grid(DR.plot.main, DR.plot.legend, ncol = 1, rel_heights = c(1, 0.05)) 
  
  save_plot(paste0("figure/EFig1_C2_validation_study_phenotype.pdf"), 
            p.final, base_height = 18, base_width = 18)
  
  

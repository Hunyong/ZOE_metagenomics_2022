### 0.1 library
  source("scripts/F.generic.R")  # typeDRNA(), sample.exclude, taxa.exclude
  library(ggplot2); library(ggrepel); library(dplyr); library(cowplot)


### 0.2 parameters
  ZOE = 1  # ZOE == 2: ZOE 2.0 main study, ZOE == 1: ZOE 2.0 pilot study
  type1 = type = "bact"
  humann = 2; bracken = TRUE
  nrm = TRUE #normalization
  threshold.pval = 0.10

# Figure 1 for DR.no == 1, pheno = "microt3cb", ZOE = 1
for (DR.no in 2:1) {
  DRNA = c("DNA", "RNA")[DR.no]
  for (ZOE in 1:2) {
  summary.fn <- paste0("output/C2_summary_bracken_cov-", type, "-", DRNA, "-ZOE", ZOE, ".rds")
  tmp.sig <- list()
  
  for (pheno in c("t3c", "microt3cb")) {
    
    binaryDNA = ifelse(DR.no == 2, "+ binaryDNA", "")
    model1 = paste0("~ batcheffect + phenotype", binaryDNA, " + agemo + race")
    pheno.name = switch(pheno, t3c = "t3c", microt3cb = "micro t3c b")
    
  # reading
    fn = paste0("output/C2_LM_bracken_cov-", type, "-", DRNA, "-", pheno, "-ZOE", ZOE, ".rds")
    res.lm <- readRDS(fn)
  
    res <- 
      tibble(taxa = res.lm$coef %>% rownames(),
             lm.coef    = res.lm$coef$phenotype,
             lm.p    = res.lm$pval$phenotype,
             lm.q    = res.lm$pval.adj$phenotype,
             sig = ifelse(!is.na(lm.q) & lm.q <= 0.05, "L", ""),
             sig.raw = ifelse(!is.na(lm.p) & lm.p <= 0.05, "L", ""),
             direction = ifelse(!is.na(lm.q) & lm.q <= 0.05, sign2(lm.coef), ""),
             direction.raw = ifelse(!is.na(lm.q), sign2(lm.coef), "")) %>% 
      mutate(sig = ifelse(sig=="", " ", sig),
             sig.raw = ifelse(sig.raw=="", " ", sig.raw),
             genus = gsub("g\\_\\_(.*)\\.s.*", "\\1", taxa),
             species = gsub(".*\\.s\\_\\_(.*)$", "\\1", taxa),
             species = gsub("\\_", " ", species))
    # res$sig %>% table
    tmp.sig[[pheno]] <- res[res$sig != " ", c("species", "sig", "direction")] %>% arrange(sig)
    tmp.sig[[pheno]] %>% print(n = 100)
    
    
    ## saving coefficients
    # Pulling Humann2 results for comparison
    if (ZOE==2) {
      fn.zoe1 = paste0("output/C2_coef_table_bracken_", pheno, "-cov_", DRNA, "_ZOE1.csv")
      if (!file.exists(fn.zoe1)) {
        cat(">>> The Humann2 or ZOE pilot results are not available. The comparison is not reflected in the tables.\n")
        compare.previous = FALSE
      } else {
        compare.previous = TRUE
        tmp.b = read.csv(fn.zoe1)
      }
    } else {
      compare.previous = FALSE
    }
    
    val = c("*** FDR", "**  raw pval", "*   direction", "     (No)")
    
    res.table <- 
      res %>% 
      arrange(lm.q) %>% 
      transmute(species = gsub("\\.", "", species),
                `log-normal` = paste0(signif(lm.coef, 2), " (", signif(lm.p, 2), ", ", signif(lm.q, 2), ")"),
                sig,
                direction,
                sig.raw,
                direction.raw,
                lm = lm.coef) %>% 
      arrange(desc(sig)) %>%
      {if (compare.previous) dplyr::select(., -sig.raw, -direction.raw) else .} %>% 
      {if (compare.previous) 
        full_join(., tmp.b %>% 
                    transmute(species, sig_ZOE1 = sig, sig_ZOE1_raw = sig.raw, direction_ZOE1_raw = direction.raw,
                              log.normal_ZOE1 = log.normal, lm_ZOE1 = lm)) 
        # else full_join(., tmp.a %>% transmute(species, sig_ZOE2 = sig))
        else .
        } %>% 
      {if (compare.previous) 
        mutate(., 
               sig_ZOE1 = ifelse(is.na(sig_ZOE1), "(Not tested)", sig_ZOE1),
               sig_ZOE1_raw = ifelse(is.na(sig_ZOE1_raw), "(Not tested)", sig_ZOE1_raw),
               direction_ZOE1_raw = ifelse(is.na(direction_ZOE1_raw), "(Not tested)", direction_ZOE1_raw)) %>% 
          mutate(validation = ifelse(sig == " ", " ",
                                     ifelse(sig_ZOE1 != " ", val[1], # FDR
                                            ifelse(sig_ZOE1_raw != " ", val[2], # raw pval
                                                   ifelse(direction == direction_ZOE1_raw, val[3], # direction
                                                          val[4])))), # none
                 validation = ifelse(sig_ZOE1 == "(Not tested)", "(Not tested)", validation),
                 order_ZOE2 = 1:n()) %>%  # ordering by first, the number of sig models, and second, the q-values (lm, disc, cont)
          arrange(validation != val[1], validation != val[2], validation != val[3], validation != val[4]) %>% 
          mutate(order_valid = 1:n()) %>% 
          arrange(order_ZOE2, validation != val[1], validation != val[2], validation != val[3], validation != val[4])
        # else mutate(., sig_ZOE2 = ifelse(is.na(sig_ZOE2), "(Not tested)", sig_ZOE2))
        else .
        }
    # res.table %>%  View
    ######### ingredients for Extended Data Table 1 #########
    write.csv(res.table, 
              paste0("output/_C2_coef_table_bracken_" , pheno, "-cov_", DRNA, "_ZOE", ZOE, ".csv"))
    if (compare.previous) {
      min.x = min(res.table$lm_ZOE1, na.rm = TRUE) * 1.1
      min.y = min(res.table$lm, na.rm = TRUE) * 1.1
      max.x = max(res.table$lm_ZOE1, na.rm = TRUE) * 1.7
      res.figure <-
        res.table %>% 
        filter(!is.na(validation)) %>% 
        transmute(species,
                  ZOE2 = `lm`,
                  pilot = `lm_ZOE1`,
                  validation = factor(validation, levels = c(val, "(Not tested)", " "), 
                                      labels = c(val, "(Not tested)", "(N/A)")),
                  sig = ifelse(!is.na(sig) & sig != " ", "FDR significant", "FDR insignificant"))
        # mutate(ZOE2 = ifelse(is.na(ZOE2), min.y, ZOE2),
        #        pilot = ifelse(is.na(pilot), min.x, pilot))
      plot.obj <-
        res.figure %>% 
        ggplot(aes(pilot, ZOE2, col = validation, shape = sig, alpha = sig)) + 
        geom_hline(yintercept = 0, col = "gray") +
        geom_vline(xintercept = 0, col = "gray") +
        guides(shape = guide_legend("significance (ZOE 2.0)"), alpha = FALSE) +
        scale_alpha_manual(values = c(0.3, 1)) + 
        xlab("ZOE 2.0 pilot") + ylab("ZOE 2.0") + xlim(c(min.x, max.x)) +
        geom_point() + ggtitle(sprintf("coefficients for both data sets, regressed against %s (%s)", pheno, DRNA)) +
        theme_bw() + scale_x_continuous(na.value = min.x) +
        geom_text_repel(data = res.figure %>% 
                          # filter(species %in% {key.bact15 %>% gsub(".*\\.s__", "", .) %>% gsub("_", " ", .)}) %>%
                          filter(sig == "FDR significant"), 
                        aes(label = species), hjust = 1, direction = "y",
                        nudge_x = 10, segment.size = 0.2, segment.alpha = 0.5, #segment.colour = "navyblue",
                        size = 3)
      saveRDS(plot.obj, paste0("output/C2_coef_plot_bracken_", 
                               pheno, "-cov_", DRNA,"_ZOE2_validation.rds"))
      ggsave(paste0("figure/C2_coef_plot_bracken_", 
                    pheno, "-cov_", DRNA,"_ZOE2_validation.png"),
             plot = plot.obj,
             width = 10, height = 10)
    }
    
  }   
  
  ### 2. putting all model results altogether
  ## Supp Figure C
    # tmp.sig.comb <- full_join(full_join(tmp.sig[[1]], tmp.sig[[2]], by = "species"), tmp.sig[[3]], by = "species")
    tmp.sig.comb <-
      Reduce(function(x, y) full_join(x, y, by = "species"), tmp.sig) %>% 
      {mutate(., n.sig.models = apply(dplyr::select(., starts_with("sig")), 1, 
                                      function(s) nchar(paste0(ifelse(is.na(s), "", s), collapse = ""))))}
    names(tmp.sig.comb)[grepl("^sig", names(tmp.sig.comb))] <- names(tmp.sig) %>% gsub("CF", "t6", .)
    names(tmp.sig.comb)[grepl("^direction", names(tmp.sig.comb))] <- paste0("direction.", names(tmp.sig) %>% gsub("CF", "t6", .))
    tmp.sig.comb[is.na(tmp.sig.comb)] <- ""
    
    # adding the overall prevalence rate
    tmp.sig.comb$prev = NA
    .initialize(test = "none", add.epsilon = FALSE, filter = "both",  humann = 2, bracken = bracken,
                screen = TRUE, prev.threshold = 0.1, avg.detect = 0.00001, detect.rel.abund = TRUE)
    
    species.nm = gsub(".*\\.s__", "", dat$taxa$bacteria) %>% gsub("\\_", " ", .)
    for (i in 1:length(tmp.sig.comb$species)) {
      index.tmp = which(species.nm == tmp.sig.comb$species[i])
      tmp.sig.comb$prev[i] = mean(dat$otu[index.tmp, , DR.no] > 0)
    }
    
    tmp.sig.comb %>% 
      dplyr::select(prev, species, microt3cb, t3c, everything() ) %>% 
      mutate(direction = paste0(direction.microt3cb, direction.t3c)) %>% 
      mutate(overlap1 = nchar(direction)) %>% 
      mutate(overlap2 = nchar(paste0(microt3cb, t3c))) %>% 
      mutate(direction = gsub("\\-+", "-", gsub("\\++", "+", direction))) %>% 
      mutate(direction = gsub("\\-", "\u2193", gsub("\\+", "\u2191", direction))) %>% 
      arrange(desc(overlap1), desc(overlap2), desc(microt3cb), desc(t3c)) %>% 
      mutate(`overall prevalence` = round(prev, 2)) %>% 
      dplyr::select(species, `overall prevalence`, direction, everything(), -overlap1, -overlap2) %T>%
      saveRDS(summary.fn) %T>% 
      # View %>% # print(n=40)
      write.csv(paste0("output/C2_coef_table_bracken_combined-cov_", DRNA, "_ZOE", ZOE, ".csv"))
    
    
    # plotly::plot_ly(res, x = ~ cont.coef, y = ~disc.coef, z =~lm.coef, col = ~lm.coef)          
  }
}




# ### 2 x 2 validation plot ##############
# #   Added on June 5, 2021. Bracken, t3c & microt3cb, DNA & RNA.
# 
# p <- list()
# i = 1
# for (DRNA in c("DNA", "RNA")) {
#   for (pheno in c("t3c", "microt3cb")) {
#     p[[i]] = 
#       readRDS(paste0("output/C2_coef_plot_bracken_", 
#                      pheno, "-cov_", DRNA,"_ZOE2_validation.rds"))
#     i = i + 1
#   }
# }
# 
# p.legend = get_legend(p[[1]] + 
#                         guides(shape = guide_legend("significance (ZOE 2.0)", nrow = 1),
#                           color = guide_legend("validation by ZOE 2.0 pilot", nrow = 1)) + 
#                         theme(legend.direction = "horizontal", legend.position = "bottom"))
# p.main = plot_grid(plotlist = 
#                      lapply(p, function(x) x + theme(legend.position="none") + 
#                               ggtitle(x$labels$title %>% gsub(".*against ", "", .))), 
#                    nrow = 2, ncol = 2)
# p.final = plot_grid(p.main, p.legend, ncol = 1, rel_heights = c(1, 0.05))
# 
# save_plot(paste0("../figure/C2_coef_plot_bracken_", 
#                  "combined-cov_ZOE2_validation.png"), 
#           p.final, base_height = 15, base_width = 15)



### DNA vs RNA ##############
## Supp Figure D.

  # ### 1. no of significant models D vs R
  #   ZOE = 2
  #   
  #   tmp.sig.comb.D = 
  #     read.csv(paste0("output/C2_coef_table_bracken_", 
  #                     "combined-cov_", "DNA", "_ZOE", ZOE, ".csv"))
  #   tmp.sig.comb.R = 
  #     read.csv(paste0("output/C2_coef_table_bracken_", 
  #                     "combined-cov_", "RNA", "_ZOE", ZOE, ".csv"))
  #   
  #   tmp.sig.comb.DR <-
  #     full_join(tmp.sig.comb.D, tmp.sig.comb.R, by = "species", suffix = c(".D", ".R")) %>% 
  #     mutate(direction.D = ifelse(is.na(direction.D), "", direction.D),
  #            direction.R = ifelse(is.na(direction.R), "", direction.R),
  #            direction = paste0(direction.D, direction.R),
  #            direction = gsub("(.){2}", "\\1", direction)) %>% 
  #     transmute(species, prev.D, prev.R, direction, n.sig.models.D,  n.sig.models.R, 
  #               microt3cb.D, t3c.D, microt3cb.R, t3c.R)
  #   
  #   write.csv(tmp.sig.comb.DR,
  #             paste0("output/C2_coef_table_bracken_", 
  #                    "combined-cov_(DnR)", "_ZOE", ZOE, ".csv"))
  #   
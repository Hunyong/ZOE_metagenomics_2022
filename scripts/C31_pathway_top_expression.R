### C31_pathway_top_expression.R
### Author: Hunyong Cho
### Identify the top expressed pathways (free the context of ECC)
### input: "output/ETable1_C2_coef_table_bracken.csv"
### output: some intermediate outputs: 
###       output/C31_pathway_composition_humann3_RNA_sigTaxa16_Marginal.rds, 
###       output/C31_pathway_top30_humann3_RNA_sigTaxa16_Marginal.csv

### 0.1 library
  library(dplyr)
  source("scripts/F_generic.R")  # typeDRNA(), sample.exclude, taxa.exclude
  
### 0.2 parameters
  
  
### 1. Top pathway expressions
  
  total = "Data-processed/data.humann3.path.marginal.DRNA.ZOE2.rds" %>% readRDS
  total.ranking = data.frame(pathway       = total$otu %>% rownames %>% gsub(" *\\(TOTAL\\)", "", .),
                             rank.overall  = total$otu[,,2] %>% apply(1, sum, na.rm = TRUE) %>% "*"(-1) %>% rank)
  
  # The final list of 16 species.
  key.bact16 = read.csv("output/ETable1_C2_coef_table_bracken.csv")$species
  
  ### path (marginal)
  .initialize(type = "path", ZOE = 2, DR.no = 2, pheno = "t3c",
              bracken = FALSE, humann = 3,
              test = "none", add.epsilon = FALSE, screen = TRUE, 
              nrm = TRUE, nrmScale = 400000, # original scale 383K
              prev.threshold = 0.1,  avg.detect = 0, detect.rel.abund = TRUE)
  
  marginal = 
    data.frame(path = dat$otu %>% rownames,
               species = NA,
               abund.mean = dat$otu[,,1] %>% apply(1, mean, na.rm = TRUE),
               expr.mean = dat$otu[,,2] %>% apply(1, mean, na.rm = TRUE)) %>% 
    arrange(desc(expr.mean))
  
  
  
  ### pathbact (joint)
  .initialize(type = "pathbact", ZOE = 2, DR.no = 2, pheno = "t3c",
              bracken = FALSE, humann = 3,
              test = "none", add.epsilon = FALSE, screen = TRUE, 
              nrm = TRUE, nrmScale = 3000000, # original scale 3.3M
              prev.threshold = 0.1, avg.detect = 0, detect.rel.abund = TRUE)
  
  tmp = 
    data.frame(path = dat$otu %>% rownames %>% gsub(" .*", "", .),
               species = dat$otu %>% rownames %>% gsub(".* ", "", .) %>% gsub(".*\\.s__", "", .) %>% gsub("_", " ", .),
               abund.mean = dat$otu[,,1] %>% apply(1, mean, na.rm = TRUE),
               expr.mean = dat$otu[,,2] %>% apply(1, mean, na.rm = TRUE)) %>% 
    mutate(species = paste0(substr(species, 1, 3), gsub("([^ ]*) ?(.*)", ". \\2", species))) %>% 
    mutate(species = gsub("^unc", "(unclassified)", species)) %>% 
    filter(path != "UNINTEGRATED") %>% 
    arrange(desc(expr.mean))
  
  spec = tmp$species %>% unique %>% sort
  
  
  pathway = tmp$path %>% unique %>% sort
  tab = matrix(0, nrow = length(pathway), ncol = length(spec) + 1, dimnames = list(pathway, c("(total)", spec)))
  
  for (i in 1:dim(tmp)[1]) 
    tab[tmp[i, "path"], tmp[i, "species"]] = tmp[i, "expr.mean"]
  tab[, "(total)"] = apply(tab, 1, sum)
  tab = tab[order(-tab[, "(total)"]), ]  # sorting by the total expression.
  
  
  key.bact16b =
    key.bact16 %>%
    {paste0(substr(., 1, 3), gsub("([^ ]*) ?(.*)", ". \\2", .))} %>% sort
  key.bact16.tab = key.bact16b[key.bact16b %in% colnames(tab)]
  
  
  ### top expressed pathways    
  tab2 = 
    tab %>% 
    as.data.frame() %>% 
    mutate(pathway = rownames(tab)) %>% 
    transmute(pathway,
              total = `(total)`,
              unclassified = `(unclassified). `,
              `unclassified (%)` = round(unclassified / total * 100, 1),
              `key species (#)`  = apply(tab[, key.bact16.tab], 1, sum),
              `key species (%)` = round(`key species (#)`/total * 100, 1),
              `total (#)` = apply(tab[, -1], 1, function(x) sum(x > 0)),
              `key species (#)` = apply(tab[, key.bact16.tab], 1,  function(x) sum(x > 0)),
              `key species names` = apply(tab[, key.bact16.tab], 1, function(x) paste0(key.bact16.tab[x > 0], collapse = ", ")),
              `top 1 species` = apply(tab[, -1], 1, function(x) spec[-x == sort(-x)[1]]),
              `top 2 species` = apply(tab[, -1], 1, function(x) ifelse(sort(-x)[2] == 0, "", spec[-x == sort(-x)[2]])),
              `top 3 species` = apply(tab[, -1], 1, function(x) ifelse(sort(-x)[3] == 0, "", spec[-x == sort(-x)[3]]))) %>% 
    transmute(`#` = 1:n(), pathway, total, 
              `% unclassified` = `unclassified (%)`, 
              `% key species` = `key species (%)`, 
              `# species`= `total (#)`,
              `key species` = `key species names`,
              `# key species` = `key species (#)`,
              `top 3 species` = paste0(`top 1 species`, ", ", `top 2 species`, ", ", `top 3 species`)) %>% 
    left_join(total.ranking)
  saveRDS(tab2, sprintf("output/C31_pathway_composition_humann3_RNA_sigTaxa16_Marginal.rds"))
  
  
  # tab2b = Top 100
  tab2b = tab2 %>% head(101)
  tab2b[101, ] = NA  
  tab2b[101, "#"] = 0
  tab2b[101, "pathway"] = paste0("Total ", dim(tab2)[1], " pathways")
  tab2b[101, "key species"] = paste0(length(key.bact16), "Key species: ", paste(key.bact16, collapse = ", "))
  
  # tab2c = sorted by sig species. Then top 30
  tab2c = 
    tab2b %>% 
    arrange(desc(`% key species`)) %>% 
    filter(1:n()<=30) %>% 
    print

  write.csv(tab2c, sprintf("output/C31_pathway_top30_humann3_RNA_sigTaxa16_Marginal.csv"), row.names = F)

  
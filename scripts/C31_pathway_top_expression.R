### C31_pathway_top_expression.R
### Author: Hunyong Cho
### Identify the top expressed pathways (free the context of ECC)


### 0.1 library
  source("scripts/F.generic.R")  # typeDRNA(), sample.exclude, taxa.exclude
  library(ggplot2); library(ggrepel)
  
### 0.2 parameters
  DR.no = 2 ;  ZOE <- 2; pheno <- "t3c"; index <- "";
  model.nm = "cov"
  bracken = FALSE   # bracken = FALSE  for pathbact !!!!
  humann = 3       
  nrm = TRUE #normalization
  threshold.pval = 0.10
  
  
### 1. Top pathway expressions
  
  total = "Data-processed/data.humann3.path.marginal.DRNA.ZOE2.rds" %>% readRDS
  total.ranking = data.frame(pathway       = total$otu %>% rownames %>% gsub(" *\\(TOTAL\\)", "", .),
                             rank.overall  = total$otu[,,2] %>% apply(1, sum, na.rm = TRUE) %>% "*"(-1) %>% rank)
  
  # The final list of 16 species.
  key.bact16 = c("g__Prevotella.s__Prevotella_salivae",
                 "g__Prevotella.s__Prevotella_oulorum",
                 "g__Prevotella.s__Prevotella_melaninogenica",
                 "g__Prevotella.s__Prevotella_sp_oral_taxon_306",
                 "g__Prevotella.s__Prevotella_veroralis",
                 "g__Streptococcus.s__Streptococcus_mutans",
                 "g__Selenomonas.s__Selenomonas_sputigena",
                 "g__Leptotrichia.s__Leptotrichia_wadei",
                 "g__Leptotrichia.s__Leptotrichia_sp_oral_taxon_847",
                 "g__Leptotrichia.s__Leptotrichia_sp_oral_taxon_223",
                 "g__Leptotrichia.s__Leptotrichia_sp_oral_taxon_498",
                 "g__Lachnoanaerobaculum.s__Lachnoanaerobaculum_saburreum",
                 "g__Lachnospiraceae.s__Lachnospiraceae_bacterium_oral_taxon_082",
                 "g__Veillonella.s__Veillonella_atypica",
                 "g__Stomatobaculum.s__Stomatobaculum_longum",
                 "g__Centipeda.s__Centipeda_periodontii"
  )
  sigTaxa.nm16  = "sigTaxa16" 
  
  key.bact = key.bact16
  sigTaxa.nm = sigTaxa.nm16
  
  ### path (marginal)
  type1 =  "path"
  .initialize(test = "none", add.epsilon = FALSE, screen = TRUE, nrmScale = 400000,
              prev.threshold = 0.1,  avg.detect = 0, detect.rel.abund = TRUE, 
              bracken = bracken, humann = humann,        ### added March 24, 2021
              screen.among.tested.DNAs = TRUE)
  
  marginal = 
    data.frame(path = dat$otu %>% rownames,
               species = NA,
               abund.mean = dat$otu[,,1] %>% apply(1, mean, na.rm = TRUE),
               expr.mean = dat$otu[,,2] %>% apply(1, mean, na.rm = TRUE)) %>% 
    arrange(desc(expr.mean))
  
  
  
  ### pathbact (joint)
  type1 =  "pathbact"
  .initialize(test = "none", add.epsilon = FALSE, screen = TRUE, nrmScale = 3000000,
              prev.threshold = 0.1, avg.detect = 0, detect.rel.abund = TRUE, 
              bracken = bracken, humann = humann,        ### added March 24, 2021
              screen.among.tested.DNAs = TRUE)
  
  
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
  
  
  key.bact.tidy =
    key.bact %>% gsub(".*\\.s__", "", .) %>% gsub("_", " ", .) %>%
    {paste0(substr(., 1, 3), gsub("([^ ]*) ?(.*)", ". \\2", .))} %>% sort
  key.bact.tidy.tab = key.bact.tidy[key.bact.tidy %in% colnames(tab)]
  
  
  ### 1. tab, tab2: top expressed pathways    
  tab2 = 
    tab %>% 
    as.data.frame() %>% 
    mutate(pathway = rownames(tab)) %>% 
    transmute(pathway,
              total = `(total)`,
              unclassified = `(unclassified). `,
              `unclassified (%)` = round(unclassified / total * 100, 1),
              `key species (#)`  = apply(tab[, key.bact.tidy.tab], 1, sum),
              `key species (%)` = round(`key species (#)`/total * 100, 1),
              `total (#)` = apply(tab[, -1], 1, function(x) sum(x > 0)),
              `key species (#)` = apply(tab[, key.bact.tidy.tab], 1,  function(x) sum(x > 0)),
              `key species names` = apply(tab[, key.bact.tidy.tab], 1, function(x) paste0(key.bact.tidy.tab[x > 0], collapse = ", ")),
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
  saveRDS(tab2, sprintf("output/C31_pathway_composition_humann%d_%s_%sMarginal.rds", humann, "RNA", sigTaxa.nm))
  
  
  # tab2b = Top 100
  tab2b = tab2 %>% head(101)
  tab2b[101, ] = NA  
  tab2b[101, "#"] = 0
  tab2b[101, "pathway"] = paste0("Total ", dim(tab2)[1], " pathways")
  tab2b[101, "key species"] = paste0(length(key.bact), "Key species: ", paste(key.bact.tidy, collapse = ", "))
  
  # tab2c = sorted by sig species. Then top 30
  tab2c = 
    tab2b %>% 
    arrange(desc(`% key species`)) %>% 
    filter(1:n()<=30) %>% 
    print

  write.csv(tab2c, sprintf("output/C31_pathway_top30_humann%d_%s_%sMarginal.csv", humann, "RNA", sigTaxa.nm), row.names = F)

### C51_heatmap.R
### Author: Bridget Lin
### Modified by: Hunyong Cho
### output: Figure 3  


### 0. library and data
### 0.1 library
  library(dplyr)
  library(pheatmap)
  library(gridExtra)
  source("scripts/F_generic.R")  # typeDRNA(), sample.exclude, taxa.exclude

### 0.2 data
  .initialize(type = "bact", ZOE = 2, DR.no = 1, pheno = "microt3cb",
              test = "LN", screen = TRUE, prev.threshold = 0.1, 
              avg.detect = 0.00001, detect.rel.abund = TRUE, nrm = TRUE,
              bracken = TRUE, humann = 2,
              screen.among.tested.DNAs = FALSE)
  
  n.taxa = dim(dat$otu)[1]
  health = dat$meta$micro_t3c_b == 0
  n.subj0 = sum(health) # 190
  n.subj1 = sum(!health) # 110
  
### 0.3 the 16 sig species (This is obtained from the C2x.R files)
  key.bact16 = read.csv("output/ETable1_C2_coef_table_bracken.csv")$species
  
### 1. Residuals for DNA
  
  # 1.1 health
  resid_D0 <- matrix(NA, n.subj0, n.taxa)
  for(i in 1:n.taxa){
    lm1 <- lm(log(dat$otu[i, health, 1] + 1) ~ dat$meta$batch.DNA[health] +dat$meta$agemo[health])
    lm.resid <- resid(lm1)
    resid_D0[, i] = lm.resid
  }
  # 1.2 disease
  resid_D1 <- matrix(NA, n.subj1, n.taxa)
  for(i in 1:n.taxa){
    lm1 <- lm(log(dat$otu[i, !health, 1] + 1) ~ dat$meta$batch.DNA[!health] +dat$meta$agemo[!health])
    lm.resid <- resid(lm1)
    resid_D1[, i] = lm.resid
  }

  
### 2. Residuals for RNA
  included = dat$meta$exclusion == "included" # for RNA
  healthR = health & included
  diseaseR = !health & included
  
  # 2.1 health
  resid_R0 <- matrix(NA, sum(healthR), n.taxa)
  for(i in 1:n.taxa){
    lm1 <- lm(log(dat$otu[i, healthR, 2] + 1) ~ dat$meta$batch.RNA[healthR] +dat$meta$agemo[healthR])
    lm.resid <- resid(lm1)
    resid_R0[, i] = lm.resid
  }
  # 2.2 disease
  resid_R1 <- matrix(NA, sum(diseaseR), n.taxa)
  for(i in 1:n.taxa){
    lm1 <- lm(log(dat$otu[i, diseaseR, 2] + 1) ~ dat$meta$batch.RNA[diseaseR] +dat$meta$agemo[diseaseR])
    lm.resid <- resid(lm1)
    resid_R1[, i] = lm.resid
  }



### 3. Correlation matrices
  
  names_microb = dat$taxa$bacteria
  cor_D0 <- cor(resid_D0[, names_microb %in% key.bact16])
  cor_D1 <- cor(resid_D1[, names_microb %in% key.bact16])
  cor_R0 <- cor(resid_R0[, names_microb %in% key.bact16 ])
  cor_R1 <- cor(resid_R1[, names_microb %in% key.bact16 ])
  
  rownames(cor_D0) <- colnames(cor_D0) <- colnames(cor_D1) <- rownames(cor_D1) <- 
    rownames(cor_R0) <- colnames(cor_R0) <- colnames(cor_R1) <- rownames(cor_R1) <- 
    names_microb[names_microb %in% key.bact16]
 
  # removing the correlations of themselves
  for (i in 1:16) {
    cor_D0[i, i] = cor_D1[i, i] = cor_R0[i, i] = cor_R1[i, i] = NA
  }
  
  # overall correlations
  cor.overall = 
    data.frame(`health in MTG`  = cor_D0 %>% apply(1, mean, na.rm = T),
               `disease in MTG` = cor_D1 %>% apply(1, mean, na.rm = T),
               `health in MTX` = cor_R0 %>% apply(1, mean, na.rm = T),
               `disease in MTX` = cor_R1 %>% apply(1, mean, na.rm = T)) %>% 
      mutate_all(function(x) round(x, 2)) 
  cor.overall %>% View

  
### 4. plots
  # Consistent ordering between disease and health
  order1 <- hclust(dist(cor_D0))$order
  cor_D0a = cor_D0[order1, order1]
  cor_D1a = cor_D1[order1, order1]
  rownames(cor_D0a) = 
    sprintf("%s (%0.2f)", colnames(cor_D0), cor.overall[, 1])[order1]
  rownames(cor_D1a) = 
    sprintf("%s (%0.2f)", colnames(cor_D0), cor.overall[, 2])[order1]
  
  order2 <- hclust(dist(cor_R0))$order
  cor_R0a = cor_R0[order2, order2]
  cor_R1a = cor_R1[order2, order2]
  rownames(cor_R0a) = 
    sprintf("%s (%0.2f)", colnames(cor_R0), cor.overall[, 3])[order2]
  rownames(cor_R1a) = 
    sprintf("%s (%0.2f)", colnames(cor_R0), cor.overall[, 4])[order2]
  
  
  graphicParams = 
    expand.grid(DR = c("MTG", "MTX"), ecc = 0:1) %>% 
    mutate(eccNm = ifelse(ecc, "any ", "no "),
           title = c("A", "C", "B", "D"),
           label = sprintf("%s, %slocalized disease experience", DR, eccNm),
           fn = sprintf("figure/Fig3_%s%s.pdf", DR, ecc))
  corList = list(cor_D0a, cor_R0a, cor_D1a, cor_R1a)

  plotList = list()
  for (i in 1:4) {
    plotList[[i]] = 
      pheatmap(corList[[i]], 
               main = graphicParams$label[i],
               breaks=seq(-1, 1, by = 0.02), cluster_rows = F, cluster_cols = F,
               treeheight_row = 0, treeheight_col = 0, fontsize = 12)[[4]]
    # dev.off()
  }
  plotAll <- grid.arrange(arrangeGrob(grobs= plotList, ncol=2, labels = c("A", "C", "B", "D")))
  ggsave("figure/Fig3_heatmap.pdf", plotAll, width = 22, height = 20)
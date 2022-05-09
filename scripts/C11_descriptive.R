###### Nonzero mean vs zero proportion ######################################################
###### C06.00.trait_screen.R       ######################################################

bracken = TRUE; humann = 2
type <- type1 <- "bact"; ZOE <- 2; DR.no = 1; 
pheno <- "t3c"
nrm = TRUE #normalization
brk = ifelse(bracken, "bracken", "humman2")
trait <- c("t3c", "micro_t3c_b")

library(tidyverse); library(magrittr); library(energy)
source("scripts/F.generic.R")  # typeDRNA(), sample.exclude, taxa.exclude
source("scripts/F.eco.R")  # typeDRNA(), sample.exclude, taxa.exclude

for (ZOE in 1:2) {
    
  # Subjects filtered / Virus not included / TPM-Normalized (=rescaled) / noninformative taxa excluded
  .initialize(test = "none", add.epsilon = F, 
              filter = "both", screen = FALSE, humann = humann, bracken = bracken)
  
  
  # 1. Demographics and phenotype distribution
  dat$meta[, trait] %>% apply(2, function(x) c(mean = mean(x, na.rm = TRUE),
                                               sd = sd(x, na.rm = TRUE)))
  
  
  # 2. dCor  
  set.seed(10)
  result <-
    sapply(trait, function(s) {
      print(s)
      if (is.null(dat$meta[[s]])) return(c(est.dCor = NA, pval = NA))
      tmp.complete <- !is.na(dat$meta[[s]])
      dist.cor(t(dat$otu[,tmp.complete, DR.no]), dat$meta[[s]][tmp.complete], pval = TRUE, R = 1000)
    }) %>% print
  saveRDS(result, paste0("output/C11_dCor_bact_DNA_ZOE", ZOE, "_", brk, ".rds"))
  write.csv(result, paste0("output/C11_dCor_bact_DNA_ZOE", ZOE, "_", brk, ".csv"))
  
  
  # 3. Shannon entropy
  ecoDat <- ecoSys(dat, DR.no = DR.no)
  # shannon <- as.matrix(microbiome::global(dat$otu[,, DR.no], index = "shannon"))
  shannon.p <-
    trait %>% 
    sapply(function(x) {
      form <- paste0("shannon  ~ ", x)
      lmout <- lm(as.formula(form), data = ecoDat)
      summary(lmout)$coef[x, c("Estimate", "Pr(>|t|)")]
    }) %>% print
}

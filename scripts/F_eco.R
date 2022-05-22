
### phyloseq conversion

DRNA2phylo <- function(data, DR.no = 1) {
  require(phyloseq)
  cat("converting", c("DNA", "RNA")[DR.no],"\n")
  taxa <- as.matrix(data$taxa); rownames(taxa) = data$taxa[,1]
  meta <- data$meta; rownames(meta) = data$meta$id
  meta <- sample_data(meta)
  rownames(meta) <- meta$id
  phyloseq(otu_table(data$otu[,,DR.no], taxa_are_rows = TRUE), 
           tax_table(taxa),
           meta)  
}

ecoSys <- function(data, DR.no = 1) { #DR.no = 1 <=> DNA
  require(microbiome)
  eco = data.frame(mean.exp = apply(data$otu[,, DR.no], 2, mean, na.rm = TRUE),
                   shannon = microbiome::alpha(data$otu[,, DR.no], index = "shannon")[[1]],
                   rich = microbiome::richness(data$otu[,, DR.no], index = "observed")[[1]],
                   domin = microbiome::dominance(data$otu[,, DR.no], index = "relative")[[1]])
  
  eco <- cbind(data$meta, eco)
  eco
}


#' @examples 
#' divergence2(dat, DR.no = 1, group = cariesfree)  #group argument requires nonstandard evaluation
divergence2 <- function(data, DR.no, group = "cariesfree") {
  require(dplyr)
  require(microbiome)
  require(rlang)
  dat.phylo <- DRNA2phylo(data, DR.no = DR.no)
  # levels <- data$meta[[deparse(substitute(group))]] %>% unique
  levels <- data$meta[[group]] %>% unique
  levels <- levels[!is.na(levels)]
  cat("Divergence by ", group, " groups with levels = ", levels, ".\n")
  # beta <- lapply(levels, function(s) eval(substitute(divergence(subset_samples(dat.phylo, group == s)))))
  
  beta <- list()
  k = 1
  for (i in levels) {
    # beta[[k]] <- eval(substitute(divergence(subset_samples(dat.phylo, group == i))))
    # dat.tmp <- subset_samples(dat.phylo, {{sym(group)}} == i)
    samp_vec <- sample_data(dat.phylo)[[group]]
    dat.tmp   <- prune_samples(samp_vec == i, dat.phylo)
    reference <- apply(abundances(dat.tmp), 1, median)
    beta[[k]] <- divergence(dat.tmp, reference, method = "bray")
    # this failed
    # grp.var <- enquo(group)
    # dat.tmp <- subset_samples(dat.phylo, !!grp.var == i)
    k = k + 1
  }
  beta <- data.frame(group = rep(levels, sapply(beta, length)), beta = unlist(beta))
  # names(beta)[1] <- deparse(substitute(group))
  names(beta)[1] <- "group"
  beta
}

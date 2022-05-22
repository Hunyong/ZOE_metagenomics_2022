# time calculator
tt <- function(s){
  if (s==1) {time.tmp <<- Sys.time() # record time
  } else if (s==2) { # calculate time
    return(data.frame(begin = time.tmp, end = Sys.time(), elapsed = Sys.time() - time.tmp))
  }
}

eps <- function(mat, min.thres = 1) {
  if (class(mat)[1] == "data.frame") {
    mat <- as.data.frame(mat)
  }
  if (min(mat, na.rm = TRUE) > 0) {
    warning("All elements are already positive.")
    return(0)
  }
  eps.1 <- min(mat[mat > 0], na.rm = TRUE)
  eps <- min(eps.1, min.thres)
  cat("epsilon = ", eps, " (min positive value was ", eps.1, ").", "\n")
  return(eps)
}

epsilon.add <- function(mat, eps = NULL, min.thres = 1) {
  if (class(mat)[1] == "data.frame") {
    mat <- as.data.frame(mat)
  }
  if (is.null(eps)) eps <- eps(mat, min.thres = min.thres)
  mat + eps
}

pick.rows <- function(dat, row.index = NULL, row = NULL, exclude = FALSE, otu.only = FALSE) {
  if (is.null(row.index) & is.null(row)) stop("Either row.index or row should be provided.")
  if (is.null(dim(dat$taxa))) dat$taxa <- as.data.frame(dat$taxa)
  
  vec <- if (!is.null(row.index)) {
    as.character(dat$taxa[,1]) %in% row.index
  } else {1:length(dat$taxa[,1]) %in% row}
  if (exclude) {vec <- !vec}
  
  dat$otu <- dat$otu[vec,,, drop=FALSE]
  dat$taxa <- dat$taxa[vec, , drop=FALSE]
  if ("reads" %in% names(dat)) dat$reads <- dat$reads[vec,,, drop=FALSE]
  
  cat("dimension after pick.cols = ", dim(dat$otu), "\n")  
  
  if (otu.only) {return(dat$otu)}
  dat
}

pick.cols <- function(dat, col.index = NULL, col = NULL, exclude = FALSE, otu.only = FALSE) {
  if (is.null(col.index) & is.null(col)) stop("Either col.index or col should be provided.")
  
  vec <- if (!is.null(col.index)) {
    as.character(dat$meta$id) %in% col.index
  } else {1:length(dat$meta$id) %in% col}
  if (exclude) {vec <- !vec}
  dat$otu <- dat$otu[,vec,, drop=FALSE]
  dat$meta <- dat$meta[vec,, drop=FALSE]
  cat("dimension after pick.cols = ", dim(dat$otu), "\n")  
  
  if (otu.only) {return(dat$otu)}
  dat
}
if (FALSE) { pick.rows(dat, "UniRef90_A0A011MIA4 g__Haemophilus.s__Haemophilus_paraphrohaemolyticus")}

#' @examples 
#' typeDRNA()
typeDRNA <- function(type = c("path", "gene", "bact", "pathbact", "genebact",
                              "bactRPK", "geneRPK", "genebactRPK"), 
                     DR.no = 1, ZOE = 1, withVirus = TRUE, normalize = FALSE,
                     bracken = FALSE,
                     humann = 2) {
  type <- match.arg(type)
  DRNA <- c("DNA", "RNA")[DR.no]
  batch <- c("batch.DNA", "batch.RNA")[DR.no]
  # phase <- ifelse(as.numeric(ZOE) == 1, "180", "ZOE2")
  phase <- ifelse(as.numeric(ZOE) == 1, "ZOE1", "ZOE2") # updated 2020/05/27
  humann3. <- if (humann == 2) "" else "humann3."
  reads <- paste0("reads.", DRNA)
  if (!bracken) {
    if (type == "path") {
      dat.file <- sprintf("Data-processed/data.%spath.marginal.DRNA.%s.rds", humann3., phase)
      Type <- "path";       Unit <- "counts"; type.name <- "pathways"
    } else if (type %in% c("gene", "geneRPK")) {
      dat.file <- sprintf("Data-processed/data.%sgeneRPK.marginal.DRNA.%s.rds", humann3., phase)
      Type <- "geneRPK";    Unit <- "RPKs"; type.name <- "gene families"
    } else if (type %in% c("bact", "bactRPK")) {
      dat.file <- sprintf("Data-processed/data.%sbactRPK.marginal.DRNA.%s.rds", humann3., phase)
      Type <- "bactRPK";    Unit <- "RPKs"; type.name <- "species"
    } else if (type %in% c("genebact", "genebactRPK")) {
      dat.file <- sprintf("Data-processed/data.%sgeneRPK.joint.DRNA.%s.rds", humann3., phase)
      Type <- "genebactRPK";    Unit <- "RPKs"; type.name <- "gene-species"
    } else if (type == "pathbact") {
      dat.file <- sprintf("Data-processed/data.%spath.joint.DRNA.%s.rds", humann3., phase)
      Type <- "pathbact";    Unit <- "counts"; type.name <- "path-species"
    } else stop("type is not correct.")
  } else {
    if (type == "bact") {
      dat.file <- sprintf("Data-processed/data.bracken2.full.DRNA.%s.rds", phase)
      Type <- "bact";       Unit <- "counts"; type.name <- "species"
    } else stop("type is not correct.")
  }
  if (withVirus) {dat.file <- gsub("Data-processed", "Data-processed-withVirus", dat.file)}
  abundance <- c("abundance", "expression")[DR.no]
  if (normalize) {
    Unit <- gsub("RPKs", "TPM", Unit)
    Unit <- gsub("^counts", "normalized counts (M/sum)", Unit)
    # Unit <- gsub("^reads", "normalized reads (M/sum)", Unit)
  }
  cat("Type: ", Type, "\n")
  cat("DRNA: ", DR.no, ", ", if (is.null(DR.no)) "DRNA not selected." else DRNA, "\n")
  cat("file to read: ", dat.file, "\n")
  list(Type = Type, DR.no = DR.no, DRNA = DRNA, batch = batch, reads = reads,
       type.name = type.name,
       dat.file = dat.file, Unit = Unit, abundance = abundance, withVirus = withVirus)
}
normalize <- function(dat, scale = 1E6) {
  dim.otu <- dim(dat$otu)
  multiplier <- apply(dat$otu, 2:3, sum, na.rm = TRUE)/scale
  # regularize all zero vector
  multiplier <- ifelse(multiplier == 0, 1, multiplier)
  # duplicate the multiplier for all taxa
  multiplier <- array(multiplier, dim.otu[c(2,3,1)])
  # reorder the dimensions to agree with the original data
  multiplier <- aperm(multiplier, c(3,1,2))
  dat$otu <- dat$otu / multiplier
  dat
}

#' @examples 
#' screen.gene(dat)
screen.gene <- function(data, DR.no, avg.detect = 2, prevalence = 0.03, screen.among.tested.DNAs = FALSE) {
  # data: each column = sample, each row = genes
  # gene.col: the column where gene names are provided. If null, ignored.
  gene.name = data$taxa[,1]
  n.genes = dim(data$otu)[1]
  n.sample = dim(data$otu)[2]
  
  if (DR.no == 2 && screen.among.tested.DNAs) {
    cat("Screening for RNAs is done only among those screened taxa for DNAs.\n")
    useF1.D = rowSums(data$otu[,,1], na.rm = TRUE) >= avg.detect * n.sample
    useF2.D = apply(data$otu[,,1], 1, function(s) mean(s > 0, na.rm = TRUE)) >= prevalence
    useF.D = (useF1.D & useF2.D)
    cat("total taxa = ", n.genes, ", those tested for DNAs = ", sum(useF.D), "\n")
  } else {
    useF.D = rep(TRUE, n.genes) # by default
  }
  
  useF1 = rowSums(data$otu[,,DR.no], na.rm = TRUE) >= avg.detect * n.sample
  useF2 = apply(data$otu[,,DR.no], 1, function(s) mean(s > 0, na.rm = TRUE)) >= prevalence
  useF = which(useF1 & useF2 & useF.D)
  cat("total taxa = ", n.genes, 
      if (DR.no) {paste0(", candidates (species tested for DNAs) = ", sum(useF.D))},
      ", after prevalence = ", sum(useF2 & useF.D), 
      ", after abundance w/o prevalence = ", sum(useF1 & useF.D), 
      ", after both = ", length(useF))
  print(paste(length(useF), "out of ", n.genes, " will be used.\n"))
  data <- pick.rows (data, row = useF, row.index = NULL, exclude = FALSE, otu.only = FALSE)
  data$screen <- list(feature = data.frame(index = useF, names = gene.name[useF]),
                      stat = c(use = length(useF), out.of = n.genes, 
                               out.of.DNA.tested = if(DR.no == 1) NA else sum(useF.D),
                               percentage = round(length(useF)/n.genes, 2)*100,
                               after.prev = sum(useF2 & useF.D), 
                               after.abund = sum(useF1 & useF.D)))
  data
}

typePheno <- function(pheno) {
  cat("pheno.var, pheno.name, lvl, n.lvl are globally assigned.")
  pheno.var <<- switch(pheno, 
                       CF = "caries", ecc = "ecc", t3c = "t3c", t6c = "t6c",
                       t3cb = c("meta_t3c_b", "micro_t3c_b"), t3 = "t3", t6 = "t6",
                       microt3b = "micro_t3_b", microt3cb = "micro_t3c_b",
                       lca1 = "ecc_type", lca2 = "ecc2_type",  # old code
                       lca1ur = "ecc_type_ur", lca2ur = "ecc2_type_ur")  # old code
  # LCA, old traits
  if (grepl("lca([0-9])(o|u)([0-9])", pheno) & is.null(pheno.var))
    pheno.var <<- # maps, e.g., lca6o5 to t6_overall_5cm
      gsub("lca([0-9])(o|u)([0-9])", "lca_t\\1_\\2\\3", pheno)
  # LCA, new traits
  if (grepl("lca\\.(t[0-9])\\.([0-9])cm", pheno) & is.null(pheno.var))
    pheno.var <<- pheno
  
  pheno.name <<- switch(pheno, CF = "caries", ecc = "ECC", t3c = "t3c", t6c = "t6c",
                        t3cb = "t3c b", t3 = "t3", t6 = "t6",
                        microt3b = "micro t3 b", microt3cb = "micro t3c b",
                        lca1 = "LCA1", lca2 = "LCA2", 
                        lca1ur = "LCA1 (untreated)", lca2ur = "LCA2 (untreated)")
  if (grepl("lca([0-9])(o|u)([0-9])", pheno) & is.null(pheno.name))
    pheno.name <<- # maps, e.g., lca6o5 to t6_overall_5cm
    gsub("lca([0-9])(o|u)([0-9])", "LCA (t\\1_\\2_\\3 levels)", pheno) %>% 
    gsub("_o_", ", overall, ", .) %>% 
    gsub("_u_", ", unrestored, ", .)
  
  if (grepl("lca\\.(t[0-9])\\.([0-9])cm", pheno) & is.null(pheno.name))
    pheno.name <<- # maps, e.g., lca6o5 to t6_overall_5cm
    gsub("lca\\.(t[0-9])\\.([0-9])cm", "LCA (\\1 \\2-class membership)", pheno)
    
  lvl       <<- switch(pheno, CF = c("caries", "cariesfree"), ecc = c("healthy", "restored", "caries"), 
                       t3c = NA, t6c = NA, t3cb = NA, t3 = NA, t6 = NA,
                       microt3b = NA, microt3cb = NA,
                       lca1 = c("healthy", "I", "II", "III", "IV"),
                       lca2 = c("healthy", "I", "II", "III", "IV", "V"),
                       lca1ur = c("healthy", "restored", "I", "II", "III"),
                       lca2ur = c("healthy", "restored", "I", "II", "III", "IV"))
  if (grepl("lca", pheno) & is.null(pheno.name)) {
    lvl <<- levels(dat$meta[[pheno.var]])
  }
  
  n.lvl     <<- switch(pheno, CF = 1, ecc = 2, t3c = 1, t6c = 1, t3cb = 1, 
                       t3 = 1, t6 = 1, microt3b = 1, microt3cb = 1,
                       lca1 = 5, lca2 = 6, lca1ur = 5, lca2ur = 6)
  if (grepl("lca", pheno) & is.null(pheno.name))
    n.lvl <<- length(lvl)
  
  
  metadata  <<- dat$meta[, pheno.var]
  if (pheno == "t3cb") metadata <<- apply(metadata, 1, sum)
  dat$meta$pheno <<- unlist(metadata) ## added Sep 9, 2020
}

.initialize <- function(type, ZOE, DR.no, pheno,
                        test = c("none", "LN"), add.epsilon = FALSE, 
                        filter = c("both", "subjects only", "taxa only", "none"),
                        screen = TRUE, pheno.as.factor = NULL, nrmScale = NULL,
                        prev.threshold = 0.1, nrm = TRUE,
                        avg.detect = NULL, detect.rel.abund = FALSE, 
                        humann, bracken, screen.among.tested.DNAs = FALSE  # screen RNAs among those species screened for DNAs.
                       ) {
  # filter = removing irrelevant subjects and taxa
  # screen = removing by prevalence and abundance screening
  filter = match.arg(filter)
  filter_taxa = (filter %in% c("both", "taxa only"))
  filter_subjects = (filter %in% c("both", "subjects only"))
  
  test = match.arg(test)
  
  library(tidyverse); library(magrittr); library(gamlss); library(coin)
  if (test == "LN") add.epsilon = TRUE
  if (is.null(pheno.as.factor)) {
    pheno.as.factor = if (test == "KW") TRUE else FALSE
  }
  
  ts   <-  typeDRNA(type = type, DR.no = DR.no, ZOE = ZOE, withVirus = FALSE, normalize = nrm,
                    bracken = bracken, humann = humann)
  Type <<- ts$Type; dat.file <<- ts$dat.file;    batch <<- ts$batch;    DRNA <<- ts$DRNA; 
  reads <<- ts$reads
  Unit <<- ts$Unit; abundance <<- ts$abundance;  type.name <<- ts$type.name
  dat  <<- readRDS(dat.file)
  
  ### filtering
  virus.taxa <<- dat$taxa[grep(virus.keyword, dat$taxa[,1]), 1] %>% as.character
  dim.before.filter.subject = dim(dat$otu)
  if (filter_subjects) {
    exclude.subj = if (DR.no==1) sample.exclude else sample.exclude.RNA
    if (length(exclude.subj) > 0)
      dat <<- dat %>%  # excluding specific taxa and samples
        pick.cols(col.index = exclude.subj, exclude = TRUE, otu.only = FALSE)
  }
  dim.after.filter.subject = dim(dat$otu)
  if (filter_taxa)
    dat <<- dat %>%
      pick.rows(row.index = taxa.exclude, exclude = TRUE, otu.only = FALSE)
  rownames(dat$otu) <<- gsub(" \\(TOTAL\\)", "", rownames(dat$otu)); dat$taxa[,1] <<- gsub(" \\(TOTAL\\)", "", dat$taxa[,1])
  dim.after.filter.taxa = dim(dat$otu)
  
  ### normalizing
  mean.library.size = apply(dat$otu[,, DR.no], 2, sum, na.rm = TRUE) %>% mean
  if (nrm & is.null(nrmScale)) {
    if (bracken) {
      if (ZOE == 2) # Last updated on Jan 22, 2021
        nrmScale = switch(DR.no, `1` =  8e+6, `2` = 11e+6) # Bracken
      if (ZOE == 1) # Last updated on Jan 22, 2021
        nrmScale = switch(DR.no, `1` =  5e+6, `2` = 3e+6) # Bracken
    } else if (humann == 3) {                                       # Humann3
      if (ZOE == 2) # Last updated on March 24, 2021
        nrmScale = switch(DR.no, `1` =  7e+6, `2` = 14e+6) # Humann3 7,212,041 / 13,928,475
      if (ZOE == 1)
        stop("Not defined.")
    } else {
        stop("Not defined.")
    }
    
    warning(sprintf("nrmScale is not specified. Set as %0.0f, while the mean library size is %0.0f", 
                    nrmScale, mean.library.size))
  }
  
  if (nrm) dat <<- normalize(dat, scale = nrmScale) # nrmScale = 1e+6
  if (screen & !nrm & detect.rel.abund) stop("normalization is not done, so the abundance screening cannot be done by the relative abundance.")
  dat$meta[, c("ST.DNA", "ST.RNA")] <<- dat$otu %>% apply(2:3, sum, na.rm = TRUE)
  
  # check congruence of ids
  cat("congruence of ids: ", all(dat$otu %>% colnames == dat$meta$id), "\n")
  
  ### Screening
  ############################################################################################
  ## 1.1 regularization: adding epsilon for log transformation, screening out
  # before screening, 1. find epsilon (min positive value). 2. Then screen out. 3. Add the epsilon
  
  ## screening first! then epsilon addition
  if (is.null(avg.detect)) {
    avg.detect = switch(Type, bactRPK = 2, path = 2, geneRPK = 2/10, 
                        genebactRPK = 2/10, pathbact = 2/10,
                        bact = 2, gene = 2/10, genebact = 2/10)
    warning(sprintf("avg.detect is not specified. avg.detect is set as %0.2f", avg.detect))
  } 
  
  # # of species ~= 200, pathways ~= 400, genes = 400K, path+species = 6K, gene+species = 500K (ZOE2, DNA)
  # the numbers for DNAs are comparable with those of RNAs
  # screening by abundance is done to remove errors from low valued data.
  # Thus the more aggregated (e.g. species), less erratic, and singletons (gene-species) are more likely to have errors,
  # but we applied scales less stringent filters to those finely defined taxa (gene-species).
  #  i.e., instead of the same threshold of 2, applied a ten times smaller threshold.

  
  if (screen)
    dat       <<- screen.gene(dat, DR.no = DR.no, avg.detect = avg.detect * nrmScale, 
                              prevalence = prev.threshold, screen.among.tested.DNAs = screen.among.tested.DNAs)
  dim.after.screen = dim(dat$otu)
  print(paste0("dimension after screening: ", paste(dim(dat$otu), collapse = " / ")))
  
  ## Epsilon addition
  if (add.epsilon) {
    dat$epsilon       <<- data.frame(eps.DNA = eps(dat$otu[,,1]), eps.RNA = eps(dat$otu[,,2]))
    dat$otu[,, DR.no] <<- epsilon.add(dat$otu[,, DR.no], eps = as.numeric(dat$epsilon[DR.no]))  # adding epsilon
  }

  ## phenotype information
  typePheno(pheno = pheno)

  
  ## summary of the data.
  dat$info <<- data.frame(filename = dat.file,
                          reference = ifelse(bracken, "Bracken", 
                                             ifelse (humann == 2, "Humann2", "Humann3")),
                          ZOE = ZOE,
                          DRNA = DRNA,
                          type = Type,
                          taxa.level = "species",
                          normalize = nrm,
                          nrmScale = if (!is.null(nrmScale)) nrmScale else NA,
                          originalScale = mean.library.size,
                          virus.included = FALSE,
                          filter_taxa = filter_taxa,
                          filter_subjects = filter_subjects,
                          screen = screen,
                          screen.prev.cutoff = prev.threshold,
                          screen.abund.cutoff = avg.detect,
                          screen.rel.abund = detect.rel.abund,
                          n.before.filter = dim.before.filter.subject[2],
                          n.after.filter = dim.after.filter.subject[2],
                          p.before.filter = dim.before.filter.subject[1],
                          p.after.filter = dim.after.filter.subject[1],
                          p.DNA.tested = if(DR.no == 1 | !screen) NA else dat$screen$stat["out.of.DNA.tested"],
                          p.after.screen = dim.after.screen[1],
                          p.after.screen.prev.only = if(screen) dat$screen$stat["after.prev"] else NA,
                          p.after.screen.abund.only = if(screen) dat$screen$stat["after.abund"] else NA,
                          epsilon.added = add.epsilon)
  
  n.taxa    <<- dim(dat$taxa)[1]
  
  # if (filter_subjects) {
  #   # removing NA samples
  #   complete <- (!is.na(dat$meta$cariesfree) & !is.na(dat$meta$cariesfree))
  #   if (sum(!complete) > 0) {
  #     dat$otu <- dat$otu[, complete, ]
  #     dat$meta <- dat$meta[complete, ]
  #   }
  # }
  # # rm(metadata)
  
  #other phenotypes
  if (ZOE == 2) {
    race <- factor(dat$meta$race2)
    levels(race) <- c("1black", "3other", "2white", "3other", "3other")
    race <- factor(race, levels = c("1black", "2white", "3other"))
    dat$meta$race <<- race
    race <<- race
  } else {
    race <<- dat$meta$race <<- dat$meta$race %>% as.factor
  }
  
  ## CDR (cellular detection rate. Before adding an epsilon.)
  cdr = colMeans(dat$otu[,, DR.no] > 
                   if(add.epsilon) dat$epsilon[1, DR.no] else 0)
  
  ## 1.2 data frame for regression
  dat.reg <<- data.frame(Unit = as.numeric(dat$otu[1, , DR.no]), 
                         ST = dat$meta[, c("ST.DNA", "ST.RNA")[DR.no]] %>% 
                           unlist %>% as.numeric,
                         id = dat$meta$id, 
                         binaryDNA = ifelse(as.numeric(dat$otu[1, , 1]), 1, 0),
                         phenotype = metadata  %>% unlist %>% 
                           {if (!pheno.as.factor & pheno == "CF") {as.numeric(.) - 1} else {.}}, 
                         batcheffect = dat$meta[[batch]],
                         cdr  = cdr,
                         race = race,
                         agemo = dat$meta$agemo)
}

### Genus-level analysis
genus.transform = function(dat, col = "bacteria", na.rm = TRUE) {
  genera = sub(" .*", "", dat$taxa[, col])
  # species = sub("[^ ]* ", "", dat$taxa)
  genera.unique = genera %>% unique %>% sort
  maps = sapply(genera.unique, function(x) (genera == x))
  n.sample = dim(dat$otu)[2]
  n.genera = length(genera.unique)
  
  # a new copy
  newdat = list(otu = array(NA, dim = c(n.genera, n.sample, 2), 
                            dimnames = c(list(genera.unique), dimnames(dat$otu)[2:3])),
                taxa = setNames(data.frame(genus = genera.unique), nm = col),
                meta = dat$meta)
  
  # fill in aggregated otus.
  for (i in 1:n.genera) {
    newdat$otu[i,,] = apply(dat$otu[maps[,i],,, drop = FALSE], 2:3, sum, na.rm = na.rm)
  }
  
  # if there are additional information other than otu, taxa, and meta, ignore them.
  # prevalence info, screening, etc based on species would not make sense.
  
  newdat
}



dist.cor <- function(x, y, check.y.only = TRUE, pval = FALSE, R = 200) {
  if (is.character(y)) {
    y <- as.factor(y)
  }
  if (is.factor(y)) {
    #lvls = levels(y)
    y = dummies::dummy(y) # categorical into dummies.
  }
  if (!check.y.only) {
    if (is.character(x)) {
      x <- as.factor(x)
    }
    if (is.factor(x)) {
      x = dummies::dummy(x) # categorical into dummies.
    }
  }
  if (pval == TRUE) {
    res <- energy::dcor.test(x, y, R = R)
    return(c(est = res$estimates["dCor"], pval = res$p.value))
  } else {
    energy::dcor(x, y)
  }
}

sign2 <- function(x) ifelse(x > 0 , "+", ifelse(x < 0, "-", "0"))

signif2 = function(x, digits = 2, as.char = FALSE) {
  rnd = round(x, digits)
  sig = signif(x, digits)
  is.sig = abs(rnd - sig) < 10^-(digits+1)
  # if (as.char) val[] <- as.character(val)
  if (as.char) {
    ifelse(abs(x) <= 10^-(digits), 
           sprintf("%.0e", x),
           sprintf("%.2f", x))
  } else {
    ifelse(is.sig, sig, rnd)
  }
}


##### CONSTANTS
  # exclusion (taxa) #For gene "U/M", for bact "U/M" and "U/I", and for path "U/C"
  taxa.exclude <- c("UNMAPPED (TOTAL)", "UNINTEGRATED (TOTAL)", "unclassified", "UniRef90_unknown (TOTAL)")
  virus.keyword <- ".*(virus|virid|viroid|phage).*"
  
  # exclusion (meta): pilot data (2016) and very low-expressed samples
  sample.exclude <- c()
  sample.exclude.RNA <- c(sample.exclude, "subj00000045", "subj000000A7", "subj0000012A")

##### GRAPHICAL CONSTANTS and CONSTANT FUNCTIONS
  basePlot <- function(textsize = 16, titlesize = 19) {
    theme(plot.title = element_text(size=titlesize), 
          text = element_text(size=textsize), legend.text=element_text(size=textsize), legend.title=element_text(size=textsize),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())}

  col.caries <- list(values=c("#F8766D", "#00BA38"), breaks=c("caries", "caries-free"), name="caries status", labels=c("caries  ", "caries-free"))
  shape.caries <- list(values=c(15, 16), breaks=c("caries", "caries-free"), name="caries status", labels=c("caries  ", "caries-free"))
  col.batch.RNA <-list(values=c("#F8766D", "#00BA38", "#619CFF"), breaks=c("170628", "170718", "180530.RNA"), name="sequencing date", labels=c("June 2017", "July 2017", "May 2018"))
  
  col.caries2 <- list(values=c("#F8766D", "#00BA38", "#F8766D", "#00BA38"), 
                      breaks=c("caries_180530", "caries-free_180530", "caries_190812", "caries-free_190812"), 
                      name="caries & batch", 
                      labels=c("caries   (2018)", "caries-free (2018)", "caries   (2019)", "caries-free (2019)"))
  shape.caries2 <- list(values=c(15, 15, 16, 16), 
                        breaks=c("caries_180530", "caries-free_180530", "caries_190812", "caries-free_190812"), 
                        name="caries & batch", 
                        labels=c("caries   (2018)", "caries-free (2018)", "caries   (2019)", "caries-free (2019)"))
  ZOE2.levels <- c("caries_180530", "caries-free_180530", "caries_190812", "caries-free_190812")
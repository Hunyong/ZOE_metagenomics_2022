####################################################################################################
########## Functions - tabulation of otu data for joint (species-genes) analysis
####################################################################################################
library(reshape2)
####################################################################################################


jAnalyze <- function(type, direction = "long", includeReads = FALSE, ...) {
  # ... arguments:  study = 170421, id, type="genefamilies", info.file=humann2 
  # long: Yijk
  # hybrid: Yjk x i (gene-bacteria combination x subject)
  # array: j x k x i (3 dimension)
  data <- reader(..., type = type)
  if (class(data) == "data.frame") { data <- list(data) } else if (class(data) != "list") {
    stop ("The data is neither a data.frame nor a list.")
  }
  id = names(data)

  # joint-analysis by individual
  data <- lapply(data, .jAnalyze.ind, type = type, direction = direction,
                 includeReads = includeReads)

  # Stack up and add id variable
  id <- rep(names(data), sapply(data, nrow))
  data <- as.data.frame(do.call("rbind", data))
  data$id <- id
  if (direction != "long") {
    stop("wide form: TBD")
    # long to wide
  }
  
  # cleansing
  rownames(data) <- NULL                   # remove rownames
  data[,1] <- gsub("\\:.*","", data[,1])  # remove description from gene names
  data[is.na(data)] <- 0                  # replace NA with 0
  return(data)
}

# marginalized analysis - 1. individual integration
# type: genefamilies, metaphlan, pathway, # margin: gene, pathway, bacteria
.jAnalyze.ind <- function(data, type = "genefamilies", direction = "long",
                          includeReads = FALSE) {
  # data: first column = gene (or pathway) + bacteria, second = count (or percentage)
  if (dim(data)[1] == 0) return(data.frame())
  if (direction != "long") {
    stop("wide form: TBD")
  }
  
  if (grepl(type, "metaphlan2.tsv")) {
    type = "methaphlan"
    name = c("bacteria", "category")
    #unit = "percent"
  } else if (grepl(type, "pathabundance.tsv.pathcoverage.tsv")){
    type = "path"
    name = c("path", "bacteria")
    #unit = ifelse (grepl(type, "pathabundance"), "count" , "coverage")
  } else if (type == "bracken") {
    type = "bracken"
    name = c("bacteria", "category")
    #unit = ifelse (grepl("cpm", type), "CPM" , "RPK")
  } else {
    type = "gene"
    name = c("gene", "bacteria")
    #unit = ifelse (grepl("cpm", type), "CPM" , "RPK")
  }
  
  if (!includeReads) {
    name[3] <- names(data)[2] # adding units
    data <- data.frame(a = data[,1], b = NA, d = data[,2])
  } else {
    name = c(name, names(data)[-1])
    data <- cbind(data.frame(a = data[,1], b = NA), d = data[, -1])
  }
  names(data) <- name
  data[,2] <- ifelse(grepl("\\|", data[,1]), gsub(".*\\|","",data[,1]), "(TOTAL)")
  data[,1] <- gsub("\\|.*", "", data[,1])
  
  # removed these lines (this causes inconsistency across subjects) on Oct 9, 2018
  #tmp <- table(data[,1])      # table of counts
  #tmp <- names(tmp)[tmp == 1] # if the counts is 1 (e.g. UNMAPPED), the (TOTAL) should remain, o/w to be removed.
  #data <- data[(data[,1] %in% tmp) | data[,2] != "(TOTAL)", ]
  
  return(data)
}
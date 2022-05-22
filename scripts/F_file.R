####################################################################################################
########## Functions - file
########## functions extracting folder names, file names
####################################################################################################

### Read me: 
  # data.dir, humann2 should be defined. ... -> See code0.0.data.R
  # For outcome object, see code0.0.data.R
    
####################################################################################################
### 2. finder: Getting address of files, given study, id, file type.
  # type: "genefamilies-cpm", "genefamilies", "metaphlan", "pathabundance", "pathcoverage"
  finder <- function (study = study.list, id, type = "genefamilies-cpm", info.file = humann2, include.path = TRUE) {
    # study.list = c("160707", "170421", "170628", "170718", "180227", "180530D", "180530R", "190812D", "190812R")  
    study.list = names(info.file)
      # 160707 added in 201802
      # 180227 added in 201809 (SAMIPRE)
      # 180530D, 180530R added in 201809
      # 190812D, 190812R added in 201908
      # study.list should be matched with the order of humann2!
    study = as.character(study)
    study = match.arg(study)
    study.no = which(study == study.list)    #force into index
    if (grepl("190812", study)) zipped = TRUE else zipped = FALSE
    if (length(study) == 0) stop(paste("No such study. Choose among ", paste(study.list, collapse = ",")))
    if (grepl(type, "genefamilies.")) {type = "genefamilies\\."}    #to prevent to choose genefamilies-cpm
    
    id.set = info.file[[study.no]]$no
    if (id[1] == "all") {id = id.set} else { id  = id.set[id.set %in% id]}
    path <- info.file[[study.no]]$path[id.set %in% id]
    file <- if (type %in%  c("bracken", "quant")) {
      path
    } else if (!zipped) {
      sapply(path, list.files, pattern = type)
    } else {
      sapply(path, function(s) {
        zip.file = gsub(".*\\/(.*)\\.zip", "\\1", s)
        files = unzip(s, list = T)$Name[-1]
        files = files[grepl(type, files)]
        files = gsub(zip.file, "", files)
        files
        })
    }

    no.exist <- (sapply(file, length) == 0)
    if (any(no.exist)) {warning(paste("There is no ", type, " file(s) for ", paste(id[no.exist], collapse=", ")))}
    if (include.path) {file = paste(path, file, sep = ifelse(zipped, "@", "/"))}
    names(file) <- id
    return(file[!no.exist])
  }
    if (FALSE) { # example
      finder(170421, id = 222, type = "genefamilies-cpm", info.file = humann2)
      finder(170421, id = 222, type = "genefamilies-cpm", info.file = humann2, include.path=FALSE)
      finder(170421, id = "all", type = "genefamilies-cpm", info.file = humann2, include.path=FALSE)
      finder(170421, id = "all", type = "genefamilies", info.file = humann2, include.path=FALSE)      
      finder(190812, id = 222, type = "genefamilies-cpm", info.file = humann2)
      finder("190812D", id = "00001FFAD6", type = "genefamilies-cpm", info.file = humann2)
      finder("190812D", id = "all", type = "genefamilies-cpm", info.file = humann2) #zip files
    }
  
  # container code - ID mapper
  .IDmapper.sub <- function(cont.code, mapping.table, mapping.file=NULL, 
                            code.col = "containerCode", id.col = "donor",
                            sep = ",") {
    if (!is.null(mapping.file)) mapping.table = read.csv(mapping.file, header=TRUE, sep = sep)
    row = which(mapping.table[,code.col] %in% cont.code)
    if (length(row) == 0) # if cont.code is not found, try using the numeric version (the Excel file automatically drops leading zeros if all elements are numeric.)
      row = which(mapping.table[,code.col] %in% as.numeric(cont.code))
    id = mapping.table[row, id.col]
    if (length(id) == 0) id = NA
    return(id)
  }
  IDmapper = Vectorize(.IDmapper.sub)
  if (FALSE) {
    IDmapper("ID00001C9E5E", mapping.file="/Users/hycho/Documents/1STAT_Research/201708 Dental/Data-pheno/Vials in BSP_TO_MICROBIOME Wed Mar 16 2016 01 46  05 PM.csv")
    IDmapper(humann2[[1]]$no, mapping.file="/Users/hycho/Documents/1STAT_Research/201708 Dental/Data-pheno/Vials in BSP_TO_MICROBIOME Wed Mar 16 2016 01 46  05 PM.csv")
    IDmapper(c("ID00001C9E5E", 0), mapping.file="/Users/hycho/Documents/1STAT_Research/201708 Dental/Data-pheno/Vials in BSP_TO_MICROBIOME Wed Mar 16 2016 01 46  05 PM.csv")
  }
  
####################################################################################################
### 3. reader: reading files, given study, id, file type.
  reader <- function (study = "170421", id, type = "genefamilies-cpm", info.file = humann2, 
                      colClasses = c("character", "numeric"), 
                      col = if (type == "bracken") c(1,6) else "all", 
                      complete = !(type == "bracken")) {
    files <- finder(study = study, id = id, type = type, info.file = info.file, 
                    include.path = !type %in% c("bracken", "quant"))
    # read.table(path, sep="\t", header=FALSE, skip=1, stringsAsFactors=FALSE) is not correct. See example below
    id <- names(files)
    if (grepl("190812", study)) zipped = TRUE else zipped = FALSE
# tmp.files <<- files
    if (length(files) > 1) {
      data <- lapply (files, .read.delim.col, header=(type %in% c("bracken", "quant")), 
                      skip = ifelse(type  %in% c("bracken", "quant"), 0, 1), 
                      colClasses = colClasses, col = col, type = type, complete = complete, zipped = zipped)
      names(data) <- id
    } else {
      data <- .read.delim.col(files, header=(type == "bracken"), 
                              skip = ifelse(type == "bracken", 0, 1),
                              colClasses = colClasses, col = col, type = type, complete = complete, zipped = zipped)
    }
    return(data)
  }
  # wrapper function of read.delim2 (extracting specific columns only)
  .read.delim.col <- function (fn, ..., col, type = NULL, complete = TRUE, zipped = FALSE) {
    if (zipped) {
      zip.file = gsub("\\@.*", "", fn)
      member.file = gsub(".*\\@", "", fn)
      fn <- unz(zip.file, member.file)
    }
    data <- read.csv(fn, ..., stringsAsFactors=FALSE, sep = "\t")
# print(data)    
    if (complete) {  # this module was added on Oct 9, 2018 to fill the gap between total and sum(items).
      require(magrittr); require(dplyr)
      if (!is.null(type) & type == "quant") {
        names(data) <- c("gene", "length",  "eff.length", "TPM", "reads")
        data = data[, c("gene", "TPM", "reads", "eff.length", "length")]
        # V1 = name, V2 = tpm, V3 = numreads
      } else {
        names(data)[1:2] <- c("V1", "V2")
      
        data$index = seq_len(nrow(data))
        tmp.marginal <- data %>% dplyr::filter(!grepl("\\|", V1))
        tmp.joint    <- data %>% 
          dplyr::filter(grepl("\\|", V1)) %>% 
          dplyr::mutate(V1 = gsub("\\|.*", "", V1)) %>% 
          dplyr::select(-index) %>%
          group_by(V1) %>% 
          summarize(V2 = sum(V2)) %>% 
          ungroup
        tmp <- dplyr::left_join(tmp.marginal, tmp.joint, by = "V1")
        tmp[is.na(tmp$V2.y), "V2.y"] <- 0  # forcing NA into 0
        tmp %<>% 
          transmute(V1 = paste0(V1, "|(GAP)"), V2 = V2.x - V2.y, index = index + .5) %>%
          filter(V2 >= 1)  #throw away small errors
        data <- rbind(data, tmp) %>% arrange(index) %>% dplyr::select(-index)
      }
    }
    
    if (col[1] == "all") {col = 1:dim(data)[2]} else {col = which((1:dim(data)[2]) %in% col)}
    if (grepl(type, "bracken")) {
      names(data) <- gsub("new\\_est\\_reads", "reads", names(data))
    } else if (!grepl(type, "quant")) {
      if (!is.null(type)) {
        if (grepl(type, "genefamilies.")) {name = c("gene.bacteria", "RPK")
        } else if (grepl(type, "genefamilies-cpm.")) {name = c("gene.bacteria", "CPM")
        } else if (grepl(type, "pathabundance.")) {name = c("path.bacteria", "count")
        } else if (grepl(type, "pathabundance.")) {name = c("path.bacteria", "proportion")
        } else if (grepl(type, "metaphlan2.")) {name = c("bacteria", "percentage")
        }
      }
      names(data) <- name
    }
    data <- data[,col]
    data
  }
  
    if (FALSE) { # example
      # read.delim vs read.table
      a <- "/Users/hycho/Documents/1STAT_Research/201708 Dental/Data/170421_UNC31-K00269_0061_AHHTHVBBXX/HUMANN2/00001AE8D7-222/00001AE8D7-222_genefamilies.tsv"
      tmp <- read.delim2(a, header=FALSE, skip=1, stringsAsFactors=FALSE)
      tmp2 <- read.table(a, sep="\t", header=FALSE, skip=1, stringsAsFactors=FALSE); rm(a)
      dim(tmp); dim(tmp2); #actual n = 308885 (linux code: sed -n '$=' '00001AE8D7-222_genefamilies.tsv')
      
      # examples
      a <- reader(170421, id = "all", type="abundance", info.file=humann2) # list of 118
      a <- reader(170421, id = c(222, 235), type="abundance", info.file=humann2) # list of 2
      a <- reader(170421, id = 222, type="genefamilies", info.file=humann2) # a data.frame (not a list)
      a <- reader(170421, id = 222, type="abundance", info.file=humann2) # a data.frame (not a list)
      a <- reader(170421, id = 222, type="abundance", info.file=humann2, col = 1)
      a <- reader(170421, id = 222, type="abundance", info.file=humann2, col = 2)
      a <- reader(170628, id = "all", type="abundance", info.file=humann2, col = 2)
      a <- reader(170628, id = "all", type="abundance", info.file=humann2, col = 1)
        # 170628 - ID 322 - pathabundance file is missing
        # 170628 - ID 420 - pathabundance file is empty
      a1 <- .read.delim.col(fn = "/Users/hycho/Documents/1Research/201708 Dental/Data-bracken/2019-RNA-OUT/10015-0000201A2B-8.bracken.out", 
                            type = "bracken", col = c(1,7), complete = FALSE, zipped = FALSE)
      a1 <- .read.delim.col(fn = "/Users/hycho/Documents/1Research/201708 Dental/Data-bracken/2019-RNA-OUT/10015-0000201A2B-8.bracken.out", 
                            type = "bracken", col = c(1:7), complete = FALSE, zipped = FALSE)
      a <- reader("2016-RNA-OUT", id = "all", type="bracken", info.file=bracken)
    }

####################################################################################################
### 9. useful functions
  if (FALSE) { # example
    list.files(recursive=FALSE)   # All files in current folder        
    list.files(recursive=TRUE)    # All files in current and daughter folders
    list.dirs(recursive=FALSE)    # folders only in current folder
    
    
  }
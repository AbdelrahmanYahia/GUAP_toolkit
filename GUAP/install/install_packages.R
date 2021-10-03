ins_norm <- function(x){
  write(paste0("\033[0;", 33, "m","Installing ",x,"...","\033[0m"), stderr())
  suppressMessages(if(!require(x, quietly=TRUE, character.only = TRUE)){install.packages(x, repos='http://cran.us.r-project.org', quiet = TRUE, verbose = FALSE, character.only = TRUE)})
  write(paste0("\033[0;", 32, "m","Done ",x,"...","\033[0m"), stderr())
}

ins_BcM <- function(x){
  write(paste0("\033[0;", 33, "m","Installing ",x,"...","\033[0m"), stderr())
  suppressMessages(if(!require(x, quietly=TRUE,character.only = TRUE)){BiocManager::install(x, update = FALSE,ask = FALSE,character.only = TRUE)})
  write(paste0("\033[0;", 32, "m","Done ",x,"...","\033[0m"), stderr())
}

packages.nrom <- c( 
  "optparse", 
  "ggplot2", 
  "gridExtra",
  "pheatmap",
  "RColorBrewer"
)

packages.bcm <- c( 
  "dada2",
  "viridis",
  "hrbrthemes",
  "tidyverse",
  "phyloseq",
  "Biostrings",
  "vegan",
  "knitr",
  "DECIPHER",
  "phangorn",
  "ComplexHeatmap",
  "microbiome",
  "eulerr",
  "DESeq2",
  "PCAtools",
  "WGCNA",
  "EnhancedVolcano",
  "DEGreport"
)

package.check.norm <- lapply(
  packages.nrom,FUN = function(x) { 
    ins_norm(x)
  }
)

write(paste0("\033[0;", 33, "m","Installing ","qiime2R","...","\033[0m"), stderr())
suppressMessages(if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools",repos='http://cran.us.r-project.org',quiet = TRUE, verbose = FALSE)})
suppressMessages(devtools::install_github("jbisanz/qiime2R",quiet = TRUE))
write(paste0("\033[0;", 32, "m","Done ","qiime2R","...","\033[0m"), stderr())
write(paste0("\033[0;", 33, "m","Installing ","microbiomeutilities","...","\033[0m"), stderr())
suppressMessages(devtools::install_github("microsud/microbiomeutilities",quiet = TRUE))
write(paste0("\033[0;", 32, "m","Done ","microbiomeutilities","...","\033[0m"), stderr())

package.check.BCM <- lapply(
  packages.bcm,FUN = function(x) { 
  ins_BcM(x)
  }
)

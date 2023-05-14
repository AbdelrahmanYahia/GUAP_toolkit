suppressMessages(library("optparse", quietly = T))
defaultW <- getOption("warn") 
options(warn = -1) 
#############################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
            help="input path", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="output file Directory [default= %default]", metavar="character"),
  make_option(c("--database"),type="character", default="EnsDb",
              help="database to use with AnnotationHub"),
  make_option("--organism", type="character", default="Homo sapeins", metavar="character",
              help="organism to downlaod data for"),
  make_option("--transcript", action="store_true",default=FALSE,
              help="data supplied is Ensemble transcripts")
); 
#############################################
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("No [input] was supplied", call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("No [output] was supplied", call.=FALSE)
}
############################################
# Load required libraries
suppressMessages(library("AnnotationHub", quietly = T))
suppressMessages(library("AnnotationDbi", quietly = T))
suppressMessages(library("ensembldb", quietly = T)) 
suppressMessages(library("data.table", quietly = T))
############################################

# Set file path, database, organism and output file
file <- opt$input
db <- opt$database
org <- opt$organism
out <- opt$output

# Read input counts file
df <- read.csv(file = file, sep = '\t', row.names = 1, check.names = FALSE) 
original_names <- names(df)
# Create an AnnotationHub instance
ah <- AnnotationHub() 

# Query AnnotationHub to get organism and database information
HsEnsDb <- query(ah, c("EnsDb", "Homo sapiens"))[[1]]


# Get annotation information from ensembldb
annotations <- genes(HsEnsDb, return.type = "data.frame") 
if (opt$transcript){
  annotations_transcripts <- transcripts(HsEnsDb, return.type = "data.frame") 
  
  # Map annotation files to input counts data frame
  annot <- annotations_transcripts %>%
    dplyr::select(tx_id, gene_id) %>%
    dplyr::filter(tx_id %in% rownames(df)) 
  
  
  annot2 <- annotations %>%
    dplyr::select(gene_id, gene_name, entrezid) %>%
    dplyr::filter(gene_id %in% annot$gene_id) %>%
    dplyr::right_join(annot, by = c("gene_id" = "gene_id")) 
  
  # Convert row names to a column in the input counts data frame
  df <- df %>%
    tibble::rownames_to_column(var = "tx_id") 
  
  # Join the two dataframes based on the common column tx_id
  merged_df <- merge(df, annot2, by = "tx_id")
  
  # Extract specific columns from annotations and join with input counts data frame
  result <- merged_df %>%
    dplyr::select(tx_id, gene_name) %>%
    dplyr::right_join(df, by = c("tx_id" = "tx_id")) 
  # Check for missing values in gene_name column and fill with gene_id if missing
  result$Gene <- ifelse(is.na(result$gene_name), result$tx_id, result$gene_name) 
}else{
  # Map annotation files to input counts data frame
  annot <- annotations %>%
    dplyr::select(gene_id, gene_name, entrezid) %>%
    dplyr::filter(gene_id %in% rownames(df)) 
  
  # Convert row names to a column in the input counts data frame
  df <- df %>%
    tibble::rownames_to_column(var = "gene_id") 
  
  # Extract specific columns from annotations and join with input counts data frame
  result <- annot %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::right_join(df, by = c("gene_id" = "gene_id")) 

  # Check for missing values in gene_name column and fill with gene_id if missing
  result$Gene <- ifelse(is.na(result$gene_name), result$gene_id, result$gene_name) 
}


# Aggregate results by gene, calculating the mean of other columns
agg_results <- aggregate(result, list(result$Gene), FUN = mean) 
# Set row names of aggregated results to Group.1 column
rownames(agg_results) <- agg_results$Group.1 
# Extract columns 5 to ncol(newdf)-1 from aggregated results and store in newdf
newdf <- agg_results[5:ncol(agg_results)-1] 
final_table <- as.data.frame(lapply(newdf, as.integer))
rownames(final_table) <- rownames(newdf)
colnames(final_table) <- colnames(newdf)
write.table(final_table, paste0(opt$output), quote=FALSE, col.names=TRUE, sep="\t", row.names = TRUE)
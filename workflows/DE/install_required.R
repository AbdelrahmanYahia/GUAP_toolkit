suppressMessages(install.packages("BiocManager", dependencies = T, verbose = F, 
                quiet = T, repos = "https://cloud.r-project.org"))


# if (!require("BiocManager", quietly = TRUE))
#     suppressMessages(install.packages("BiocManager")) # nolint


write(paste0("\033[0;", 33, "m", "installing ",
             "optparse", "\033[0m"), stderr())

tryCatch({
    suppressMessages(install.packages("optparse", dependencies = T, verbose = F, 
                quiet = T, repos = "https://cloud.r-project.org"))
  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error","\033[0m",err), stderr())
  }
)

tryCatch({
  write(paste0("\033[0;", 33, "m", "installing ", 
    "DESeq2", "\033[0m"), stderr())

    suppressMessages(BiocManager::install("DESeq2",ask = F, quiet = TRUE))

  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error","\033[0m",err), stderr())
  }
)

tryCatch({
    write(paste0("\033[0;", 33, "m", "installing ",
            "pheatmap", "\033[0m"), stderr())
    suppressMessages(install.packages("pheatmap", dependencies = T, 
                verbose = F, quiet = T, repos = "https://cloud.r-project.org"))


  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error","\033[0m",err), stderr())
  }
)


tryCatch({
    write(paste0("\033[0;", 33, "m", "installing ",
            "PCAtools", "\033[0m"), stderr())
    suppressMessages(BiocManager::install("PCAtools", ask = F, quiet = TRUE))


  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error","\033[0m",err), stderr())
  }
)

tryCatch({
    write(paste0("\033[0;", 33, "m", "installing ", "WGCNA", "\033[0m"), stderr())
    suppressMessages(BiocManager::install("WGCNA", ask = F, quiet = TRUE))

  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error","\033[0m",err), stderr())
  }
)

tryCatch({
write(paste0("\033[0;", 33, "m", "installing ", "EnhancedVolcano",
                "\033[0m"), stderr())
    suppressMessages(BiocManager::install("EnhancedVolcano", ask = F, quiet = TRUE))
  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error","\033[0m",err), stderr())
  }
)

tryCatch({
write(paste0("\033[0;", 33, "m", "installing ",
            "DEGreport", "\033[0m"), stderr())
    suppressMessages(BiocManager::install("DEGreport", ask = F, quiet = TRUE))


  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error","\033[0m",err), stderr())
  }
)

tryCatch({
write(paste0("\033[0;", 33, "m", "installing ",
                "RColorBrewer", "\033[0m"), stderr())

    suppressMessages(install.packages("RColorBrewer", dependencies = T, verbose = F,
                quiet = T, repos = "https://cloud.r-project.org"))


  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error","\033[0m",err), stderr())
  }
)

tryCatch({
write(paste0("\033[0;", 33, "m", "installing ", 
            "tidyverse", "\033[0m"), stderr())
    suppressMessages(install.packages("tidyverse", dependencies = T, verbose = F,
                quiet = T, repos = "https://cloud.r-project.org"))


  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error","\033[0m",err), stderr())
  }
)

tryCatch({
write(paste0("\033[0;", 33, "m", "installing ", "hrbrthemes",
                 "\033[0m"), stderr())
    suppressMessages(install.packages("hrbrthemes", dependencies = T, verbose = F,
                    quiet = T, repos = "https://cloud.r-project.org"))



  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error","\033[0m",err), stderr())
  }
)

tryCatch({
write(paste0("\033[0;", 33, "m", "installing ", "viridis", "\033[0m"), stderr())
    suppressMessages(install.packages("viridis", dependencies = T, verbose = F,
                quiet = T, repos = "https://cloud.r-project.org"))

  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error","\033[0m",err), stderr())
  }
)
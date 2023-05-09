options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

check_present <- function(pkg){
    if (!pkg %in% installed.packages()[, "Package"])  {
        return (FALSE)
    }else{
        return (TRUE)
    }
}

try_install_package <- function(pkg){
  status <- list(
    out = "",
    error = "",
    warning = "",
    pass = F
  )
    if (check_present(pkg)){
      status$pass <- T
      status$out <- paste0(pkg, "already installed ")
      return(status)
    }else{
      tryCatch(
        expr = {
          out <- capture.output(install.packages(pkg, dependencies = TRUE,verbose = F,quiet = T),type = c("output", "message"),split = F)
          status$out <- out
          if (!check_present(pkg)){
            stop(paste0("Error while installing: ", pkg))
          }
          status$pass <- T
          return(status)
        },
        error = function(error) {
          status$pass <- F
          status$error <- as.character(conditionMessage(error))
          return(status)
        },
        warning = function(warning) {
          status$warning <- as.character(conditionMessage(warning))
          status$pass <- F
          return(status)
        }
      )
    }
}

check_available <- function(pkg){
  av <- suppressMessages(BiocManager::available(pkg))
  return(pkg %in% av)
}

try_install_BiocManager_package <- function(pkg){
  status <- list(
    out = "",
    error = "",
    warning = "",
    pass = F
  )
  if (check_present(pkg)){
    status$pass <- T
    status$out <- paste0(pkg, "already installed ")
    return(status)
  }else{
    tryCatch(
      expr = {
        out <- capture.output(suppressMessages(BiocManager::install(pkg, dependencies = TRUE, quiet = TRUE)),type = c("output", "message"),split = F)
        status$out <- out
        if (!check_present(pkg)){
          stop(paste0("Error while installing: ", pkg))
        }
        status$pass <- T
        return(status)
      },
      error = function(error) {
        status$pass <- F
        status$error <- as.character(conditionMessage(error))
        return(status)
      },
      warning = function(warning) {
        status$warning <- as.character(conditionMessage(warning))
        status$pass <- F
        return(status)
      }
    )
  }
}

return_and_print <- function(m){
  message(m)
  return(m)
}

install_package <- function(pkg){
  run_info <- list(
    pkg = pkg,
    pass = F,
    log = c()
  )
  tryCatch(
    {
      if (!check_present(pkg)) {
          run_info$log <- c(run_info$log, return_and_print(paste0("Installing '", pkg, "' using install.package()")))
          tryCatch(
            expr = {
              out <- try_install_package(pkg)
              run_info$log <- c(run_info$log, out)
              if (!out$pass){
                run_info$log <- c(run_info$log, return_and_print(paste0("Installing '", pkg, "' using install.package() failed")))
                run_info$log <- c(run_info$log, return_and_print(paste0("trying to install '", pkg, "' using BiocManager::install()")))
                out_bioc <- try_install_BiocManager_package(pkg)
                run_info$log <- c(run_info$log, out_bioc)
                if (!out_bioc$pass){
                  if(check_available(pkg)){
                    run_info$log <- c(run_info$log, return_and_print(paste0("Installing '", pkg, "'failed and not found in BiocManager check the log")))
                    run_info$pass <- F
                    return(run_info)
                  }else{
                    run_info$log <- c(run_info$log, return_and_print(paste0("Installing '", pkg, "'failed and not found in BiocManager check the log")))
                    run_info$pass <- F
                    return(run_info)
                  }
                }else{
                  run_info$log <- c(run_info$log, return_and_print(paste0("Installing '", pkg, "'Success")))
                  run_info$pass <- T
                  return(run_info)
                }
              }else{
                run_info$log <- c(run_info$log, return_and_print(paste0("Installing '", pkg, "'Success")))
                run_info$pass <- T
                return(run_info)
              }
            },
            error= function(error) {
              if (!check_present(pkg)) {
                run_info$log <- c(run_info$log, return_and_print(paste0("Installing '", pkg, "'failed with error ", error)))
                run_info$pass <- F
                return(run_info)
              }else{
                run_info$log <- c(run_info$log, return_and_print(paste0("Installing '", pkg, "'seem to success but there was an error ", error)))
                run_info$pass <- T
                return(run_info)
              }
            },
            warning = function(warning) {
              if (!check_present(pkg)) {
                run_info$log <- c(run_info$log, return_and_print(paste0("Installing '", pkg, "'failed with error ", warning)))
                run_info$pass <- F
                return(run_info)
              }else{
                run_info$log <- c(run_info$log, return_and_print(paste0("Installing '", pkg, "'seem to success but there was an error ", warning)))
                run_info$pass <- T
                return(run_info)
          }
        })
      }else{
        run_info$log <- c(run_info$log, return_and_print(paste0( pkg, " already installed")))
        run_info$pass <- T
        return(run_info)
    }
  })
}

install_packagesList <- function(pkgs){
    failed_packages <- c()
    log <- c()
    for (pkg in pkgs){
      one <- install_package(pkg)
      log <- c(log, one$log)
      if (!one$pass){
        failed_packages <- c(failed_packages, one$pkg)
        message(paste0("\033[0;", 31, "m", pkg," Failed","\033[0m"))
      }
    }
    if (length(failed_packages) > 0) {
        message("\033[0;", 31, "m", "Failed installations:","\033[0m")
        message(failed_packages)
    } else {
        message("All packages installed successfully.")
    }
    return(log)
}

packages <- c("DESeq2", "ggplot2", "pheatmap", "RColorBrewer", "PCAtools", "WGCNA", 
              "EnhancedVolcano", "tidyverse", "hrbrthemes", "viridis", "DEGreport", "tools")

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
  "DEGreport",
  "optparse", 
  "ggplot2", 
  "gridExtra",
  "pheatmap",
  "RColorBrewer"
)

# Open the file for redirection
output_file <- "output.log"
file_conn <- file(output_file, open = "wt")

full_log <- install_packagesList(packages)
full_log <- install_packagesList(packages.bcm)

# Redirect output to the file
sink(file_conn, type = "output")
sink(file_conn, type = "message")
print(full_log)
sink(type = "output")
sink(type = "message")
close(file_conn)


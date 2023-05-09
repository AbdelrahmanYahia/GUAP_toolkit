
install_pkg <- function(pkg){
  if (!check_present(pkg)) {
    if (requireNamespace("BiocManager", quietly = TRUE)) {
      message(paste0("Installing '", pkg, "' using BiocManager::install()"))
      tryCatch(
        expr = {
          BiocManager::install(pkg, dependencies = TRUE)
          if (!check_present(pkg)){
            stop(paste0("Error in BiocManager while installing: ", pkg))
          }
        },
        error = function(e) {
          message(paste0("Error while installing ", pkg))
          message(paste0("Installing '", pkg, "' using install.packages()"))
          tryCatch(
            expr = {
              install.packages(pkg, dependencies = TRUE)
              if (!check_present(pkg)){
                stop(paste0("Error while installing: ", pkg))
              }
            },
            error = function(e) {
              failed_packages <- c(failed_packages, pkg)
              warning(paste0("Error while installing: ", pkg))
            }
          )
        }
      )
    } else {
      warning(paste0("Installing 'BiocManager' seemed to fail, will use normal installion instead"))
      message(paste0("Installing '", pkg, "' using install.packages()"))
      tryCatch(
        expr = {
          install.packages(pkg, dependencies = TRUE)
          if (!check_present(pkg)){
            stop(paste0("Error while installing: ", pkg))
          }
        },
        error = function(e) {
          failed_packages <- c(failed_packages, pkg)
          warning(paste0("Error while installing: ", pkg))
        }
      )
    }
  } else {
    message(paste0("Package '", pkg, "' is already installed."))
  }
  return(failed_packages)
}

stop_if_failed <- function(pkg){
  if (!pkg %in% installed.packages()[, "Package"]) {
    stop(paste0("Package '", pkg, "' installed successfully."))
  }
}


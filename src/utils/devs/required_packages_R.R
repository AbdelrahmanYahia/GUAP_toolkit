# List of packages to install
options(repos = c(CRAN = "https://cloud.r-project.org"))

packages <- c("DESeq2", "ggplot2", "pheatmap", "RColorBrewer", "PCAtools", "WGCNA", 
              "EnhancedVolcano", "tidyverse", "hrbrthemes", "viridis", "DEGreport", "tools")

# Create a vector to store failed installations
failed_packages <- c()

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}





########################################
if (!"Deriv" %in% installed.packages()[, "Package"]) {
  stop("Package Deriv not installed successfully.")
}

#__________________#
install_or_fail <- function(package_name){ 

   tryCatch({install.packages(package_name, dependencies = TRUE) 
         library(package_name)}, 
         error = function(e){ print(e) }, 
         warning = function(w){
           catch <-
             grepl("download of package .* failed", w$message) ||
             grepl("(dependenc|package).*(is|are) not available", w$message) ||
             grepl("installation of package.*had non-zero exit status", w$message) ||
             grepl("installation of one or more packages failed", w$message)
           if(catch){ print(w$message)
             stop(paste("installation failed for:",package_name ))}}
         )
}

 #_________________#
# Start writing to an output file
sink('analysis-output.txt')

set.seed(12345)
x <-rnorm(10,10,1)
y <-rnorm(10,11,1)
# Do some stuff here
cat(sprintf("x has %d elements:\n", length(x)))
print(x)
cat("y =", y, "\n")

cat("=============================\n")
cat("T-test between x and y\n")
cat("=============================\n")
t.test(x,y)

# Stop writing to the file
sink()


# Append to the file
sink('analysis-output.txt', append=TRUE)
cat("Some more stuff here...\n")
sink()

#############################################




# Iterate over each package and install it
for (pkg in packages) {
  # Check if the package is already installed
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Try installing using BiocManager::install()
    if (requireNamespace("BiocManager", quietly = TRUE)) {
      message(paste0("Installing '", pkg, "' using BiocManager::install()"))
      tryCatch(
        expr = {
          BiocManager::install(pkg, dependencies = TRUE)
        },
        error = function(e) {
          failed_packages <- c(failed_packages, pkg)  # Store failed package
          message(paste0("Error '", e, "' while installing ", pkg))
          message(paste0("Installing '", pkg, "' using install.packages()"))
          install.packages(pkg, dependencies = TRUE)
        }
      )
    } else {
      failed_packages <- c(failed_packages, pkg)  # Store failed package
      message(paste0("Installing '", pkg, "' using install.packages()"))
      install.packages(pkg, dependencies = TRUE)
    }
  } else {
    message(paste0("Package '", pkg, "' is already installed."))
  }
}

# Print or store failed installations
if (length(failed_packages) > 0) {
  message("Failed installations:")
  message(failed_packages)
} else {
  message("All packages installed successfully.")
}
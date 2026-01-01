# Install BayesPrismExt package

setwd("c:/Users/amirh/Desktop")

cat("Building BayesPrismExt package...\n")

# Try to install the package
tryCatch({
  # First check if devtools is available
  if (!require("devtools", quietly=TRUE)) {
    cat("Installing devtools...\n")
    install.packages("devtools", repos="http://cran.r-project.org", quiet=TRUE)
  }
  
  library(devtools, quietly=TRUE)
  
  cat("Installing from local source...\n")
  devtools::install("BayesPrism", quiet=FALSE, upgrade="never")
  
  cat("\n✓ Package built and installed successfully!\n\n")
  
  cat("Testing package load...\n")
  library(BayesPrismExt)
  cat("✓ Package loaded successfully!\n\n")
  
  cat("Checking available classes:\n")
  tryCatch({
    cl <- getClassDef("BayesPrism")
    if (!is.null(cl)) cat("✓ BayesPrism class found\n")
  }, error=function(e) cat("Note: Could not retrieve class definition\n"))
  
  cat("\nPackage installation test complete!\n")
  
}, error = function(e) {
  cat("\n✗ Error during installation:\n")
  print(e)
  quit(status=1)
})

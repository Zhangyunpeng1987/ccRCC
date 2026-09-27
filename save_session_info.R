# Run this script after loading the packages used by the analysis.
# It records the current R session and the versions of currently loaded packages.

output_dir <- getwd()

session_file <- file.path(output_dir, "sessionInfo_current.txt")
package_file <- file.path(output_dir, "loaded_package_versions_current.csv")

sink(session_file)
cat("Session information recorded after loading the analysis packages.\n\n")
print(sessionInfo())
sink()

loaded_packages <- loadedNamespaces()
package_versions <- data.frame(
  package = sort(loaded_packages),
  version = vapply(sort(loaded_packages), function(pkg) {
    as.character(packageVersion(pkg))
  }, character(1)),
  stringsAsFactors = FALSE
)

write.csv(package_versions, package_file, row.names = FALSE)

message("Saved session information to: ", normalizePath(session_file, winslash = "/", mustWork = FALSE))
message("Saved loaded package versions to: ", normalizePath(package_file, winslash = "/", mustWork = FALSE))

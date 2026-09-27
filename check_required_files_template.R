# Copy paths_template.R to paths.R and edit paths.R before running this check.
# This script checks expected input paths only; it does not run the analysis.

config_file <- if (file.exists("paths.R")) "paths.R" else "paths_template.R"
source(config_file)

check_exists <- function(path, expected_type) {
  exists <- if (expected_type == "directory") {
    dir.exists(path)
  } else {
    file.exists(path)
  }

  data.frame(
    expected_type = expected_type,
    path = path,
    status = if (exists) "FOUND" else "MISSING",
    stringsAsFactors = FALSE
  )
}

checks <- list(
  check_exists(data_root, "directory"),
  check_exists(scrna_root, "directory"),
  check_exists(scatac_root, "directory"),
  check_exists(tcga_root, "directory"),
  check_exists(emtab_root, "directory"),
  do.call(rbind, lapply(scrna_sample_dirs, check_exists, expected_type = "directory")),
  do.call(rbind, lapply(scatac_sample_dirs, check_exists, expected_type = "directory")),
  do.call(rbind, lapply(scatac_peak_files, check_exists, expected_type = "file")),
  do.call(rbind, lapply(scatac_metadata_files, check_exists, expected_type = "file")),
  do.call(rbind, lapply(scatac_fragment_files, check_exists, expected_type = "file")),
  check_exists(tcga_metadata_file, "file"),
  check_exists(tcga_clinical_file, "file"),
  check_exists(tcga_gdc_download_dir, "directory"),
  check_exists(emtab_expression_file, "file"),
  check_exists(emtab_clinical_file, "file")
)

input_check_report <- do.call(rbind, checks)
write.csv(input_check_report, "input_check_report.csv", row.names = FALSE)

missing_items <- input_check_report[input_check_report$status == "MISSING", ]

cat("Input check completed.\n")
cat("Configuration file used: ", config_file, "\n", sep = "")
cat("Total checked paths: ", nrow(input_check_report), "\n", sep = "")
cat("Found paths: ", sum(input_check_report$status == "FOUND"), "\n", sep = "")
cat("Missing paths: ", nrow(missing_items), "\n", sep = "")
cat("Report saved to: input_check_report.csv\n")

if (nrow(missing_items) > 0) {
  cat("\nMissing paths should be reviewed before rerunning the analysis scripts.\n")
}

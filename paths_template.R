# Copy this file to paths.R and edit the paths for your local environment.
# This template is provided for rerun preparation; the main analysis scripts
# still contain section-level local paths that should be updated before rerun.

project_root <- "path/to/project"
data_root <- file.path(project_root, "data")
output_root <- file.path(project_root, "outputs")

scrna_root <- file.path(data_root, "GSE207493", "scRNA")
scatac_root <- file.path(data_root, "GSE207493", "scATAC")
tcga_root <- file.path(data_root, "TCGA_KIRC")
emtab_root <- file.path(data_root, "E-MTAB-1980")

intermediate_root <- file.path(output_root, "intermediate_objects")
figure_root <- file.path(output_root, "figures")
network_root <- file.path(output_root, "regulatory_networks")
survival_root <- file.path(output_root, "survival_modeling")
multinichenet_root <- file.path(output_root, "multinichenet")
infercnv_root <- file.path(output_root, "infercnv")
cicero_root <- file.path(output_root, "cicero")

sample_ids <- c(
  "RCC81", "RCC84", "RCC86", "RCC87", "RCC94", "RCC96", "RCC99",
  "RCC100", "RCC101", "RCC103", "RCC104", "RCC106", "RCC112",
  "RCC113", "RCC114", "RCC115", "RCC116", "RCC119", "RCC120"
)

scrna_sample_dirs <- file.path(scrna_root, sample_ids)
scatac_sample_dirs <- file.path(scatac_root, sample_ids)

scatac_peak_files <- file.path(scatac_sample_dirs, "peaks.bed")
scatac_metadata_files <- file.path(scatac_sample_dirs, "singlecell.csv")
scatac_fragment_files <- file.path(scatac_sample_dirs, "fragments.tsv.gz")

tcga_metadata_file <- file.path(tcga_root, "metadata.cart.json")
tcga_clinical_file <- file.path(tcga_root, "clinical.tsv")
tcga_gdc_download_dir <- file.path(tcga_root, "gdc_download")

emtab_expression_file <- file.path(emtab_root, "ccRCC_exp_log_quantile_normalized.txt")
emtab_clinical_file <- file.path(emtab_root, "clinical_metadata.xlsx")

required_output_dirs <- c(
  output_root,
  intermediate_root,
  figure_root,
  network_root,
  survival_root,
  multinichenet_root,
  infercnv_root,
  cicero_root
)

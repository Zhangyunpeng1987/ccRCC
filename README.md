# Single-Cell Multi-Omics Dissection of Malignant Evolutionary Mechanisms and Construction of a Prognostic Model for Clear Cell Renal Cell Carcinoma

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)

## Research summary

Clear cell renal cell carcinoma (ccRCC) is characterized by pronounced heterogeneity across WHO histological grades, complicating systematic characterization and clinical management. In this study, we comprehensively dissected the dynamic tumor microenvironment (TME) and uncovered underlying grade-transition mechanisms by integrating single-cell RNA and single-cell ATAC sequencing (scRNA-seq and scATAC-seq) data.

Epigenetic alterations were found to consistently precede metabolic reprogramming and invasive adaptation within malignant cells. We identified a core-node-based prognostic signature (CBG) using machine learning ensembles—highlighting key antagonistic survival determinants, namely SLC11A1 and SH3YL1—whose robust performance in predicting patient outcomes was validated across independent clinical cohorts.

Intercellular communication analysis revealed a temporal progression from early inflammation, through vascular remodeling, to an immunosuppression-dominant state. Additionally, key transcription factors, IRF7 and ZNF683, were found to drive the pseudotime trajectory of CD8⁺ T cell exhaustion, while monocyte differentiation toward M1 and M2 macrophages was orchestrated by NFIC/IL1B and CEBPD/GLI2. These multi-omics insights collectively decipher the dynamic immunomodulatory mechanisms of ccRCC progression, providing a robust framework for subtype-specific prognostication and precision therapeutic targeting.

## Repository Structure

The code and processed data are organized into sequential steps corresponding to the analysis workflow in the manuscript.

```text
ccRCC/
├── README.md                 # Project Overview & Instructions
├── LICENSE                   # MIT License
├── Code_summary/             # Analysis Source Code
│   ├── ATAC_code      # Code for all scATAC-seq analysis results
│   ├── RNA_code             # Code for all scRNA-seq analysis results

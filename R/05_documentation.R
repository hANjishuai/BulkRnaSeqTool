#' BulkRNAseqTool: Comprehensive Pipeline for Bulk RNA-Seq Data Analysis
#'
#' @description 
#' The BulkRNAseqTool package provides an end-to-end solution for analyzing bulk 
#' RNA-seq data. It includes functionality for raw count processing, quality control, 
#' differential expression analysis, functional enrichment, and visualization.
#'
#' @details
#' The core of the package is the `BulkRNASeqAnalysis` R6 class, which manages the 
#' entire analysis workflow. Users can either run a predefined pipeline of tasks 
#' or execute individual analysis steps as needed.
#'
#' Key features include:
#' \itemize{
#'   \item Raw count matrix processing and normalization
#'   \item Quality control and sample filtering
#'   \item Differential expression analysis with DESeq2 or edgeR
#'   \item Visualization (PCA, heatmaps, volcano plots, MA plots)
#'   \item Gene Set Enrichment Analysis (GSEA)
#'   \item GO and KEGG pathway enrichment
#'   \item Automated report generation
#' }
#'
#' @section Getting Started:
#' To begin using BulkRNAseqTool:
#' 
#' \preformatted{
#' # Install the package
#' remotes::install_github("yourusername/BulkRNAseqTool")
#' 
#' # Load the package
#' library(BulkRNAseqTool)
#' 
#' # Initialize the analysis engine
#' analysis <- BulkRNASeqAnalysis$new(output_root = "my_analysis")
#' 
#' # Add tasks and run pipeline
#' analysis$add_task("process_counts", process_counts_matrix_task,
#'                  params = list(counts_path = "data/raw_counts.txt"))
#' analysis$run_pipeline()
#' }
#'
#' @section Vignettes:
#' Detailed tutorials are available in the package vignettes:
#' \itemize{
#'   \item \code{vignette("Getting-Started", package = "BulkRNAseqTool")}
#'   \item \code{vignette("Advanced-Analysis", package = "BulkRNAseqTool")}
#' }
#'
#' @keywords internal
"_PACKAGE"

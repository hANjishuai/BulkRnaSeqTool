#!/usr/bin/env Rscript
# BulkRNAseqTool - Full Analysis Example
# Usage: Rscript full_analysis.R

# Load package(change directory where BulkRNAseqTool is)
library(BulkRnaSeqTool)

# Initialize analysis engine
analysis <- BulkRNASeqAnalysis$new(output_root = "results")

# Add tasks to pipeline
analysis$add_task(
  "process_counts",
  BulkRnaSeqTool::process_counts_matrix_task,
  params = list(
    counts_path = "./inst/extdata/sample_data/raw_counts.txt",
    gene_mapping_path = "./inst/extdata/db/g2s_vm36m_gencode.txt"
  )
)
analysis$add_task(
  "filter_genes",
  BulkRnaSeqTool::filter_low_expression_task,
  params = list(min_count = 10, min_samples = 3)
)

analysis$add_task(
  "filter_samples",
  BulkRnaSeqTool::filter_samples_task,
  params = list(sample_pattern = ".*")
)

analysis$add_task(
  "prepare_dge",
  BulkRnaSeqTool::prepare_dge_analysis_task,
  params = list(ctrl_pattern = "C", paired = FALSE)
)

analysis$add_task(
  "qc_visualization",
  BulkRnaSeqTool::perform_qc_visualization_task,
  params = list(output_dir = "01_check")
)

analysis$add_task(
  "dge",
  BulkRnaSeqTool::perform_dge_analysis_task,
  params = list(reference = "Ctrl", alpha = 0.85, fc = 0.1)
)

analysis$add_task(
  "enrichment",
  BulkRnaSeqTool::perform_go_kegg_enrichment_task,
  params = list(
    organism = "mmu", 
    OrgDb = "org.Mm.eg.db",
    alpha = 0.85,
    fc = 0.1)
)

analysis$add_task(
  "volcano",
  BulkRnaSeqTool::generate_volcano_plots_task,
  params = list(alpha = 0.85, fc = 0.1)
)

analysis$add_task(
  "gsea",
  BulkRnaSeqTool::perform_gsea_analysis_task,
  params = list(
    pathway_db = "./inst/extdata/db/selected_pathway.xlsx",
    species = "Mus musculus"
  )
)

analysis$add_task(
  "MAplot",
  BulkRnaSeqTool::perform_MAplot_task,
  params = list(output_dir = "02_DEG")
)

# Execute full pipeline
analysis$run_pipeline()

# 检查输出
cat("\n===== 测试结果摘要 =====\n")
cat("处理基因数量:", nrow(analysis$counts), "\n")
cat("差异分析结果组数:", length(analysis$results$dge_analysis), "\n")
cat("生成的火山图:", length(analysis$results$volcano), "\n")
message("warnings:", warnings())

# Save analysis object
saveRDS(analysis, file.path(analysis$.__enclos_env__$private$.output_root, "analysis_object.rds"))

message("\nAnalysis completed successfully!")
message("Results saved to: ", analysis$.__enclos_env__$private$.output_root)

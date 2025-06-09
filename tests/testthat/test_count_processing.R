# tests/testthat/test_count_processing.R
test_that("计数矩阵处理正确", {
  analyzer <- BulkRnaSeqTool::BulkRNASeqAnalysis$new()
  counts_path <- system.file("extdata/sample_data/raw_counts.txt", package = "BulkRnaSeqTool")
  gene_path <- system.file("extdata/db/g2s_vm36m_gencode.txt", package = "BulkRnaSeqTool")
  
  result <- analyzer$process_counts_matrix(counts_path, gene_path)
  
  expect_s3_class(analyzer$counts, "data.table")
  expect_true("Gene_symbol" %in% names(analyzer$counts))
  expect_gt(nrow(analyzer$counts), 1000)
})

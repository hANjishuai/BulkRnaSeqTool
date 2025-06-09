# tests/testthat/test_dge_workflow.R
test_that("完整的差异表达分析流程", {
  analyzer <- BulkRnaSeqTool::BulkRNASeqAnalysis$new()
  
  # 设置测试数据路径
  counts_path <- system.file("extdata/sample_data/raw_counts.txt", package = "BulkRnaSeqTool")
  gene_path <- system.file("extdata/db/g2s_vm36m_gencode.txt", package = "BulkRnaSeqTool")
  
  # 执行核心流程
  analyzer$process_counts_matrix(counts_path, gene_path)
  analyzer$filter_low_expression(min_count = 5, min_samples = 2)
  analyzer$prepare_dge_analysis(ctrl_pattern = "Ctrl")
  analyzer$perform_dge_analysis()
  
  # 验证结果
  expect_s4_class(analyzer$dds, "DESeqDataSet")
  expect_true("results" %in% names(analyzer))
  expect_true("normalized_counts" %in% names(analyzer))
})

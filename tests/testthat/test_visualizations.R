# tests/testthat/test_visualizations.R
test_that("可视化函数生成有效输出", {
  # 创建测试数据
  test_data <- data.frame(
    log2FoldChange = rnorm(1000),
    padj = runif(1000, 0, 0.1),
    significance = sample(c("up", "down", "ns"), 1000, replace = TRUE)
  )
  
  # 测试火山图
  vplot <- BulkRnaSeqTool:::generate_volcano_plots(test_data, "Test Group", 0.05, 1)
  expect_s3_class(vplot, "ggplot")
  
  # 测试是否生成PDF文件
  test_file <- tempfile(fileext = ".pdf")
  ggsave(test_file, vplot)
  expect_gt(file.size(test_file), 1000)  # 文件应大于1KB
})

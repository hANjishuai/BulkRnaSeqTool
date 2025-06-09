# tests/testthat/test_utility_functions.R
test_that("私有方法正常工作", {
  # 测试目录创建
  test_dir <- tempdir()
  BulkRnaSeqTool:::setup_directories(test_dir)
  expect_true(dir.exists(file.path(test_dir, "00_data")))
  
  # 测试日志记录
  BulkRnaSeqTool:::save_log(test_dir, "test_task", "success", list(param = "value"))
  log_file <- file.path(test_dir, "logs", "analysis_log.csv")
  expect_true(file.exists(log_file))
})

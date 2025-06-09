# tests/testthat/test_class_initialization.R
test_that("类初始化正确", {
  analyzer <- BulkRnaSeqTool::BulkRNASeqAnalysis$new()
  expect_s3_class(analyzer, "BulkRNASeqAnalysis")
  expect_true(dir.exists("results"))
  expect_type(analyzer$params, "list")
})

test_that("任务添加功能正常", {
  analyzer <- BulkRnaSeqTool::BulkRNASeqAnalysis$new()
  analyzer$add_task("test_task", function(self) { return("success") })
  expect_length(analyzer$tasks, 1)
  expect_equal(analyzer$tasks$test_task$func(NULL), "success")
})

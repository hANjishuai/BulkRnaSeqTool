# 分析任务函数 ------------------------------------------------------------

#' @title 处理RNA-seq计数矩阵 - 任务函数
process_counts_matrix_task <- function(self, ...) {
  self$process_counts_matrix(...)
}

#' @title 过滤低表达基因 - 任务函数
filter_low_expression_task <- function(self, ...) {
  self$filter_low_expression(...)
}

#' @title 样本筛选 - 任务函数
filter_samples_task <- function(self, ...) {
  self$filter_samples(...)
}

#' @title 准备差异分析数据 - 任务函数
prepare_dge_analysis_task <- function(self, ...) {
  self$prepare_dge_analysis(...)
}

#' @title 数据质控可视化 - 任务函数
perform_qc_visualization_task <- function(self, ...) {
  self$perform_qc_visualization(...)
}

#' @title 执行差异表达分析 - 任务函数
perform_dge_analysis_task <- function(self, ...) {
  self$perform_dge_analysis(...)
}

#' @title 生成MA图 - 任务函数
perform_MAplot_task <- function(self, ...) {
  self$perform_MAplot(...)
}

#' @title 生成火山图 - 任务函数
generate_volcano_plots_task <- function(self, ...) {
  self$generate_volcano_plots(...)
}

#' @title 执行GSEA分析 - 任务函数
perform_gsea_analysis_task <- function(self, ...) {
  self$perform_gsea_analysis(...)
}

#' @title 执行GO/KEGG富集分析 - 任务函数
perform_go_kegg_enrichment_task <- function(self, ...) {
  self$perform_go_kegg_enrichment(...)
}

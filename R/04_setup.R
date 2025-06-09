#' @import R6
#' @import data.table
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import DESeq2
#' @import edgeR
#' @import msigdbr
#' @import clusterProfiler
#' @import limma
#' @importFrom FactoMineR PCA
#' @importFrom factoextra fviz_pca_ind
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom GseaVis gseaNb
#' @importFrom writexl write_xlsx
#' @importFrom readxl read_xlsx
#' @importFrom jsonlite toJSON
#' @importFrom parallel detectCores
#' @importFrom future plan multisession availableCores
#' @importFrom stringr str_split
#' @importFrom stats dist
#' @importFrom methods as
NULL

# 在加载包时执行的环境设置
.onAttach <- function(libname, pkgname) {
    .libPaths("/home/jifanghan/R/4.4.1/library/")
    packageStartupMessage("欢迎使用BulkRNASeqAnalysis工具箱！")
    packageStartupMessage("版本: 0.1.0 | 最后更新: ", date())
    
    # 检查必要的目录结构
    if (!dir.exists("results")) {
        dir.create("results", recursive = TRUE)
        packageStartupMessage("创建默认输出目录: results/")
    }
    
    # 设置并行计算（如果可用）
    if (requireNamespace("future", quietly = TRUE) && 
        requireNamespace("parallel", quietly = TRUE)) {
        cores <- parallel::detectCores()
        if (cores > 1) {
            future::plan(future::multisession, workers = cores/2 - 1)
            packageStartupMessage("设置并行计算，使用 ", cores/2 - 1, " 个核心")
        }
    }
}

#' Bulk RNA-seq Analysis Pipeline
#'
#' @description
#' An R6 class providing a comprehensive pipeline for bulk RNA-seq data analysis.
#' The class manages data, parameters, tasks, and results throughout the analysis workflow.
#'
#' @details
#' The `BulkRNASeqAnalysis` class encapsulates the entire RNA-seq analysis workflow,
#' from raw count processing to differential expression and functional enrichment.
#' Users can configure analysis parameters, add tasks to a pipeline, and execute
#' the workflow with consistent logging and error handling.
#'
#' @field counts Raw count matrix
#' @field meta_data Sample metadata
#' @field dds DESeq2 dataset object
#' @field results List of analysis results
#' @field params List of analysis parameters
#' @field tasks List of analysis tasks
#' @field normalized_counts Normalized count matrix
#'
#' @section Public Methods:
#' \describe{
#'   \item{\code{initialize(output_root = "results")}}{Initialize the analysis engine}
#'   \item{\code{add_task(task_name, func, params = list())}}{Add an analysis task}
#'   \item{\code{run_pipeline(task_names = NULL)}}{Execute the analysis pipeline}
#'   \item{\code{process_counts_matrix(...)}}{Process raw count matrix}
#'   \item{\code{filter_low_expression(...)}}{Filter low-expression genes}
#'   \item{\code{filter_samples(...)}}{Filter samples based on patterns}
#'   \item{\code{prepare_dge_analysis(...)}}{Prepare for differential expression analysis}
#'   \item{\code{perform_qc_visualization(...)}}{Perform quality control visualization}
#'   \item{\code{perform_dge_analysis(...)}}{Perform differential expression analysis}
#'   \item{\code{generate_volcano_plots(...)}}{Generate volcano plots}
#'   \item{\code{perform_gsea_analysis(...)}}{Perform GSEA analysis}
#'   \item{\code{perform_go_kegg_enrichment(...)}}{Perform GO/KEGG enrichment}
#'   \item{\code{perform_MAplot(...)}}{Generate MA plots}
#' }
#'
#' @import R6
#' @import DESeq2
#' @export
BulkRNASeqAnalysis <- R6::R6Class(
    "BulkRNASeqAnalysis",
    public = list(
        counts = NULL,
        meta_data = NULL,
        dds = NULL,
        results = list(),
        params = list(),
        tasks = list(),
        normalized_counts = NULL,
    
        #' @description Initialize a new BulkRNASeqAnalysis object
        #' @param output_root Root directory for analysis outputs (default: "results")
        initialize = function(output_root = "results") {
            message("初始化Bulk RNA-seq分析引擎 [", Sys.time(), "]")
            private$.output_root <- output_root
            private$setup_directories()
            private$load_config()
        },
    
        #' @description Add an analysis task to the pipeline
        #' @param task_name Unique name for the task
        #' @param func Task function to execute
        #' @param params List of parameters for the task
        add_task = function(task_name, func, params = list()) {
            self$tasks[[task_name]] <- list(
              func = func, 
              params = params
              )
            message("添加任务: ", task_name)
        },
    
        #' @description Execute the analysis pipeline
        #' @param task_names Optional vector of task names to run (default: all tasks)
        run_pipeline = function(task_names = NULL) {
            if (is.null(task_names)) {
              task_names <- names(self$tasks)
            }

            message("\n===== 启动分析流程 =====")
            message("将执行以下任务: ", paste(task_names, collapse = ", "))

            for (task_name in task_names) {
              if (!task_name %in% names(self$tasks)) {
                warning("跳过未知任务: ", task_name)
                next
              }

            task <- self$tasks[[task_name]]
            message("\n[执行任务] ", task_name)

            # 合并默认参数和用户参数
            full_params <- if (!is.null(self$params[[task_name]])) {
              modifyList(self$params[[task_name]], task$params)
            } else {
              task$params
            }

            # 执行任务函数
            # 智能参数传递
            if (".__self__" %in% names(attributes(task$func)) && 
                is(attr(task$func, ".__self__"), "R6")) {
                # 如果是R6方法，不传递self参数
                result <- do.call(task$func, args = full_params)
            } else {
                # 普通函数传递self参数
                result <- do.call(task$func, args = c(list(self = self), full_params))
            }  

            # 保存结果
            self$results[[task_name]] <- result
            message("\n✓ 任务完成: ", task_name)
        }
        message("\n===== 分析流程完成 =====")
        },
    
        #' Process RNA-seq Count Matrix
        #' 
        #' @description
        #' Load and preprocess raw count matrix, including gene ID conversion and
        #' sample name cleaning.
        #'
        #' @param counts_path Path to raw count file
        #' @param gene_mapping_path Path to gene ID mapping file
        #' @param output_dir Output directory for processed data
        #' @param expected_cols Expected columns in count matrix
        #' @param suffix_pattern Pattern to remove from sample names
        #' 
        #' @return Processed count matrix
        process_counts_matrix = function(
            counts_path = "data/raw_counts.txt",
            gene_mapping_path = "data/gene_mapping.txt",
            output_dir = "00_data",
            expected_cols = c("Chr", "Start", "End", "Strand", "Length"),
            suffix_pattern = "_sorted\\.bam$"
            ) {
            message("处理原始计数矩阵...")
            # 检查并创建输出目录
            full_output_dir <- file.path(private$.output_root, output_dir)
            if (!dir.exists(full_output_dir)) {
                dir.create(full_output_dir, recursive = TRUE)
                message("创建输出目录: ", full_output_dir)
            }
            # 1. 读取原始计数数据
            if (!file.exists(counts_path)) {
                stop("计数文件不存在: ", counts_path)
            }
            counts <- data.table::fread(counts_path, header = TRUE)

            # 2. 验证前导列
            actual_cols <- names(counts)[2:6]
            if (!identical(actual_cols, expected_cols)) {
                warning("前导列不匹配!\n预期: ", paste(expected_cols, collapse = ", "),
                        "\n实际: ", paste(actual_cols, collapse = ", "))
                # 保留实际列名继续处理
                expected_cols <- actual_cols
            }

            # 3. 移除前导列
            counts <- counts[, !..expected_cols]

            # 4. 处理样本列名
            sample_cols <- names(counts)[-1]
            new_colnames <- basename(sample_cols)
            new_colnames <- gsub(suffix_pattern, "", new_colnames)
            data.table::setnames(counts, sample_cols, new_colnames)

            # 5. 添加基因符号映射
            if (!file.exists(gene_mapping_path)) {
                stop("基因映射文件不存在: ", gene_mapping_path)
            }
            gene_id2symbol <- data.table::fread(
                gene_mapping_path,
                col.names = c("Geneid", "Gene_symbol"),
                header = FALSE
            )

            # 处理重复的Geneid
            if (any(duplicated(gene_id2symbol$Geneid))) {
                warning("发现重复的Geneid，仅保留每个Geneid的第一个映射")
                gene_id2symbol <- gene_id2symbol[!duplicated(Geneid)]
            }

            counts[, Gene_symbol := gene_id2symbol$Gene_symbol[match(Geneid, gene_id2symbol$Geneid)]]

            # 6. 重新排列列顺序
            data.table::setcolorder(counts, c("Geneid", "Gene_symbol", new_colnames))

            # 7. 保存带符号的矩阵
            symbols_path <- file.path(full_output_dir, "counts_symbols.txt")
            data.table::fwrite(counts, symbols_path)
            message("保存带基因符号的计数矩阵: ", symbols_path)

            # 8. 按基因符号聚合计数
            # 移除Geneid列
            counts[, Geneid := NULL]

            # 聚合计数
            counts_agg <- counts[, lapply(.SD, sum), by = Gene_symbol, .SDcols = new_colnames]

            # 9. 保存聚合后的矩阵
            clean_path <- file.path(full_output_dir, "counts_clean.txt")
            data.table::fwrite(counts_agg, clean_path)
            message("保存聚合后的计数矩阵: ", clean_path)

            # 10. 存储处理后的矩阵
            self$counts <- counts_agg

            # 11. 返回处理信息
            result <- list(
              symbols_matrix = counts,
              aggregated_matrix = counts_agg,
              parameters = list(
                counts_path = counts_path,
                gene_mapping_path = gene_mapping_path,
                output_dir = full_output_dir,
                samples_processed = length(new_colnames)
              )
            )

            # 保存日志
            private$save_log("process_counts_matrix", "success", list(counts_path = counts_path))

            return(result)
        },
    
        #' 质控过滤
        #' @description 过滤低表达基因
        #' @param min_count 最小计数阈值
        #' @param min_samples 最小样本数阈值
        #' @param output_dir 输出目录
        filter_low_expression = function(
            min_count = 10,
            min_samples = 3,
            output_dir = "00_data") {
            message("过滤低表达基因...")
            if (is.null(self$counts)) {
              stop("请先处理计数矩阵")
            }

            full_output_dir <- file.path(private$.output_root, output_dir)
            if (!dir.exists(full_output_dir)) {
              dir.create(full_output_dir, recursive = TRUE)
            }
            # 计算满足条件的基因
            count_data <- as.data.frame(self$counts)
            rownames(count_data) <- count_data$Gene_symbol
            count_data$Gene_symbol <- NULL
            keep_genes <- rowSums(count_data > min_count) >= min_samples
            counts_filt <- count_data[keep_genes, ]

            # 保存过滤后的矩阵
            data.table::fwrite(as.data.table(counts_filt, keep.rownames = "Gene_symbol"), 
            file.path(full_output_dir, "counts_filt.txt"))
            message("过滤后保留 ", sum(keep_genes), "/", nrow(count_data), " 个基因")

            # 更新计数数据
            self$counts <- as.data.table(counts_filt, keep.rownames = "Gene_symbol")

            # 保存日志
            private$save_log("filter_low_expression", "success", 
                list(min_count = min_count, min_samples = min_samples))
            return(counts_filt)
        },

        #' 样本筛选
        #' @description 根据样本命名规则筛选样本
        #' @param sample_pattern 样本名匹配模式
        #' @param exclude_pattern 排除样本模式
        #' @param output_dir 输出目录
        filter_samples = function(
            sample_pattern = ".*",
            exclude_pattern = NULL,
            output_dir = "00_data"
            ) {
            message("筛选样本...")
            
            if (is.null(self$counts)) {
                stop("请先处理计数矩阵")
            }

            full_output_dir <- file.path(private$.output_root, output_dir)
            if (!dir.exists(full_output_dir)) {
                dir.create(full_output_dir, recursive = TRUE)
            }

            sample_cols <- grep(sample_pattern, names(self$counts[,-"Gene_symbol"]), value = TRUE)
            message("\n正在提取包含", sample_pattern, "的样本！")

            if (!is.null(exclude_pattern)) {
                message("\n正在剔除包含", exclude_pattern, "的样本！")
                sample_cols <- sample_cols[!grepl(exclude_pattern, sample_cols)]
            }

            # 筛选样本
            counts_filt <- self$counts[, .SD, .SDcols = c("Gene_symbol", sample_cols)]

            # 保存结果
            output_file <- paste0(
                "counts_filtered_", 
                gsub("[^[:alnum:]]", "_", sample_pattern),
                ifelse(!is.null(exclude_pattern), paste0("_exclude_", exclude_pattern), ""),
                ".txt"
            )

            fwrite(counts_filt, file.path(full_output_dir, output_file))
            message("筛选后保留 ", length(sample_cols), " 个样本")

            # 更新计数数据
            self$counts <- counts_filt

            # 保存日志
            private$save_log("filter_samples", "success", 
                            list(sample_pattern = sample_pattern, 
                            exclude_pattern = exclude_pattern))

            return(counts_filt)
        },
    
        #' 准备差异分析数据
        #' @description 创建样本元数据并准备DESeq2输入
        #' @param ctrl_pattern 对照组的特征
        #' @param split_pattern 命名中的分隔符号
        #' @param paired 是否为成对（逻辑值）
        #' @param output_dir 输出目录
        prepare_dge_analysis = function(
            ctrl_pattern = "C",
            split_pattern = "",
            paired = FALSE,
            output_dir = "02_DEG"
            ) {
            message("准备差异分析数据...")

            if (is.null(self$counts)) {
                stop("请先处理计数矩阵")
            }

            full_output_dir <- file.path(private$.output_root, output_dir)
            if (!dir.exists(full_output_dir)) {
              dir.create(full_output_dir, recursive = TRUE)
            }

            # 提取样本名
            sample_names <- names(self$counts)[-1]  # 排除Gene_symbol列

            # 创建元数据
            meta_data <- data.table(
                sample_id = sample_names,
                Group = ifelse(grepl(ctrl_pattern, sample_names), "Ctrl",
                    paste0("Treat-", toupper(stringr::str_split(
                        sample_names,
                        split_pattern, 
                        simplify = TRUE)[, 1])
                    )
                ),
                pair_id = if (paired) {
                    paste0("Sample-", stringr::str_split(
                        sample_names, 
                        split_pattern, 
                        simplify = TRUE)[, 3]
                    )
                } else {
                    NA
                }
            )
            # 保存元数据
            fwrite(meta_data, file.path(full_output_dir, "metaData.csv"))
            self$meta_data <- meta_data

            # 准备DESeq2输入
            count_data <- as.data.frame(self$counts)
            rownames(count_data) <- count_data$Gene_symbol
            count_data$Gene_symbol <- NULL
      
            dds_data <- list(
                counts = count_data,
                meta_data = meta_data
            )
            # 保存数据
            save(dds_data, file = file.path(full_output_dir, "dds_data.Rdata"))
            message("差异分析数据准备完成，包含 ", nrow(meta_data), " 个样本")

            # 保存日志
            private$save_log("prepare_dge_analysis", "success",
                             list(ctrl_pattern = ctrl_pattern, paired = paired))
            return(dds_data)
        },

        #' 数据质控可视化
        #' @description 执行数据质控并生成可视化报告
        #' @param output_dir 输出目录
        perform_qc_visualization = function(
            output_dir = "01_check"
            ) {
            message("执行数据质控可视化...")

            if (is.null(self$counts) || is.null(self$meta_data)) {
              stop("请先准备计数矩阵和元数据")
            }

            full_output_dir <- file.path(private$.output_root, output_dir)
            if (!dir.exists(full_output_dir)) {
                dir.create(full_output_dir, recursive = TRUE)
            }

            # 获取计数数据和元数据
            count_data <- as.data.frame(self$counts)
            rownames(count_data) <- count_data$Gene_symbol
            count_data$Gene_symbol <- NULL
            meta_data <- self$meta_data

            # 样本排序
            ordered_samples <- meta_data[order(Group), sample_id]
            counts_ordered <- count_data[, ordered_samples]

            # 归一化处理
            cpm_data <- edgeR::cpm(counts_ordered)
            log_cpm <- log2(cpm_data + 1)

             # 样本间距离热图
            sample_dists <- dist(t(log_cpm))
            dist_matrix <- as.matrix(sample_dists)

            p1 <- pheatmap(
                dist_matrix,
                fontsize = 7,
                clustering_distance_rows = sample_dists,
                clustering_distance_cols = sample_dists,
                angle_col = "45",
                show_rownames = TRUE,
                show_colnames = TRUE,
                main = "样本间距离热图"
            )

            # PCA分析
            pca_res <- PCA(t(log_cpm), graph = FALSE)
            p2 <- fviz_pca_ind(
                pca_res,
                title = "主成分分析",
                col.ind = meta_data$Group,
                geom = "text",
                repel = TRUE,
                mean.point = FALSE,
                legend.title = "分组",
                palette = "Dark2"
            )

            # 保存结果
            pdf_file <- file.path(full_output_dir, "qc_report.pdf")
            pdf(pdf_file, width = 12, height = 10)
            print(p1)
            print(p2)
            dev.off()

            message("生成质控报告: ", pdf_file)
            # 保存日志
            private$save_log("perform_qc_visualization", "success", list())
            return(pdf_file)
            # 方法实现保持不变
        },

        #' 执行差异表达分析
        #' @description 使用DESeq2进行差异表达分析
        #' @param reference 指定对照
        #' @param alpha 矫正p值阈值
        #' @param fc 差异倍数阈值
        #' @param output_dir 输出目录
        perform_dge_analysis = function(
            reference = "Ctrl",
            alpha = 0.05,
            fc = 0.5,
            output_dir = "02_DEG"
            ) {
            message("执行差异表达分析...")

            if (is.null(self$counts) || is.null(self$meta_data)) {
                stop("请先准备计数矩阵和元数据")
            }

            full_output_dir <- file.path(private$.output_root, output_dir)
            if (!dir.exists(full_output_dir)) {
                dir.create(full_output_dir, recursive = TRUE)
            }

            # 获取计数数据和元数据
            count_data <- as.data.frame(self$counts)
            count_data <- tibble::column_to_rownames(
                count_data,
                "Gene_symbol")
            count_data$Gene_symbol <- NULL
            meta_data <- self$meta_data

            # 创建基础数据框
            colData = data.frame(
                group = factor(meta_data$Group),
                row.names = meta_data$sample_id
            )

            # 检查是否为配对设计
            if (!all(is.na(meta_data$pair_id))) {
                colData$pair = factor(meta_data$pair_id)
                # 准备DESeq2对象（配对设计）
                dds <- DESeqDataSetFromMatrix(
                        countData = as.matrix(count_data),
                        colData = colData,
                        design = ~ pair + group
                )
                } else {
                # 准备DESeq2对象（非配对设计）
                dds <- DESeqDataSetFromMatrix(
                        countData = as.matrix(count_data),
                        colData = colData,
                        design = ~ group
                  )
                }

            # 在比较之前，进行样本分布检查
            vsdata <- vst(dds, blind = FALSE)
            if (!all(is.na(meta_data$pair_id))) {
              assay(vsdata) <- limma::removeBatchEffect(assay(vsdata), vsdata$pair)
            } 

            p1 <- plotPCA(vsdata, intgroup = "group")

            # 设置对照组
            dds$group <- relevel(dds$group, ref = reference)

            # 执行差异分析
            dds <- DESeq(dds, quiet = TRUE)

            # 保存质控报告
            qc_file1 <- file.path(full_output_dir, "dds_qc_report.pdf")
            pdf(qc_file1, width = 12, height = 10)
            print(p1)
            dev.off()

            qc_file2 <- file.path(full_output_dir, "dds_qc_report2.pdf")
            pdf(qc_file2, width = 12, height = 10)
            par(mar = c(8, 5, 2, 2))
            boxplot(log10(SummarizedExperiment::assays(dds)[["cooks"]]), range = 0, las = 2)
            dev.off()

            # 获取标准化表达矩阵
            if (nrow(meta_data) < 30) {
                normalized_counts <- SummarizedExperiment::assay(rlog(dds))
            } else {
                normalized_counts <- SummarizedExperiment::assay(vst(dds, blind = FALSE))
            }

            # 保存标准化矩阵
            norm_file <- file.path(full_output_dir, "normalized_counts.csv")
            fwrite(as.data.table(normalized_counts, keep.rownames = "Gene"), norm_file)
            self$normalized_counts <- normalized_counts

            # 对每个处理组执行差异分析
            groups <- setdiff(levels(dds$group), reference)
            diff_results <- list()

            for (grp in groups) {
                message("\n正在进行差异分析:", grp, " VS ", reference)
                grp_dir <- file.path(full_output_dir, grp)
                dir.create(grp_dir, showWarnings = FALSE, recursive = TRUE)

                # 提取差异结果
                res <- results(dds, contrast = c("group", grp, reference), alpha = alpha, cooksCutoff = FALSE)
                res_df <- as.data.table(res, keep.rownames = "Gene")

                # 添加差异表达标记
                res_df[, significance := fcase(
                    padj < alpha & log2FoldChange > fc, "up",
                    padj < alpha & log2FoldChange < -fc, "down",
                    default = "ns"
                )]
            # 保存结果
            result_file <- file.path(grp_dir, "dge_results.csv")
            fwrite(res_df, result_file)

            # 保存子集表达矩阵
            grp_samples <- meta_data[Group %in% c(reference, grp), sample_id]
            subset_counts <- normalized_counts[, grp_samples]
            subset_file <- file.path(grp_dir, "normalized_subset.csv")
            fwrite(as.data.table(subset_counts, keep.rownames = "Gene"), subset_file)

            diff_results[[grp]] <- list(
                results = res_df,
                normalized_subset = subset_counts
            )
            }
            message("差异分析完成，结果保存在: ", full_output_dir)
            self$dds <- dds

            # 保存日志
            private$save_log("perform_dge_analysis", "success", 
                list(reference = reference, alpha = alpha, fc = fc))
            return(diff_results)
        },
    
        #' 绘制火山图
        #' @description 为每个比较组绘制差异表达火山图
        #' @param dge_dir 差异分析结果目录
        #' @param alpha padj阈值
        #' @param fc logfc阈值
        #' @param output_dir 输出目录
        generate_volcano_plots = function(
            dge_dir = "02_DEG",
            alpha = 0.05,
            fc = 0.5,
            output_dir = "02_DEG") {
            message("生成火山图...")
            
            full_dge_dir <- file.path(private$.output_root, dge_dir)
            full_output_dir <- file.path(private$.output_root, output_dir)
      
            # 获取所有分组目录
            group_dirs <- list.dirs(full_dge_dir, recursive = FALSE, full.names = TRUE)
            message("\n正在对以下分组制作火山图：", paste(basename(group_dirs), collapse = " "))
      
            plot_files <- list()

            for (grp_dir in group_dirs) {
                grp <- basename(grp_dir)
                res_file <- file.path(grp_dir, "dge_results.csv")

                if (file.exists(res_file)) {
                    dge_res <- fread(res_file)

                    # 计算统计信息
                    up_count <- sum(dge_res$significance == "up", na.rm = TRUE)
                    down_count <- sum(dge_res$significance == "down", na.rm = TRUE)

                    # 创建标题
                    plot_title <- paste(grp, "vs Control")
                    plot_subtitle <- sprintf("UP: %d\nDOWN: %d\n|log2FC| > %0.2f\npadj < %0.2f", 
                                          up_count, down_count, fc, alpha)

                    # 创建火山图
                    volcano_plot <- ggplot(dge_res, aes(x = log2FoldChange, y = -log10(padj))) +
                        geom_point(aes(color = significance), alpha = 0.6, size = 2) +
                        scale_color_manual(values = c("down" = "blue", "ns" = "gray", "up" = "red")) +
                        geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "darkgray") +
                        geom_vline(xintercept = c(-fc, fc), linetype = "dashed", color = "darkgray") +
                        labs(title = plot_title, subtitle = plot_subtitle, 
                             x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
                        theme_bw() +
                        theme(legend.position = "none",
                              plot.title = element_text(face = "bold", size = 14),
                              plot.subtitle = element_text(size = 10))

                    # 保存图像
                    plot_file <- file.path(grp_dir, paste0(grp, "_volcano.pdf"))
                    ggsave(plot_file, volcano_plot, width = 8, height = 6)
                    plot_files[[grp]] <- plot_file
                }
            }
            message("火山图生成完成")
            # 保存日志
            private$save_log("generate_volcano_plots", "success", 
                                list(alpha = alpha, fc = fc))
            return(plot_files)
        },
    
        #' 执行GSEA分析
        #' @description 对差异表达结果进行基因集富集分析
        #' @param dge_dir 差异分析结果目录
        #' @param pathway_db 自选通路数据库路径
        #' @param output_dir 输出目录
        #' @param species 物种
        #' @param gs_cat 条目分类
        #' @param pvaluecut gsea的阈值
        #' @param gs_subcat 条目亚分类
        perform_gsea_analysis = function(
            dge_dir = "02_DEG",
            pathway_db = "db/selected_pathway.xlsx",
            species = "Homo sapiens",
            gs_cat = "C5",
            pvaluecut = 0.25,
            gs_subcat = NULL,
            output_dir = "03_Enrichment") {
            
            message("执行GSEA分析...")
                  
            full_dge_dir <- file.path(private$.output_root, dge_dir)
            full_output_dir <- file.path(private$.output_root, output_dir)

            # 检查通路数据库
            if (!file.exists(pathway_db)) {
                stop("通路数据库文件不存在: ", pathway_db)
            }
            pathway_list <- read_xlsx(pathway_db)$ID

            # 获取所有分组目录
            group_dirs <- list.dirs(full_dge_dir, recursive = FALSE, full.names = TRUE)
            gsea_results <- list()

            for (grp_dir in group_dirs) {
                grp <- basename(grp_dir)
                res_file <- file.path(grp_dir, "dge_results.csv")

                if (file.exists(res_file)) {
                    dge_res <- fread(res_file)

                    # 准备排序基因列表
                    gene_list <- dge_res[!is.na(log2FoldChange), .(Gene, log2FoldChange)]
                    gene_list <- setNames(gene_list$log2FoldChange,  gene_list$Gene)
                    gene_list <- sort(gene_list, decreasing = TRUE) 

                    # 准备基因集
                    go_df <- msigdbr(species = species, category = gs_cat)
                    go_df <- go_df[, c("gs_name", "gene_symbol", "gs_exact_source")] %>% data.table()
                    go_sub <- go_df[gs_exact_source %in% pathway_list]
                    go_sub <- unique(go_sub[, .(gs_name, gene_symbol)])

                    # 执行GSEA    
                    set.seed(666) 
                    gsea_res <- GSEA( 
                      geneList = gene_list,   
                      TERM2GENE = go_sub,
                      verbose = FALSE,
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = pvaluecut
                    )

                    # 保存结果
                    gsea_dir <- file.path(full_output_dir, grp)
                    dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)
                    result_file <- file.path(gsea_dir, "gsea_results.csv")
                    fwrite(as.data.table(gsea_res@result), result_file)

                    # 绘制通路图
                    plot_files <- c()
                    if (nrow(gsea_res@result) > 0) {
                        message("\n分组 ", grp, " 发现 ", nrow(gsea_res@result), " 条显著通路")
                        for (term_id in gsea_res@result$ID) {
                            tryCatch({
                                gsea_plot <- gseaNb(
                                  object = gsea_res,
                                  geneSetID = term_id,
                                  subPlot = 2,
                                  addPval = TRUE,
                                  pvalX = 0.7,
                                  pvalY = 0.8
                                )
                                plot_file <- file.path(gsea_dir, paste0(term_id, "_gsea.pdf"))
                                ggsave(plot_file, gsea_plot, width = 8.5, height = 4.7)
                                plot_files <- c(plot_files, plot_file)
                            }, error = function(e) {
                                message("无法绘制通路 ", term_id, ": ", e$message)
                            })
                        }
                    }

                    gsea_results[[grp]] <- list(
                        result = gsea_res,
                        plots = plot_files
                    )
                }
            }
            message("GSEA分析完成")
            # 保存日志
            private$save_log("perform_gsea_analysis", "success", 
                list(species = species, gs_cat = gs_cat))
            return(gsea_results)
            # 方法实现保持不变
        },
    
        #' 执行GO和KEGG富集分析
        #' @description 对差异表达结果进行基因本体论(GO)和KEGG通路富集分析
        #' @param dge_dir 差异分析结果目录
        #' @param alpha 显著性阈值（调整后p值）
        #' @param fc 差异表达倍数变化阈值（绝对值log2FC）
        #' @param reference 参考组名称
        #' @param dge_file 差异分析结果文件名
        #' @param qv 富集分析q值阈值
        #' @param organism KEGG物种缩写
        #' @param OrgDb OrgDb数据库对象
        #' @param ont GO子本体论类型
        perform_go_kegg_enrichment = function(
            dge_dir = "02_DEG",
            alpha = 0.05,
            fc = 1,
            reference = "Ctrl",
            dge_file = "dge_results.csv",
            qv = 0.05,
            organism = "hsa",
            OrgDb = "org.Hs.eg.db",
            ont = "BP") {

            message("执行GO和KEGG富集分析...")
            full_dge_dir <- file.path(private$.output_root, dge_dir)
            enrichment_results <- list()

            # 辅助函数：执行富集分析并保存结果和图表
            perform_enrichment <- function(
                enrich_func, 
                genes,
                file_prefix, 
                plot_title,
                color_var = "p.adjust",
                ...) {

                # 执行富集分析
                enrich_result <- enrich_func(gene = genes, ...)

                if (is.null(enrich_result) || nrow(enrich_result) == 0) {
                  message("未找到显著富集的通路: ", plot_title)
                  return(NULL)
                }

                # 保存结果
                csv_file <- paste0(file_prefix, ".csv")
                fwrite(as.data.table(enrich_result@result), csv_file)

                # 绘制点图
                plot_file <- paste0(file_prefix, ".pdf")
                dot_plot <- dotplot(
                    enrich_result,
                    x = "GeneRatio",
                    color = color_var,
                    showCategory = 20,
                    title = plot_title
                ) + theme(axis.text.y = element_text(size = 5))

                # 保存图表
                ggsave(plot_file, dot_plot, width = 6.5, height = 6.6)

                list(result = enrich_result, plot = plot_file)
            }

            # 获取所有分组目录
            group_dirs <- list.dirs(full_dge_dir, recursive = FALSE, full.names = TRUE)

            for (grp_dir in group_dirs) {
                grp <- basename(grp_dir)
                res_file <- file.path(grp_dir, dge_file)

                if (!file.exists(res_file)) {
                    warning("文件不存在: ", res_file)
                    next
                }

                message("\n处理分组: ", grp, " vs ", reference)
                dge_res <- fread(res_file)

                # 筛选差异表达基因
                deg_genes <- dge_res[padj < alpha & abs(log2FoldChange) >= fc, Gene]

                if (length(deg_genes) == 0) {
                  warning("在分组 ", grp, " 中未找到满足条件的差异表达基因")
                  next
                }

                # 基因ID转换
                gene_symbol <- tryCatch({
                  bitr(
                    geneID = deg_genes,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = OrgDb
                  )
                }, error = function(e) {
                    stop("基因ID转换失败: ", e$message)
                })

                # 设置文件前缀
                base_prefix <- file.path(grp_dir, paste0(grp, "_VS_", reference))

                # 执行GO富集分析
                go_result <- perform_enrichment(
                    enrich_func = enrichGO,
                    genes = gene_symbol$ENTREZID,
                    file_prefix = paste0(base_prefix, "_go"),
                    plot_title = paste0(grp, " vs ", reference, " - GO Enrichment"),
                    OrgDb = OrgDb,
                    ont = ont,
                    qvalueCutoff = qv,
                    readable = TRUE
                )

                # 执行KEGG富集分析
                kegg_result <- perform_enrichment(
                    enrich_func = enrichKEGG,
                    genes = gene_symbol$ENTREZID,
                    file_prefix = paste0(base_prefix, "_kegg"),
                    plot_title = paste0(grp, " vs ", reference, " - KEGG Pathway Enrichment"),
                    color_var = "qvalue",
                    organism = organism,
                    qvalueCutoff = qv
                )

                enrichment_results[[grp]] <- list(
                    go = go_result,
                    kegg = kegg_result
                )
            }

            message("\nGO和KEGG富集分析完成")
            # 保存日志
            private$save_log("perform_go_kegg_enrichment", "success", 
                            list(organism = organism, ont = ont))
            return(enrichment_results)
            # 方法实现保持不变
        },
    
        #' Generate MA Plots
        #' 
        #' @description
        #' Create MA plots for each comparison group showing log fold change vs mean expression.
        #'
        #' @param output_dir Output directory for plots
        #' 
        #' @return List of plot file paths
        perform_MAplot = function(
            output_dir = "02_DEG"
            ) {
            message("生成MA图...")

            if (is.null(self$dds)) {
              stop("请先执行差异表达分析")
            }

            full_output_dir <- file.path(private$.output_root, output_dir)
            if (!dir.exists(full_output_dir)) {
              dir.create(full_output_dir, recursive = TRUE)
            }

            # 获取所有比较组
            results_names <- resultsNames(self$dds)[-1]  # 去掉Intercept

            plot_files <- list()
            for (res_name in results_names) {
                # 提取结果名称中的组名
                grp <- gsub("group_", "", gsub("_vs_.*", "", res_name))

                # 创建分组目录
                grp_dir <- file.path(full_output_dir, grp)
                dir.create(grp_dir, showWarnings = FALSE, recursive = TRUE)

                # 原始MA图
                ma_file <- file.path(grp_dir, paste0(grp, "_MAplot.pdf"))
                pdf(ma_file, width = 6, height = 8)
                plotMA(results(self$dds, name = res_name), main = paste("MA Plot:", res_name))
                dev.off()

                # 收缩后的MA图
                res_shrink <- lfcShrink(self$dds, coef = res_name, type = "apeglm")
                ma_shrink_file <- file.path(grp_dir, paste0(grp, "_shrink_MAplot.pdf"))
                pdf(ma_shrink_file, width = 6, height = 8)
                plotMA(res_shrink, ylim = c(-10, 10), main = paste("Shrink MA Plot:", res_name))
                dev.off()

                plot_files[[res_name]] <- list(
                  raw_ma = ma_file,
                  shrink_ma = ma_shrink_file
                )
            }

        message("MA图生成完成")
        private$save_log("perform_MAplot", "success", list())
        return(plot_files)
        }
    ),

    private = list(
        # 私有字段
        .output_root = NULL,

        # 设置目录结构
        setup_directories = function() {
            # 调用私有方法文件中的实现
            setup_directories(private$.output_root)
        },

        # 加载默认配置
        load_config = function() {
            # 调用私有方法文件中的实现
            self$params <- load_config()
        },

        # 保存分析日志
        save_log = function(task_name, status, params) {
            # 调用私有方法文件中的实现
            save_log(private$.output_root, task_name, status, params)
        }
    )
)

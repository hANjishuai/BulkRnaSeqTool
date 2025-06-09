#' Setup Analysis Directories
#'
#' @description
#' Create the standard directory structure for RNA-seq analysis outputs.
#'
#' @param output_root Root output directory
#' 
#' @keywords internal
setup_directories <- function(output_root) {
    dirs <- c(
        "00_data", "01_check", "02_DEG", "03_Enrichment",
        "logs", "results", "reports"
    )
    
    for (dir in dirs) {
        full_dir <- file.path(output_root, dir)
        if (!dir.exists(full_dir)) {
            dir.create(full_dir, recursive = TRUE)
            message("创建目录: ", full_dir)
        }
    }
}

#' Load Default Configuration
#'
#' @description
#' Initialize the default analysis parameters for the RNA-seq pipeline.
#'
#' @return List of default parameters
#' 
#' @keywords internal
load_config <- function() {
    params <- list(
        process_counts_matrix = list(
            counts_path = "data/raw_counts.txt",
            gene_mapping_path = "data/gene_mapping.txt",
            output_dir = "00_data"
        ),
        filter_low_expression = list(
            min_count = 10,
            min_samples = 3,
            output_dir = "00_data"
        ),
        filter_samples = list(output_dir = "00_data"),
        prepare_dge_analysis = list(
            ctrl_pattern = "C",
            split_pattern = "",
            paired = FALSE,
            output_dir = "02_DEG"
        ),
        perform_qc_visualization = list(output_dir = "01_check"),
        perform_dge_analysis = list(
            reference = "Ctrl",
            alpha = 0.05,
            fc = 0.5,
            output_dir = "02_DEG"
        ),
        perform_MAplot = list(output_dir = "02_DEG"),
        generate_volcano_plots = list(
            dge_dir = "02_DEG",
            alpha = 0.05,
            fc = 0.5,
            output_dir = "02_DEG"
        ),
        perform_gsea_analysis = list(
            dge_dir = "02_DEG",
            pathway_db = "db/selected_pathway.xlsx",
            species = "Homo sapiens",
            gs_cat = "C5",
            pvaluecut = 0.25,
            output_dir = "03_Enrichment"
        ),
        perform_go_kegg_enrichment = list(
            dge_dir = "02_DEG",
            alpha = 0.05,
            fc = 1,
            reference = "Ctrl",
            dge_file = "dge_results.csv",
            qv = 0.05,
            organism = "hsa",
            OrgDb = "org.Hs.eg.db",
            ont = "BP"
        )
    )
    message("加载默认分析参数")
    return(params)
}

#' Save Analysis Log
#'
#' @description
#' Record analysis steps, status, and parameters in a timestamped log file.
#'
#' @param output_root Root output directory
#' @param task_name Name of the analysis task
#' @param status Task status (e.g., "success", "error")
#' @param params Parameters used for the task
#' 
#' @keywords internal
save_log <- function(output_root, task_name, status, params) {
    log_dir <- file.path(output_root, "logs")
    if (!dir.exists(log_dir)) {
        dir.create(log_dir, recursive = TRUE)
    }

    log_entry <- data.table::data.table(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        task = task_name,
        status = status,
        params = jsonlite::toJSON(params, auto_unbox = TRUE)
    )

    log_file <- file.path(log_dir, paste0("analysis_log_", Sys.Date(), ".csv"))

    if (file.exists(log_file)) {
        data.table::fwrite(log_entry, log_file, append = TRUE)
    } else {
        data.table::fwrite(log_entry, log_file)
    }
}

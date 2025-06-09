# tests/testthat/helper_generate_test_data.R
set.seed(123)
counts <- matrix(rnbinom(200*6, mu=100, size=1/0.5), nrow=200, ncol=6)
colnames(counts) <- c("Control1", "Control2", "Control3", 
                     "Treatment1", "Treatment2", "Treatment3")
gene_ids <- paste0("ENSG", 10000 + 1:200)
write.table(cbind(Geneid=gene_ids, counts), 
            file="/home/jifanghan/R_development1/BulkRnaSeqTool/inst/extdata/sample_data/test_counts.txt", 
            sep="\t", row.names=FALSE)

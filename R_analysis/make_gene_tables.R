# Collapse the per-sample RSEM output into the gene-level tables in data/.
#
# The RSEM files (<run>.genes.results and <run>.isoforms.results for each run,
# about 600 MB) are not part of the repository. This script documents how the
# tables that analysis.R reads were derived from them; it only needs to be run
# again if the quantification changes.
#
# Run from this directory:  Rscript make_gene_tables.R

suppressPackageStartupMessages(library(tximport))

runs <- read.csv("SraRunTable.txt", header = TRUE)$Run
dir.create("data", showWarnings = FALSE)

write_gene_table <- function(x, name){
  con <- gzfile(file.path("data", paste0(name, ".tsv.gz")), "w")
  write.table(data.frame(gene_id = rownames(x), x, check.names = FALSE),
              con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
}

# TPM as reported by RSEM at the gene level.
tpm <- sapply(runs, function(run){
  quant <- read.csv(file.path("rsem", paste0(run, ".genes.results")), sep = "\t", header = TRUE)
  setNames(quant$TPM, quant$gene_id)
})
write_gene_table(tpm, "gene_tpm")

# Expected counts and average transcript lengths per gene, summarised from the
# isoform-level results with tximport. DESeq2 uses both.
iso <- read.delim(file.path("rsem", paste0(runs[1], ".isoforms.results")))
tx2gene <- iso[, c("transcript_id", "gene_id")]
files <- setNames(file.path("rsem", paste0(runs, ".isoforms.results")), runs)
txi <- tximport(files, type = "rsem", tx2gene = tx2gene)
stopifnot(identical(colnames(txi$counts), runs), txi$countsFromAbundance == "no")

write_gene_table(round(txi$counts, 2), "gene_counts")
write_gene_table(round(txi$length, 3), "gene_length")

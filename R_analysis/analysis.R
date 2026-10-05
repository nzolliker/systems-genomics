# Bulk RNA-seq analysis: human forebrain and hindbrain before and after birth
#
# Data: 25 samples from Cardoso-Moreira et al. 2019 (ENA/SRA PRJEB26969),
# quantified with RSEM. The gene-level tables in data/ are derived from the
# RSEM output with make_gene_tables.R.
#
# Run from this directory:  Rscript analysis.R
# Figures are written to figures/, tables to results/.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(pbapply)
  library(DESeq2)
  library(gplots)
  library(viridis)
  library(dendextend)
  library(reshape2)
  library(msigdbr)
  library(fgsea)
  library(gridExtra)
})

set.seed(1)
fig_dir <- "figures"
res_dir <- "results"
dir.create(fig_dir, showWarnings = FALSE)
dir.create(res_dir, showWarnings = FALSE)

stage_levels <- c("16 wpc", "18 wpc", "19 wpc", "neonate", "infant", "toddler")
stage_colors <- c("16 wpc"  = "lightblue",
                  "18 wpc"  = "#0000FF",
                  "19 wpc"  = "darkblue",
                  "neonate" = "lightgoldenrod",
                  "infant"  = "darkorange",
                  "toddler" = "#FF4500")
birth_colors <- c("prenatal" = "orange", "postnatal" = "blue")

theme_big <- theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 15),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))


# ---- 1. Expression (RSEM TPM) and sample metadata ---------------------------

read_gene_table <- function(name)
  as.matrix(read.delim(file.path("data", paste0(name, ".tsv.gz")), row.names = 1, check.names = FALSE))

expr <- read_gene_table("gene_tpm")

meta <- read.csv("SraRunTable.txt", header = TRUE) %>%
  transmute(Run,
            stage = factor(sub(" week post conception,late embryo", " wpc", Developmental_Stage),
                           levels = stage_levels),
            brain_part = factor(Experimental_Factor._organism_part..exp.,
                                levels = c("hindbrain", "forebrain")),
            sex = Experimental_Factor._sex..exp.,
            birth = factor(ifelse(stage %in% c("16 wpc", "18 wpc", "19 wpc"), "prenatal", "postnatal"),
                           levels = c("prenatal", "postnatal")))
rownames(meta) <- meta$Run
stopifnot(!anyNA(meta$stage), !anyNA(meta$brain_part), setequal(meta$Run, colnames(expr)))
expr <- expr[, meta$Run]

print(table(meta$stage, meta$brain_part))


# ---- 2. Gene annotation ------------------------------------------------------

# Stable IDs are taken from the quantification itself; Ensembl BioMart is only
# queried for gene symbols and the result is cached, so reruns do not need the
# network.
meta_genes <- data.frame(ensembl_gene_id_version = rownames(expr),
                         ensembl_gene_id = sub("\\..*$", "", rownames(expr)))

annot_file <- file.path(res_dir, "gene_annotation.tsv")
if (!file.exists(annot_file)) {
  biomart_query <- paste0(
    '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>',
    '<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" datasetConfigVersion="0.6">',
    '<Dataset name="hsapiens_gene_ensembl" interface="default">',
    '<Attribute name="ensembl_gene_id"/><Attribute name="hgnc_symbol"/><Attribute name="chromosome_name"/>',
    '</Dataset></Query>')
  download.file(paste0("https://www.ensembl.org/biomart/martservice?query=",
                       URLencode(biomart_query, reserved = TRUE)),
                annot_file, quiet = TRUE)
}
annot <- read.delim(annot_file, header = FALSE, colClasses = "character", quote = "",
                    col.names = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name")) %>%
  mutate(hgnc_symbol = na_if(hgnc_symbol, "")) %>%
  arrange(is.na(hgnc_symbol)) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)
meta_genes <- meta_genes %>%
  left_join(annot, by = "ensembl_gene_id")
cat("Genes with HGNC symbol:  ", sum(!is.na(meta_genes$hgnc_symbol)), "of", nrow(meta_genes), "\n")


# ---- 3. Filter expressed genes -----------------------------------------------

# Keep genes detected in at least half of the samples or with mean TPM >= 1.
meta_genes$expressed <- rowMeans(expr > 0) >= 0.5 | rowMeans(expr) >= 1

png(file.path(fig_dir, "01_expression_filter.png"), width = 1600, height = 700, res = 150)
layout(matrix(1:2, nrow = 1))
hist(rowSums(expr > 0), breaks = 25, main = "All genes",
     xlab = "Number of samples with TPM > 0")
hist(rowSums(expr[meta_genes$expressed, ] > 0), breaks = 25, main = "Expressed genes",
     xlab = "Number of samples with TPM > 0")
dev.off()

write.table(unique(meta_genes$ensembl_gene_id[meta_genes$expressed]),
            file = file.path(res_dir, "genes_expressed.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)


# ---- 4. Sample similarity: hierarchical clustering ---------------------------

corr_pearson <- cor(log1p(expr[meta_genes$expressed, ]))
corr_spearman <- cor(expr[meta_genes$expressed, ], method = "spearman")
hcl_pearson <- hclust(as.dist(1 - corr_pearson))
hcl_spearman <- hclust(as.dist(1 - corr_spearman))

# Leaf labels have to follow the leaf order of the dendrogram, not the order
# of the samples in the metadata table.
label_dendrogram <- function(hcl, labels, colors){
  dend <- as.dendrogram(hcl)
  ord <- order.dendrogram(dend)
  dend %>%
    set("labels", as.character(labels)[ord]) %>%
    set("labels_colors", colors[ord]) %>%
    set("labels_cex", 1.1)
}

png(file.path(fig_dir, "02_sample_dendrograms.png"), width = 2000, height = 1700, res = 150)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
par(mar = c(7, 4, 3, 1))
plot(label_dendrogram(hcl_pearson, meta$stage, birth_colors[as.character(meta$birth)]),
     main = "Pearson correlation - developmental stage")
legend("topright", legend = names(birth_colors), fill = birth_colors, bty = "n")
plot(label_dendrogram(hcl_spearman, meta$stage, birth_colors[as.character(meta$birth)]),
     main = "Spearman correlation - developmental stage")
part_colors <- c("hindbrain" = "orange", "forebrain" = "blue")
plot(label_dendrogram(hcl_pearson, meta$brain_part, part_colors[as.character(meta$brain_part)]),
     main = "Pearson correlation - brain part")
legend("topright", legend = names(part_colors), fill = part_colors, bty = "n")
plot(label_dendrogram(hcl_spearman, meta$brain_part, part_colors[as.character(meta$brain_part)]),
     main = "Spearman correlation - brain part")
dev.off()


# ---- 5. PCA and highly variable genes ----------------------------------------

plot_pca <- function(pca, title){
  percent_variance <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)
  ggplot(data.frame(pca$x, meta)) +
    geom_point(aes(x = PC1, y = PC2, color = stage, shape = brain_part), size = 5) +
    labs(title = title,
         x = paste0("PC1 (", percent_variance[1], "% variance)"),
         y = paste0("PC2 (", percent_variance[2], "% variance)"),
         color = "Developmental stage",
         shape = "Brain part") +
    scale_color_manual(values = stage_colors) +
    theme_big
}

pca <- prcomp(log1p(t(expr[meta_genes$expressed, ])), center = TRUE, scale. = TRUE)
ggsave(file.path(fig_dir, "03_pca_expressed.png"),
       plot_pca(pca, "PCA - all expressed genes"),
       width = 10, height = 6.5, dpi = 150, bg = "white")

# Genes whose squared coefficient of variation exceeds the fitted
# mean-dependent technical trend (Brennecke et al. 2013).
estimate_variability <- function(expr){
  means <- apply(expr, 1, mean)
  vars <- apply(expr, 1, var)
  cv2 <- vars / means^2

  minMeanForFit <- unname(median(means[which(cv2 > 0.3)]))
  useForFit <- means >= minMeanForFit
  fit <- glm.fit(x = cbind(a0 = 1, a1tilde = 1/means[useForFit]),
                 y = cv2[useForFit],
                 family = Gamma(link = "identity"))
  a0 <- unname(fit$coefficients["a0"])
  a1 <- unname(fit$coefficients["a1tilde"])

  df <- ncol(expr) - 1
  afit <- a1/means + a0
  varFitRatio <- vars/(afit*means^2)
  pval <- pchisq(varFitRatio*df, df = df, lower.tail = FALSE)

  data.frame(mean = means,
             var = vars,
             cv2 = cv2,
             useForFit = useForFit,
             pval = pval,
             padj = p.adjust(pval, method = "BH"),
             row.names = rownames(expr))
}

var_genes <- estimate_variability(expr[meta_genes$expressed, ])
meta_genes$highvar <- meta_genes$ensembl_gene_id_version %in%
  rownames(var_genes)[which(var_genes$padj < 0.01)]

pca_highvar <- prcomp(log1p(t(expr[meta_genes$highvar, ])), center = TRUE, scale. = TRUE)
ggsave(file.path(fig_dir, "04_pca_highvar.png"),
       plot_pca(pca_highvar, "PCA - highly variable genes"),
       width = 10, height = 6.5, dpi = 150, bg = "white")

cat("Genes in quantification: ", nrow(expr), "\n")
cat("Expressed genes:         ", sum(meta_genes$expressed), "\n")
cat("Highly variable genes:   ", sum(meta_genes$highvar), "\n")


# ---- 6. Differential expression across stages: linear model (ANCOVA) --------

# Per gene, an F-test compares log(TPM + 1) ~ cond + covariates against the
# model without cond. The fold change is the ratio between the highest and the
# lowest group mean, so it is always >= 1.
DE_test <- function(expr, cond, covar = NULL, padj_method = "holm"){
  pval_fc <- data.frame(t(pbapply(expr, 1, function(e){
    dat <- data.frame(y = log1p(e), cond = cond)
    if (! is.null(covar))
      dat <- data.frame(dat, covar)

    m1 <- lm(y ~ ., data = dat)
    m0 <- lm(y ~ . - cond, data = dat)
    test <- anova(m1, m0)
    pval <- test$Pr[2]

    avgs <- tapply(log1p(e), cond, mean)
    fc <- exp(max(avgs) - min(avgs))

    c(pval = unname(pval), fc = unname(fc))
  })), row.names = rownames(expr))
  padj <- p.adjust(pval_fc$pval, method = padj_method)
  data.frame(pval_fc, padj = padj)[, c("pval", "padj", "fc")]
}

res_DE <- DE_test(expr = expr[meta_genes$highvar, ],
                  cond = meta$stage,
                  covar = meta %>% dplyr::select(brain_part)) %>%
  rownames_to_column("gene") %>%
  left_join(meta_genes %>% dplyr::select(gene = ensembl_gene_id_version, hgnc_symbol), by = "gene") %>%
  mutate(DE = padj < 0.1 & fc > 2,
         label = ifelse(DE, hgnc_symbol, NA))

p_volcano_ancova <- ggplot(res_DE, aes(x = log(fc), y = -log10(padj), col = DE, label = label)) +
  geom_point() +
  geom_text_repel(max.overlaps = 20, show.legend = FALSE) +
  geom_vline(xintercept = log(2), col = "#303030", linetype = "dotted") +
  geom_hline(yintercept = -log10(0.1), col = "#303030", linetype = "dotted") +
  scale_color_manual(values = c("FALSE" = "#909090", "TRUE" = "red")) +
  labs(title = "Stage effect in highly variable genes (linear model)",
       x = "log(max/min fold change between stages)",
       y = "-log10(adjusted p-value)") +
  theme_big
ggsave(file.path(fig_dir, "05_volcano_ancova.png"), p_volcano_ancova,
       width = 10, height = 6.5, dpi = 150, bg = "white")

cat("ANCOVA DEGs (padj < 0.1, fc > 2): ", sum(res_DE$DE), "\n")


# ---- 7. Differential expression across stages: DESeq2 (LRT) ------------------

# Expected counts and average transcript lengths per gene, as summarised by
# tximport from the RSEM isoform-level results.
txi <- list(counts = read_gene_table("gene_counts")[, meta$Run],
            length = read_gene_table("gene_length")[, meta$Run],
            countsFromAbundance = "no")

dds <- DESeqDataSetFromTximport(txi,
                                colData = meta,
                                design = ~ brain_part + stage)
stopifnot(identical(colnames(dds), meta$Run))
dds <- dds[intersect(rownames(expr)[meta_genes$expressed], rownames(dds)), ]

# Likelihood-ratio test for the stage effect, with brain part as covariate,
# on the same highly variable genes as the linear model.
dds_stage <- dds[intersect(rownames(expr)[meta_genes$highvar], rownames(dds)), ]
dds_stage <- DESeq(dds_stage, test = "LRT", reduced = ~ brain_part)
res_DESeq2 <- results(dds_stage)

# The p-value tests all stages jointly; the reported log2 fold change is the
# last stage coefficient (toddler vs. 16 wpc).
res_df <- as.data.frame(res_DESeq2) %>%
  rownames_to_column("gene") %>%
  left_join(meta_genes %>% dplyr::select(gene = ensembl_gene_id_version, hgnc_symbol), by = "gene") %>%
  mutate(DE = !is.na(padj) & padj < 0.1)

p_volcano_deseq2 <- ggplot(res_df %>% filter(!is.na(padj)),
                           aes(x = log2FoldChange, y = -log10(padj), color = DE)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "#909090", "TRUE" = "red")) +
  labs(title = "Stage effect in highly variable genes (DESeq2 LRT)",
       x = "log2 fold change (toddler vs. 16 wpc)",
       y = "-log10(adjusted p-value)",
       color = "padj < 0.1") +
  theme_big
ggsave(file.path(fig_dir, "06_volcano_deseq2.png"), p_volcano_deseq2,
       width = 10, height = 6.5, dpi = 150, bg = "white")

cat("DESeq2 genes tested:      ", nrow(res_df), "\n")
cat("DESeq2 DEGs (padj < 0.1): ", sum(res_df$DE), "\n")

# Agreement of the two tests on the genes tested by both.
cmp <- inner_join(res_df %>% transmute(gene,
                                       pval_deseq2 = pvalue,
                                       deseq2 = p.adjust(pvalue, method = "bonferroni") < 0.1),
                  res_DE %>% transmute(gene, pval_ancova = pval, ancova = padj < 0.1),
                  by = "gene")
cat("Spearman correlation of p-values (DESeq2 vs. ANCOVA): ",
    round(cor(cmp$pval_deseq2, cmp$pval_ancova, method = "spearman", use = "complete.obs"), 3), "\n")
print(table(DESeq2 = cmp$deseq2, ANCOVA = cmp$ancova, useNA = "ifany"))

write.table(res_df %>% arrange(padj),
            file = file.path(res_dir, "deseq2_stage_lrt.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)


# ---- 8. Grouping DEGs by expression profile across stages --------------------

avg_expr <- sapply(stage_levels, function(s)
  rowMeans(expr[, meta$stage == s, drop = FALSE]))

DEG_list <- res_df$gene[res_df$DE]
avg_expr_DEG <- avg_expr[DEG_list, , drop = FALSE]
corr_DEG <- cor(t(avg_expr_DEG), method = "spearman")
hcl_DEG <- hclust(as.dist(1 - corr_DEG), method = "complete")

n_clusters <- 5
cl_DEG <- cutree(hcl_DEG, k = n_clusters)
cluster_colors <- setNames(scales::hue_pal()(n_clusters), 1:n_clusters)
print(table(cl_DEG))

png(file.path(fig_dir, "07_deg_correlation_heatmap.png"), width = 1500, height = 1500, res = 150)
heatmap.2(corr_DEG, Rowv = as.dendrogram(hcl_DEG), Colv = as.dendrogram(hcl_DEG),
          trace = "none", scale = "none", labRow = NA, labCol = NA, col = viridis,
          ColSideColors = cluster_colors[cl_DEG],
          main = "Spearman correlation between DEGs")
dev.off()

# Per-gene scaled expression across stages, one panel per cluster.
avg_expr_DEG_list <- tapply(names(cl_DEG), cl_DEG, function(x) avg_expr[x, , drop = FALSE])
scaled_expr_DEG_list <- lapply(avg_expr_DEG_list, function(x) t(scale(t(x))))

png(file.path(fig_dir, "08_cluster_boxplots.png"), width = 2000, height = 1300, res = 150)
layout(matrix(1:6, nrow = 2, byrow = TRUE))
par(mar = c(7, 3, 3, 1))
for (cl in 1:n_clusters) {
  boxplot(scaled_expr_DEG_list[[cl]],
          main = paste0("Cluster ", cl, " (", nrow(scaled_expr_DEG_list[[cl]]), " genes)"),
          col = cluster_colors[cl], las = 2)
}
dev.off()

avg_expr_cluster <- aggregate(avg_expr_DEG, by = list(Cluster = cl_DEG), FUN = mean)
avg_expr_melted <- melt(avg_expr_cluster, id.vars = "Cluster")
p_clusters <- ggplot(avg_expr_melted, aes(x = variable, y = value, fill = as.factor(Cluster))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Average expression of DEG clusters across development",
       x = "Developmental stage", y = "Average expression (TPM)", fill = "Cluster") +
  theme_big
ggsave(file.path(fig_dir, "09_cluster_profiles.png"), p_clusters,
       width = 10, height = 6.5, dpi = 150, bg = "white")

# Gene lists per cluster, e.g. for functional annotation with DAVID
# (background: results/genes_expressed.txt).
for (i in 1:n_clusters){
  write.table(meta_genes$ensembl_gene_id[meta_genes$ensembl_gene_id_version %in% names(which(cl_DEG == i))],
              file = file.path(res_dir, paste0("genes_cluster", i, ".txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}


# ---- 9. GSEA: forebrain vs. hindbrain ----------------------------------------

# Likelihood-ratio test for the brain part effect on all expressed genes, with
# stage as covariate. hindbrain is the reference level, so positive scores mean
# higher expression in forebrain.
dds_gsea <- dds
design(dds_gsea) <- ~ stage + brain_part
dds_gsea <- DESeq(dds_gsea, test = "LRT", reduced = ~ stage)
res_part <- as.data.frame(results(dds_gsea)) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(pvalue), !is.na(log2FoldChange))

scores <- setNames(sign(res_part$log2FoldChange) * (-log10(pmax(res_part$pvalue, .Machine$double.xmin))),
                   sub("\\..*$", "", res_part$gene))
scores <- scores[!duplicated(names(scores))]
scores_ordered <- sort(scores, decreasing = TRUE)

# MSigDB C8: cell type signature gene sets. msigdbr >= 10 calls this argument
# `collection` instead of `category`.
genesets_celltype <- msigdbr(species = "Homo sapiens", category = "C8")
genesets_celltype_list <- split(genesets_celltype$ensembl_gene, genesets_celltype$gs_name)

fgsea_celltype <- fgsea(pathways = genesets_celltype_list,
                        stats = scores_ordered,
                        minSize = 15,
                        maxSize = 500,
                        nproc = 1)
fgsea_celltype <- fgsea_celltype[order(fgsea_celltype$NES, decreasing = TRUE), ]
print(head(fgsea_celltype[, 1:7], 10))
print(tail(fgsea_celltype[, 1:7], 10))

data.table::fwrite(fgsea_celltype, file.path(res_dir, "gsea_celltype_forebrain_vs_hindbrain.tsv"),
                   sep = "\t", sep2 = c("", ";", ""))

plot_enrichment_grid <- function(gene_sets){
  plots <- lapply(gene_sets, function(gs)
    plotEnrichment(genesets_celltype_list[[gs]], scores_ordered) +
      labs(title = gs) +
      theme(plot.title = element_text(size = 8)))
  arrangeGrob(grobs = plots, ncol = 3)
}

ggsave(file.path(fig_dir, "10_gsea_forebrain.png"),
       plot_enrichment_grid(head(fgsea_celltype$pathway, 9)),
       width = 12, height = 9, dpi = 150, bg = "white")
ggsave(file.path(fig_dir, "11_gsea_hindbrain.png"),
       plot_enrichment_grid(rev(tail(fgsea_celltype$pathway, 9))),
       width = 12, height = 9, dpi = 150, bg = "white")

writeLines(capture.output(sessionInfo()), file.path(res_dir, "session_info.txt"))

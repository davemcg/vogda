# ============================================================================
# Differential Transcript & Splice Junction Usage: youth vs adult, controlling for study
# Salmon -> SummarizedExperiment -> edgeR (version 4) pipeline
# Integrated with Splice Junction Differential Testing (EiaD Splice Diff)
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(data.table)
  library(SummarizedExperiment)
  library(HDF5Array)
  library(txdbmaker)
  library(GenomicFeatures)
  library(GenomicRanges)
})

# ---------------------------------------------------------
# 1. Set import paths and metadata
# ---------------------------------------------------------
salmon_base <- "/Volumes/OGVFB_BG/EiaD_2026/salmon_quant/"
salmon_dirs <- list.dirs(salmon_base, recursive = FALSE)

sj_h5_path <- "~/data/vogda/eiad/se_h5/se_sj_human.20261007.hdf5/"
gtf_path   <- "/Users/mcgaugheyd/data/vogda/eiad/references/gencode.v48.primary_assembly.annotation.gtf.gz"

nmeta <- data.table::fread('paper/2026/analysis/data/2026_05_28_EiaD.meta.tsv.gz')

nmeta_filter <- nmeta %>%
  filter(sample_accession %in% basename(salmon_dirs)) %>%
  filter(Organism == 'Homo sapiens', Tissue == 'Retina') %>%
  filter(is.na(Perturbation) | Perturbation == 'MGS 1')

nmeta_filter <- nmeta_filter %>%
  mutate(Age = case_when(Source %in% c("Organoid", "ESC") ~ 'Fetal', TRUE ~ Age))

nmeta_filter$study_accession <- factor(nmeta_filter$study_accession)
nmeta_filter$Age <- factor(nmeta_filter$Age, levels = c("Adult", "Fetal"))

# ---------------------------------------------------------
# 2. Read in transcript counts / bootstraps using catchSalmon
# ---------------------------------------------------------
y <- catchSalmon(paste0(salmon_base, nmeta_filter$sample_accession))

# ---------------------------------------------------------
# 3. Extract Transcript Annotations & Build Transcript SE
# ---------------------------------------------------------
anno <- data.table::fread(paste0(salmon_base, nmeta_filter$sample_accession[1], "/quant.sf"))
anno$Gene       <- sapply(anno$Name, function(x) strsplit(x, '\\|')[[1]][6])
anno$GeneID     <- sapply(anno$Name, function(x) strsplit(x, '\\|')[[1]][2])
anno$Transcript <- sapply(anno$Name, function(x) strsplit(x, '\\|')[[1]][1])

se <- SummarizedExperiment(
  assays = list(counts = y$counts / y$annotation$Overdispersion,
                raw_counts = y$counts),
  rowData = DataFrame(
    TXNAME = anno$Transcript,
    GENEID = anno$GeneID,
    Gene   = anno$Gene,
    Length = anno$Length,
    EffectiveLength = anno$EffectiveLength
  ),
  colData = nmeta_filter
)

# Transcript fraction calculation
counts_mat <- assays(se)$raw_counts
eff_len    <- rowData(se)$EffectiveLength

length_adjusted_counts <- counts_mat / eff_len
gene_totals            <- rowsum(length_adjusted_counts, group = rowData(se)$GENEID)
gene_totals_aligned    <- gene_totals[rowData(se)$GENEID, ]

transcript_fractions <- length_adjusted_counts / gene_totals_aligned
transcript_fractions[is.nan(transcript_fractions)] <- 0

assay(se, "tx_fraction") <- transcript_fractions
stopifnot(!anyNA(rowData(se)$GENEID))

# ---------------------------------------------------------
# 4. Transcript Filtering (Expressed & Multi-isoform genes)
# ---------------------------------------------------------
keep_expr <- filterByExpr(se, group = se$Age)
se_filter <- se[keep_expr, ]

tx_per_gene    <- table(rowData(se_filter)$GENEID)
multi_tx_genes <- names(tx_per_gene)[tx_per_gene >= 2]
se_filter      <- se_filter[rowData(se_filter)$GENEID %in% multi_tx_genes, ]

message("Retained ", nrow(se_filter), " transcripts across ",
        length(unique(rowData(se_filter)$GENEID)), " multi-isoform genes")

# ---------------------------------------------------------
# 5. Normalize & Design Matrix (Transcript Analysis)
# ---------------------------------------------------------
dge    <- normLibSizes(se_filter)
design <- model.matrix(~ study_accession + Age, data = colData(se_filter))
colnames(design) <- make.names(colnames(design))
coef_of_interest  <- grep("^Age", colnames(design), value = TRUE)
stopifnot(length(coef_of_interest) == 1L)
message("Testing coefficient: ", coef_of_interest)

qr_d <- qr(design)
stopifnot(qr_d$rank == ncol(design))

# ---------------------------------------------------------
# 6. Transcript Dispersion & QL Fit / DTU Analysis
# ---------------------------------------------------------
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design, robust = TRUE)

sp <- diffSpliceDGE(
  fit,
  coef   = coef_of_interest,
  geneid = "GENEID",
  exonid = "TXNAME"
)

topsplice_simes <- topSpliceDGE(sp, test = "Simes", number = Inf)
topsplice_gene  <- topSpliceDGE(sp, test = "gene",  number = Inf)
topsplice_exon  <- topSpliceDGE(sp, test = "exon",  number = Inf)

# ---------------------------------------------------------
# 7. Transcript-level & Gene-level DE
# ---------------------------------------------------------
qlf_tx <- glmQLFTest(fit, coef = coef_of_interest)
tx_de  <- topTags(qlf_tx, n = Inf, sort.by = "PValue")$table

all_gene_lookup <- rowData(se)$GENEID
gene_counts     <- rowsum(assay(se, 'raw_counts'), group = all_gene_lookup, reorder = FALSE, na.rm = FALSE)
gene_counts     <- gene_counts[!is.na(rownames(gene_counts)), ]

se_g <- SummarizedExperiment(
  assays = list(counts = gene_counts),
  rowData = DataFrame(
    GeneID = rownames(gene_counts),
    Gene   = unique(data.frame(GeneID = rowData(se)$GENEID,
                               Gene   = rowData(se)$Gene))[match(rownames(gene_counts),
                                                                 unique(data.frame(GeneID = rowData(se)$GENEID,
                                                                                   Gene   = rowData(se)$Gene))$GeneID), "Gene"]
  ),
  colData = colData(se)
)

se_g   <- se_g[filterByExpr(se_g, group = se_g$Age), ]
dge_g  <- normLibSizes(se_g)
dge_g  <- estimateDisp(dge_g, design)
fit_g  <- glmQLFit(dge_g, design, robust = TRUE)
qlf_g  <- glmQLFTest(fit_g, coef = coef_of_interest)
gene_de <- topTags(qlf_g, n = Inf, sort.by = "PValue")$table

# ============================================================================
# 8. Splice Junction Data Loading, Filtering & Annotation
# ============================================================================
message("Loading and processing Splice Junction HDF5 SummarizedExperiment...")
sj <- HDF5Array::loadHDF5SummarizedExperiment(sj_h5_path)

sj_mat <- assay(sj, 'counts')
annot_vec <- rowData(sj)$annot

# Cutoff 1: Annotated Junctions (>=5 unique reads in >=10% samples)
ann_keep <- annot_vec == 1L & (rowSums(sj_mat >= 5L) >= 0.10 * ncol(sj_mat))

# Cutoff 2: Novel Junctions (>=5 unique reads in >=25% samples AND seen in >3 studies)
study_seen <- vapply(
  split(seq_len(ncol(sj_mat)), colData(sj)$study_accession),
  function(ix) rowSums(sj_mat[, ix, drop = FALSE] >= 3L) > 0L,
  logical(nrow(sj_mat))
)
n_studies_seen <- rowSums(study_seen)

nov_keep <- annot_vec == 0L &
  (rowSums(sj_mat >= 5L) >= 0.25 * ncol(sj_mat)) &
  (n_studies_seen > 3L)

keep_sj_rows <- ann_keep | nov_keep
sj_filter    <- sj[keep_sj_rows, ]

message("Retained ", sum(annot_vec[keep_sj_rows] == 1L), " annotated + ",
        sum(annot_vec[keep_sj_rows] == 0L), " novel junctions")

# Annotate Junctions with Overlapping Transcripts via GTF
if (file.exists(gtf_path)) {
  message("Annotating splice junctions with overlapping transcript IDs from GTF...")
  txdb_v48 <- txdbmaker::makeTxDbFromGFF(file = gtf_path, format = "gtf")
  tx_introns <- GenomicFeatures::intronsByTranscript(txdb_v48, use.names = TRUE)
  flat_introns <- unlist(tx_introns)
  flat_introns$transcript_id <- names(flat_introns)

  sj_gr <- GRanges(
    seqnames = rowData(sj_filter)$chr,
    ranges   = IRanges(start = rowData(sj_filter)$start, end = rowData(sj_filter)$end),
    strand   = rowData(sj_filter)$strand
  )

  olaps <- GenomicRanges::findOverlaps(sj_gr, flat_introns, ignore.strand = FALSE)
  tx_list <- tapply(
    flat_introns$transcript_id[subjectHits(olaps)],
    queryHits(olaps),
    function(x) paste(unique(x), collapse = ",")
  )

  rowData(sj_filter)$TRANSCRIPT_ID <- NA_character_
  rowData(sj_filter)$TRANSCRIPT_ID[as.numeric(names(tx_list))] <- tx_list
}

# ---------------------------------------------------------
# 9. Splice Junction Subset, DGE Modeling & Differential Testing
# ---------------------------------------------------------
# Filter sample metadata matching the transcript selection logic
sj_sub <- sj_filter[, sj_filter$Tissue == "Retina" &
                      !grepl("MGS 2|MGS 3|MGS 4", sj_filter$Perturbation)]
sj_sub$Age[is.na(sj_sub$Age)] <- "Fetal"

# Remove junctions without gene annotation
sj_sub <- sj_sub[!is.na(rowData(sj_sub)$GENEID) & rowData(sj_sub)$GENEID != "", ]

sj_sub$Age             <- factor(sj_sub$Age, levels = c("Adult", "Fetal"))
sj_sub$study_accession <- factor(tolower(sj_sub$study_accession))

dge_sj <- edgeR::DGEList(
  counts = assay(sj_sub, "counts"),
  genes  = as.data.frame(rowData(sj_sub))
)

design_sj <- model.matrix(~ study_accession + Age, data = colData(sj_sub))
colnames(design_sj) <- make.names(colnames(design_sj))

keep_sj_expr <- edgeR::filterByExpr(dge_sj, design = design_sj)
dge_sj       <- dge_sj[keep_sj_expr, , keep.lib.sizes = FALSE]
dge_sj       <- edgeR::normLibSizes(dge_sj)

dge_sj <- edgeR::estimateDisp(dge_sj, design_sj)
fit_sj <- edgeR::glmQLFit(dge_sj, design_sj, robust = TRUE)

# Differential Splice Junction Testing
sp_sj <- edgeR::diffSpliceDGE(
  fit_sj,
  coef   = coef_of_interest,
  geneid = "GENEID",
  exonid = "jid"
)

sj_topsplice_exon  <- edgeR::topSpliceDGE(sp_sj, test = "exon",  number = Inf)
sj_topsplice_gene  <- edgeR::topSpliceDGE(sp_sj, test = "gene",  number = Inf)
sj_topsplice_simes <- edgeR::topSpliceDGE(sp_sj, test = "Simes", number = Inf)

# ---------------------------------------------------------
# 10. Splice Junction Group-Level Summary Metrics
# ---------------------------------------------------------
group_vec_sj <- colData(sj_sub)$Age

cpm_mat_sj <- edgeR::cpm(dge_sj)
cpm_by_group <- vapply(
  split(seq_along(group_vec_sj), group_vec_sj),
  function(idx) rowMeans(cpm_mat_sj[, idx, drop = FALSE]),
  numeric(nrow(cpm_mat_sj))
) |>
  as.data.frame() |>
  rename_with(~ paste0("mean_cpm_", .x))

raw_by_group <- vapply(
  split(seq_along(group_vec_sj), group_vec_sj),
  function(idx) rowMeans(dge_sj$counts[, idx, drop = FALSE]),
  numeric(nrow(dge_sj$counts))
) |>
  as.data.frame() |>
  rename_with(~ paste0("mean_raw_", .x))

sj_count_summary <- bind_cols(
  jid = dge_sj$genes$jid,
  cpm_by_group,
  raw_by_group
)

sj_topsplice_exon <- sj_topsplice_exon %>%
  left_join(sj_count_summary, by = "jid")

# ---------------------------------------------------------
# 11. Output Generation & Saving
# ---------------------------------------------------------
saveRDS(se, file = 'paper/2026/analysis/data/eiad_retina_adult_fetal.se.20261014.rds')

fdr_threshold <- 0.05

pan_diff_table <- topsplice_gene |>
  select(GENEID, Gene, FDR.ts_gene = FDR) |>
  inner_join(topsplice_simes |> select(GENEID, FDR.ts_simes = FDR), by = 'GENEID') |>
  inner_join(topsplice_exon  |> select(GENEID, TXNAME, FDR.ts_exon = FDR, logFC.ts_exon = logFC), by = 'GENEID') |>
  left_join(gene_de |> select(GENEID = GeneID, logFC.gene_de = logFC, FDR.gene_de = FDR), by = 'GENEID') |>
  left_join(tx_de   |> select(TXNAME, logFC.gtx_de = logFC, FDR.tx_de = FDR), by = 'TXNAME') |>
  mutate(
    Gene_DTU_Sig = (FDR.ts_simes < fdr_threshold) | (FDR.ts_gene < fdr_threshold),
    Category     = case_when(
      FDR.gene_de < fdr_threshold & !Gene_DTU_Sig ~ "DE only (No DTU)",
      FDR.gene_de < fdr_threshold & Gene_DTU_Sig  ~ "Both (DE + DTU)",
      FDR.gene_de >= fdr_threshold & Gene_DTU_Sig ~ "DTU only (No DE)",
      TRUE ~ "None"
    )
  ) |>
  as_tibble()

out <- list(
  pan_diff_table     = pan_diff_table,
  gene_de            = gene_de,
  tx_de              = tx_de,
  topsplice_exon     = topsplice_exon,
  topsplice_gene     = topsplice_gene,
  topsplice_simes    = topsplice_simes,
  sj_topsplice_exon  = sj_topsplice_exon,
  sj_topsplice_gene  = sj_topsplice_gene,
  sj_topsplice_simes = sj_topsplice_simes,
  sj_count_summary   = sj_count_summary
)

saveRDS(out, file = 'paper/2026/analysis/data/eiad_retina_adult_vs_fetal_tx_diff.20261014.rds')
message("Pipeline completed successfully!")

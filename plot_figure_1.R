################################################################################
# plot_figure_1.R
# Expression ~ doubling time: Figure 1, Supplemental Figure 1 and
# Supplemental Table 2.
#
# Runs the doubling-time (DT) analysis for three datasets (CCLE cell lines,
# PDX samples and HCC) and writes ONLY:
#   - Figure 1 (figure_1.pdf / .png) and its standalone panels
#   - Supplemental Figure 1 (supp_figure_1.pdf / .jpg) and its panels
#   - Supplemental Table 2 (sup_tab2.xlsx)
# all directly into OUT_ROOT.
#
# The per-dataset analysis lives in one function, run_fig1_pipeline(cfg); what
# differs between datasets (loader, terminology, symbol harmonisation) is
# declared in the DATASETS list below.
#
# Gene symbols: PDX and HCC symbols are harmonised to approved HGNC symbols
# right after loading (harmonise_symbols, requires HGNChelper), with
# multi-mapping symbols resolved against the CCLE gene set. So RUN_DATASETS
# must list "ccle" first.
#
# Correlation methods:
#   cor_all : pan-cancer pooled Spearman
#   Fisher  : direction-penalised Fisher combination of per-group p-values
#   Deconf  : group-mean-subtracted Spearman
################################################################################

################################################################################
# WORKING DIRECTORY — defaults to the folder this script is in
#
# Every input file is read from, and every output written to, a folder inside
# WORK_DIR. Expected layout:
#   WORK_DIR/input_files/ccle/       Models_DT.csv, DT_Expression_CCLE_Data.csv
#   WORK_DIR/input_files/hcc/        GSE54236_qnorm.tsv,
#                                    GSE54236_clinical_metadata_compact.tsv
#   WORK_DIR/input_files/pdx/        logTPM_matrix_aver.csv,
#                                    PDM_Samples_com_gc.csv
#   WORK_DIR/input_files/pdx/growth_curve_points/
#                                    all_markers.csv, gc.png, points.csv
#   WORK_DIR/figure_1/       outputs (created by the script)
################################################################################

# Leave WORK_DIR as NULL to use the folder this script is saved in. To use a
# different folder, set it here, e.g. WORK_DIR <- "/path/to/project".
WORK_DIR <- NULL

if (is.null(WORK_DIR) || !nzchar(WORK_DIR)) {
  WORK_DIR <- local({
    script_dir <- NULL
    # Rscript plot_figure_x.R  ->  the path is in the --file= argument
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    # (Rscript writes spaces in that path as "~+~", so they are restored first)
    if (length(file_arg) > 0)
      script_dir <- dirname(normalizePath(
        gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE)))
    # source("plot_figure_x.R") (incl. RStudio's Source button)  ->  ofile
    if (is.null(script_dir)) {
      for (i in rev(seq_len(sys.nframe()))) {
        f <- sys.frame(i)$ofile
        if (!is.null(f)) { script_dir <- dirname(normalizePath(f)); break }
      }
    }
    # Lines run from an open file in RStudio  ->  the editor's file path
    if (is.null(script_dir) && requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable()) {
      p <- rstudioapi::getSourceEditorContext()$path
      if (!is.null(p) && nzchar(p)) script_dir <- dirname(normalizePath(p))
    }
    # Anything else: the current R working directory
    if (is.null(script_dir)) {
      script_dir <- getwd()
      message("Could not detect the script's folder; using the current ",
              "working directory instead.")
    }
    script_dir
  })
}
message("WORK_DIR: ", WORK_DIR)

if (!dir.exists(WORK_DIR))
  stop(sprintf("WORK_DIR does not exist: %s — set WORK_DIR at the top of the script.",
               WORK_DIR))

# ── Install missing packages ───────────────────────────────────────────────────
# Installs only the packages that are not already installed; packages that are
# present are left as they are (no updates). CRAN packages via
# install.packages(), Bioconductor packages via BiocManager::install().
CRAN_PKGS <- c("ggplot2", "ggrepel", "patchwork", "msigdbr", "openxlsx",
               "igraph", "ggraph", "ggpubr", "latex2exp", "stringr",
               "tidyr", "scales", "png", "HGNChelper")
BIOC_PKGS <- c("clusterProfiler", "GOSemSim", "GO.db", "org.Hs.eg.db")

is_installed <- function(pkgs)
  vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)

missing_cran <- CRAN_PKGS[!is_installed(CRAN_PKGS)]
if (length(missing_cran) > 0) {
  message("Installing CRAN package(s): ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
}

missing_bioc <- BIOC_PKGS[!is_installed(BIOC_PKGS)]
if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  message("Installing Bioconductor package(s): ",
          paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

still_missing <- c(CRAN_PKGS, BIOC_PKGS)[!is_installed(c(CRAN_PKGS, BIOC_PKGS))]
if (length(still_missing) > 0)
  stop("These packages could not be installed: ",
       paste(still_missing, collapse = ", "),
       ". Install them manually and re-run.")
rm(missing_cran, missing_bioc, still_missing)

loadNamespace("tidyr")            # used as tidyr:: in the HCC loader
# ── Libraries ──────────────────────────────────────────────────────────────────
library(ggplot2);  library(ggrepel);  library(patchwork)
library(clusterProfiler); library(msigdbr); library(openxlsx)
library(igraph);   library(ggraph)
library(ggpubr)                  # stat_compare_means(), stat_cor()
library(latex2exp)               # TeX() axis labels
library(stringr)                 # str_wrap() for GO term names in the GSEA table

MIN_LINEAGE <- 10    # minimum samples/cell lines per group


################################################################################
# MAIN OUTPUT DIRECTORY
#
# Every output is written into OUT_ROOT, a folder inside WORK_DIR. Change the
# folder name to write the outputs somewhere else instead of overwriting them.
################################################################################

OUT_ROOT <- file.path(WORK_DIR, "figure_1")
if (file.exists(file.path(OUT_ROOT, "figure_1.pdf")))
  warning(sprintf("OUT_ROOT already holds a figure_1.pdf, which will be overwritten: %s",
                  OUT_ROOT))
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

################################################################################
# DATASET CONFIGURATIONS — edit paths here
################################################################################

DATASETS <- list(

  # ── 1. CCLE cell lines ──────────────────────────────────────────────────────
  ccle = list(
    tag        = "ccle",
    nice_name  = "CCLE cell lines",

    # Terminology (log messages)
    unit_plural    = "cell lines",
    group_singular = "lineage",
    group_plural   = "lineages",

    # CCLE is the reference gene set, so its symbols are left untouched
    harmonise_symbols = FALSE,

    # Loader must return list(expr_mat, dt_vec, lineage_vec)
    loader = function() {
      Models   <- read.csv(file.path(WORK_DIR, "input_files", "ccle",
                                     "Models_DT.csv"),
                           header = TRUE, row.names = 1)
      rna_data <- read.csv(file.path(WORK_DIR, "input_files", "ccle",
                                     "DT_Expression_CCLE_Data.csv"),
                           row = 2)[, -c(1, 3)]
      rna_data <- subset(rna_data, complete.cases(Avg_2DT))

      if (!"OncotreeLineage" %in% colnames(Models))
        stop("Models_DT.csv has no OncotreeLineage column; found: ",
             paste(colnames(Models), collapse = ", "))

      common_cells <- intersect(rownames(rna_data), rownames(Models))
      message(sprintf("    %d cell lines in RNA and Models (of %d RNA, %d Models).",
                      length(common_cells), nrow(rna_data), nrow(Models)))
      if (length(common_cells) == 0)
        stop("No shared row names between the expression CSV and Models_DT.csv ",
             "— check that both are keyed on the same cell-line ID.")

      rna_aligned    <- rna_data[common_cells, ]
      # drop = FALSE: keeps a data frame even if Models_DT.csv is single-column
      models_aligned <- Models[common_cells, , drop = FALSE]
      gene_cols      <- setdiff(colnames(rna_aligned), "Avg_2DT")

      list(expr_mat    = rna_aligned[, gene_cols],
           dt_vec      = rna_aligned[["Avg_2DT"]],
           lineage_vec = as.character(models_aligned$OncotreeLineage),
           # Cell-line IDs, in the same row order as dt_vec / lineage_vec
           id_vec      = rownames(rna_aligned))
    }
  ),

  # ── 2. PDX samples (two CSVs) ───────────────────────────────────────────────
  pdx = list(
    tag        = "pdx",
    nice_name  = "PDX samples",

    unit_plural    = "samples",
    group_singular = "tumour type",
    group_plural   = "tumour types",

    # Correct PDX symbols to approved HGNC symbols, resolving multi-mapping
    # suggestions against the CCLE gene set. Requires ccle to run first.
    harmonise_symbols = TRUE,
    symbol_ref_from   = "ccle",

    # Second gene kept for the figure_1 scatter panel (MKI67 is always kept)
    focus_gene = "RAN",

    loader = function() {
      RNA_CSV  <- file.path(WORK_DIR, "input_files", "pdx",
                            "logTPM_matrix_aver.csv")
      META_CSV <- file.path(WORK_DIR, "input_files", "pdx",
                            "PDM_Samples_com_gc.csv")

      rna_data <- as.data.frame(read.csv(RNA_CSV, row.names = 1,
                                         check.names = FALSE))
      message(sprintf("    RNA  : %d samples x %d genes",
                      nrow(rna_data), ncol(rna_data)))

      meta_raw <- read.csv(META_CSV, row.names = 1, check.names = FALSE)
      message(sprintf("    Meta : %d samples, columns: %s",
                      nrow(meta_raw), paste(colnames(meta_raw), collapse = ", ")))

      missing_cols <- setdiff(c("dt_days", "OncoTreeCode"), colnames(meta_raw))
      if (length(missing_cols) > 0)
        stop(sprintf("META_CSV is missing required column(s): %s",
                     paste(missing_cols, collapse = ", ")))

      common_samples <- intersect(rownames(rna_data), rownames(meta_raw))
      if (length(common_samples) == 0)
        stop("No shared sample IDs between RNA_CSV and META_CSV row names.")
      message(sprintf("    %d samples shared between RNA and metadata.",
                      length(common_samples)))

      rna_aligned  <- rna_data[common_samples, , drop = FALSE]
      meta_aligned <- meta_raw[common_samples, , drop = FALSE]

      dt_vec <- as.numeric(meta_aligned[["dt_days"]])   # days
      message(sprintf("    DT range: %.1f - %.1f  |  median: %.1f",
                      min(dt_vec, na.rm = TRUE), max(dt_vec, na.rm = TRUE),
                      median(dt_vec, na.rm = TRUE)))

      list(expr_mat    = rna_aligned[, colnames(rna_aligned), drop = FALSE],
           dt_vec      = dt_vec,
           lineage_vec = as.character(meta_aligned[["OncoTreeCode"]]),
           # Sample IDs, in the same row order as dt_vec / lineage_vec
           id_vec      = rownames(meta_aligned))
    }
  ),

  # ── 3. HCC — a SINGLE cancer type ───────────────────────────────────────────
  # Everything that compares groups is skipped automatically for a one-group
  # dataset (see SINGLE-GROUP BEHAVIOUR at the top of run_fig1_pipeline):
  #   - panel a is one violin, with no Kruskal-Wallis p / eta2 in the subtitle
  #   - the Fisher meta-analysis is undefined and is not run
  #   - the deconfounded run is mathematically identical to the pooled run
  #     (subtracting one group mean does not change Spearman ranks), so the
  #     pooled result is reused rather than recomputed
  # The pooled Spearman, its volcano and the gene scatters all run exactly as
  # they do for CCLE/PDX.
  hcc = list(
    tag        = "hcc",
    nice_name  = "HCC",

    unit_plural    = "samples",
    group_singular = "tumour type",
    group_plural   = "tumour types",

    # Minimum n for a group to be analysed. MIN_LINEAGE (10) is the default;
    # set it lower here if the single cancer type has fewer samples than that,
    # otherwise no group is eligible and the pipeline falls back to all samples
    # with a warning.
    min_group_n = MIN_LINEAGE,

    # Correct symbols to approved HGNC symbols against the CCLE gene set.
    # Requires "ccle" to appear before "hcc" in RUN_DATASETS.
    harmonise_symbols = TRUE,
    symbol_ref_from   = "ccle",

    # Expression TSV (genes in rows, samples in columns; transposed on load)
    # and a clinical metadata TSV keyed on the same sample IDs. Returns
    # list(expr_mat, dt_vec, lineage_vec, id_vec), like the other loaders.
    loader = function() {

      RNA_CSV  <- file.path(WORK_DIR, "input_files", "hcc",
                            "GSE54236_qnorm.tsv")
      META_CSV <- file.path(WORK_DIR, "input_files", "hcc",
                            "GSE54236_clinical_metadata_compact.tsv")

      DT_COL      <- "doubling.time..days."   # DT column in META_CSV (days)
      GROUP_COL   <- NA      # cancer-type column; NA = none, use CANCER_TYPE
      CANCER_TYPE <- "HCC"   # label used when GROUP_COL is missing/blank

      rna_data <- as.data.frame(t(read.table(RNA_CSV, sep = "\t",
                                             header = TRUE, row.names = 1,
                                             check.names = FALSE)))
      
      message(sprintf("    RNA  : %d samples x %d genes",
                      nrow(rna_data), ncol(rna_data)))

      meta_raw <- read.table(META_CSV, sep = "\t", header = TRUE, row = 1)[,c("Sample_title","doubling.time..days.")]
      
      # Parse Sample_title into sample type / patient / biopsy. Rows 1-81 are
      # "<type> <sample>", rows 82-161 "<type> <patient> <sample>"; the ranges
      # assume the file's row order. The sample field is then split on "_".
      meta_raw <- rbind(tidyr::separate(meta_raw[c(1:81),],col = (Sample_title), sep = " ", into = c("type", "sample"), remove = FALSE),
                        tidyr::separate(meta_raw[c(82:161),],col = (Sample_title), sep = " ", into = c("type", "patient","sample"), remove = FALSE)[,-c(3)])[,-c(2)]
      
      meta_raw <- tidyr::separate(meta_raw,col = (sample), sep = "_", into = c("type", "patient","biopsy"), remove = FALSE)
      
      # Tumour samples only ("T"); biopsies 047, 127 and 107 are excluded (the
      # same exclusion is applied in plot_figure_3.R and predict_hcc_LOOCV.R).
      meta_raw <- subset(meta_raw, (type == "T") & !(biopsy %in% c("047","127","107")))
      
      message(sprintf("    Meta : %d samples, columns: %s",
                      nrow(meta_raw), paste(colnames(meta_raw), collapse = ", ")))

      if (!DT_COL %in% colnames(meta_raw))
        stop(sprintf("META_CSV has no '%s' column; found: %s",
                     DT_COL, paste(colnames(meta_raw), collapse = ", ")))

      common_samples <- intersect(rownames(rna_data), rownames(meta_raw))
      if (length(common_samples) == 0)
        stop("No shared sample IDs between RNA_CSV and META_CSV row names.")
      message(sprintf("    %d samples shared between RNA and metadata.",
                      length(common_samples)))

      rna_aligned  <- rna_data[common_samples, , drop = FALSE]
      meta_aligned <- meta_raw[common_samples, , drop = FALSE]

      dt_vec <- as.numeric(meta_aligned[[DT_COL]])
      message(sprintf("    DT range: %.1f - %.1f  |  median: %.1f",
                      min(dt_vec, na.rm = TRUE), max(dt_vec, na.rm = TRUE),
                      median(dt_vec, na.rm = TRUE)))

      # One cancer type: use the column if it exists, otherwise label every
      # sample with CANCER_TYPE so the group machinery still has something to
      # key on (and panel a / the scatters get a sensible legend).
      lineage_vec <- if (!is.na(GROUP_COL) && GROUP_COL %in% colnames(meta_aligned))
        as.character(meta_aligned[[GROUP_COL]]) else
          rep(CANCER_TYPE, length(common_samples))
      lineage_vec[is.na(lineage_vec) | lineage_vec == ""] <- CANCER_TYPE

      message(sprintf("    Cancer type(s): %s",
                      paste(sprintf("%s (n=%d)", names(table(lineage_vec)),
                                    as.integer(table(lineage_vec))),
                            collapse = ", ")))

      list(expr_mat    = rna_aligned,
           dt_vec      = dt_vec,
           lineage_vec = lineage_vec,
           id_vec      = rownames(meta_aligned))
    }
  )
)

# Which datasets to run this session (subset to re-run just one).
# Keep "ccle" first: pdx and hcc resolve multi-mapping gene symbols against it.
RUN_DATASETS <- c("ccle", "pdx", "hcc")

################################################################################
# SHARED HELPERS — used for every dataset
################################################################################

# ── harmonise_symbols() ───────────────────────────────────────────────────────
# Maps gene symbols onto approved HGNC symbols using HGNChelper:
#   1. replace unapproved symbols with the suggested approved one
#   2. multi-mapping suggestions ("A /// B") resolve to whichever candidate is
#      present in `reference` (falling back to the first candidate)
#   3. if a resolution collides with another gene, the losing rows revert to
#      their original symbols, so the output stays one-to-one
# Returns the input unchanged (with a message) if HGNChelper isn't installed.
harmonise_symbols <- function(genes, reference = NULL, label = "Gene symbols") {

  if (!requireNamespace("HGNChelper", quietly = TRUE)) {
    message(paste("    NOTE: HGNChelper not installed — symbols left as-is.",
                  "Run install.packages(\"HGNChelper\")."))
    return(genes)
  }
  genes <- as.character(genes)
  chk   <- suppressWarnings(HGNChelper::checkGeneSymbols(genes))

  updated      <- genes
  fix          <- !chk$Approved & !is.na(chk$Suggested.Symbol)
  updated[fix] <- chk$Suggested.Symbol[fix]

  multi    <- grepl("///", updated, fixed = TRUE)
  resolved <- updated
  if (any(multi)) {
    resolved[multi] <- vapply(updated[multi], function(g) {
      cand <- trimws(strsplit(g, "///", fixed = TRUE)[[1]])
      hit  <- if (is.null(reference)) character(0) else cand[cand %in% reference]
      if (length(hit) > 0) hit[1] else cand[1]
    }, character(1), USE.NAMES = FALSE)
  }
  resolved <- trimws(resolved)

  final <- resolved
  for (g in unique(resolved[duplicated(resolved)])) {
    i    <- which(resolved == g)
    keep <- i[genes[i] == g]
    keep <- if (length(keep) > 0) keep[1] else i[1]
    drop <- setdiff(i, keep)
    final[drop] <- genes[drop]
  }

  message(sprintf(
    "    %s: %d updated, %d multi-mapping resolved, %d reverted on collision.",
    label, sum(fix), sum(multi), sum(final != resolved)))

  final
}

# Gene sets of datasets already loaded, so a later dataset can resolve
# multi-mapping symbols against an earlier one (PDX and HCC against CCLE).
.GENE_REF <- new.env(parent = emptyenv())

# ── format_p_value() ──────────────────────────────────────────────────────────
# Formats p-values as plotmath expression strings, e.g.
#   "italic(P) == 1.23 %*% 10^-04"
# Use with parse = TRUE, i.e. annotate("text", parse = TRUE, label = ...) or
# geom_text(aes(label = ...), parse = TRUE). Plain labs(subtitle=) will NOT
# parse it and would show the raw string.
format_p_value <- function(p, threshold = 0.001, scientific = TRUE) {
  to_sci <- function(x, digits = 2) {
    sci <- formatC(x, format = "e", digits = digits)
    gsub("e", " %*% 10^", sci)
  }
  vapply(p, function(pi) {
    if (is.na(pi)) return(NA_character_)   # guard against NA
    if (!scientific) {
      if (pi >= threshold) {
        paste0("italic(P) == ", signif(pi, 3))
      } else {
        paste0("italic(P) == ", pi)
      }
    } else {
      if (pi >= threshold) {
        paste0("italic(P) == ", to_sci(pi))
      } else {
        paste0("italic(P) == ", to_sci(pi, digits = 0))
      }
    }
  }, character(1))
}

# ── fisher_dir_p() ────────────────────────────────────────────────────────────
# Direction penalty used by the Fisher meta-analysis in run_fig1_pipeline().
#
# The two-sided p is kept when the lineage agrees with the pooled sign and
# replaced by 1 - p when it disagrees. Note that this makes a flat lineage
# that leans the wrong way (p ~ 1) count as strong support, and an exact p = 1
# with the opposite sign becomes 1e-15 after flooring.
FISHER_METHOD <- "p if agree, 1 - p if disagree"
fisher_dir_p <- function(p, agree) if (agree) p else 1 - p

# ── MSigDB GOBP gene sets — fetched once, reused for every dataset ─────────────
.MSIG_CACHE <- new.env(parent = emptyenv())

get_msig_gobp <- function() {
  if (!is.null(.MSIG_CACHE$gobp) && !is.null(.MSIG_CACHE$gobp_go_map))
    return(.MSIG_CACHE$gobp)
  message("\n>>> Fetching GOBP gene sets (MSigDB)...")
  m <- tryCatch(
    # msigdbr >= 10 argument names
    msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"),
    error = function(e)
      # msigdbr < 10 argument names
      msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"))
  # Gene-set name -> GO ID map, kept from the same download for panel h
  # (msigdbr carries the GO ID in gs_exact_source; GOBP_* names are not GO IDs)
  src_col <- intersect(c("gs_exact_source", "gs_source_id"), colnames(m))[1]
  go_map  <- unique(as.data.frame(m)[, c("gs_name", src_col)])
  names(go_map) <- c("gs_name", "go_id")
  .MSIG_CACHE$gobp_go_map <- go_map

  m <- m[, c("gs_name", "gene_symbol")]
  .MSIG_CACHE$gobp <- m
  message(sprintf("    %d gene-set/gene pairs loaded.", nrow(m)))
  m
}

# Gene-set name -> GO ID, from the same cached download as get_msig_gobp()
get_msig_gobp_go_map <- function() {
  if (is.null(.MSIG_CACHE$gobp_go_map)) get_msig_gobp()
  .MSIG_CACHE$gobp_go_map
}

# ── run_spearman() ─────────────────────────────────────────────────────────────
# Spearman correlation of every gene with DT (cor.test, exact = FALSE), with BH
# FDR across genes. Genes with fewer than min_obs finite pairs are skipped.
# Returns one row per gene, sorted by p-value.
run_spearman <- function(expr_matrix, dt_vector, label = "all", min_obs = 10,
                         unit_plural = "samples") {
  expr_matrix <- as.matrix(expr_matrix)
  storage.mode(expr_matrix) <- "numeric"
  finite_dt   <- is.finite(dt_vector)
  keep        <- apply(expr_matrix, 2, function(col)
    sum(finite_dt & is.finite(col)) >= min_obs)
  expr_matrix <- expr_matrix[, keep, drop = FALSE]
  message(sprintf("    %-40s | %d %s | %d genes",
                  label, nrow(expr_matrix), unit_plural, ncol(expr_matrix)))
  res <- sapply(colnames(expr_matrix), function(g) {
    ok <- finite_dt & is.finite(expr_matrix[, g])
    t  <- cor.test(expr_matrix[ok, g], dt_vector[ok],
                   method = "spearman", exact = FALSE)
    c(rho = unname(t$estimate), p_value = t$p.value)
  })
  df         <- as.data.frame(t(res))
  df$gene    <- rownames(df)
  df$p_adj   <- p.adjust(df$p_value, method = "BH")
  df$n_cells <- nrow(expr_matrix)
  df$group   <- label
  df[order(df$p_value), c("group","gene","n_cells","rho","p_value","p_adj")]
}

################################################################################
# PIPELINE — runs the whole Figure 1 analysis for one dataset
################################################################################

# SINGLE-GROUP BEHAVIOUR
# A dataset with only one cancer type has no between-group structure, so three
# things change automatically (each is announced in the log):
#   panel a  no Kruskal-Wallis p, no eta2, no between/within variance line
#   Fisher   not defined with one group -> fisher_df is NULL and every Fisher
#            output (volcano) is skipped
#   Deconf   subtracting a single group mean is a constant shift, which leaves
#            Spearman ranks untouched -> the pooled result is reused verbatim
# Everything else — pooled Spearman, volcano, gene scatters — runs unchanged.
run_fig1_pipeline <- function(cfg) {

  message(sprintf("\n\n%s\n### DATASET: %s\n%s",
                  strrep("#", 80), cfg$nice_name, strrep("#", 80)))

  UNITS  <- cfg$unit_plural
  GRP    <- cfg$group_singular
  GRPS   <- cfg$group_plural

  ##############################################################################
  # 1. LOAD & ALIGN DATA
  ##############################################################################

  message(">>> Loading data...")
  dat         <- cfg$loader()
  expr_mat    <- dat$expr_mat
  dt_vec      <- dat$dt_vec
  lineage_vec <- dat$lineage_vec
  # Sample / cell-line IDs. Optional: a loader that doesn't supply id_vec falls
  # back to the expression matrix row names, then to positional labels, so the
  # supplemental tables always have a first column to key on.
  id_vec      <- dat$id_vec
  if (is.null(id_vec)) {
    id_vec <- rownames(expr_mat)
    if (is.null(id_vec)) id_vec <- sprintf("%s_%d", cfg$tag, seq_along(dt_vec))
  }
  id_vec <- as.character(id_vec)
  if (length(id_vec) != length(dt_vec))
    stop(sprintf("id_vec has %d entries but dt_vec has %d — loader '%s' returned ",
                 length(id_vec), length(dt_vec), cfg$tag),
         "IDs that don't line up with the DT vector.")
  rm(dat)

  # ── Harmonise gene symbols before anything downstream ───────────────────────
  # Applied here so the corrected symbols flow into every correlation, volcano
  # and scatter. Datasets are processed in RUN_DATASETS order, so the reference
  # dataset (CCLE) must come first for a later one to resolve multi-mapping
  # symbols against it. Columns that end up with a duplicate symbol are dropped
  # (first occurrence kept).
  if (isTRUE(cfg$harmonise_symbols)) {
    ref <- if (is.null(cfg$symbol_ref_from)) NULL else .GENE_REF[[cfg$symbol_ref_from]]
    if (!is.null(cfg$symbol_ref_from) && is.null(ref))
      message(sprintf(paste("    NOTE: no '%s' gene set cached yet — multi-mapping",
                            "symbols will take their first candidate."),
                      cfg$symbol_ref_from))

    old_names <- colnames(expr_mat)
    new_names <- harmonise_symbols(old_names, ref,
                                   label = sprintf("%s symbols", cfg$tag))
    new_names <- as.character(new_names)         # drop the attribute
    dup <- duplicated(new_names)


    if (any(dup)) {
      message(sprintf("    Dropping %d duplicate column(s) after harmonisation.",
                      sum(dup)))
      expr_mat  <- expr_mat[, !dup, drop = FALSE]
      new_names <- new_names[!dup]
    }
    colnames(expr_mat) <- new_names
  }
  .GENE_REF[[cfg$tag]] <- colnames(expr_mat)

  # DT table for the figure panels (sample ID, DT, group)
  pa_raw <- data.frame(sample_id = id_vec,
                       dt = dt_vec, lineage = as.character(lineage_vec),
                       stringsAsFactors = FALSE)

  ##############################################################################
  # 2. CORRELATIONS
  ##############################################################################

  # ── 2a: cor_all — pooled ────────────────────────────────────────────────────
  message(sprintf("\n>>> [1/3] Correlation — all %s pooled...", UNITS))
  cor_all <- run_spearman(expr_mat, dt_vec, label = "All_Pooled",
                          unit_plural = UNITS)

  # ── 2b: per group + Fisher meta p-value ─────────────────────────────────────
  # Spearman within each group that has >= MIN_GRP samples (genes with zero
  # variance in that group are dropped first).
  message(sprintf("\n>>> [2/3] Correlation — per %s + Fisher meta p-value...", GRP))
  lineage_counts    <- table(lineage_vec)
  MIN_GRP           <- if (is.null(cfg$min_group_n)) MIN_LINEAGE else cfg$min_group_n
  eligible_lineages <- sort(names(lineage_counts[lineage_counts >= MIN_GRP]))
  message(sprintf("    %d %s with >= %d %s",
                  length(eligible_lineages), GRPS, MIN_GRP, UNITS))

  # Small single-group datasets can end up with nothing eligible. Rather than
  # failing, treat every sample as one group and say so.
  if (length(eligible_lineages) == 0) {
    warning(sprintf(paste("No %s reaches the %d-%s minimum — treating all %s as",
                          "a single group. Set min_group_n in the dataset config",
                          "to silence this."),
                    GRPS, MIN_GRP, sub("s$", "", UNITS), UNITS))
    eligible_lineages <- sort(unique(as.character(lineage_vec)))
  }

  # Fisher and the deconfounded run only mean something with >= 2 groups.
  multi_group <- length(eligible_lineages) >= 2

  cor_per_lineage <- list()

  for (lin in eligible_lineages) {
    idx      <- which(lineage_vec == lin)
    expr_sub <- expr_mat[idx, , drop = FALSE]
    dt_sub   <- dt_vec[idx]
    gene_var <- apply(expr_sub, 2, var, na.rm = TRUE)
    expr_sub <- expr_sub[, !is.na(gene_var) & gene_var > 0, drop = FALSE]
    cor_lin  <- run_spearman(expr_sub, dt_sub, label = lin, unit_plural = UNITS)
    cor_per_lineage[[lin]] <- cor_lin
  }

  # Direction-penalised Fisher combination — needs at least two groups.
  # For each gene: the reference sign is the pooled rho (or the median group rho
  # if the gene is missing from the pooled run); each group's p is adjusted by
  # fisher_dir_p(), floored to [1e-15, 1 - 1e-15], and combined with Fisher's
  # method (-2 * sum(log p) ~ chi-squared, 2k df). Genes tested in fewer than two
  # groups get NA and are dropped.
  n_elig_cells <- sum(lineage_vec %in% eligible_lineages)

  if (!multi_group) {
    message(sprintf(paste("    Single %s — Fisher meta-analysis skipped",
                          "(nothing to combine across %s)."), GRP, GRPS))
    fisher_df <- NULL
  } else {
  message(sprintf("    Computing Fisher meta p-values (%s)...", FISHER_METHOD))
  all_genes_tested <- unique(unlist(lapply(cor_per_lineage, `[[`, "gene")))
  ref_rho_vec <- setNames(cor_all$rho, cor_all$gene)
  get_ref_rho <- function(g) {
    r <- ref_rho_vec[g]
    if (!is.na(r)) return(r)
    rhos <- sapply(cor_per_lineage, function(df) {
      idx <- match(g, df$gene); if (is.na(idx)) NA_real_ else df$rho[idx]
    })
    median(rhos, na.rm = TRUE)
  }

  fisher_p_raw <- sapply(all_genes_tested, function(g) {
    ref <- get_ref_rho(g); if (is.na(ref)) return(NA_real_)
    ps  <- sapply(cor_per_lineage, function(df) {
      idx <- match(g, df$gene); if (is.na(idx)) return(NA_real_)
      fisher_dir_p(df$p_value[idx], sign(df$rho[idx]) == sign(ref))
    })
    ps <- ps[!is.na(ps)]; ps <- pmax(pmin(ps, 1 - 1e-15), 1e-15)
    if (length(ps) < 2) return(NA_real_)
    pchisq(-2 * sum(log(ps)), df = 2 * length(ps), lower.tail = FALSE)
  })
  fisher_p_adj <- p.adjust(fisher_p_raw, method = "BH")

  fisher_df <- data.frame(
    group   = "Fisher",
    gene    = names(fisher_p_raw),
    n_cells = n_elig_cells,
    rho     = sapply(names(fisher_p_raw), get_ref_rho),
    p_value = fisher_p_raw,
    p_adj   = fisher_p_adj,
    stringsAsFactors = FALSE)
  fisher_df <- fisher_df[!is.na(fisher_df$p_value), ]
  fisher_df <- fisher_df[order(fisher_df$p_value), ]
  }   # end multi_group Fisher block

  # ── 2c: deconfounded — group-mean-subtracted ────────────────────────────────
  # Within each eligible group, subtract the group mean from every gene and from
  # DT, then run one Spearman over all samples together.
  message(sprintf("\n>>> [3/3] Correlation — %s-mean-subtracted (deconfounded)...", GRP))

  if (!multi_group) {
    # Removing a single group's mean subtracts the same constant from every
    # sample, which leaves the Spearman ranks — and therefore rho, p and FDR —
    # exactly as they are in the pooled run. Reuse it rather than recompute.
    message(sprintf(paste("    Single %s — deconfounded is identical to the",
                          "pooled run; reusing it."), GRP))
    cor_dc          <- cor_all
    cor_dc$group    <- "Deconfounded (single group: same as pooled)"
  } else {
  elig_idx <- which(lineage_vec %in% eligible_lineages)
  expr_dc  <- as.matrix(expr_mat[elig_idx, ])
  storage.mode(expr_dc) <- "numeric"
  dt_dc    <- dt_vec[elig_idx]
  lin_dc   <- lineage_vec[elig_idx]

  for (lin in eligible_lineages) {
    idx_lin <- which(lin_dc == lin)
    expr_dc[idx_lin, ] <- sweep(
      expr_dc[idx_lin, , drop = FALSE], 2,
      colMeans(expr_dc[idx_lin, , drop = FALSE], na.rm = TRUE), "-")
    dt_dc[idx_lin] <- dt_dc[idx_lin] - mean(dt_dc[idx_lin], na.rm = TRUE)
  }
  cor_dc <- run_spearman(as.data.frame(expr_dc), dt_dc, label = "Deconfounded",
                         unit_plural = UNITS)

  rm(expr_dc); gc(verbose = FALSE)
  }   # end multi_group deconfounded block


  ##############################################################################
  # 3. EXPRESSION OF MKI67 + FOCUS GENE for the figure scatter panels
  ##############################################################################

  # Expression for the figure_1 scatter panels — just these genes, so the whole
  # matrix does not have to be kept in fig1_results.
  newfig_genes <- unique(c("MKI67", cfg$focus_gene))
  newfig_expr  <- expr_mat[, intersect(newfig_genes, colnames(expr_mat)),
                           drop = FALSE]
  if (ncol(newfig_expr) < length(newfig_genes))
    message(sprintf("    NOTE: %s missing from expression matrix.",
                    paste(setdiff(newfig_genes, colnames(newfig_expr)),
                          collapse = ", ")))

  message(sprintf("\n>>> [%s] outputs complete.", cfg$nice_name))


  invisible(list(
    # Data the figure panels are rebuilt from
    panel_data = list(dt_data      = pa_raw,
                      scatter_expr = newfig_expr,
                      scatter_dt   = dt_vec),

    cor_all     = cor_all,
    fisher_df   = fisher_df,
    cor_dc      = cor_dc))
}

################################################################################
# ANALYSIS — run the pipeline once per dataset
#
# This is the expensive part. Everything below the FIGURES banner is layout only
# and can be re-run against `fig1_results` without touching the analysis.
################################################################################

fig1_results <- list()

for (ds in RUN_DATASETS) {
  cfg <- DATASETS[[ds]]
  if (is.null(cfg)) {
    warning(sprintf("Unknown dataset '%s' — skipped.", ds)); next
  }
  res <- tryCatch(run_fig1_pipeline(cfg),
                  error = function(e) {
                    message(sprintf("\n!!! DATASET '%s' FAILED: %s", ds, e$message))
                    NULL
                  })
  fig1_results[[ds]] <- res
  gc(verbose = FALSE)
}

message("\n\n>>> ALL DATASETS DONE.")
for (ds in names(fig1_results)) {
  r <- fig1_results[[ds]]
  message(sprintf("    %-6s : %s", ds,
                  if (is.null(r)) "FAILED (see log above)" else "OK"))
}


################################################################################
################################################################################
##                                                                            ##
##  FIGURE 1 — ONE SELF-CONTAINED BLOCK PER PANEL                             ##
##                                                                            ##
##  Every panel below does its own data prep and builds its own ggplot with    ##
##  all values written out in full. Editing one panel cannot affect any        ##
##  other. To re-run a panel, select its block and run it — the analysis       ##
##  results it reads from `fig1_results` are already in memory.                ##
##                                                                            ##
##  The only functions used are the ones that cannot reasonably be inlined:    ##
##    format_p_value()           p-value -> plotmath string                   ##
##                                                                            ##
##  Figure 1 has nine panels (a-i): see the comment above the assembly for    ##
##  the layout.                                                               ##
##                                                                            ##
##  The standalone panel PDFs are named after the object that draws them,      ##
##  e.g. panel_volcano_pooled_hcc.pdf, not by layout letter.                   ##
##                                                                            ##
##  One assembled figure comes out of this section:                            ##
##    <FIG_DIR>/figure_1.pdf (+ .png)       the full figure (panels a-i)       ##
##                                                                            ##
################################################################################
################################################################################

FIG_DIR <- OUT_ROOT

# figure_1 includes HCC panels and a three-way shared programme, so all
# three datasets have to be in fig1_results. Failing here with a clear message
# beats a confusing NULL error twenty panels later.
for (need in c("ccle", "pdx", "hcc"))
  if (is.null(fig1_results[[need]]))
    stop(sprintf(paste("figure_1 needs the '%s' results. Add it to RUN_DATASETS",
                       "and re-run the analysis section."), need))
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

################################################################################
# PANEL A — example PDX growth curve with its exponential fit
#
# One digitised growth curve (gc3492) drawn over the original chart image,
# with the fitted y = a * exp(b * x) curve and its statistics.
#
# Everything here is prefixed pa_ so it cannot overwrite any other object
# (a, b, p, s, d, labs, ... would otherwise clash with the rest of the script).
################################################################################

# ── Inputs ────────────────────────────────────────────────────────────────────
PA_DIR       <- file.path(WORK_DIR, "input_files", "pdx", "growth_curve_points")
PA_MARKERS   <- file.path(PA_DIR, "all_markers.csv")
PA_IMAGE_PNG <- file.path(PA_DIR, "gc.png")
PA_IMAGE_ID  <- "111316-319-R_P0_gc3492.jpg"   # which digitised image to analyse

if (!requireNamespace("png", quietly = TRUE))
  stop("Panel a needs the 'png' package: install.packages(\"png\")")

# ── Data ──────────────────────────────────────────────────────────────────────
# All digitised marker points; keep only the image being analysed
pa_all_markers <- read.csv(PA_MARKERS)
pa_d <- subset(pa_all_markers, image == PA_IMAGE_ID)
if (nrow(pa_d) < 3)
  stop(sprintf("Panel a: only %d points for image %s in %s.",
               nrow(pa_d), PA_IMAGE_ID, PA_MARKERS))

# ── Fit y = a * exp(b * x) ────────────────────────────────────────────────────
# Gaussian GLM with log link: E[y] = exp(b0 + b1*x), errors on the
# original y scale. Same estimates as nls(y ~ a * exp(b * x)).
pa_fit <- glm(y ~ x, family = gaussian(link = "log"), data = pa_d)
pa_s   <- summary(pa_fit)

pa_a    <- unname(exp(coef(pa_fit)[1]))             # value of y at x = 0
pa_b    <- unname(coef(pa_fit)[2])                  # continuous growth rate per day
pa_p    <- pa_s$coefficients["x", "Pr(>|t|)"]       # t-test p-value for b
pa_res  <- pa_d$y - fitted(pa_fit)                  # residuals on the y scale
pa_r2   <- 1 - sum(pa_res^2) / sum((pa_d$y - mean(pa_d$y))^2)   # pseudo-R^2
pa_dt   <- log(2) / pa_b                            # doubling time, days

message(sprintf("Panel a: %s | n = %d | a = %.4g, b = %.3g/day | R2 = %.2f | DT = %.2f days",
                PA_IMAGE_ID, nrow(pa_d), pa_a, pa_b, pa_r2, pa_dt))

# Fitted curve evaluated on a fine x grid, for a smooth line.
# (Not named grid: grid() is also a base R function.)
pa_fit_curve   <- data.frame(x = seq(min(pa_d$x), max(pa_d$x), length.out = 200))
pa_fit_curve$y <- predict(pa_fit, newdata = pa_fit_curve, type = "response")

# ── Background image ──────────────────────────────────────────────────────────
pa_img   <- png::readPNG(PA_IMAGE_PNG)
pa_img_h <- dim(pa_img)[1]   # height in pixels
pa_img_w <- dim(pa_img)[2]   # width in pixels

# ── Plot ──────────────────────────────────────────────────────────────────────
# All plotting parameters are set inside the ggplot call below. The image
# edges (xmin/xmax/ymin/ymax, repeated in coord_cartesian) are the data values
# at the four edges of gc.png, calibrated from its tick marks (day 2 = col 178,
# day 92 = col 1310; 0 = row 980, 2600 = row 28). Specific to this image:
# recalibrate for other images.
panel_a_gc <- ggplot(pa_d, aes(x, y)) +
  # 1. Original chart as the background, stretched to the calibrated edges
  annotation_raster(pa_img, xmin = -12.152, xmax = 100.428,
                            ymin = -278.571, ymax = 2676.471) +
  # 2. Semi-transparent box behind the statistics
  annotate("rect", xmin = 3, xmax = 60, ymin = 1500, ymax = 2560,
           fill = "white", alpha = 0.4,
           colour = "grey40", linewidth = 0.4) +
  # 3. Data points and fitted exponential curve
  geom_point(size = 1.5, colour = "blue") +
  geom_line(data = pa_fit_curve, colour = "blue", linewidth = 1) +
  # 4. Statistics text, one plotmath line each; parse = TRUE renders them.
  #    Numbers are quoted so trailing zeros survive parsing; %% is a literal %
  #    inside sprintf, so %%.%% becomes the plotmath centred dot %.%.
  geom_text(data = data.frame(
              label = c(
                sprintf('italic(y) == "%.4g" %%.%% italic(e)^{"%.3g" * italic(x)}',
                        pa_a, pa_b),
                sprintf('italic(R)^2 == "%.2f"', pa_r2),
                format_p_value(pa_p),
                sprintf('italic(DT) == "%.2f"~days', pa_dt)),
              x = 6,                          # left edge of the text block (days)
              y = 2480 - (0:3) * 250),        # one line every 250 mm^3, top down
            aes(x = x, y = y, label = label),
            hjust = 0, vjust = 1, size = 6, parse = TRUE,
            inherit.aes = FALSE) +
  # Panel exactly covers the image, no padding
  coord_cartesian(xlim = c(-12.152, 100.428),
                  ylim = c(-278.571, 2676.471), expand = FALSE) +
  # No ggplot axes: the image supplies its own
  theme_void() +
  theme(aspect.ratio = pa_img_h / pa_img_w, plot.margin = margin(0, 0, 0, 0))

# Standalone copy at the original image's pixel size (width * dpi = img_w)
ggsave(file.path(FIG_DIR, "panel_a_growth_curve.png"), panel_a_gc,
       width = pa_img_w / 300, height = pa_img_h / 300, dpi = 300)

################################################################################
# PANEL B — PDX mitotic density vs doubling time
#
# Doubling times are the PDX values already used by the other PDX panels
# (fig1_results$pdx$panel_data$dt_data: sample_id, dt in days, lineage), so
# only the mitosis table is read here. Samples are matched by ID: row names of
# points.csv against dt_data$sample_id. dt_data holds the samples shared by the
# PDX RNA and metadata tables, so this panel uses the same sample set as the
# rest of the PDX panels.
#
# Everything here is prefixed pb_ so it cannot overwrite any other object.
################################################################################

PB_MITOSIS_CSV <- file.path(WORK_DIR, "input_files", "pdx", "growth_curve_points",
                            "points.csv")
PB_MITOSIS_COL <- "mitotic_density_per_mm2"

pb_mitosis <- read.csv(PB_MITOSIS_CSV, row.names = 1, check.names = FALSE)
if (!PB_MITOSIS_COL %in% colnames(pb_mitosis))
  stop(sprintf("Panel b: column '%s' not found in %s (columns: %s).",
               PB_MITOSIS_COL, PB_MITOSIS_CSV,
               paste(colnames(pb_mitosis), collapse = ", ")))

pb_dt   <- fig1_results$pdx$panel_data$dt_data
pb_ids  <- intersect(pb_dt$sample_id, rownames(pb_mitosis))
pb_df   <- data.frame(sample_id = pb_ids,
                      dt_days   = pb_dt$dt[match(pb_ids, pb_dt$sample_id)],
                      mitotic_density_per_mm2 =
                        as.numeric(pb_mitosis[pb_ids, PB_MITOSIS_COL]))

# Both axes are log10, so keep strictly positive, finite values only
pb_ok <- is.finite(pb_df$dt_days) & pb_df$dt_days > 0 &
         is.finite(pb_df$mitotic_density_per_mm2) & pb_df$mitotic_density_per_mm2 > 0
message(sprintf(paste("Panel b: %d mitosis samples, %d PDX DT samples, %d shared,",
                      "%d plotted (%d dropped: missing or non-positive values)"),
                nrow(pb_mitosis), nrow(pb_dt), length(pb_ids),
                sum(pb_ok), sum(!pb_ok)))
pb_df <- pb_df[pb_ok, ]
if (nrow(pb_df) < 3)
  stop("Panel b: fewer than 3 samples shared between points.csv and the PDX DT data.")

panel_b_mitosis <- ggplot(data = pb_df,
                          aes(x = log10(dt_days), y = log10(mitotic_density_per_mm2))) +
  geom_point(size = 1, alpha = 0.5, colour = "#808080") +
  geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
  stat_cor(aes(label = paste(after_stat(r.label),
                             after_stat(format_p_value(p)),
                             sep = "~`,`~")),
           r.accuracy = 0.001, method = "spearman", parse = TRUE,
           size = 7,
           label.y = max(log10(pb_df$mitotic_density_per_mm2)) +
             0.10 * diff(range(log10(pb_df$mitotic_density_per_mm2)))) +
  xlab(TeX("$log_{10}(DT)\\ (days)$")) +
  ylab(TeX("$log_{10}(MD)\\ (mm^{-2})$")) +
  theme(legend.position  = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(colour = "black"),
        axis.text        = element_text(size = 16,  face = "bold"),
        axis.title       = element_text(size = 16, face = "bold"))

ggsave(file.path(FIG_DIR, "panel_mitosis_vs_dt_pdx.pdf"), panel_b_mitosis,
       width = 5, height = 5)

################################################################################
# PANEL C — PDX doubling-time distribution by tumour type
#
# Tumour types with >= 20 samples, ordered by median DT; Kruskal-Wallis p.
################################################################################

dtl_pdx_raw <- fig1_results$pdx$panel_data$dt_data
dtl_pdx_df  <- dtl_pdx_raw[is.finite(dtl_pdx_raw$dt) & dtl_pdx_raw$dt > 0 &
                 !is.na(dtl_pdx_raw$lineage) & dtl_pdx_raw$lineage != "", ]
dtl_pdx_df$lineage <- as.character(dtl_pdx_df$lineage)

dtl_pdx_n   <- table(dtl_pdx_df$lineage)
dtl_pdx_df  <- dtl_pdx_df[dtl_pdx_df$lineage %in% names(dtl_pdx_n[dtl_pdx_n >= 20]), ]

dtl_pdx_med        <- tapply(dtl_pdx_df$dt, dtl_pdx_df$lineage, median, na.rm = TRUE)
dtl_pdx_df$lineage <- factor(dtl_pdx_df$lineage, levels = names(sort(dtl_pdx_med)))
dtl_pdx_df$log_dt  <- log10(dtl_pdx_df$dt)

dt_by_lineage_pdx <- ggplot(dtl_pdx_df, aes(x = lineage, y = log_dt, fill = lineage)) +
  geom_jitter(aes(colour = lineage), width = 0.22, alpha = 1.0,
              size = 3, shape = 16, na.rm = TRUE) +
  geom_violin(alpha = 0.55, linewidth = 0.25, scale = "width",
              bw = "nrd0", na.rm = TRUE) +
  geom_boxplot(width = 0.18, fill = "white", outlier.shape = NA,
               linewidth = 0.70, colour = "black", na.rm = TRUE) +
  stat_compare_means(aes(label = after_stat(format_p_value(..p..))),
                     parse = TRUE, na.rm = TRUE, size = 5,
                     vjust = 0.4, hjust = 0.5,
                     label.x = (nlevels(dtl_pdx_df$lineage) + 1) / 4,
                     label.y = 2.1) +
  # x labels: tumour type with its n, e.g. "Breast\n(n=34)"
  scale_x_discrete(labels = function(l)
    paste0(l, "\n(n=", table(dtl_pdx_df$lineage)[l], ")")) +
  scale_fill_manual(values = scales::hue_pal(l = 60, c = 70)(nlevels(dtl_pdx_df$lineage)),
                    guide = "none") +
  scale_colour_manual(values = scales::hue_pal(l = 60, c = 70)(nlevels(dtl_pdx_df$lineage)),
                      guide = "none") +
  xlab(NULL) +
  ylab(TeX("$log_{10}(DT)\\ (days)$")) +
  # Visible range only: coord_cartesian keeps every sample in the violins,
  # boxplots and the Kruskal-Wallis test (scale limits would drop them).
  scale_y_continuous(breaks = seq(0, 2, 0.5)) +
  coord_cartesian(ylim = c(0, 2.3)) +
  theme(legend.position  = "none",
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(colour = "black"),
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 14,
                                        face = "bold"),
        axis.text.y      = element_text(size = 14, face = "bold"),
        axis.title       = element_text(size = 14, face = "bold"))
ggsave(file.path(FIG_DIR, "panel_dt_by_lineage_pdx.pdf"), dt_by_lineage_pdx, width = 12, height = 4)

################################################################################
# PANEL F — CCLE pooled volcano
#
# Pooled Spearman rho vs -log10(p). Red / blue: FDR < 0.1 and |rho| >= 0.1;
# green: proliferation markers; purple: CCNB1 and KIF18A. The same rules are
# used by every volcano in this script.
################################################################################

vol_ccle_df   <- fig1_results$ccle$cor_all
vol_ccle_df$x <- vol_ccle_df$rho
vol_ccle_df$y <- -log10(vol_ccle_df$p_value)
vol_ccle_df   <- vol_ccle_df[is.finite(vol_ccle_df$x) & is.finite(vol_ccle_df$y), ]

vol_ccle_df$significance <- "not_significant"
vol_ccle_df$significance[vol_ccle_df$p_adj < 0.1 & vol_ccle_df$x >=  0.1] <- "high_corr"
vol_ccle_df$significance[vol_ccle_df$p_adj < 0.1 & vol_ccle_df$x <= -0.1] <- "low_corr"
vol_ccle_df$significance[vol_ccle_df$gene %in% c("MKI67","PCNA","MCM2","CDKN1A","CDKN1B")] <-
  "key_gene"
vol_ccle_df$significance[vol_ccle_df$gene %in% c("CCNB1","KIF18A")] <-
  "key_gene_2"
vol_ccle_df$significance <- factor(vol_ccle_df$significance,
  levels = c("not_significant", "low_corr", "high_corr", "key_gene", "key_gene_2"))

# Labels: the 5 key genes (green), CCNB1 + KIF18A (purple), plus the 3 strongest
# positive and 3 strongest negative
vol_ccle_key  <- intersect(c("MKI67","PCNA","MCM2","CDKN1A","CDKN1B","CCNB1","KIF18A"), vol_ccle_df$gene)
vol_ccle_cand <- vol_ccle_df[!vol_ccle_df$gene %in% vol_ccle_key, ]
vol_ccle_lab  <- vol_ccle_df[vol_ccle_df$gene %in% c(vol_ccle_key,
                                head(vol_ccle_cand$gene[order(-vol_ccle_cand$x)], 3),
                                head(vol_ccle_cand$gene[order( vol_ccle_cand$x)], 3)), ]

# grey first, key genes last, so highlights are not buried
vol_ccle_df <- vol_ccle_df[order(match(as.character(vol_ccle_df$significance),
                         c("not_significant","low_corr","high_corr","key_gene","key_gene_2"))), ]

volcano_pooled_ccle <- ggplot(vol_ccle_df, aes(x = x, y = y, colour = significance)) +
  geom_point(aes(alpha = significance), size = 0.7) +
  geom_text_repel(data = vol_ccle_lab, aes(label = gene, colour = significance),
                  fontface = ifelse(vol_ccle_lab$gene %in% vol_ccle_key, "bold", "plain"),
                  size = 4, box.padding = 0.5, max.overlaps = Inf,
                  point.size = 2, segment.size = 0.2, segment.alpha = 0.5,
                  show.legend = FALSE) +
  scale_colour_manual(values = c("not_significant" = "#808080",
                                 "low_corr"        = "#005AB5",
                                 "high_corr"       = "#DC3220",
                                 "key_gene"        = "#2CA02C",
                                 "key_gene_2"      = "#8E44AD")) +
  scale_alpha_manual(values = c("not_significant" = 0.1, "low_corr" = 1,
                                "high_corr" = 1, "key_gene" = 1,
                                "key_gene_2" = 1), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  xlab(TeX("$\\rho$")) +
  ylab(TeX("$-log_{10}(P_{value})$")) +
  ggtitle("CCLE") +
  theme(legend.position  = "none",
        plot.title       = element_text(size = 14, face = "bold",
                                        hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(colour = "black"),
        axis.text        = element_text(size = 14, face = "bold"),
        axis.title       = element_text(size = 16, face = "bold"))
ggsave(file.path(FIG_DIR, "panel_volcano_pooled_ccle.pdf"), volcano_pooled_ccle, width =  4, height = 4)

################################################################################
# PANEL D — PDX pooled volcano (same rules as panel f)
################################################################################

vol_pdx_df   <- fig1_results$pdx$cor_all
vol_pdx_df$x <- vol_pdx_df$rho
vol_pdx_df$y <- -log10(vol_pdx_df$p_value)
vol_pdx_df   <- vol_pdx_df[is.finite(vol_pdx_df$x) & is.finite(vol_pdx_df$y), ]

vol_pdx_df$significance <- "not_significant"
vol_pdx_df$significance[vol_pdx_df$p_adj < 0.1 & vol_pdx_df$x >=  0.1] <- "high_corr"
vol_pdx_df$significance[vol_pdx_df$p_adj < 0.1 & vol_pdx_df$x <= -0.1] <- "low_corr"
vol_pdx_df$significance[vol_pdx_df$gene %in% c("MKI67","PCNA","MCM2","CDKN1A","CDKN1B")] <-
  "key_gene"
vol_pdx_df$significance[vol_pdx_df$gene %in% c("CCNB1","KIF18A")] <-
  "key_gene_2"
vol_pdx_df$significance <- factor(vol_pdx_df$significance,
  levels = c("not_significant", "low_corr", "high_corr", "key_gene", "key_gene_2"))

vol_pdx_key  <- intersect(c("MKI67","PCNA","MCM2","CDKN1A","CDKN1B","CCNB1","KIF18A"), vol_pdx_df$gene)
vol_pdx_cand <- vol_pdx_df[!vol_pdx_df$gene %in% vol_pdx_key, ]
vol_pdx_lab  <- vol_pdx_df[vol_pdx_df$gene %in% c(vol_pdx_key,
                                head(vol_pdx_cand$gene[order(-vol_pdx_cand$x)], 3),
                                head(vol_pdx_cand$gene[order( vol_pdx_cand$x)], 3)), ]

vol_pdx_df <- vol_pdx_df[order(match(as.character(vol_pdx_df$significance),
                         c("not_significant","low_corr","high_corr","key_gene","key_gene_2"))), ]

volcano_pooled_pdx <- ggplot(vol_pdx_df, aes(x = x, y = y, colour = significance)) +
  geom_point(aes(alpha = significance), size = 0.7) +
  geom_text_repel(data = vol_pdx_lab, aes(label = gene, colour = significance),
                  fontface = ifelse(vol_pdx_lab$gene %in% vol_pdx_key, "bold", "plain"),
                  size = 4, box.padding = 0.5, max.overlaps = Inf,
                  point.size = 2, segment.size = 0.2, segment.alpha = 0.5,
                  show.legend = FALSE) +
  scale_colour_manual(values = c("not_significant" = "#808080",
                                 "low_corr"        = "#005AB5",
                                 "high_corr"       = "#DC3220",
                                 "key_gene"        = "#2CA02C",
                                 "key_gene_2"      = "#8E44AD")) +
  scale_alpha_manual(values = c("not_significant" = 0.1, "low_corr" = 1,
                                "high_corr" = 1, "key_gene" = 1,
                                "key_gene_2" = 1), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  xlab(TeX("$\\rho$")) +
  ylab(TeX("$-log_{10}(P_{value})$")) +
  ggtitle("PDX") +
  theme(legend.position  = "none",
        plot.title       = element_text(size = 14, face = "bold",
                                        hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(colour = "black"),
        axis.text        = element_text(size = 14, face = "bold"),
        axis.title       = element_text(size = 16, face = "bold"))
ggsave(file.path(FIG_DIR, "panel_volcano_pooled_pdx.pdf"), volcano_pooled_pdx,
       width = 4, height = 4)

################################################################################
# PANEL G — HCC pooled volcano (same rules as panel f)
################################################################################

vol_hcc_df   <- fig1_results$hcc$cor_all
vol_hcc_df$x <- vol_hcc_df$rho
vol_hcc_df$y <- -log10(vol_hcc_df$p_value)
vol_hcc_df   <- vol_hcc_df[is.finite(vol_hcc_df$x) & is.finite(vol_hcc_df$y), ]

vol_hcc_df$significance <- "not_significant"
vol_hcc_df$significance[vol_hcc_df$p_adj < 0.1 & vol_hcc_df$x >=  0.1] <- "high_corr"
vol_hcc_df$significance[vol_hcc_df$p_adj < 0.1 & vol_hcc_df$x <= -0.1] <- "low_corr"
vol_hcc_df$significance[vol_hcc_df$gene %in% c("MKI67","PCNA","MCM2","CDKN1A","CDKN1B")] <-
  "key_gene"
vol_hcc_df$significance[vol_hcc_df$gene %in% c("CCNB1","KIF18A")] <-
  "key_gene_2"
vol_hcc_df$significance <- factor(vol_hcc_df$significance,
  levels = c("not_significant", "low_corr", "high_corr", "key_gene", "key_gene_2"))

vol_hcc_key  <- intersect(c("MKI67","PCNA","MCM2","CDKN1A","CDKN1B","CCNB1","KIF18A"), vol_hcc_df$gene)
vol_hcc_cand <- vol_hcc_df[!vol_hcc_df$gene %in% vol_hcc_key, ]
vol_hcc_lab  <- vol_hcc_df[vol_hcc_df$gene %in% c(vol_hcc_key,
                                head(vol_hcc_cand$gene[order(-vol_hcc_cand$x)], 3),
                                head(vol_hcc_cand$gene[order( vol_hcc_cand$x)], 3)), ]

vol_hcc_df <- vol_hcc_df[order(match(as.character(vol_hcc_df$significance),
                         c("not_significant","low_corr","high_corr","key_gene","key_gene_2"))), ]

volcano_pooled_hcc <- ggplot(vol_hcc_df, aes(x = x, y = y, colour = significance)) +
  geom_point(aes(alpha = significance), size = 0.7) +
  geom_text_repel(data = vol_hcc_lab, aes(label = gene, colour = significance),
                  fontface = ifelse(vol_hcc_lab$gene %in% vol_hcc_key, "bold", "plain"),
                  size = 4, box.padding = 0.5, max.overlaps = Inf,
                  point.size = 2, segment.size = 0.2, segment.alpha = 0.5,
                  show.legend = FALSE) +
  scale_colour_manual(values = c("not_significant" = "#808080",
                                 "low_corr"        = "#005AB5",
                                 "high_corr"       = "#DC3220",
                                 "key_gene"        = "#2CA02C",
                                 "key_gene_2"      = "#8E44AD")) +
  scale_alpha_manual(values = c("not_significant" = 0.1, "low_corr" = 1,
                                "high_corr" = 1, "key_gene" = 1,
                                "key_gene_2" = 1), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  xlab(TeX("$\\rho$")) +
  ylab(TeX("$-log_{10}(P_{value})$")) +
  ggtitle("HCC") +
  theme(legend.position  = "none",
        plot.title       = element_text(size = 14, face = "bold",
                                        hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(colour = "black"),
        axis.text        = element_text(size = 14, face = "bold"),
        axis.title       = element_text(size = 16, face = "bold"))
ggsave(file.path(FIG_DIR, "panel_volcano_pooled_hcc.pdf"), volcano_pooled_hcc,
       width = 4, height = 4)

################################################################################
# PANEL E — PDX, RAN and MKI67 vs doubling time (one panel, two facets)
################################################################################

mki67_pdx_df <- data.frame(
  x = as.numeric(fig1_results$pdx$panel_data$scatter_expr[, "MKI67"]),
  y = log10(fig1_results$pdx$panel_data$scatter_dt))
mki67_pdx_df <- mki67_pdx_df[is.finite(mki67_pdx_df$x) & is.finite(mki67_pdx_df$y), ]
mki67_pdx_df$gene <- "MKI67"

ran_pdx_df <- data.frame(
  x = as.numeric(fig1_results$pdx$panel_data$scatter_expr[, "RAN"]),
  y = log10(fig1_results$pdx$panel_data$scatter_dt))
ran_pdx_df <- ran_pdx_df[is.finite(ran_pdx_df$x) & is.finite(ran_pdx_df$y), ]
ran_pdx_df$gene <- "RAN"

ran_mki67_pdx_df <- rbind(ran_pdx_df,mki67_pdx_df)

scatter_ran_mki67_pdx <- ggplot(ran_mki67_pdx_df, aes(x = x, y = y)) +
  geom_point(alpha = 0.4, colour = "#808080", size = 0.9) +
  geom_smooth(method = "lm", colour = "#005AB5", se = FALSE, na.rm = TRUE) +
  stat_cor(aes(label = ..r.label..),
           method = "spearman", r.accuracy = 0.001, parse = TRUE, size = 4.5,
           label.x.npc = 0.03, label.y.npc = 0.97,
           geom = ggrepel::GeomTextRepel, bg.color = "white", bg.r = 0.25,
           force = 0, seed = 1, segment.color = NA) +
  stat_cor(aes(label = after_stat(format_p_value(..p..))),
           method = "spearman", parse = TRUE, size = 4.5,
           label.x.npc = 0.03, label.y.npc = 0.80,
           geom = ggrepel::GeomTextRepel, bg.color = "white", bg.r = 0.25,
           force = 0, seed = 1, segment.color = NA) +
  facet_wrap(~ gene, nrow = 2, scales = "free") +
  xlab(TeX("$log_{2}(RNA_{TPM}+1)$")) +
  ylab(TeX("$log_{10}(DT)\\ (days)$")) +
  theme(legend.position  = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(colour = "black"),
        axis.text        = element_text(size = 12, face = "bold"),
        axis.title       = element_text(size = 12, face = "bold"),
        strip.background = element_blank(),
        strip.text       = element_text(size = 12,  face = "bold"))
ggsave(file.path(FIG_DIR, "scatter_ran_mki67_pdx.pdf"), scatter_ran_mki67_pdx,
       width = 2, height = 4)


################################################################################
################################################################################
##                                                                            ##
##  GSEA ON GOBP — figure_1 panels h and i                                     ##
##                                                                            ##
##  Input: each dataset's pan-cancer pooled Spearman table,                    ##
##  fig1_results[[ds]]$cor_all, already in memory with harmonised symbols.    ##
##                                                                            ##
##    panel h  gsea_emap_combined — enrichment map of the PDX GSEA terms        ##
##             (fast-growth side), clustered by gene overlap + GO semantics     ##
##    panel i  gsea_table — GSEA results for CCLE | PDX | HCC                   ##
##                                                                            ##
##  (The log messages from this section call these "Panel i" / "Panel j".)   ##
##                                                                            ##
##  Everything here is prefixed gsea_ so it cannot overwrite any other object. ##
##                                                                            ##
################################################################################
################################################################################

library(GOSemSim); library(GO.db); library(org.Hs.eg.db)


################################################################################
# GSEA — one run per dataset
#
# Ranked vector: rho minus its median, sorted DECREASING, so the top of the list
# is the most POSITIVE rho and the bottom the most NEGATIVE. NES < 0 therefore
# means the set is high in FAST growers.
#
# minGSSize/maxGSSize set explicitly (GSEA's defaults are 10 and 500).
# eps = 0 stops p-values being floored at 1e-10. pvalueCutoff = 1 returns
# everything. Gene sets are the cached MSigDB GOBP sets (get_msig_gobp).
#
# Cached per dataset, so re-running the panels below does not re-run GSEA.
################################################################################

GSEA_MIN_GS <- 5
GSEA_MAX_GS <- 2000

.GSEA_CACHE <- if (exists(".GSEA_CACHE")) .GSEA_CACHE else new.env(parent = emptyenv())

gsea_gobp_run <- function(ds, force = FALSE) {
  if (!force && !is.null(.GSEA_CACHE[[ds]])) return(.GSEA_CACHE[[ds]])
  message(sprintf("\n>>> GSEA GOBP — %s", ds))

  cor_df <- fig1_results[[ds]]$cor_all
  gl <- cor_df$rho
  names(gl) <- cor_df$gene
  gl <- gl[is.finite(gl) & !is.na(names(gl)) & names(gl) != ""]

  # Duplicate symbols: keep the strongest absolute rho for each gene
  gl <- gl[order(abs(gl), decreasing = TRUE)]
  gl <- gl[!duplicated(names(gl))]

  gl <- sort(gl, decreasing = TRUE)          # REQUIRED by GSEA()
  gl <- gl - median(gl)
  message(sprintf("    %d genes | rho range [%.3f, %.3f] | %d tied values",
                  length(gl), min(gl), max(gl), length(gl) - length(unique(gl))))

  set.seed(42)
  res <- GSEA(geneList      = gl,
              TERM2GENE     = as.data.frame(get_msig_gobp()),
              minGSSize     = GSEA_MIN_GS,
              maxGSSize     = GSEA_MAX_GS,
              pAdjustMethod = "BH",
              pvalueCutoff  = 1,
              eps           = 0,
              seed          = TRUE,
              nPermSimple   = 1000,
              verbose       = FALSE)

  df <- as.data.frame(res)
  message(sprintf("    %d sets tested | FDR<0.05: %d (%d up, %d down)",
                  nrow(df), sum(df$p.adjust < 0.05),
                  sum(df$p.adjust < 0.05 & df$NES > 0),
                  sum(df$p.adjust < 0.05 & df$NES < 0)))

  # core_enrichment = the leading-edge genes, the GSEA equivalent of geneID
  df$term      <- gsub("^GOBP_", "", df$ID)
  df$direction <- ifelse(df$NES < 0, "fast-growth (neg rho)", "slow-growth (pos rho)")
  df$n_leading <- lengths(strsplit(df$core_enrichment, "/", fixed = TRUE))
  df <- df[order(df$p.adjust, -abs(df$NES)),
           c("ID", "term", "setSize", "NES", "pvalue", "p.adjust",
             "n_leading", "direction", "core_enrichment")]
  rownames(df) <- NULL

  out <- list(df = df)
  .GSEA_CACHE[[ds]] <- out
  out
}

gsea_all <- lapply(setNames(nm = c("ccle", "pdx", "hcc")), gsea_gobp_run)


################################################################################
# PANEL H — PDX GSEA enrichment map (combined similarity, two-level layout)
#
# Terms are clustered by BOTH gene overlap and ontology semantics:
#   S = W_SEM * S_semantic + (1 - W_SEM) * S_jaccard
# Both matrices are rank-normalised to [0,1] first (nonzero entries only),
# because raw Wang values (typically 0.3-0.9) would otherwise swamp raw Jaccard
# values (0-0.3).
#
# All parameters for this panel, in one place:
################################################################################

gsea_PADJ_CUT  <- 0.1        # which GSEA terms to cluster
gsea_NES_SIDE  <- "neg"      # "neg" = fast-growth only, "pos", or "both"
gsea_MAX_TERMS <- 500        # semantic matrix is O(n^2); keep this sane
gsea_W_SEM     <- 0.5        # 1 = ontology only, 0 = gene overlap only
gsea_K_CLUST   <- 11         # number of Ward clusters (cutree k)
gsea_USE_LEADING_EDGE <- TRUE  # TRUE = Jaccard on leading-edge genes;
                               # FALSE = on full gene sets
gsea_MIN_CLUST <- 3          # drop clusters with fewer terms than this
gsea_EDGE_Q    <- 0.90       # keep the top (1 - EDGE_Q) fraction of edges
gsea_GAP        <- 0.35      # clear space between cluster discs
gsea_DISC_SCALE <- 1.00      # cluster disc size; lower = more white space
gsea_META_KEEP  <- 0.50      # fraction of between-cluster links drawn

# Hand-written cluster names, keyed by cluster number. Ward cluster numbers
# follow the order terms enter cutree (most significant term first), so they
# are stable for the same input but WILL move if the data, filters, k or W_SEM
# change. The check table printed below lists each cluster's top terms next to
# its manual name — read it after every change. Set gsea_USE_MANUAL_LABELS to
# FALSE to name each cluster by its most significant term instead.
gsea_USE_MANUAL_LABELS <- TRUE
gsea_manual_labels <- c(
  "1"  = "Ribonucleoprotein and chromatin complex assembly",
  "2"  = "M phase: chromosome segregation and spindle",
  "3"  = "Translation and mRNA processing",
  "4"  = "DNA replication and repair",
  "5"  = "Cell cycle regulation and checkpoints",
  "6"  = "rRNA processing and Pol I/III transcription",
  "7"  = "Nucleocytoplasmic and mitochondrial transport",
  "8"  = "Signalling and sensory perception (heterogeneous)",
  "9"  = "Oxidative phosphorylation and nucleotide biosynthesis",
  "10" = "Cell fusion and differentiation",
  "11" = "Telomere maintenance and DNA biosynthesis")

# Keeps the top (1 - q) fraction of edges. If more than q of all pairs are
# exactly zero, the plain quantile is 0 and would keep every nonzero edge, so
# the quantile is then taken over the nonzero pairs only.
gsea_edge_cut <- function(M, q) {
  v   <- M[upper.tri(M)]
  cut <- quantile(v, q, names = FALSE)
  if (cut == 0 && any(v > 0)) {
    cut <- quantile(v[v > 0], q, names = FALSE)
    message(sprintf(paste("    edge threshold: >%.0f%% of pairs are zero,",
                          "quantile taken over nonzero pairs (%.3f)"), 100 * q, cut))
  }
  A <- M; diag(A) <- 0; A[A < cut] <- 0
  A
}

gsea_rank01 <- function(M) {
  ut <- upper.tri(M)
  v  <- M[ut]
  r  <- numeric(length(v))
  nz <- v > 0
  if (any(nz)) r[nz] <- rank(v[nz], ties.method = "average") / sum(nz)
  R  <- matrix(0, nrow(M), ncol(M), dimnames = dimnames(M))
  R[ut] <- r; R <- R + t(R); diag(R) <- 1
  R
}

# ── Term selection ────────────────────────────────────────────────────────────
gsea_sig <- gsea_all$pdx$df
gsea_sig <- gsea_sig[gsea_sig$p.adjust < gsea_PADJ_CUT, ]
if (gsea_NES_SIDE == "neg") gsea_sig <- gsea_sig[gsea_sig$NES < 0, ]
if (gsea_NES_SIDE == "pos") gsea_sig <- gsea_sig[gsea_sig$NES > 0, ]
gsea_sig <- gsea_sig[order(gsea_sig$p.adjust), ]
if (nrow(gsea_sig) > gsea_MAX_TERMS) gsea_sig <- gsea_sig[seq_len(gsea_MAX_TERMS), ]
message(sprintf("Panel i: clustering %d PDX GSEA terms", nrow(gsea_sig)))

# ── GO accessions from the MSigDB names ───────────────────────────────────────
# Taken from the cached MSigDB download (see get_msig_gobp_go_map()).
gsea_map <- get_msig_gobp_go_map()

gsea_sig$go_id <- gsea_map$go_id[match(gsea_sig$ID, gsea_map$gs_name)]
gsea_sig$go_id[!grepl("^GO:[0-9]{7}$", gsea_sig$go_id)] <- NA
message(sprintf("Panel i: %d/%d terms mapped to GO IDs",
                sum(!is.na(gsea_sig$go_id)), nrow(gsea_sig)))
gsea_sig <- gsea_sig[!is.na(gsea_sig$go_id), ]
gsea_sig <- gsea_sig[!duplicated(gsea_sig$go_id), ]
rownames(gsea_sig) <- gsea_sig$ID
if (nrow(gsea_sig) < gsea_K_CLUST)
  stop(sprintf("Panel i: only %d terms to cluster, fewer than gsea_K_CLUST = %d.",
               nrow(gsea_sig), gsea_K_CLUST))

# ── Jaccard on genes ──────────────────────────────────────────────────────────
gsea_msig_gobp <- as.data.frame(get_msig_gobp())
gsea_gene_lists <- if (gsea_USE_LEADING_EDGE) {
  strsplit(gsea_sig$core_enrichment, "/", fixed = TRUE)
} else {
  lapply(gsea_sig$ID, function(x)
    unique(gsea_msig_gobp$gene_symbol[gsea_msig_gobp$gs_name == x]))
}
names(gsea_gene_lists) <- gsea_sig$ID

gsea_n <- length(gsea_gene_lists)
gsea_S_jac <- matrix(0, gsea_n, gsea_n,
                     dimnames = list(gsea_sig$ID, gsea_sig$ID))
for (i in seq_len(gsea_n)) for (j in i:gsea_n) {
  a <- gsea_gene_lists[[i]]; b <- gsea_gene_lists[[j]]
  v <- length(intersect(a, b)) / length(union(a, b))
  gsea_S_jac[i, j] <- v; gsea_S_jac[j, i] <- v
}
message(sprintf("Panel i: Jaccard median off-diagonal %.3f, max %.3f",
                median(gsea_S_jac[upper.tri(gsea_S_jac)]),
                max(gsea_S_jac[upper.tri(gsea_S_jac)])))

# ── Semantic similarity (Wang, GO DAG topology only) ──────────────────────────
gsea_sem_bp <- godata("org.Hs.eg.db", ont = "BP", computeIC = FALSE)
gsea_S_sem  <- mgoSim(gsea_sig$go_id, gsea_sig$go_id, semData = gsea_sem_bp,
                      measure = "Wang", combine = NULL)
dimnames(gsea_S_sem) <- list(gsea_sig$ID, gsea_sig$ID)
gsea_S_sem[is.na(gsea_S_sem)] <- 0
gsea_S_sem <- gsea_S_sem[rownames(gsea_S_jac), colnames(gsea_S_jac)]
message(sprintf("Panel i: Wang median off-diagonal %.3f, max %.3f",
                median(gsea_S_sem[upper.tri(gsea_S_sem)]),
                max(gsea_S_sem[upper.tri(gsea_S_sem)])))

# ── Combined similarity ───────────────────────────────────────────────────────
gsea_S <- gsea_W_SEM * gsea_rank01(gsea_S_sem) +
  (1 - gsea_W_SEM) * gsea_rank01(gsea_S_jac)
diag(gsea_S) <- 1
message(sprintf("Panel i: combined S, %.1f%% of pairs exactly zero",
                100 * mean(gsea_S[upper.tri(gsea_S)] == 0)))

# ── Ward clusters ─────────────────────────────────────────────────────────────
gsea_hc <- hclust(as.dist(1 - gsea_S), method = "ward.D2")
gsea_sig$cluster_active <- cutree(gsea_hc, k = gsea_K_CLUST)

# ── Small-cluster filter ──────────────────────────────────────────────────────
gsea_keep <- gsea_sig$cluster_active %in%
  names(which(table(gsea_sig$cluster_active) >= gsea_MIN_CLUST))
if (any(!gsea_keep))
  message(sprintf("Panel i: dropping %d terms in clusters of <%d:\n  %s",
                  sum(!gsea_keep), gsea_MIN_CLUST,
                  paste(gsea_sig$term[!gsea_keep], collapse = "\n  ")))
gsea_sig     <- gsea_sig[gsea_keep, ]
gsea_Sf      <- gsea_S[gsea_sig$ID, gsea_sig$ID]

# Automatic name = most significant term in the cluster
gsea_auto_lab <- tapply(seq_len(nrow(gsea_sig)), gsea_sig$cluster_active,
                        function(i) gsub("_", " ", gsea_sig$term[i][which.min(gsea_sig$p.adjust[i])]))

# ── Graph ─────────────────────────────────────────────────────────────────────
gsea_g <- igraph::graph_from_adjacency_matrix(gsea_edge_cut(gsea_Sf, gsea_EDGE_Q),
                                              mode = "undirected", weighted = TRUE)
gsea_vi <- match(igraph::V(gsea_g)$name, gsea_sig$ID)
igraph::V(gsea_g)$cluster    <- gsea_sig$cluster_active[gsea_vi]
igraph::V(gsea_g)$neg_logFDR <- -log10(pmax(gsea_sig$p.adjust[gsea_vi], 1e-300))
igraph::V(gsea_g)$NES        <- gsea_sig$NES[gsea_vi]
igraph::V(gsea_g)$n_leading  <- gsea_sig$n_leading[gsea_vi]

gsea_iso <- which(igraph::degree(gsea_g) == 0)
message(sprintf("Panel i: %d of %d terms have no surviving edges and are dropped",
                length(gsea_iso), igraph::vcount(gsea_g)))
if (length(gsea_iso) > 0) gsea_g <- igraph::delete_vertices(gsea_g, gsea_iso)

gsea_el <- igraph::as_edgelist(gsea_g, names = FALSE)
igraph::E(gsea_g)$same_cluster <-
  igraph::V(gsea_g)$cluster[gsea_el[, 1]] == igraph::V(gsea_g)$cluster[gsea_el[, 2]]
message(sprintf("Panel i: graph %d nodes, %d edges (%d within-cluster), %d clusters",
                igraph::vcount(gsea_g), igraph::ecount(gsea_g),
                sum(igraph::E(gsea_g)$same_cluster),
                length(unique(igraph::V(gsea_g)$cluster))))

# ── Cluster names, with the check table ───────────────────────────────────────
gsea_cl_ids <- sort(unique(igraph::V(gsea_g)$cluster))
gsea_top3 <- sapply(gsea_cl_ids, function(k) {
  s <- gsea_sig[gsea_sig$ID %in% igraph::V(gsea_g)$name & gsea_sig$cluster_active == k, ]
  paste(head(gsub("_", " ", s$term[order(s$p.adjust)]), 3), collapse = " | ")
})
gsea_lab_raw <- if (gsea_USE_MANUAL_LABELS) {
  l <- gsea_manual_labels[as.character(gsea_cl_ids)]
  l[is.na(l)] <- gsea_auto_lab[as.character(gsea_cl_ids)][is.na(l)]
  l
} else gsea_auto_lab[as.character(gsea_cl_ids)]
names(gsea_lab_raw) <- gsea_cl_ids

message("Panel i: cluster names — CHECK the manual name matches the top terms")
print(data.frame(cluster   = gsea_cl_ids,
                 n_terms   = as.vector(table(igraph::V(gsea_g)$cluster)[as.character(gsea_cl_ids)]),
                 label     = unname(gsea_lab_raw),
                 top_terms = unname(gsea_top3)),
      right = FALSE)

gsea_cl_lab <- sapply(gsea_lab_raw, function(x) paste(strwrap(x, 28), collapse = "\n"))
names(gsea_cl_lab) <- gsea_cl_ids

# ── Two-level layout ──────────────────────────────────────────────────────────
#   1. cluster centres placed by a layout on the cluster-cluster similarity
#   2. centres scaled apart until no two cluster discs can overlap
#   3. each cluster's terms laid out inside its own disc
set.seed(99)
gsea_cl_n <- table(igraph::V(gsea_g)$cluster)[as.character(gsea_cl_ids)]

gsea_Sc <- matrix(0, length(gsea_cl_ids), length(gsea_cl_ids),
                  dimnames = list(gsea_cl_ids, gsea_cl_ids))
for (a in seq_along(gsea_cl_ids)) for (b in seq_along(gsea_cl_ids)) {
  ia <- igraph::V(gsea_g)$name[igraph::V(gsea_g)$cluster == gsea_cl_ids[a]]
  ib <- igraph::V(gsea_g)$name[igraph::V(gsea_g)$cluster == gsea_cl_ids[b]]
  gsea_Sc[a, b] <- mean(gsea_Sf[ia, ib])
}
diag(gsea_Sc) <- 0
gsea_g_meta <- igraph::graph_from_adjacency_matrix(gsea_Sc, mode = "undirected",
                                                   weighted = TRUE)

gsea_cen <- igraph::layout_with_fr(gsea_g_meta, weights = igraph::E(gsea_g_meta)$weight)
gsea_cen <- sweep(gsea_cen, 2, colMeans(gsea_cen), "-")
gsea_cen <- gsea_cen / max(abs(gsea_cen))
rownames(gsea_cen) <- as.character(gsea_cl_ids)

gsea_rad <- pmax(gsea_DISC_SCALE * 0.45 * sqrt(as.numeric(gsea_cl_n) / max(gsea_cl_n)), 0.10)
names(gsea_rad) <- as.character(gsea_cl_ids)

gsea_need <- 1
for (a in seq_along(gsea_cl_ids)) for (b in seq_along(gsea_cl_ids)) if (a < b) {
  d <- sqrt(sum((gsea_cen[a, ] - gsea_cen[b, ])^2))
  if (d > 0) gsea_need <- max(gsea_need, (gsea_rad[a] + gsea_rad[b] + gsea_GAP) / d)
}
gsea_cen <- gsea_cen * gsea_need
message(sprintf("Panel i: centres scaled by %.2f to clear %d discs",
                gsea_need, length(gsea_rad)))

gsea_xy_all <- matrix(NA_real_, igraph::vcount(gsea_g), 2,
                      dimnames = list(igraph::V(gsea_g)$name, c("x", "y")))
for (k in as.character(gsea_cl_ids)) {
  idx <- which(as.character(igraph::V(gsea_g)$cluster) == k)
  sub <- igraph::induced_subgraph(gsea_g, idx)
  xy  <- if (igraph::vcount(sub) == 1) matrix(c(0, 0), 1, 2) else
         igraph::layout_with_fr(sub, weights = igraph::E(sub)$weight)
  xy  <- sweep(xy, 2, colMeans(xy), "-")
  r   <- max(sqrt(rowSums(xy^2)))
  if (r > 0) xy <- xy / r
  gsea_xy_all[idx, 1] <- xy[, 1] * gsea_rad[[k]] + gsea_cen[k, 1]
  gsea_xy_all[idx, 2] <- xy[, 2] * gsea_rad[[k]] + gsea_cen[k, 2]
}
gsea_lo <- create_layout(gsea_g, layout = "manual",
                         x = gsea_xy_all[, 1], y = gsea_xy_all[, 2])

# Between-cluster ribbons, strongest gsea_META_KEEP fraction
gsea_mw <- igraph::as_data_frame(gsea_g_meta, what = "edges")
gsea_mw <- gsea_mw[gsea_mw$weight >= quantile(gsea_mw$weight, 1 - gsea_META_KEEP), ]
gsea_meta_seg <- data.frame(
  x    = gsea_cen[gsea_mw$from, 1], y    = gsea_cen[gsea_mw$from, 2],
  xend = gsea_cen[gsea_mw$to,   1], yend = gsea_cen[gsea_mw$to,   2],
  w    = gsea_mw$weight)

gsea_cent <- data.frame(x   = gsea_cen[, 1],
                        y   = gsea_cen[, 2] + gsea_rad[as.character(gsea_cl_ids)] + 0.16,
                        lbl = gsea_cl_lab[as.character(gsea_cl_ids)],
                        stringsAsFactors = FALSE)


gsea_emap_combined <- ggraph(gsea_lo) +
  geom_curve(data = gsea_meta_seg,
             aes(x = x, y = y, xend = xend, yend = yend,
                 linewidth = w, alpha = w),
             curvature = 0.15, colour = "#5b7fa6", lineend = "round",
             inherit.aes = FALSE) +
  geom_edge_link0(aes(alpha = weight, edge_width = weight, filter = same_cluster),
                  colour = "grey55") +
  geom_node_point(aes(fill = factor(cluster), size = n_leading),
                  shape = 21, colour = "grey25", stroke = 0.4) +
  annotate("label", x = gsea_cent$x, y = gsea_cent$y, label = gsea_cent$lbl,
           size = 4, fontface = "bold", lineheight = 0.9,
           fill = scales::alpha("white", 0.92), colour = "grey10",
           linewidth = 0.3, label.r = unit(0.15, "lines")) +
  scale_linewidth_continuous(range = c(0.6, 4), guide = "none") +
  scale_alpha_continuous(range = c(0.18, 0.55), guide = "none") +
  # colour tied to cluster id, not plotting order
  scale_fill_manual(values = setNames(scales::hue_pal(l = 62, c = 70)(length(gsea_cl_ids)),
                                      as.character(gsea_cl_ids)),
                    guide = "none") +
  scale_edge_alpha_continuous(range = c(0.06, 0.40), guide = "none") +
  scale_edge_width_continuous(range = c(0.2, 1.4),  guide = "none") +
  scale_size_continuous(name = "Leading-edge genes", range = c(2, 10)) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0.010, add = 0.03)) +
  scale_y_continuous(expand = expansion(mult = 0.010, add = 0.04)) +
  theme_graph(base_family = "sans") +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "panel_gsea_emap_combined.pdf"), gsea_emap_combined,
       width = 6, height = 6, device = cairo_pdf)


################################################################################
# PANEL I — GSEA GOBP table: CCLE | PDX | HCC
#
# Rows are the union of
# each column's top gsea_TAB_TOP terms by BH FDR, on the side set by
# gsea_TAB_SIDE. Cell values come from that column's FULL GSEA table, so a term
# that is top in one column still shows the others' actual FDR; "n/a" means the
# term was not tested there (outside the 5-2000 size window in that dataset).
# Cell text is the FDR, cell fill is the NES. Terms in more columns' top lists
# sort first, then by best FDR.
################################################################################

gsea_TAB_TOP  <- 10
gsea_TAB_SIDE <- "neg"      # "neg" = fast-growth terms, "pos", or "both"

gsea_tab_cols    <- lapply(gsea_all, `[[`, "df")
gsea_tab_headers <- c(ccle = "CCLE", pdx = "PDX", hcc = "HCC")[names(gsea_tab_cols)]

gsea_tab_tops <- lapply(gsea_tab_cols, function(z) {
  z <- z[is.finite(z$p.adjust), ]
  if (gsea_TAB_SIDE == "neg") z <- z[z$NES < 0, ]
  if (gsea_TAB_SIDE == "pos") z <- z[z$NES > 0, ]
  head(z$ID[order(z$p.adjust, -abs(z$NES))], gsea_TAB_TOP)
})
gsea_tab_ids  <- unique(unlist(gsea_tab_tops, use.names = FALSE))
gsea_tab_ntop <- vapply(gsea_tab_ids, function(id)
  sum(vapply(gsea_tab_tops, function(t) id %in% t, logical(1))), integer(1))

gsea_tab_p   <- sapply(gsea_tab_cols, function(z) z$p.adjust[match(gsea_tab_ids, z$ID)])
gsea_tab_nes <- sapply(gsea_tab_cols, function(z) z$NES[match(gsea_tab_ids, z$ID)])
gsea_tab_p   <- matrix(as.numeric(gsea_tab_p), nrow = length(gsea_tab_ids),
                       dimnames = list(NULL, names(gsea_tab_cols)))
gsea_tab_nes <- matrix(as.numeric(gsea_tab_nes), nrow = length(gsea_tab_ids),
                       dimnames = list(NULL, names(gsea_tab_cols)))

gsea_tab_best <- apply(gsea_tab_p, 1, function(v) { v <- v[is.finite(v)]
                                         if (length(v)) min(v) else Inf })
gsea_tab_ord  <- order(-gsea_tab_ntop, gsea_tab_best)
gsea_tab_ids  <- gsea_tab_ids[gsea_tab_ord]
gsea_tab_p    <- gsea_tab_p[gsea_tab_ord, , drop = FALSE]
gsea_tab_nes  <- gsea_tab_nes[gsea_tab_ord, , drop = FALSE]
gsea_tab_ntop <- gsea_tab_ntop[gsea_tab_ord]

gsea_tab_nrow <- length(gsea_tab_ids); gsea_tab_ncol <- ncol(gsea_tab_p)
message(sprintf("Panel j: %d rows x %d cols", gsea_tab_nrow, gsea_tab_ncol))

# GOBP_DNA_REPLICATION -> "DNA Replication", wrapped rather than cut
gsea_tab_lab <- sub("^GOBP_", "", gsea_tab_ids)
gsea_tab_lab <- tools::toTitleCase(tolower(gsub("_", " ", gsea_tab_lab)))
gsea_tab_lab <- gsub("Rrna", "rRNA", gsea_tab_lab)
gsea_tab_lab <- gsub("Dna", " DNA", gsea_tab_lab)
gsea_tab_lab <- gsub("Rna", " RNA", gsea_tab_lab)
gsea_tab_desc <- gsea_tab_lab
gsea_tab_lab  <- stringr::str_wrap(gsea_tab_lab, width = 30)

gsea_tab_cw <- 0.5 / gsea_tab_ncol
gsea_tab_cx <- 0.5 + (seq_len(gsea_tab_ncol) - 0.5) * gsea_tab_cw
gsea_tab_y  <- rev(seq_len(gsea_tab_nrow))

gsea_tab_names <- data.frame(x = 0, y = gsea_tab_y, label = gsea_tab_lab)

gsea_tab_cells <- do.call(rbind, lapply(seq_len(gsea_tab_ncol), function(j) {
  v   <- gsea_tab_p[, j]
  txt <- rep("\"n/a\"", length(v))
  ok  <- !is.na(v)
  if (any(ok)) txt[ok] <- gsub("e", " %*% 10^",
                               formatC(pmax(v[ok], 1e-300), format = "e", digits = 1))
  data.frame(x = gsea_tab_cx[j], y = gsea_tab_y, nes = gsea_tab_nes[, j], label = txt)
}))
gsea_tab_hdr <- data.frame(x     = c(0, gsea_tab_cx),
                           y     = gsea_tab_nrow + 1,
                           label = c("GO biological process", unname(gsea_tab_headers)),
                           hjust = c(0, rep(0.5, gsea_tab_ncol)))

gsea_table <- ggplot() +
  geom_rect(data = data.frame(y = gsea_tab_y[seq(1, gsea_tab_nrow, by = 2)]),
            aes(xmin = -0.03, xmax = 0.48, ymin = y - 0.5, ymax = y + 0.5),
            fill = "grey94") +
  geom_tile(data = gsea_tab_cells, aes(x = x, y = y, fill = nes),
            width = 0.5 / gsea_tab_ncol * 0.94, height = 0.88,
            colour = "white", linewidth = 0.4) +
  geom_text(data = gsea_tab_names, aes(x = x, y = y, label = label),
            hjust = 0, size = 3.0, lineheight = 0.85) +
  geom_text(data = gsea_tab_cells, aes(x = x, y = y, label = label),
            parse = TRUE, size = 3.0,
            # white text on the strong ends of the NES scale, dark elsewhere
            colour = ifelse(!is.na(gsea_tab_cells$nes) &
                              abs(gsea_tab_cells$nes) /
                                max(abs(gsea_tab_cells$nes), na.rm = TRUE) > 0.55,
                            "white", "grey10")) +
  geom_text(data = gsea_tab_hdr, aes(x = x, y = y, label = label, hjust = hjust),
            fontface = "bold", size = 3.2) +
  geom_hline(yintercept = gsea_tab_nrow + 0.5, linewidth = 0.5, colour = "black") +
  scale_fill_gradient2(low = "#005AB5", mid = "#f7f7f7", high = "#DC3220",
                       midpoint = 0,
                       # symmetric around 0, set by the largest |NES| shown
                       limits = c(-1, 1) * max(abs(gsea_tab_cells$nes), na.rm = TRUE),
                       na.value = "grey92", name = "NES",
                       guide = guide_colourbar(barwidth = 6, barheight = 0.6,
                                               title.position = "top",
                                               title.hjust = 0.5)) +
  scale_x_continuous(limits = c(-0.05, 1.05), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, gsea_tab_nrow + 1.6), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title    = element_text(size = 9, face = "bold"),
        legend.text     = element_text(size = 8),
        plot.margin     = margin(t = 2, r = 2, b = 2, l = 2))
ggsave(file.path(FIG_DIR, "panel_gsea_table.pdf"), gsea_table,
       width = 5, height = 9, device = cairo_pdf)


# Terms and clusters behind panel h (Table 2.3) and the table as drawn in
# panel i (Table 2.4)
gsea_terms_out <- gsea_sig[order(gsea_sig$cluster_active, gsea_sig$p.adjust),
                           c("ID", "term", "go_id", "setSize", "NES", "p.adjust",
                             "n_leading", "cluster_active", "core_enrichment")]
gsea_terms_out$in_map        <- gsea_terms_out$ID %in% igraph::V(gsea_g)$name
gsea_terms_out$cluster_label <- unname(gsea_lab_raw[as.character(gsea_terms_out$cluster_active)])
gsea_tab_out <- data.frame(term = gsea_tab_desc, ID = gsea_tab_ids,
                           n_top_lists = gsea_tab_ntop, check.names = FALSE)
for (j in seq_len(gsea_tab_ncol)) {
  gsea_tab_out[[paste(gsea_tab_headers[j], "FDR")]] <- gsea_tab_p[, j]
  gsea_tab_out[[paste(gsea_tab_headers[j], "NES")]] <- gsea_tab_nes[, j]
}


################################################################################
# ASSEMBLE AND SAVE
################################################################################

# Panel labels: tagged a-i, in the order the panels are added below.
#
#   a  PDX growth curve with exponential fit (panel_a_gc)
#   b  PDX mitotic density vs doubling time (panel_b_mitosis)
#   c  doubling time by tumour type, PDX (dt_by_lineage_pdx)
#   d  PDX  pooled volcano (volcano_pooled_pdx)
#   e  PDX RAN and MKI67 vs doubling time, two facets (scatter_ran_mki67_pdx)
#   f  CCLE pooled volcano (volcano_pooled_ccle)
#   g  HCC  pooled volcano (volcano_pooled_hcc)
#   h  GSEA enrichment map, PDX only (gsea_emap_combined)
#   i  GSEA GOBP table, CCLE | PDX | HCC (gsea_table)
#

FIG1_W           <- 15            # figure width, inches
FIG1_IN_PER_UNIT <- 17.3 / 4.5    # inches per layout-height unit (~3.84 in)

# Every panel is its own area in the design. patchwork fills the areas in
# alphabetical order with the panels in the order they are added; the design
# has no area H, so areas A-G hold panels a-g and areas I and J hold panels h
# and i. tag_levels = "a" letters them a-i. The grid is 21 columns wide (each
# column is 1/21 of the figure width):
#   row 1    a 9 | b 12
#   row 2    c 21
#   rows 3-4 d 6 | e 3 | f 6 | g 6   (all four span both rows)
#   row 5    h 15 | i 6
fig1 <- wrap_elements(full = panel_a_gc) +            # a
  wrap_elements(full = panel_b_mitosis) +                 # b
  wrap_elements(full = dt_by_lineage_pdx) +               # c
  wrap_elements(full = volcano_pooled_pdx) +              # d
  wrap_elements(full = scatter_ran_mki67_pdx) +           # e
  wrap_elements(full = volcano_pooled_ccle) +             # f
  wrap_elements(full = volcano_pooled_hcc) +              # g
  wrap_elements(full = gsea_emap_combined) +              # h
  wrap_elements(full = gsea_table) +                      # i
  plot_layout(
    design = paste(c("AAAAAAAAABBBBBBBBBBBB",
                     "CCCCCCCCCCCCCCCCCCCCC",
                     "DDDDDDEEEFFFFFFGGGGGG",
                     "DDDDDDEEEFFFFFFGGGGGG",
                     "IIIIIIIIIIIIIIIJJJJJJ"),
                   collapse = "\n"),
    ## one entry per design row
    heights = c(1.1,    # a | b
                0.9,    # c
                0.5,    # d | e | f | g (upper half)
                0.5,    # d | e | f | g (lower half)
                1.8)) + # h | i
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

# Height follows the sum of `heights` x FIG1_IN_PER_UNIT
# (5.1 units -> ~19.6 in with the defaults above).
# Heights are read back from the layout, so they are only written once above.
fig1_row_h <- fig1$patches$layout$heights
FIG1_H     <- sum(fig1_row_h) * FIG1_IN_PER_UNIT

# Panel a has a fixed aspect ratio (that of gc.png). The design gives it 9 of
# 21 columns in row 1. This reports how many columns the image would need to
# fill its slot exactly at row 1's height; if it differs, move the A/B boundary
# in the first design row.
message(sprintf("Panel a: fills its slot exactly at %.1f of 21 columns (design uses 9)",
                fig1_row_h[1] * FIG1_IN_PER_UNIT * pa_img_w / pa_img_h / (FIG1_W / 21)))
message(sprintf("figure_1: %.1f x %.1f in", FIG1_W, FIG1_H))
ggsave(file.path(FIG_DIR, "Figure_1.pdf"), fig1,
       width = FIG1_W, height = FIG1_H, limitsize = FALSE)
ggsave(file.path(FIG_DIR, "Figure_1.png"), fig1,
       width = FIG1_W, height = FIG1_H, dpi = 300,
       limitsize = FALSE)
message(sprintf("\n>>> figure_1 -> %s", file.path(FIG_DIR, "Figure_1.pdf")))


################################################################################
################################################################################
##                                                                            ##
##  SUPPLEMENTAL FIGURE 1 — ONE SELF-CONTAINED BLOCK PER PANEL                ##
##                                                                            ##
##  Same style and structure as figure_1 above.                                ##
##                                                                            ##
##    a, b   DT distributions, PDX (a) then CCLE (b), with NO minimum-group-  ##
##           size filter, so ALL cancer types appear.                          ##
##    c, d   volcanoes for the Fisher meta-analysis, PDX (c) then CCLE (d).    ##
##    e, f   ... and the deconfounded (group-mean-subtracted) Spearman,        ##
##           PDX (e) then CCLE (f).                                            ##
##                                                                            ##
##  Panel object names do not follow the printed letters (s1a = CCLE DT,      ##
##  s1b = PDX DT, s1c/s1e = CCLE, s1d/s1f = PDX); the assembly adds them      ##
##  PDX-first.                                                                ##
##                                                                            ##
##  So the three methods appear as: pooled -> Figure 1 d, f, g; Fisher ->     ##
##  supp c/d; deconfounded -> supp e/f. Thresholds, colours, key genes and the ##
##  "3 strongest each way" labelling rule match figure_1 exactly.              ##
##                                                                            ##
##  Object names are prefixed s1a_ ... s1f_ so nothing here can clash with,    ##
##  or be clobbered by, the figure_1 blocks. Re-running a block only needs     ##
##  `fig1_results`, which is already in memory.                                ##
##                                                                            ##
##  Layout       Panel sizes (inches, as used for the standalone PDFs)         ##
##    AAAAAA       a, b  width scales with group count x 4.5                   ##
##    BBBBBB       c-f   4.00 x 4.00                                           ##
##    CCDDEEFF                                                                ##
##                                                                            ##
################################################################################
################################################################################

SUPP_DIR <- FIG_DIR          # set to file.path(FIG_DIR, "supplemental") to split
dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)

# Minimum group size. 1 = keep every cancer type, including singletons.
SUPP_MIN_N <- 1

# Groups with < 3 points get no violin (a density needs >= 2 distinct values);
# their jitter points and box still show.
SUPP_MIN_VIOLIN_N <- 3

# Per-panel width in inches: fixed margin + a strip per group, so the x labels
# stay legible however many cancer types survive the filter.
supp_panel_width <- function(n_groups, per_group = 0.42, base = 2.2,
                             min_w = 12, max_w = 40) {
  min(max(base + per_group * n_groups, min_w), max_w)
}

################################################################################
# SUPP. PANEL B — CCLE doubling-time distribution, ALL lineages (object s1a)
################################################################################

s1a_raw <- fig1_results$ccle$panel_data$dt_data  # columns: sample_id, dt, lineage
s1a_df  <- s1a_raw[is.finite(s1a_raw$dt) & s1a_raw$dt > 0 &
                     !is.na(s1a_raw$lineage) & s1a_raw$lineage != "", ]
s1a_df$lineage <- as.character(s1a_df$lineage)

# No size filter here — this is the whole point of the supplemental panel.
s1a_n  <- table(s1a_df$lineage)
s1a_df <- s1a_df[s1a_df$lineage %in% names(s1a_n[s1a_n >= SUPP_MIN_N]), ]

s1a_med        <- tapply(s1a_df$dt, s1a_df$lineage, median, na.rm = TRUE)
s1a_df$lineage <- factor(s1a_df$lineage, levels = names(sort(s1a_med)))
s1a_df$log_dt  <- log10(s1a_df$dt)

s1a_keep  <- table(s1a_df$lineage)
s1a_xlabs <- setNames(paste0(levels(s1a_df$lineage), "\n(n=",
                             s1a_keep[levels(s1a_df$lineage)], ")"),
                      levels(s1a_df$lineage))
s1a_pal   <- setNames(scales::hue_pal(l = 60, c = 70)(nlevels(s1a_df$lineage)),
                      levels(s1a_df$lineage))

# Violins only where a density is meaningful
s1a_vio <- s1a_df[s1a_df$lineage %in%
                    names(s1a_keep[s1a_keep >= SUPP_MIN_VIOLIN_N]), ]

# Kruskal–Wallis label position. Computed so it clears the tallest group;
# replace with a number in data coords (log10 DT) to pin it manually.
s1a_p_y <- max(s1a_df$log_dt, na.rm = TRUE) +
  0.10 * diff(range(s1a_df$log_dt, na.rm = TRUE))

message(sprintf(">>> Supp. panel a: %d cell lines, %d lineages (all)",
                nrow(s1a_df), nlevels(s1a_df$lineage)))

panel_s1a <- ggplot(s1a_df, aes(x = lineage, y = log_dt, fill = lineage)) +
  geom_jitter(aes(colour = lineage), width = 0.22, alpha = 1.0,
              size = 2, shape = 16, na.rm = TRUE) +
  geom_violin(data = s1a_vio, alpha = 0.55, linewidth = 0.25, scale = "width",
              bw = "nrd0", na.rm = TRUE) +
  geom_boxplot(width = 0.18, fill = "white", outlier.shape = NA,
               linewidth = 0.70, colour = "black", na.rm = TRUE) +
  stat_compare_means(aes(label = after_stat(format_p_value(..p..))),
                     parse = TRUE, na.rm = TRUE, size = 5,
                     vjust = 0.4, hjust = 0.5,
                     label.x = (nlevels(s1a_df$lineage) + 1) / 4,
                     label.y = s1a_p_y) +
  scale_x_discrete(labels = s1a_xlabs) +
  scale_fill_manual(values = s1a_pal, guide = "none") +
  scale_colour_manual(values = s1a_pal, guide = "none") +
  coord_cartesian(clip = "off") +
  xlab(" ") +
  ylab(TeX("$log_{10}(Average\\ DT)\\ (hours)$")) +
  theme(legend.position  = "none",
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(colour = "black"),
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 9,
                                        face = "bold"),
        axis.text.y      = element_text(size = 10, face = "bold"),
        axis.title       = element_text(size = 10, face = "bold"))

s1a_w <- supp_panel_width(nlevels(s1a_df$lineage))
# (the standalone file names follow the object names, not the printed letters)
ggsave(file.path(SUPP_DIR, "panel_a_supp_figure_1.pdf"), panel_s1a,
       width = s1a_w, height = 4.5, device = cairo_pdf, limitsize = FALSE)

################################################################################
# SUPP. PANEL A — PDX doubling-time distribution, ALL tumour types (object s1b)
################################################################################

s1b_raw <- fig1_results$pdx$panel_data$dt_data
s1b_df  <- s1b_raw[is.finite(s1b_raw$dt) & s1b_raw$dt > 0 &
                     !is.na(s1b_raw$lineage) & s1b_raw$lineage != "", ]
s1b_df$lineage <- as.character(s1b_df$lineage)

s1b_n  <- table(s1b_df$lineage)
s1b_df <- s1b_df[s1b_df$lineage %in% names(s1b_n[s1b_n >= SUPP_MIN_N]), ]

s1b_med        <- tapply(s1b_df$dt, s1b_df$lineage, median, na.rm = TRUE)
s1b_df$lineage <- factor(s1b_df$lineage, levels = names(sort(s1b_med)))
s1b_df$log_dt  <- log10(s1b_df$dt)

s1b_keep  <- table(s1b_df$lineage)
s1b_xlabs <- setNames(paste0(levels(s1b_df$lineage), "\n(n=",
                             s1b_keep[levels(s1b_df$lineage)], ")"),
                      levels(s1b_df$lineage))
s1b_pal   <- setNames(scales::hue_pal(l = 60, c = 70)(nlevels(s1b_df$lineage)),
                      levels(s1b_df$lineage))

s1b_vio <- s1b_df[s1b_df$lineage %in%
                    names(s1b_keep[s1b_keep >= SUPP_MIN_VIOLIN_N]), ]

s1b_p_y <- max(s1b_df$log_dt, na.rm = TRUE) +
  0.10 * diff(range(s1b_df$log_dt, na.rm = TRUE))

message(sprintf(">>> Supp. panel b: %d samples, %d tumour types (all)",
                nrow(s1b_df), nlevels(s1b_df$lineage)))

panel_s1b <- ggplot(s1b_df, aes(x = lineage, y = log_dt, fill = lineage)) +
  geom_jitter(aes(colour = lineage), width = 0.22, alpha = 1.0,
              size = 2, shape = 16, na.rm = TRUE) +
  geom_violin(data = s1b_vio, alpha = 0.55, linewidth = 0.25, scale = "width",
              bw = "nrd0", na.rm = TRUE) +
  geom_boxplot(width = 0.18, fill = "white", outlier.shape = NA,
               linewidth = 0.70, colour = "black", na.rm = TRUE) +
  stat_compare_means(aes(label = after_stat(format_p_value(..p..))),
                     parse = TRUE, na.rm = TRUE, size = 5,
                     vjust = 0.4, hjust = 0.5,
                     label.x = (nlevels(s1b_df$lineage) + 1) / 4,
                     label.y = s1b_p_y) +
  scale_x_discrete(labels = s1b_xlabs) +
  scale_fill_manual(values = s1b_pal, guide = "none") +
  scale_colour_manual(values = s1b_pal, guide = "none") +
  coord_cartesian(clip = "off") +
  xlab(" ") +
  ylab(TeX("$log_{10}(Average\\ DT)\\ (days)$")) +
  theme(legend.position  = "none",
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(colour = "black"),
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 9,
                                        face = "bold"),
        axis.text.y      = element_text(size = 10, face = "bold"),
        axis.title       = element_text(size = 10, face = "bold"))

s1b_w <- supp_panel_width(nlevels(s1b_df$lineage))
ggsave(file.path(SUPP_DIR, "panel_b_supp_figure_1.pdf"), panel_s1b,
       width = s1b_w, height = 4.5, device = cairo_pdf, limitsize = FALSE)

################################################################################
# SUPP. PANEL D — CCLE Fisher meta-analysis volcano (object s1c)
################################################################################

s1c_df   <- fig1_results$ccle$fisher_df
s1c_df$x <- s1c_df$rho
s1c_df$y <- -log10(s1c_df$p_value)
s1c_df   <- s1c_df[is.finite(s1c_df$x) & is.finite(s1c_df$y), ]

s1c_df$significance <- "not_significant"
s1c_df$significance[s1c_df$p_adj < 0.1 & s1c_df$x >=  0.1] <- "high_corr"
s1c_df$significance[s1c_df$p_adj < 0.1 & s1c_df$x <= -0.1] <- "low_corr"
s1c_df$significance[s1c_df$gene %in% c("MKI67","PCNA","MCM2","CDKN1A","CDKN1B")] <-
  "key_gene"
s1c_df$significance[s1c_df$gene %in% c("CCNB1","KIF18A")] <-
  "key_gene_2"
s1c_df$significance <- factor(s1c_df$significance,
  levels = c("not_significant", "low_corr", "high_corr", "key_gene", "key_gene_2"))

s1c_key  <- intersect(c("MKI67","PCNA","MCM2","CDKN1A","CDKN1B","CCNB1","KIF18A"), s1c_df$gene)
s1c_cand <- s1c_df[!s1c_df$gene %in% s1c_key, ]
s1c_lab  <- s1c_df[s1c_df$gene %in% c(s1c_key,
                                      head(s1c_cand$gene[order(-s1c_cand$x)], 3),
                                      head(s1c_cand$gene[order( s1c_cand$x)], 3)), ]

s1c_df <- s1c_df[order(match(as.character(s1c_df$significance),
                         c("not_significant","low_corr","high_corr","key_gene","key_gene_2"))), ]

panel_s1c <- ggplot(s1c_df, aes(x = x, y = y, colour = significance)) +
  geom_point(aes(alpha = significance), size = 0.7) +
  geom_text_repel(data = s1c_lab, aes(label = gene, colour = significance),
                  fontface = ifelse(s1c_lab$gene %in% s1c_key, "bold", "plain"),
                  size = 4, box.padding = 0.5, max.overlaps = Inf,
                  point.size = 2, segment.size = 0.2, segment.alpha = 0.5,
                  show.legend = FALSE) +
  scale_colour_manual(values = c("not_significant" = "#808080",
                                 "low_corr"        = "#005AB5",
                                 "high_corr"       = "#DC3220",
                                 "key_gene"        = "#2CA02C",
                                 "key_gene_2"      = "#8E44AD")) +
  scale_alpha_manual(values = c("not_significant" = 0.1, "low_corr" = 1,
                                "high_corr" = 1, "key_gene" = 1,
                                "key_gene_2" = 1), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  ggtitle("CCLE — Fisher") +          # delete this line for an untitled panel
  xlab(TeX("$\\rho$")) +
  ylab(TeX("$-log_{10}(P_{value})$")) +
  theme(legend.position  = "none",
        plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(colour = "black"),
        axis.text        = element_text(size = 14, face = "bold"),
        axis.title       = element_text(size = 18, face = "bold"))
ggsave(file.path(SUPP_DIR, "panel_c_supp_figure_1.pdf"), panel_s1c,
       width = 4, height = 4, device = cairo_pdf)

################################################################################
# SUPP. PANEL C — PDX Fisher meta-analysis volcano (object s1d)
################################################################################

s1d_df   <- fig1_results$pdx$fisher_df
s1d_df$x <- s1d_df$rho
s1d_df$y <- -log10(s1d_df$p_value)
s1d_df   <- s1d_df[is.finite(s1d_df$x) & is.finite(s1d_df$y), ]

s1d_df$significance <- "not_significant"
s1d_df$significance[s1d_df$p_adj < 0.1 & s1d_df$x >=  0.1] <- "high_corr"
s1d_df$significance[s1d_df$p_adj < 0.1 & s1d_df$x <= -0.1] <- "low_corr"
s1d_df$significance[s1d_df$gene %in% c("MKI67","PCNA","MCM2","CDKN1A","CDKN1B")] <-
  "key_gene"
s1d_df$significance[s1d_df$gene %in% c("CCNB1","KIF18A")] <-
  "key_gene_2"
s1d_df$significance <- factor(s1d_df$significance,
  levels = c("not_significant", "low_corr", "high_corr", "key_gene", "key_gene_2"))

s1d_key  <- intersect(c("MKI67","PCNA","MCM2","CDKN1A","CDKN1B","CCNB1","KIF18A"), s1d_df$gene)
s1d_cand <- s1d_df[!s1d_df$gene %in% s1d_key, ]
s1d_lab  <- s1d_df[s1d_df$gene %in% c(s1d_key,
                                      head(s1d_cand$gene[order(-s1d_cand$x)], 3),
                                      head(s1d_cand$gene[order( s1d_cand$x)], 3)), ]

s1d_df <- s1d_df[order(match(as.character(s1d_df$significance),
                         c("not_significant","low_corr","high_corr","key_gene","key_gene_2"))), ]

panel_s1d <- ggplot(s1d_df, aes(x = x, y = y, colour = significance)) +
  geom_point(aes(alpha = significance), size = 0.7) +
  geom_text_repel(data = s1d_lab, aes(label = gene, colour = significance),
                  fontface = ifelse(s1d_lab$gene %in% s1d_key, "bold", "plain"),
                  size = 4, box.padding = 0.5, max.overlaps = Inf,
                  point.size = 2, segment.size = 0.2, segment.alpha = 0.5,
                  show.legend = FALSE) +
  scale_colour_manual(values = c("not_significant" = "#808080",
                                 "low_corr"        = "#005AB5",
                                 "high_corr"       = "#DC3220",
                                 "key_gene"        = "#2CA02C",
                                 "key_gene_2"      = "#8E44AD")) +
  scale_alpha_manual(values = c("not_significant" = 0.1, "low_corr" = 1,
                                "high_corr" = 1, "key_gene" = 1,
                                "key_gene_2" = 1), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  ggtitle("PDX — Fisher") +
  xlab(TeX("$\\rho$")) +
  ylab(TeX("$-log_{10}(P_{value})$")) +
  theme(legend.position  = "none",
        plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(colour = "black"),
        axis.text        = element_text(size = 14, face = "bold"),
        axis.title       = element_text(size = 18, face = "bold"))
ggsave(file.path(SUPP_DIR, "panel_d_supp_figure_1.pdf"), panel_s1d,
       width = 4, height = 4, device = cairo_pdf)

################################################################################
# SUPP. PANEL F — CCLE deconfounded volcano (object s1e)
################################################################################

s1e_df   <- fig1_results$ccle$cor_dc
s1e_df$x <- s1e_df$rho
s1e_df$y <- -log10(s1e_df$p_value)
s1e_df   <- s1e_df[is.finite(s1e_df$x) & is.finite(s1e_df$y), ]

s1e_df$significance <- "not_significant"
s1e_df$significance[s1e_df$p_adj < 0.1 & s1e_df$x >=  0.1] <- "high_corr"
s1e_df$significance[s1e_df$p_adj < 0.1 & s1e_df$x <= -0.1] <- "low_corr"
s1e_df$significance[s1e_df$gene %in% c("MKI67","PCNA","MCM2","CDKN1A","CDKN1B")] <-
  "key_gene"
s1e_df$significance[s1e_df$gene %in% c("CCNB1","KIF18A")] <-
  "key_gene_2"
s1e_df$significance <- factor(s1e_df$significance,
  levels = c("not_significant", "low_corr", "high_corr", "key_gene", "key_gene_2"))

s1e_key  <- intersect(c("MKI67","PCNA","MCM2","CDKN1A","CDKN1B","CCNB1","KIF18A"), s1e_df$gene)
s1e_cand <- s1e_df[!s1e_df$gene %in% s1e_key, ]
s1e_lab  <- s1e_df[s1e_df$gene %in% c(s1e_key,
                                      head(s1e_cand$gene[order(-s1e_cand$x)], 3),
                                      head(s1e_cand$gene[order( s1e_cand$x)], 3)), ]

s1e_df <- s1e_df[order(match(as.character(s1e_df$significance),
                         c("not_significant","low_corr","high_corr","key_gene","key_gene_2"))), ]

panel_s1e <- ggplot(s1e_df, aes(x = x, y = y, colour = significance)) +
  geom_point(aes(alpha = significance), size = 0.7) +
  geom_text_repel(data = s1e_lab, aes(label = gene, colour = significance),
                  fontface = ifelse(s1e_lab$gene %in% s1e_key, "bold", "plain"),
                  size = 4, box.padding = 0.5, max.overlaps = Inf,
                  point.size = 2, segment.size = 0.2, segment.alpha = 0.5,
                  show.legend = FALSE) +
  scale_colour_manual(values = c("not_significant" = "#808080",
                                 "low_corr"        = "#005AB5",
                                 "high_corr"       = "#DC3220",
                                 "key_gene"        = "#2CA02C",
                                 "key_gene_2"      = "#8E44AD")) +
  scale_alpha_manual(values = c("not_significant" = 0.1, "low_corr" = 1,
                                "high_corr" = 1, "key_gene" = 1,
                                "key_gene_2" = 1), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  ggtitle("CCLE — Deconfounded") +
  xlab(TeX("$\\rho$")) +
  ylab(TeX("$-log_{10}(P_{value})$")) +
  theme(legend.position  = "none",
        plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(colour = "black"),
        axis.text        = element_text(size = 14, face = "bold"),
        axis.title       = element_text(size = 18, face = "bold"))
ggsave(file.path(SUPP_DIR, "panel_e_supp_figure_1.pdf"), panel_s1e,
       width = 4, height = 4, device = cairo_pdf)

################################################################################
# SUPP. PANEL E — PDX deconfounded volcano (object s1f)
################################################################################

s1f_df   <- fig1_results$pdx$cor_dc
s1f_df$x <- s1f_df$rho
s1f_df$y <- -log10(s1f_df$p_value)
s1f_df   <- s1f_df[is.finite(s1f_df$x) & is.finite(s1f_df$y), ]

s1f_df$significance <- "not_significant"
s1f_df$significance[s1f_df$p_adj < 0.1 & s1f_df$x >=  0.1] <- "high_corr"
s1f_df$significance[s1f_df$p_adj < 0.1 & s1f_df$x <= -0.1] <- "low_corr"
s1f_df$significance[s1f_df$gene %in% c("MKI67","PCNA","MCM2","CDKN1A","CDKN1B")] <-
  "key_gene"
s1f_df$significance[s1f_df$gene %in% c("CCNB1","KIF18A")] <-
  "key_gene_2"
s1f_df$significance <- factor(s1f_df$significance,
  levels = c("not_significant", "low_corr", "high_corr", "key_gene", "key_gene_2"))

s1f_key  <- intersect(c("MKI67","PCNA","MCM2","CDKN1A","CDKN1B","CCNB1","KIF18A"), s1f_df$gene)
s1f_cand <- s1f_df[!s1f_df$gene %in% s1f_key, ]
s1f_lab  <- s1f_df[s1f_df$gene %in% c(s1f_key,
                                      head(s1f_cand$gene[order(-s1f_cand$x)], 3),
                                      head(s1f_cand$gene[order( s1f_cand$x)], 3)), ]

s1f_df <- s1f_df[order(match(as.character(s1f_df$significance),
                         c("not_significant","low_corr","high_corr","key_gene","key_gene_2"))), ]

panel_s1f <- ggplot(s1f_df, aes(x = x, y = y, colour = significance)) +
  geom_point(aes(alpha = significance), size = 0.7) +
  geom_text_repel(data = s1f_lab, aes(label = gene, colour = significance),
                  fontface = ifelse(s1f_lab$gene %in% s1f_key, "bold", "plain"),
                  size = 4, box.padding = 0.5, max.overlaps = Inf,
                  point.size = 2, segment.size = 0.2, segment.alpha = 0.5,
                  show.legend = FALSE) +
  scale_colour_manual(values = c("not_significant" = "#808080",
                                 "low_corr"        = "#005AB5",
                                 "high_corr"       = "#DC3220",
                                 "key_gene"        = "#2CA02C",
                                 "key_gene_2"      = "#8E44AD")) +
  scale_alpha_manual(values = c("not_significant" = 0.1, "low_corr" = 1,
                                "high_corr" = 1, "key_gene" = 1,
                                "key_gene_2" = 1), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  ggtitle("PDX — Deconfounded") +
  xlab(TeX("$\\rho$")) +
  ylab(TeX("$-log_{10}(P_{value})$")) +
  theme(legend.position  = "none",
        plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(colour = "black"),
        axis.text        = element_text(size = 14, face = "bold"),
        axis.title       = element_text(size = 18, face = "bold"))
ggsave(file.path(SUPP_DIR, "panel_f_supp_figure_1.pdf"), panel_s1f,
       width = 4, height = 4, device = cairo_pdf)

################################################################################
# ASSEMBLE AND SAVE — SUPPLEMENTAL FIGURE 1
################################################################################

# Volcanoes sit in one row (c d e f). For a 2 x 2 block instead, swap the last
# two design rows for "CCCCCCCCCCCCDDDDDDDDDDDD" / "EEEEEEEEEEEEFFFFFFFFFFFF"
# and use heights = c(1, 1, 1.1, 1.1).
# PDX first in each pair: tags follow the order panels are added, so the PDX
# panels are a, c, e and the CCLE panels b, d, f.
supp_fig1 <- panel_s1b +   # a  PDX  DT by tumour type
  panel_s1a +              # b  CCLE DT by lineage
  panel_s1d +              # c  PDX  Fisher volcano
  panel_s1c +              # d  CCLE Fisher volcano
  panel_s1f +              # e  PDX  deconfounded volcano
  panel_s1e +              # f  CCLE deconfounded volcano
  plot_layout(
    design = paste(c("AAAAAAAAAAAAAAAAAAAAAAAA",
                     "BBBBBBBBBBBBBBBBBBBBBBBB",
                     "CCCCCCDDDDDDEEEEEEFFFFFF"),
                   collapse = "\n"),
    heights = c(1, 1, 1.1)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(size = 15,
                                                        face = "bold"))) &
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1))

supp_fig1_w <- max(s1a_w, s1b_w)

ggsave(file.path(SUPP_DIR, "supp_figure_1.pdf"), supp_fig1,
       width = supp_fig1_w, height = 15, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(SUPP_DIR, "supp_figure_1.png"), supp_fig1,
       width = supp_fig1_w, height = 15, dpi = 300, limitsize = FALSE)
message(sprintf("\n>>> supp_figure_1 -> %s",
                file.path(SUPP_DIR, "supp_figure_1.pdf")))


################################################################################
################################################################################
##                                                                            ##
##  SUPPLEMENTAL TABLE 2 — source data for figure 1 and supplemental          ##
##  figure 1                                                                   ##
##                                                                            ##
##  One sheet per panel (or group of panels), in figure order, built from the  ##
##  exact objects the panels are drawn from — so the tables and the figures    ##
##  can never drift. Run this after the figure_1, GSEA and supplemental        ##
##  figure blocks above.                                                       ##
##                                                                            ##
##  Supplemental table 2 (sup_tab2.xlsx) — genes and gene sets                 ##
##    Table 2.1  fig_1 d,f,g,   Gene-level correlations merged on gene: pooled ##
##               supp c-f       (PDX, CCLE, HCC), Fisher and deconfounded      ##
##                              (PDX, CCLE)                                    ##
##    Table 2.2  fig_1 h,i      Full GSEA GOBP results, PDX, CCLE and HCC      ##
##    Table 2.3  fig_1 h        PDX GSEA enrichment map: terms and clusters    ##
##    Table 2.4  fig_1 i        GSEA table exactly as drawn                    ##
##                                                                            ##
################################################################################
################################################################################

# ── Table 2.1 — gene-level correlations (figure_1 d/f/g, supp c-f) ────────────
# One row per gene, one column block per analysis. Column names are
# <method>_<statistic>_<dataset>, e.g. pooled_rho_pdx. The gene sets differ
# between analyses (pooled uses every sample; Fisher and deconfounded are
# restricted to groups of >= MIN_LINEAGE; HCC has its own gene set), so a gene
# absent from an analysis gets NA in that block rather than being dropped.
supp_volcano_block <- function(df, method, dataset) {
  out <- data.frame(
    gene     = as.character(df$gene),
    rho      = df$rho,
    p_value  = df$p_value,
    p_adj_BH = df$p_adj,
    stringsAsFactors = FALSE)
  names(out)[-1] <- sprintf("%s_%s_%s", method, names(out)[-1], dataset)
  out[!duplicated(out$gene), ]
}

# Order: method-major, PDX before CCLE (as in the figures), HCC last.
sup_tab_2_1_parts <- list(
  supp_volcano_block(vol_pdx_df,  "pooled",       "pdx"),    # figure_1 d
  supp_volcano_block(vol_ccle_df, "pooled",       "ccle"),   # figure_1 f
  supp_volcano_block(vol_hcc_df,  "pooled",       "hcc"),    # figure_1 g
  supp_volcano_block(s1d_df,      "fisher",       "pdx"),    # supp c
  supp_volcano_block(s1c_df,      "fisher",       "ccle"),   # supp d
  supp_volcano_block(s1f_df,      "deconfounded", "pdx"),    # supp e
  supp_volcano_block(s1e_df,      "deconfounded", "ccle"))   # supp f

sup_tab_2_1 <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE),
                      sup_tab_2_1_parts)
sup_tab_2_1 <- sup_tab_2_1[order(sup_tab_2_1$gene), ]

message(sprintf("    Table 2.1: %d genes across 7 analyses (%s).",
                nrow(sup_tab_2_1),
                paste(sprintf("%s n=%d",
                              c("pooled PDX", "pooled CCLE", "pooled HCC",
                                "Fisher PDX", "Fisher CCLE",
                                "deconf PDX", "deconf CCLE"),
                              vapply(sup_tab_2_1_parts, nrow, integer(1))),
                      collapse = ", ")))

# ── Table 2.2 — full GSEA GOBP results (behind figure_1 h and i) ──────────────
# Every gene set tested, for each dataset, stacked and tagged by dataset.
# Ranked by pooled Spearman rho (median-centred); NES < 0 = enriched in fast
# growers. core_enrichment = the leading-edge genes.
sup_tab_2_2 <- do.call(rbind, lapply(c("pdx", "ccle", "hcc"),
  function(ds) {
    z <- gsea_all[[ds]]$df
    data.frame(dataset         = toupper(ds),
               ID              = z$ID,
               Description     = z$term,
               setSize         = z$setSize,
               NES             = z$NES,
               pvalue          = z$pvalue,
               p_adjust        = z$p.adjust,
               n_leading_edge  = z$n_leading,
               direction       = z$direction,
               core_enrichment = z$core_enrichment,
               stringsAsFactors = FALSE)[order(z$p.adjust, -abs(z$NES)), ]
  }))
rownames(sup_tab_2_2) <- NULL

# ── Table 2.3 — PDX GSEA enrichment map (figure_1 h) ──────────────────────────
# Only the terms drawn in the panel (in_map), with their cluster and its label.
sup_tab_2_3 <- gsea_terms_out[gsea_terms_out$in_map,
                              c("cluster_active", "cluster_label", "ID", "term",
                                "go_id", "setSize", "NES", "p.adjust",
                                "n_leading", "core_enrichment")]
names(sup_tab_2_3) <- c("cluster", "cluster_label", "ID", "Description",
                        "GO_ID", "setSize", "NES", "p_adjust",
                        "n_leading_edge", "core_enrichment")
sup_tab_2_3$cluster_size <- as.integer(table(sup_tab_2_3$cluster)[as.character(sup_tab_2_3$cluster)])
sup_tab_2_3 <- sup_tab_2_3[order(-sup_tab_2_3$cluster_size, sup_tab_2_3$cluster,
                                 sup_tab_2_3$p_adjust),
                           c("cluster", "cluster_label", "cluster_size", "ID",
                             "Description", "GO_ID", "setSize", "NES", "p_adjust",
                             "n_leading_edge", "core_enrichment")]
rownames(sup_tab_2_3) <- NULL

message(sprintf("    Table 2.3: %d terms in %d clusters (all drawn in figure_1 i).",
                nrow(sup_tab_2_3), length(unique(sup_tab_2_3$cluster))))

# ── Table 2.4 — GSEA table as drawn (figure_1 i) ──────────────────────────────
# Same rows, order and values as the panel: FDR and NES per dataset.
sup_tab_2_4 <- gsea_tab_out
rownames(sup_tab_2_4) <- NULL

## ---- Contents sheets, written by hand --------------------------------------
## One row per sheet. Edit the text freely; keep the vectors the same length
## or data.frame() will recycle and silently mislabel rows.
sup_tab_2_contents <- data.frame(
  sheet = c(
    "Table 2.1",
    "Table 2.2",
    "Table 2.3",
    "Table 2.4"),
  description = c(
    "Correlation of gene expression with doubling time: pan-cancer pooled Spearman for PDX, CCLE and HCC, and Fisher meta-analysis across groups and deconfounded (group-mean-subtracted) Spearman for PDX and CCLE. Columns are named <method>_<statistic>_<dataset>; NA means the gene was not tested by that analysis.",
    "Gene set enrichment analysis (GSEA) of GO Biological Process gene sets for PDX, CCLE and HCC, with genes ranked by pan-cancer pooled Spearman rho with doubling time. NES < 0 means the set is enriched among genes higher in faster-growing samples. core_enrichment lists the leading-edge genes.",
    "PDX GSEA terms shown in the enrichment map (fast-growth side, FDR < 0.1), with the cluster each belongs to; clusters combine GO semantic similarity and leading-edge gene overlap.",
    "GSEA terms shown in the Fig. 1i table, with FDR and NES in each dataset."),
  stringsAsFactors = FALSE)

# -- Workbooks ----------------------------------------------------------------
# One workbook per supplemental table: a Contents sheet, then "Table <t>.<i>"
# sheets taken from the objects sup_tab_<t>_<i> built above. The contents are
# hand-written, so check they still match the tables before anything is written.
supp_write_table <- function(t, n, contents, file) {
  stopifnot(identical(contents$sheet, sprintf("Table %d.%d", t, seq_len(n))))
  wb <- createWorkbook()
  addWorksheet(wb, "Contents")
  writeData(wb, sheet = "Contents", x = contents)
  for (i in seq_len(n)) {
    sh <- sprintf("Table %d.%d", t, i)
    addWorksheet(wb, sh)
    writeData(wb, sheet = sh, x = get(sprintf("sup_tab_%d_%d", t, i)))
  }
  saveWorkbook(wb, file, overwrite = TRUE)
  message(sprintf(">>> supplemental table %d -> %s  (%s)", t, file,
                  paste(sprintf("%d rows",
                                vapply(mget(sprintf("sup_tab_%d_%d", t, seq_len(n)),
                                            envir = globalenv()),
                                       nrow, integer(1))),
                        collapse = ", ")))
}

supp_write_table(2, 4, sup_tab_2_contents, file.path(SUPP_DIR, "sup_tab2.xlsx"))

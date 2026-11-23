################################################################################
# plot_figure_2.R
# Figure 2, Supplementary Figure 2 and Supplementary Table 3: cross-validated
# doubling-time prediction models on PDX and CCLE, the CCLE <-> PDX domain swap,
# and the figures.
#
# Run after plot_figure_1.R (panel 2f reads its Supplemental Table 2).
#
# WHAT THIS SCRIPT WRITES
#   <out_dir>/pipeline_objects.rds          one per dataset (read by the
#                                           figures-only script)
#   <NEW_FIG_DIR>/domain_swap.rds           7.2
#   <NEW_FIG_DIR>/figure2_data.rds          7.3
#   <NEW_FIG_DIR>/panel_f_gene_corr.rds     7.4c
#   <NEW_FIG_DIR>/panel_ab_ensemble.rds     7.4d
#   <NEW_FIG_DIR>/Figure_2a-f_<description>.pdf, Figure_2.pdf / .png
#   <NEW_FIG_DIR>/Sup_Figure_2_panela-c.pdf, Sup_Figure_2.pdf / .png
#   <NEW_FIG_DIR>/sup_tab3.xlsx              (+ sup_tab3_*.csv.gz if a table
#                                            is too long for Excel)
#
# LAYOUT
#   Section 1  Libraries
#   Section 2  Global parameters
#   Section 3  DATASETS + output paths          <- edit paths here
#   Section 4  format_p_value                   <- the only custom plotting fn
#   Section 5  Pipeline functions (modelling)
#   Section 6  PIPELINE RUN (loop over datasets)
#   Section 7  FIGURES  <<<<<<  every plot lives below this line
#
# WHY THE SPLIT
# The pipeline saves everything each dataset produced to
#   <out_dir>/pipeline_objects.rds
# Section 7 loads those files, so you can re-run any figure WITHOUT re-running
# the modelling. Open a fresh R session and run:
#
#     source("plot_figure_2.R", echo = FALSE)    # full run, once
#
# then afterwards, to redo figures only, run plot_figure_2_from_cache.R — or,
# in this file, Sections 1-4, then 7.0 (which loads the .rds files), then
# whichever panel block you want.
#
# EACH PANEL IS SELF-CONTAINED
# Every panel in Section 7 is a standalone block: it pulls its own data, sets
# its own colours and literal values, builds its ggplot, and writes its own
# PDF. No plotting helper functions, nothing inherited from another panel. To
# change a panel, edit only that block. The blocks deliberately repeat code —
# that is the point.
################################################################################

################################################################################
# WORKING DIRECTORY — defaults to the folder this script is in
#
# Every input is read from, and every output written to, a folder inside
# WORK_DIR. This is the same WORK_DIR as plot_figure_1.R, and the input files
# shared with it (PDX expression/DT, CCLE expression) are the same files.
# Expected layout:
#   WORK_DIR/input_files/ccle/       DT_Expression_CCLE_Data.csv,
#                                    Models_DT.csv
#   WORK_DIR/input_files/pdx/        logTPM_matrix_aver.csv,
#                                    PDM_Samples_com_gc.csv
#   WORK_DIR/figure_1/               sup_tab2.xlsx — written by
#                                    plot_figure_1.R; Table 2.1 feeds panel f
#   WORK_DIR/figure_2/               Figure 2, Sup. Figure 2, sup_tab3.xlsx
#                                    and the cached .rds files
#   WORK_DIR/figure_2/pdx/           pipeline_objects.rds for PDX
#   WORK_DIR/figure_2/ccle/          pipeline_objects.rds for CCLE
# (outputs are created by the script)
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
# present are left as they are (no updates). All are on CRAN.
CRAN_PKGS <- c("glmnet", "ggplot2", "openxlsx", "pls", "ranger", "e1071",
               "xgboost", "future", "future.apply", "parallelly", "patchwork",
               "ggpubr", "ggrepel", "latex2exp", "readxl")

.is_installed <- function(pkgs)
  vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)

.missing_pkgs <- CRAN_PKGS[!.is_installed(CRAN_PKGS)]
if (length(.missing_pkgs) > 0) {
  message("Installing CRAN package(s): ", paste(.missing_pkgs, collapse = ", "))
  install.packages(.missing_pkgs)
}
.missing_pkgs <- CRAN_PKGS[!.is_installed(CRAN_PKGS)]
if (length(.missing_pkgs) > 0)
  stop("These packages could not be installed: ",
       paste(.missing_pkgs, collapse = ", "),
       ". Install them manually and re-run.")
rm(.missing_pkgs, .is_installed)

# ── Section 1. Libraries ──────────────────────────────────────────────────────
library(glmnet)
library(ggplot2)
library(openxlsx)
library(pls)
library(ranger)
library(e1071)
library(xgboost)
library(future)
library(future.apply)

suppressPackageStartupMessages({
  PATCHWORK_OK <- requireNamespace("patchwork", quietly = TRUE)
  GGPUBR_OK    <- requireNamespace("ggpubr",    quietly = TRUE)
  LATEX2EXP_OK <- requireNamespace("latex2exp", quietly = TRUE)
  READXL_OK    <- requireNamespace("readxl",    quietly = TRUE)
})
if (PATCHWORK_OK) library(patchwork)
if (GGPUBR_OK)    library(ggpubr)
if (LATEX2EXP_OK) library(latex2exp)

for (.pkg in c(patchwork = PATCHWORK_OK, ggpubr = GGPUBR_OK,
               latex2exp = LATEX2EXP_OK)[c(!PATCHWORK_OK, !GGPUBR_OK, !LATEX2EXP_OK)])
  invisible(NULL)
if (!PATCHWORK_OK) message("INFO: install.packages('patchwork') — needed to assemble multi-panel figures")
if (!GGPUBR_OK)    message("INFO: install.packages('ggpubr') — needed for stat_cor on the scatter panels")
if (!LATEX2EXP_OK) message("INFO: install.packages('latex2exp') — needed for TeX() axis labels")
if (!READXL_OK)    message("INFO: install.packages('readxl') — needed to read Supplemental Table 2 behind panel f")

################################################################################
# Section 2. GLOBAL PARAMETERS
################################################################################

GLOBAL_SEED     <- 123
N_FOLDS         <- 5
TOP_N_GENES_VEC <- c(10, 20, 50, 100, 200, 500)

## Parallel workers: one per available core, leaving one free so the machine
## stays responsive. parallelly::availableCores() also respects limits set by
## a cluster scheduler (e.g. a Slurm job's CPU allocation). Set a number here
## (e.g. N_WORKERS <- 4) to override. Results do not depend on the number of
## workers (future.seed = TRUE).
N_WORKERS <- max(1L, parallelly::availableCores() - 1L)
plan(multisession, workers = N_WORKERS)
message(sprintf("Parallel workers: %d (of %d available cores)",
                N_WORKERS, parallelly::availableCores()))
options(future.globals.maxSize = 4 * 1024^3)

## Save each dataset's objects so Section 7 can run without the pipeline.
## These files are large (they include pred_df). Set FALSE to skip.
SAVE_PIPELINE_OBJECTS <- TRUE

################################################################################
# Section 3. DATASETS AND OUTPUT PATHS   — edit here
#
#   name     : label; used as the key everywhere in Section 7
#   type     : "csv"      -> expression CSV + phenotype CSV
#              "ccle"     -> Models_DT.csv (lineage) + combined CCLE CSV
#   out_dir  : per-dataset outputs (supplementary figures, xlsx, log, .rds)
################################################################################

DATASETS <- list(

  list(
    name        = "PDX_full_sig_3",
    type        = "csv",
    rna_csv     = file.path(WORK_DIR, "input_files", "pdx",
                            "logTPM_matrix_aver.csv"),
    dt_csv      = file.path(WORK_DIR, "input_files", "pdx",
                            "PDM_Samples_com_gc.csv"),
    dt_col      = "dt_days",
    lineage_col = "OncoTreeCode",
    out_dir     = file.path(WORK_DIR, "figure_2", "pdx")
  ),

  list(
    name        = "CCLE_norm",
    type        = "ccle",
    rna_csv     = file.path(WORK_DIR, "input_files", "ccle",
                            "DT_Expression_CCLE_Data.csv"),
    models_csv  = file.path(WORK_DIR, "input_files", "ccle",
                            "Models_DT.csv"),
    dt_col      = "Avg_2DT",
    lineage_col = "OncotreeLineage",
    out_dir     = file.path(WORK_DIR, "figure_2", "ccle")
  )
)

## Where Figure 2, Supplementary Figure 2, their standalone panels, the
## cached .rds files and sup_tab3.xlsx are written.
NEW_FIG_DIR <- file.path(WORK_DIR, "figure_2")

################################################################################
## PANEL f — gene-level rho, CCLE vs PDX
##
## This panel does NOT come from the modelling pipeline. It reads the per-gene
## pooled pan-cancer rho for CCLE and PDX from Supplemental Table 2.1
## (sup_tab2.xlsx, written by plot_figure_1.R) and compares them gene by
## gene. Nothing here touches PIPE, pred_df or any pipeline_objects.rds, which
## is why panel f can be built from a cold session in seconds — see 7.4c.
##
##   min_overlap : the highlighted set is the smallest top-|rho| depth at which
##                 the CCLE and PDX lists share at least this many genes
##   max_depth   : upper bound on that search; the number of genes in the table
################################################################################
GENE_CORR <- list(
  xlsx         = file.path(WORK_DIR, "figure_1", "sup_tab2.xlsx"),
  sheet        = "Table 2.1",
  gene_col     = "gene",
  ccle_rho_col = "pooled_rho_ccle",
  pdx_rho_col  = "pooled_rho_pdx",
  min_overlap = 100,
  max_depth   = 18507)

.out_dirs <- vapply(DATASETS, function(d) d$out_dir, character(1))
if (anyDuplicated(.out_dirs))
  stop("Duplicate out_dir in DATASETS — each dataset needs its own folder.")
rm(.out_dirs)

################################################################################
# Section 4. format_p_value
#
# The only custom function used by the plotting code. Returns a plotmath
# expression string; consume it with parse = TRUE.
#
# NOTE: cor.test(method = "spearman") underflows to exactly p = 0 on strong
# correlations with large n. Unguarded, that renders as "0 %*% 10^+00"
# (P = 0 x 10^0), which is a storage artefact rather than a result. The guard
# inside turns it into "P < 10^-308", the limit of double precision.
################################################################################

format_p_value <- function(p, threshold = 0.001, scientific = TRUE) {
  to_sci <- function(x, digits = 2) {
    sci <- formatC(x, format = "e", digits = digits)
    gsub("e", " %*% 10^", sci)
  }
  vapply(p, function(pi) {
    if (is.na(pi)) return(NA_character_)   # guard against NA
    ## UNDERFLOW GUARD — leave this ON.
    ## A double cannot hold anything below about 1e-308, so a p-value past that
    ## comes back as exactly 0 and an unguarded formatter renders it
    ## "italic(P) == 0 %*% 10^+00", i.e. "P = 0 x 10^0". That is a storage
    ## artefact, not a result: nothing was measured as zero.
    ##
    ## Panel f is where this bites — about 18,000 genes at rho ~ 0.6 gives
    ## t ~ 101 on ~18,000 df, which puts the true p near 10^-1769. The exponent
    ## can be recovered on the log scale (pt(..., log.p = TRUE)); see the note
    ## in 7.10. It is not reported here because a t-approximation quoted to
    ## 1,700 decimal places is a formal extrapolation, not a measurement.
    ##
    ## The floor is the true limit of double precision rather than R's
    ## conventional 2.2e-16, because several panels here legitimately report
    ## values far below 2.2e-16 (panel a shows ~10^-294). A floor of 2.2e-16
    ## sitting next to a measured 10^-294 would read as the LESS significant of
    ## the two, which is backwards. 10^-308 keeps every panel on one scale and
    ## is a true statement: the value is smaller than anything representable.
    if (pi <= 0) return("italic(P) < 10^-308")
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

################################################################################
# Section 5. PIPELINE FUNCTIONS  (modelling — not plotting)
################################################################################

load_dataset <- function(ds) {

  message(">>> Loading data...")

  if (ds$type == "csv") {

    rna_raw <- read.csv(ds$rna_csv, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(rna_raw) <- rna_raw[[1]];  rna_raw <- rna_raw[, -1, drop = FALSE]

    dt_raw  <- read.csv(ds$dt_csv,  check.names = FALSE, stringsAsFactors = FALSE)
    rownames(dt_raw) <- dt_raw[[1]];   dt_raw  <- dt_raw[, -1, drop = FALSE]

    common_samples <- intersect(rownames(rna_raw), rownames(dt_raw))
    if (length(common_samples) == 0)
      stop("No common sample IDs between RNA CSV and DT CSV.")

    rna_aligned <- rna_raw[common_samples, , drop = FALSE]
    dt_aligned  <- dt_raw[common_samples,  , drop = FALSE]

    expr_mat    <- as.data.frame(rna_aligned)
    dt_vec      <- as.numeric(dt_aligned[[ds$dt_col]])
    lineage_vec <- dt_aligned[[ds$lineage_col]]

    message(sprintf("    %d common samples.", length(common_samples)))

  } else if (ds$type == "ccle") {

    ## row.names = 2 then drop the first and third remaining columns.
    rna_data <- read.csv(ds$rna_csv, row.names = 2)[, -c(1, 3)]
    rna_data <- rna_data[!is.na(rna_data[[ds$dt_col]]), , drop = FALSE]

    Models       <- read.csv(ds$models_csv, header = TRUE, row.names = 1)
    if (!ds$lineage_col %in% colnames(Models))
      stop("Models file has no ", ds$lineage_col, " column; found: ",
           paste(colnames(Models), collapse = ", "))
    common_cells <- intersect(rownames(rna_data), rownames(Models))
    if (length(common_cells) == 0)
      stop("No common cell lines between CCLE CSV and Models_DT.csv.")

    rna_aligned    <- rna_data[common_cells, , drop = FALSE]
    models_aligned <- Models[common_cells, , drop = FALSE]

    gene_cols   <- setdiff(colnames(rna_aligned), ds$dt_col)
    expr_mat    <- as.data.frame(rna_aligned[, gene_cols, drop = FALSE])
    dt_vec      <- rna_aligned[[ds$dt_col]]
    lineage_vec <- as.character(models_aligned[[ds$lineage_col]])

    message(sprintf("    %d common cell lines.", length(common_cells)))

  } else {
    stop("Unknown dataset type: ", ds$type)
  }

  message(sprintf("    %d genes.", ncol(expr_mat)))

  ok_all      <- is.finite(dt_vec) & !is.na(lineage_vec)
  expr_mat    <- expr_mat[ok_all, , drop = FALSE]
  dt_vec      <- dt_vec[ok_all]
  lineage_vec <- lineage_vec[ok_all]

  message(sprintf("    %d usable samples after filtering.", length(dt_vec)))

  list(expr_mat = expr_mat, dt_vec = dt_vec, lineage_vec = lineage_vec)
}

################################################################################
# Section 5 (continued). HELPER FUNCTIONS
################################################################################

## run_spearman: Spearman rho of every gene with DT, computed as the Pearson
## correlation of ranks (fast for many genes), with t-distribution p-values and
## BH FDR. Genes with < min_obs finite values or zero variance are skipped, and
## only samples with every remaining gene finite are used.
run_spearman <- function(expr_matrix, dt_vector, label = "all", min_obs = 10) {
  expr_matrix <- as.matrix(expr_matrix); storage.mode(expr_matrix) <- "numeric"
  finite_dt   <- is.finite(dt_vector)
  keep <- apply(expr_matrix, 2, function(col)
    sum(finite_dt & is.finite(col)) >= min_obs && var(col, na.rm = TRUE) > 0)
  expr_matrix <- expr_matrix[, keep, drop = FALSE]
  message(sprintf("    %-30s | %d samples | %d genes",
                  label, nrow(expr_matrix), ncol(expr_matrix)))
  ok_rows     <- finite_dt & apply(expr_matrix, 1, function(r) all(is.finite(r)))
  expr_ok     <- expr_matrix[ok_rows, , drop = FALSE]
  dt_ok       <- dt_vector[ok_rows]
  n           <- nrow(expr_ok)
  ranked_expr <- apply(expr_ok, 2, rank, ties.method = "average")
  ranked_dt   <- rank(dt_ok, ties.method = "average")
  rho_vec     <- as.numeric(cor(ranked_expr, ranked_dt, method = "pearson"))
  t_stat      <- rho_vec * sqrt((n - 2) / pmax(1 - rho_vec^2, 1e-15))
  p_vec       <- 2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)
  df <- data.frame(group = label, gene = colnames(expr_matrix),
                   n_cells = n, rho = rho_vec, p_value = p_vec,
                   stringsAsFactors = FALSE)
  df$p_adj <- p.adjust(df$p_value, method = "BH")
  df[order(df$p_value), c("group","gene","n_cells","rho","p_value","p_adj")]
}

## select_top_genes: the top_n genes by p-value among negative-rho, positive-rho
## or all genes.
select_top_genes <- function(cor_df, direction = c("negative","positive","both"),
                             top_n = 50) {
  direction <- match.arg(direction)
  df <- cor_df[is.finite(cor_df$rho) & is.finite(cor_df$p_value), ]
  if (direction == "negative") { df <- df[df$rho < 0, ]; df <- df[order(df$p_value), ] }
  if (direction == "positive") { df <- df[df$rho > 0, ]; df <- df[order(df$p_value), ] }
  if (direction == "both")     { df <- df[order(df$p_value), ] }
  head(df$gene, top_n)
}

## prepare_xy: training / validation matrices for one fold. Missing values get
## the TRAINING median, genes flat in training are dropped, and both sets are
## z-scored with the TRAINING mean and SD.
prepare_xy <- function(expr_train, expr_valid, dt_train, dt_valid, genes) {
  genes   <- intersect(genes, colnames(expr_train))
  x_train <- as.matrix(expr_train[, genes, drop=FALSE]); storage.mode(x_train) <- "numeric"
  x_valid <- as.matrix(expr_valid[, genes, drop=FALSE]); storage.mode(x_valid) <- "numeric"
  med_train <- apply(x_train, 2, median, na.rm = TRUE)
  for (j in seq_along(genes)) {
    x_train[!is.finite(x_train[, j]), j] <- med_train[j]
    x_valid[!is.finite(x_valid[, j]), j] <- med_train[j]
  }
  sds  <- apply(x_train, 2, sd, na.rm = TRUE); keep <- is.finite(sds) & sds > 0
  x_train <- x_train[, keep, drop=FALSE]; x_valid <- x_valid[, keep, drop=FALSE]
  genes   <- genes[keep]
  mu  <- colMeans(x_train, na.rm = TRUE); sdv <- apply(x_train, 2, sd, na.rm = TRUE)
  x_train_scaled <- sweep(sweep(x_train, 2, mu, "-"), 2, sdv, "/")
  x_valid_scaled <- sweep(sweep(x_valid, 2, mu, "-"), 2, sdv, "/")
  list(x_train = x_train_scaled, x_valid = x_valid_scaled,
       y_train = dt_train, y_valid = dt_valid, genes = genes)
}

## calc_metrics: n, RMSE, MAE, R^2, Pearson and Spearman of predicted vs observed
calc_metrics <- function(y_true, y_pred) {
  ok <- is.finite(y_true) & is.finite(y_pred)
  y_true <- y_true[ok]; y_pred <- y_pred[ok]
  data.frame(n = length(y_true),
             RMSE     = sqrt(mean((y_true - y_pred)^2)),
             MAE      = mean(abs(y_true - y_pred)),
             R2       = 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2),
             Pearson  = suppressWarnings(cor(y_true, y_pred, method = "pearson")),
             Spearman = suppressWarnings(cor(y_true, y_pred, method = "spearman")))
}

## fit_models_for_gene_set: fit all 13 models on one fold's training samples and
## gene set, and return train/validation metrics, the validation predictions
## and the coefficients / importances. Models: Linear; Ridge; Lasso; Elastic
## Net (alpha 0.25, 0.5, 0.75); PLS (1, 2, 5 components); RandomForest;
## SVR linear and radial; XGBoost. Seeds are offset by fold and gene count.
fit_models_for_gene_set <- function(expr_train, expr_valid, dt_train, dt_valid,
                                    genes, ranking_method, gene_direction,
                                    fold, top_n_genes) {
  dat        <- prepare_xy(expr_train, expr_valid, dt_train, dt_valid, genes)
  x_train    <- dat$x_train; x_valid <- dat$x_valid
  y_train    <- dat$y_train; y_valid <- dat$y_valid; genes_used <- dat$genes
  if (ncol(x_train) < 2) {
    warning(sprintf("Too few usable genes for fold %d / top %d / %s / %s",
                    fold, top_n_genes, ranking_method, gene_direction))
    return(NULL)
  }
  results <- list(); predictions <- list(); coef_records <- list()

  add_model_result <- function(model_name, pred_train, pred_valid) {
    results[[model_name]] <<- rbind(
      cbind(dataset = "train",      calc_metrics(y_train, pred_train)),
      cbind(dataset = "validation", calc_metrics(y_valid, pred_valid)))
    predictions[[model_name]] <<- data.frame(
      fold = fold, top_n_genes = top_n_genes, model = model_name,
      ranking_method = ranking_method, gene_direction = gene_direction,
      cell_line = rownames(expr_valid),          # sample ID (PDX samples too)
      observed_DT = y_valid, predicted_DT = pred_valid, stringsAsFactors = FALSE)
  }
  add_coef_record <- function(model_name, values, type = "coefficient") {
    if (is.null(values)) return(NULL)
    coef_records[[paste(model_name, type, sep="_")]] <<- data.frame(
      fold = fold, top_n_genes = top_n_genes,
      ranking_method = ranking_method, gene_direction = gene_direction,
      model = model_name, feature = names(values),
      value = as.numeric(values), type = type, stringsAsFactors = FALSE)
  }

  # Linear
  x_lm_train <- cbind(intercept = 1, x_train); x_lm_valid <- cbind(intercept = 1, x_valid)
  lm_fit <- lm.fit(x = x_lm_train, y = y_train)
  add_model_result("Linear",
                   as.numeric(x_lm_train %*% lm_fit$coefficients),
                   as.numeric(x_lm_valid %*% lm_fit$coefficients))
  add_coef_record("Linear", lm_fit$coefficients)

  # Ridge
  cv_ridge <- cv.glmnet(x_train, y_train, alpha=0, nfolds=5, standardize=FALSE)
  add_model_result("Ridge",
                   as.numeric(predict(cv_ridge, x_train, s="lambda.min")),
                   as.numeric(predict(cv_ridge, x_valid, s="lambda.min")))
  add_coef_record("Ridge", as.matrix(coef(cv_ridge, s="lambda.min"))[,1])

  # Lasso
  cv_lasso <- cv.glmnet(x_train, y_train, alpha=1, nfolds=5, standardize=FALSE)
  add_model_result("Lasso",
                   as.numeric(predict(cv_lasso, x_train, s="lambda.min")),
                   as.numeric(predict(cv_lasso, x_valid, s="lambda.min")))
  add_coef_record("Lasso", as.matrix(coef(cv_lasso, s="lambda.min"))[,1])

  # Elastic Net
  for (a in c(0.25, 0.5, 0.75)) {
    mn <- paste0("ElasticNet_alpha_", a)
    cv_enet <- cv.glmnet(x_train, y_train, alpha=a, nfolds=5, standardize=FALSE)
    add_model_result(mn,
                     as.numeric(predict(cv_enet, x_train, s="lambda.min")),
                     as.numeric(predict(cv_enet, x_valid, s="lambda.min")))
    add_coef_record(mn, as.matrix(coef(cv_enet, s="lambda.min"))[,1])
  }

  # PLS
  for (nc in c(1, 2, 5)) {
    nc_use <- min(nc, ncol(x_train), nrow(x_train)-1)
    mn <- paste0("PLS_", nc, "_components")
    tryCatch({
      pls_fit <- plsr(y_train ~ x_train, ncomp=nc_use, validation="none")
      add_model_result(mn,
                       as.numeric(predict(pls_fit, x_train, ncomp=nc_use)),
                       as.numeric(predict(pls_fit, x_valid, ncomp=nc_use)))
      pls_coef <- as.numeric(coef(pls_fit, ncomp=nc_use))
      names(pls_coef) <- colnames(x_train)
      add_coef_record(mn, pls_coef)
    }, error = function(e)
      message(sprintf("    [PLS %d comp] error: %s", nc, conditionMessage(e))))
  }

  # Random Forest — matrix interface avoids formula issues with gene names
  set.seed(1000 + fold + top_n_genes)
  rf_fit <- ranger(x = x_train, y = y_train,
                   num.trees = 1000,
                   mtry = max(1, floor(sqrt(ncol(x_train)))),
                   min.node.size = 5, importance = "impurity",
                   seed = 1000 + fold + top_n_genes)
  add_model_result("RandomForest",
                   predict(rf_fit, data = x_train)$predictions,
                   predict(rf_fit, data = x_valid)$predictions)
  add_coef_record("RandomForest", rf_fit$variable.importance, "importance")

  # SVR — sanitize column names for formula interface
  safe_cn      <- make.names(colnames(x_train), unique = TRUE)
  x_tr_svr     <- x_train; colnames(x_tr_svr) <- safe_cn
  x_va_svr     <- x_valid; colnames(x_va_svr)  <- safe_cn
  svr_train_df <- data.frame(y = y_train, x_tr_svr, check.names = FALSE)
  svr_valid_df <- data.frame(x_va_svr, check.names = FALSE)

  set.seed(2000 + fold + top_n_genes)
  tune_linear <- tune.svm(y ~ ., data=svr_train_df, kernel="linear",
                           cost=c(0.1,1,10), epsilon=0.1,
                           tunecontrol=tune.control(cross=3))
  add_model_result("SVR_linear",
                   as.numeric(predict(tune_linear$best.model, svr_train_df)),
                   as.numeric(predict(tune_linear$best.model, svr_valid_df)))

  set.seed(3000 + fold + top_n_genes)
  tune_radial <- tune.svm(y ~ ., data=svr_train_df, kernel="radial",
                           cost=c(0.1,1,10),
                           gamma=c(0.01, 0.1, 1/ncol(x_train)),
                           epsilon=0.1, tunecontrol=tune.control(cross=3))
  add_model_result("SVR_radial",
                   as.numeric(predict(tune_radial$best.model, svr_train_df)),
                   as.numeric(predict(tune_radial$best.model, svr_valid_df)))

  # XGBoost
  dtrain <- xgb.DMatrix(data=x_train, label=y_train)
  dvalid <- xgb.DMatrix(data=x_valid, label=y_valid)
  xgb_params <- list(objective="reg:squarederror", eval_metric="rmse",
                     eta=0.05, max_depth=2, min_child_weight=5,
                     subsample=0.8, colsample_bytree=0.8, lambda=5, alpha=1)
  set.seed(4000 + fold + top_n_genes)
  xgb_cv <- xgb.cv(params=xgb_params, data=dtrain, nrounds=500, nfold=5,
                   early_stopping_rounds=20, verbose=0)
  best_nrounds <- xgb_cv$best_iteration
  if (is.null(best_nrounds) || length(best_nrounds) == 0 || best_nrounds < 1) best_nrounds <- 50L
  xgb_fit <- xgb.train(params=xgb_params, data=dtrain, nrounds=best_nrounds, verbose=0)
  add_model_result("XGBoost",
                   as.numeric(predict(xgb_fit, dtrain)),
                   as.numeric(predict(xgb_fit, dvalid)))
  xgb_imp <- xgb.importance(feature_names=colnames(x_train), model=xgb_fit)
  if (!is.null(xgb_imp) && nrow(xgb_imp) > 0) {
    xgb_gain <- xgb_imp$Gain; names(xgb_gain) <- xgb_imp$Feature
    add_coef_record("XGBoost", xgb_gain, "importance_gain")
  }

  res_df <- do.call(rbind, results)
  res_df$model          <- rep(names(results), each=2)
  res_df$fold           <- fold; res_df$top_n_genes <- top_n_genes
  res_df$ranking_method <- ranking_method; res_df$gene_direction <- gene_direction
  res_df$n_genes_requested <- length(genes); res_df$n_genes_used <- length(genes_used)
  res_df$genes_used <- paste(genes_used, collapse=";")
  coef_df <- if (length(coef_records) > 0) do.call(rbind, coef_records) else data.frame()
  list(metrics=res_df, predictions=do.call(rbind, predictions),
       coefficients=coef_df, genes_used=genes_used)
}
################################################################################
# Section 6. PIPELINE RUN
#
# Modelling only — no plotting. Each dataset's objects are collected into PIPE
# and written to <out_dir>/pipeline_objects.rds so Section 7 can rebuild any
# figure without re-running this.
#
# The loop body stays at TOP LEVEL (not wrapped in a function) because
# future_lapply() resolves its global exports against the calling environment.
################################################################################

## Per-dataset results, keyed by dataset name. Declared before the snapshot
## below so it survives the per-dataset cleanup.
PIPE <- list()

BASE_OBJECTS <- c(ls(), "BASE_OBJECTS", "ds_idx", "ds", "RUN_STATUS")

RUN_STATUS <- data.frame(
  dataset = vapply(DATASETS, function(d) d$name,    character(1)),
  out_dir = vapply(DATASETS, function(d) d$out_dir, character(1)),
  status  = NA_character_,
  stringsAsFactors = FALSE)

for (ds_idx in seq_along(DATASETS)) {

ds <- DATASETS[[ds_idx]]

cat(sprintf("\n\n================================================================\n"))
cat(sprintf(" PIPELINE %d / %d : %s\n", ds_idx, length(DATASETS), ds$name))
cat(sprintf("================================================================\n"))

set.seed(GLOBAL_SEED)

OUT_DIR <- ds$out_dir
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


RUN_STATUS$status[ds_idx] <- tryCatch({

message(sprintf("Log started: %s", Sys.time()))
message(sprintf("R version: %s", R.version$version.string))
message(sprintf("Dataset: %s  (type: %s)", ds$name, ds$type))
message(sprintf("OUT_DIR: %s", OUT_DIR))

################################################################################
# 6.1  LOAD DATA
################################################################################

.loaded     <- load_dataset(ds)
expr_mat    <- .loaded$expr_mat
dt_vec      <- .loaded$dt_vec
lineage_vec <- .loaded$lineage_vec
rm(.loaded)

################################################################################
# 6.2  CREATE 5-FOLD SPLITS (random, seeded by GLOBAL_SEED above)
################################################################################

message(">>> Preparing 5-fold cross-validation...")
n       <- length(dt_vec)
fold_id <- sample(rep(1:N_FOLDS, length.out = n))
message(sprintf("    Total usable samples: %d", n))
message("    Fold sizes:"); print(table(fold_id))

################################################################################
# 6.3  FOLD LOOP — gene ranking (cor_all only)
#
# Genes are ranked on each fold's TRAINING samples only, so the validation
# samples never influence gene selection.
#
# Only the cor_all ranking is computed: it is the one every Figure 2 /
# Sup Figure 2 / Sup Table 3 output uses. Nothing in this loop draws random
# numbers, so the RNG state at the parallel dispatch below depends only on
# GLOBAL_SEED and the fold split.
################################################################################

message("\n>>> Running per-fold cor_all rankings...")

all_metrics <- list(); all_predictions <- list()

## Each fold's ranking is written to a CSV and read back (this fixes the gene
## order exactly as in the published analysis). The CSVs go to a temporary
## folder that R deletes at the end of the session.
RANK_TMP <- file.path(tempdir(), paste0("rankings_", ds$name))
dir.create(RANK_TMP, recursive = TRUE, showWarnings = FALSE)

for (fold in 1:N_FOLDS) {
  message(sprintf("\n>>> Fold %d / %d", fold, N_FOLDS))

  train_idx  <- which(fold_id != fold)
  expr_train <- expr_mat[train_idx, , drop=FALSE]
  dt_train   <- dt_vec[train_idx]

  cor_all_train <- run_spearman(expr_train, dt_train,
                                label = paste0("cor_all_train_fold_", fold))
  write.csv(cor_all_train,
            file.path(RANK_TMP, sprintf("fold_%d_ranking_cor_all_train.csv", fold)),
            row.names=FALSE)
}

################################################################################
# 6.4  FLAT PARALLEL DISPATCH
#
# WHY THE JOB TABLE HAS MORE ROWS THAN ARE FITTED:
# future_lapply(future.seed = TRUE) gives job i its own random stream, chosen
# by its POSITION i in the table, and cv.glmnet (Ridge) draws its internal
# folds from that stream. The table lists every ranking x direction x gene
# count combination, but only the cor_all / "both" rows (run = TRUE) are fitted;
# the others return NULL immediately. Keeping them holds every fitted job at
# the position — and therefore the random stream — that produced the published
# results.
################################################################################

FIG_RANKING   <- "cor_all"   # must match f2_ranking in 7.3
FIG_DIRECTION <- "both"      # must match f2_direction in 7.3

fold_rankings <- lapply(seq_len(N_FOLDS), function(fold)
  read.csv(file.path(RANK_TMP, sprintf("fold_%d_ranking_cor_all_train.csv", fold))))

all_jobs <- do.call(rbind, lapply(seq_len(N_FOLDS), function(fold) {
  cor_all <- fold_rankings[[fold]]
  do.call(rbind, lapply(TOP_N_GENES_VEC, function(top_n) {
    do.call(rbind, lapply(c("cor_all", "Fisher", "Deconfounded"), function(rname) {
      do.call(rbind, lapply(c("negative", "positive", "both"), function(dir) {
        run   <- rname == FIG_RANKING && dir == FIG_DIRECTION
        genes <- character(0)
        if (rname == "cor_all") {
          g <- select_top_genes(cor_all, direction = dir, top_n = top_n)
          if (length(g) == 0) return(NULL)   # no genes -> no row
          if (run) genes <- g
        }
        data.frame(fold=fold, ranking_method=rname, gene_direction=dir,
                   top_n_genes=top_n, run=run,
                   gene_list=I(list(genes)), stringsAsFactors=FALSE)
      }))
    }))
  }))
}))
message(sprintf("    Job grid: %d rows, %d fitted (%s / %s)",
                nrow(all_jobs), sum(all_jobs$run), FIG_RANKING, FIG_DIRECTION))

message("\n>>> Dispatching ML jobs in parallel...")

all_task_results <- future_lapply(
  seq_len(nrow(all_jobs)),
  function(i) {
    jb <- all_jobs[i, ]
    if (!jb$run) return(NULL)                # not fitted; keeps the seed positions
    fold <- jb$fold; genes <- jb$gene_list[[1]]
    train_idx  <- which(fold_id != fold); valid_idx  <- which(fold_id == fold)
    expr_train <- expr_mat[train_idx, , drop=FALSE]
    expr_valid <- expr_mat[valid_idx, , drop=FALSE]
    dt_train   <- dt_vec[train_idx]; dt_valid <- dt_vec[valid_idx]
    fit_res <- fit_models_for_gene_set(
      expr_train, expr_valid, dt_train, dt_valid, genes,
      jb$ranking_method, jb$gene_direction, fold, jb$top_n_genes)
    if (is.null(fit_res)) return(NULL)
    out_name <- paste0("fold_", fold, "_", jb$ranking_method,
                       "_", jb$gene_direction, "_top", jb$top_n_genes)
    list(out_name=out_name, metrics=fit_res$metrics,
         predictions=fit_res$predictions)
  },
  future.seed = TRUE,
  future.packages = c("glmnet","pls","ranger","e1071","xgboost")
)

all_task_results <- all_task_results[!sapply(all_task_results, is.null)]
for (res in all_task_results) {
  all_metrics[[res$out_name]]     <- res$metrics
  all_predictions[[res$out_name]] <- res$predictions
}
rm(all_task_results, fold_rankings)
unlink(RANK_TMP, recursive = TRUE)

################################################################################
# 6.5  COLLECT RESULTS
################################################################################

message("\n>>> Collecting cross-validation results...")

metrics_df <- do.call(rbind, all_metrics)
pred_df    <- do.call(rbind, all_predictions)

################################################################################
# 6.6  SUMMARIZE PERFORMANCE ACROSS FOLDS  (summary_metrics_flat -> Figures 2a/2b)
################################################################################

message("\n>>> Summarizing validation performance across folds...")

valid_metrics <- metrics_df[metrics_df$dataset == "validation", ]

topn_metrics <- valid_metrics[!is.na(valid_metrics$top_n_genes), ]

summary_metrics <- aggregate(
  cbind(RMSE, MAE, R2, Pearson, Spearman) ~
    top_n_genes + ranking_method + gene_direction + model,
  data = topn_metrics,
  FUN  = function(x) c(mean=mean(x, na.rm=TRUE), sd=sd(x, na.rm=TRUE)))

summary_metrics_flat <- data.frame(
  top_n_genes    = summary_metrics$top_n_genes,
  ranking_method = summary_metrics$ranking_method,
  gene_direction = summary_metrics$gene_direction,
  model          = summary_metrics$model,
  RMSE_mean      = summary_metrics$RMSE[,"mean"],    RMSE_sd      = summary_metrics$RMSE[,"sd"],
  MAE_mean       = summary_metrics$MAE[,"mean"],     MAE_sd       = summary_metrics$MAE[,"sd"],
  R2_mean        = summary_metrics$R2[,"mean"],      R2_sd        = summary_metrics$R2[,"sd"],
  Pearson_mean   = summary_metrics$Pearson[,"mean"], Pearson_sd   = summary_metrics$Pearson[,"sd"],
  Spearman_mean  = summary_metrics$Spearman[,"mean"],Spearman_sd  = summary_metrics$Spearman[,"sd"],
  stringsAsFactors = FALSE)
summary_metrics_flat <- summary_metrics_flat[order(summary_metrics_flat$RMSE_mean), ]


################################################################################
# Collect this dataset's objects for Section 7
################################################################################

.objs <- setdiff(ls(), BASE_OBJECTS)
PIPE[[ds$name]] <- mget(.objs, envir = globalenv())
message(sprintf("\nCollected %d objects for the figures section.", length(.objs)))

if (isTRUE(SAVE_PIPELINE_OBJECTS)) {
  .rds <- file.path(OUT_DIR, "pipeline_objects.rds")
  saveRDS(PIPE[[ds$name]], .rds)
  message(sprintf("Saved %s (%.0f MB) — Section 7 can run without the pipeline.",
                  basename(.rds), file.size(.rds) / 1024^2))
}

message(sprintf("Pipeline results in: %s", OUT_DIR))
message(sprintf("Log ended: %s", Sys.time()))

"OK"

}, error = function(e) {
  message("\n!!! PIPELINE FAILED for dataset '", ds$name, "': ", conditionMessage(e))
  paste0("FAILED: ", conditionMessage(e))
})

cat(sprintf(" -> %s\n", RUN_STATUS$status[ds_idx]))

rm(list = setdiff(ls(), BASE_OBJECTS))
gc(verbose = FALSE)

}

plan(sequential)

cat("\n================================================================\n")
cat(" PIPELINE SUMMARY\n")
cat("================================================================\n")
for (i in seq_len(nrow(RUN_STATUS)))
  cat(sprintf("  %-20s %-10s %s\n", RUN_STATUS$dataset[i],
              RUN_STATUS$status[i], RUN_STATUS$out_dir[i]))
if (any(grepl("^FAILED", RUN_STATUS$status)))
  cat("\nWARNING: at least one dataset failed — see the messages above.\n")

################################################################################
################################################################################
##                                                                            ##
##   Section 7.  F I G U R E S                                                ##
##                                                                            ##
##   Everything below is plotting. To work down here you still need the SETUP  ##
##   above in memory: Sections 1-5 (libraries, parameters, DATASETS paths,     ##
##   format_p_value, pipeline functions). Section 6 — the pipeline itself —    ##
##   can be skipped once each dataset has a pipeline_objects.rds.              ##
##                                                                            ##
##   In RStudio: select from the top of the file down to just before           ##
##   "# Section 6. PIPELINE RUN" and run that, then come here.                 ##
##                                                                            ##
##   7.0   Load pipeline objects (large .rds from the pipeline)              ##
##   7.2   Domain swap: each dataset's trio predicts the other                ##
##   7.3   Build + save ALL Figure 2 data  ->  figure2_data.rds               ##
##   ---------------------------------------------------------------------   ##
##   7.4   LOAD figure2_data.rds   <<<  START HERE to redraw Figure 2 only    ##
##   7.4b  check that each cached slot holds the expected dataset            ##
##   7.4c  build panel f's data from Supplemental Table 2.1 (PIPE-free)       ##
##   7.5   Figure 2e  PDX ensemble scatter                                   ##
##   7.6   Figure 2d  the same, by cancer type                               ##
##   7.7   Figure 2a  PDX individual models                                  ##
##   7.8   Figure 2b  CCLE individual models                                 ##
##   7.9   Figure 2c  PDX ensembles                                          ##
##   7.10  Figure 2f  gene-level rho, CCLE vs PDX                            ##
##   7.11  assemble Figure 2 (panels a-f)                                     ##
##   ---------------------------------------------------------------------   ##
##   SUPPLEMENTARY FIGURE 2 — domain swap (a, b) and CCLE ensembles (c)     ##
##   7.12  swap panel a   CCLE measured vs PDX-trained ensemble               ##
##   7.13  swap panel b   PDX measured vs CCLE-trained ensemble               ##
##   7.14  assemble the domain-swap supplementary figure                      ##
##   ---------------------------------------------------------------------   ##
##   7.15  supplementary table                                                ##
##                                                                            ##
##   All DATA decisions live in 7.3. All PLOTTING decisions live in the panel  ##
##   blocks 7.5-7.13, one per panel, each reading only FIG2. Edit one panel    ##
##   without touching any other. format_p_value() is the only custom function  ##
##   the plotting code uses.                                                   ##
##                                                                            ##
################################################################################
################################################################################

## ---------------------------------------------------------------------------
## Setup check. Everything below needs these; without them you get a bare
## "object not found" error with no hint about what to do.
## ---------------------------------------------------------------------------
.need <- c("DATASETS", "NEW_FIG_DIR", "format_p_value")
.miss <- .need[!vapply(.need, exists, logical(1))]
if (length(.miss) > 0)
  stop("Section 7 cannot run yet — missing: ", paste(.miss, collapse = ", "), "\n",
       "  These are defined in Sections 1-5. Run the file from the top down to\n",
       "  just before '# Section 6. PIPELINE RUN', then come back to Section 7.\n",
       "  Section 6 itself can be skipped if pipeline_objects.rds already exist.")
rm(.need, .miss)

################################################################################
# 7.0  LOAD PIPELINE OBJECTS
#
# Uses whatever is already in memory from Section 6; otherwise reads the .rds
# files. Run this block first in a fresh session, then jump to any panel.
################################################################################

if (!exists("PIPE")) PIPE <- list()

for (.i in seq_along(DATASETS)) {
  .nm <- DATASETS[[.i]]$name
  if (!is.null(PIPE[[.nm]])) next
  .rds <- file.path(DATASETS[[.i]]$out_dir, "pipeline_objects.rds")
  if (file.exists(.rds)) {
    cat(sprintf("Loading %s from %s\n", .nm, .rds))
    PIPE[[.nm]] <- readRDS(.rds)
  } else {
    cat(sprintf("NOTE: no saved objects for '%s' (%s). Figures needing it will be skipped.\n",
                .nm, .rds))
  }
}
rm(.i, .nm)
cat(sprintf("Datasets available to the figures: %s\n\n",
            paste(names(PIPE), collapse = ", ")))

dir.create(NEW_FIG_DIR, recursive = TRUE, showWarnings = FALSE)



################################################################################
# 7.2  DOMAIN SWAP — train the trio on one dataset, predict the other
#
#      CCLE-trained trio  ->  predicts PDX
#      PDX-trained trio   ->  predicts CCLE
#
# The CV pipeline never keeps a fitted model (fit_models_for_gene_set returns
# predictions, not model objects), so the three members are REFIT here on the
# full source dataset using the same hyperparameters the pipeline uses, then
# applied to the target.
#
# Results are cached to domain_swap.rds — refitting is slow (tune.svm sweeps a
# cost x gamma grid for every gene-set size in both directions). Set
# SWAP_RECOMPUTE <- TRUE to force a rebuild.
#
# ---------------------------------------------------------------------------
# THREE THINGS THAT MATTER HERE, PLEASE READ
#
# 1. GENE OVERLAP. The two expression matrices come from different files and
#    need not share genes. Genes are ranked on the SOURCE but only among genes
#    present in BOTH. The shared-gene count is printed — if it is small, every
#    number below is limited by that, not by the models.
#
# 2. SCALING is a real methodological choice, and it changes the answer.
#      "per_dataset" (default) z-scores each dataset with ITS OWN mean/sd.
#          Standard for cross-platform transfer: it removes per-gene offset and
#          scale differences between the two processing pipelines.
#      "source" applies the SOURCE mean/sd to the target, which is what
#          prepare_xy() does for validation folds. Correct within one dataset;
#          across platforms it bakes in the source's absolute expression scale.
#    This is not cosmetic — Ridge and SVR both take linear combinations, so the
#    choice moves the Spearman. Try both before you trust either.
#
# 3. DOUBLING TIME UNITS DIFFER between the datasets (CCLE Avg_2DT vs PDX
#    dt_days). A model trained on one predicts on that one's scale, so the
#    target predictions are NOT calibrated in absolute terms. Only the
#    rank-based Spearman below is meaningful. Do not compute RMSE or R^2 on
#    these without rescaling first.
# ---------------------------------------------------------------------------
################################################################################

SWAP_TRIO      <- c("RandomForest", "Ridge", "SVR_radial")
SWAP_SIZES     <- c(10, 20, 50, 100, 200, 500)  # gene-set sizes to sweep
SWAP_GENE_DIR  <- "both"        # direction used when ranking genes on the source
SWAP_SCALING   <- "per_dataset" # "per_dataset" or "source"  — see note 2 above
SWAP_RDS       <- file.path(NEW_FIG_DIR, "domain_swap.rds")
SWAP_RECOMPUTE <- FALSE

swap_pairs <- list(
  list(source = "CCLE_norm",      target = "PDX_full_sig_3", label = "CCLE -> PDX"),
  list(source = "PDX_full_sig_3", target = "CCLE_norm",      label = "PDX -> CCLE"))

if (!SWAP_RECOMPUTE && file.exists(SWAP_RDS)) {

  SWAP <- readRDS(SWAP_RDS)
  cat(sprintf("\n7.2  Domain swap loaded from %s (set SWAP_RECOMPUTE <- TRUE to rebuild)\n",
              basename(SWAP_RDS)))

} else {

cat("\n7.2  Domain swap — refitting the trio on each dataset. This is slow.\n")
SWAP <- list()

for (sw in swap_pairs) {

SWAP[[sw$label]] <- tryCatch({

  src <- PIPE[[sw$source]]
  tgt <- PIPE[[sw$target]]
  if (is.null(src)) stop("no pipeline objects for source '", sw$source, "'")
  if (is.null(tgt)) stop("no pipeline objects for target '", sw$target, "'")

  cat(sprintf("\n     %s\n", sw$label))

  ex_s <- src$expr_mat; y_s <- as.numeric(src$dt_vec)
  ex_t <- tgt$expr_mat; y_t <- as.numeric(tgt$dt_vec)

  ## ---- 1. genes present in BOTH datasets ------------------------------------
  shared <- intersect(colnames(ex_s), colnames(ex_t))
  cat(sprintf("       source genes %d | target genes %d | shared %d\n",
              ncol(ex_s), ncol(ex_t), length(shared)))
  if (length(shared) < 10)
    stop("only ", length(shared), " genes shared between the two datasets")

  ## ---- 2. rank genes on the SOURCE, shared genes only -----------------------
  ## Uses the pipeline's own ranking so gene selection matches the CV runs.
  cor_src <- run_spearman(ex_s[, shared, drop = FALSE], y_s,
                          label = sprintf("%s source ranking", sw$label))

  per_size <- list()
  per_pred <- list()

  for (sz in SWAP_SIZES) {

    genes <- intersect(select_top_genes(cor_src, direction = SWAP_GENE_DIR,
                                        top_n = sz), shared)
    if (length(genes) < 2) { cat(sprintf("       %4d genes: too few, skipped\n", sz)); next }

    xs <- as.matrix(ex_s[, genes, drop = FALSE]); storage.mode(xs) <- "numeric"
    xt <- as.matrix(ex_t[, genes, drop = FALSE]); storage.mode(xt) <- "numeric"

    ## non-finite -> that dataset's own column median
    for (j in seq_len(ncol(xs))) {
      m <- median(xs[, j], na.rm = TRUE); if (!is.finite(m)) m <- 0
      xs[!is.finite(xs[, j]), j] <- m
    }
    for (j in seq_len(ncol(xt))) {
      m <- median(xt[, j], na.rm = TRUE); if (!is.finite(m)) m <- 0
      xt[!is.finite(xt[, j]), j] <- m
    }

    ## drop genes with no variance in EITHER dataset — they carry no signal and
    ## would divide by zero when scaling
    sd_s <- apply(xs, 2, sd); sd_t <- apply(xt, 2, sd)
    keep <- is.finite(sd_s) & sd_s > 0 & is.finite(sd_t) & sd_t > 0
    if (sum(keep) < 2) { cat(sprintf("       %4d genes: <2 usable, skipped\n", sz)); next }
    xs <- xs[, keep, drop = FALSE]; xt <- xt[, keep, drop = FALSE]

    ## ---- 3. scaling ---------------------------------------------------------
    if (identical(SWAP_SCALING, "per_dataset")) {
      xs <- sweep(sweep(xs, 2, colMeans(xs), "-"), 2, apply(xs, 2, sd), "/")
      xt <- sweep(sweep(xt, 2, colMeans(xt), "-"), 2, apply(xt, 2, sd), "/")
    } else {
      mu  <- colMeans(xs); sdv <- apply(xs, 2, sd)
      xs <- sweep(sweep(xs, 2, mu, "-"), 2, sdv, "/")
      xt <- sweep(sweep(xt, 2, mu, "-"), 2, sdv, "/")
    }
    xs[!is.finite(xs)] <- 0; xt[!is.finite(xt)] <- 0

    ## ---- 4. refit the trio on the full source, predict the target -----------
    ## Hyperparameters copied from fit_models_for_gene_set.
    set.seed(GLOBAL_SEED)
    rf_fit <- ranger(x = xs, y = y_s, num.trees = 1000,
                     mtry = max(1, floor(sqrt(ncol(xs)))),
                     min.node.size = 5, seed = GLOBAL_SEED)
    p_rf <- as.numeric(predict(rf_fit, data = xt)$predictions)

    set.seed(GLOBAL_SEED)
    cv_ridge <- cv.glmnet(xs, y_s, alpha = 0, nfolds = 5, standardize = FALSE)
    p_ridge  <- as.numeric(predict(cv_ridge, xt, s = "lambda.min"))

    safe_cn <- make.names(colnames(xs), unique = TRUE)
    xs_svr  <- xs; colnames(xs_svr) <- safe_cn
    xt_svr  <- xt; colnames(xt_svr) <- safe_cn
    set.seed(GLOBAL_SEED)
    tune_radial <- tune.svm(y ~ ., data = data.frame(y = y_s, xs_svr, check.names = FALSE),
                            kernel = "radial", cost = c(0.1, 1, 10),
                            gamma = c(0.01, 0.1, 1 / ncol(xs)), epsilon = 0.1,
                            tunecontrol = tune.control(cross = 3))
    p_svr <- as.numeric(predict(tune_radial$best.model,
                                data.frame(xt_svr, check.names = FALSE)))

    ## Ensemble = unweighted mean of the three members' PREDICTIONS.
    p_ens <- (p_rf + p_ridge + p_svr) / 3

    ## ---- 5. rank agreement with the target's measured DT --------------------
    rho <- vapply(list(RandomForest = p_rf, Ridge = p_ridge,
                       SVR_radial = p_svr, Ensemble = p_ens),
                  function(p) {
                    ok <- is.finite(p) & is.finite(y_t)
                    if (sum(ok) < 3) return(NA_real_)
                    suppressWarnings(cor(p[ok], y_t[ok], method = "spearman"))
                  }, numeric(1))

    cat(sprintf("       %4d genes (%d used): RF %.3f | Ridge %.3f | SVR %.3f | Ensemble %.3f\n",
                sz, ncol(xs), rho[["RandomForest"]], rho[["Ridge"]],
                rho[["SVR_radial"]], rho[["Ensemble"]]))

    per_size[[as.character(sz)]] <- data.frame(
      top_n_genes = sz, genes_used = ncol(xs),
      model = names(rho), Spearman = as.numeric(rho),
      stringsAsFactors = FALSE)

    per_pred[[as.character(sz)]] <- data.frame(
      top_n_genes  = sz,
      sample       = rownames(ex_t),
      observed     = y_t,
      RandomForest = p_rf, Ridge = p_ridge, SVR_radial = p_svr,
      Ensemble     = p_ens,
      stringsAsFactors = FALSE)
  }

  if (length(per_size) == 0) stop("no gene-set size produced a usable model")

  list(source         = sw$source,
       target         = sw$target,
       label          = sw$label,
       scaling        = SWAP_SCALING,
       n_shared_genes = length(shared),
       summary        = do.call(rbind, per_size),
       predictions    = do.call(rbind, per_pred))

}, error = function(e) {
  cat(sprintf("     SWAP %s FAILED: %s\n", sw$label, conditionMessage(e)))
  NULL
})

}

saveRDS(SWAP, SWAP_RDS)
cat(sprintf("\n     Saved %s\n", basename(SWAP_RDS)))

}

## Quick readout of what is available to the next two panels.
for (.k in names(SWAP)) {
  if (is.null(SWAP[[.k]])) { cat(sprintf("     %s : unavailable\n", .k)); next }
  .s <- SWAP[[.k]]$summary
  .e <- .s[.s$model == "Ensemble" & .s$top_n_genes == 100, ]
  cat(sprintf("     %-14s shared genes %5d | ensemble Spearman at 100 genes: %s\n",
              .k, SWAP[[.k]]$n_shared_genes,
              if (nrow(.e) == 1) sprintf("%.3f", .e$Spearman) else "n/a"))
}

## Drop this block's working variables. They use generic names (genes, keep,
## rho, xs, xt ...) that other figure blocks also use, and the fitted models
## are large. Only SWAP survives.
rm(list = intersect(ls(), c(
  ".k", ".s", ".e", "sw", "sz", "src", "tgt", "ex_s", "ex_t", "y_s", "y_t",
  "shared", "cor_src", "genes", "xs", "xt", "sd_s", "sd_t", "keep", "mu", "sdv",
  "j", "m", "rf_fit", "cv_ridge", "tune_radial", "safe_cn", "xs_svr", "xt_svr",
  "p_rf", "p_ridge", "p_svr", "p_ens", "rho", "per_size", "per_pred")))
invisible(gc(verbose = FALSE))

################################################################################
# 7.3  BUILD AND SAVE THE FIGURE 2 DATA
#
# Extracts exactly what each panel needs and writes it to ONE small file:
#     <NEW_FIG_DIR>/figure2_data.rds
#
# Everything expensive or bulky happens here — the big pipeline objects, the
# per-fold ensemble averaging, the swap results. What comes out is a few
# hundred KB.
#
# EVERY DATA DECISION LIVES IN THIS SECTION: which model, how many genes, which
# ranking, the minimum samples per cancer type. Change one here, re-run this
# block, then re-run the plots. The panels below hold only plotting choices —
# colours, sizes, themes.
#
# Requires: 7.0 (PIPE) and 7.2 (SWAP).
################################################################################

FIG2_RDS <- file.path(NEW_FIG_DIR, "figure2_data.rds")

cat("\n7.3  Building Figure 2 data\n")

FIG2 <- list()

## ---- shared data choices ---------------------------------------------------
f2_model     <- "SVR_radial"   # model shown in the scatter panels A and B
f2_top_n     <- 100            # gene-set size for the scatter panels
f2_ranking   <- "cor_all"      # ranking method used throughout
f2_direction <- "both"         # gene direction used throughout
f2_min_n     <- 10             # minimum samples for a cancer type to get a sub-panel
f2_swap_n    <- 100            # gene-set size, domain-swap supplementary figure

f2_combos <- list(
  "RandomForest + Ridge"              = c("RandomForest", "Ridge"),
  "RandomForest + SVR_radial"         = c("RandomForest", "SVR_radial"),
  "Ridge + SVR_radial"                = c("Ridge", "SVR_radial"),
  "RandomForest + Ridge + SVR_radial" = c("RandomForest", "Ridge", "SVR_radial"))

## DATA SLOTS. The FIG2 slot letters are internal keys, not panel letters.
## They are printed as:
##   $C -> Figure 2a   $D -> 2b   $E -> 2c   $B -> 2d   $A -> 2e   $F -> 2f

## ---- SLOT A (printed as Figure 2e) : PDX measured vs predicted ----------
f2_panelA_ds <- "PDX_full_sig_3"     # dataset shown in Figure 2e (and 2d)
FIG2$A <- tryCatch({
  pr <- PIPE[[f2_panelA_ds]]$pred_df
  if (is.null(pr)) stop("no pred_df for ", f2_panelA_ds, " — run 7.0")
  pr <- pr[pr$model == f2_model &
           as.numeric(pr$top_n_genes) == f2_top_n &
           pr$ranking_method == f2_ranking &
           pr$gene_direction == f2_direction, ]
  if (nrow(pr) == 0) stop("no predictions matched the filter")
  ## each sample is validated in exactly one fold, so pooling gives one point
  pr <- pr[!duplicated(pr$cell_line), ]
  d  <- data.frame(x = as.numeric(pr$observed_DT),
                   y = as.numeric(pr$predicted_DT))
  d  <- d[is.finite(d$x) & is.finite(d$y), ]
  attr(d, "dataset") <- f2_panelA_ds
  cat(sprintf("     A: %d samples from %s (%s, %d genes)\n",
              nrow(d), f2_panelA_ds, f2_model, f2_top_n))
  d
}, error = function(e) { cat("     A FAILED:", conditionMessage(e), "\n"); NULL })

## ---- SLOT B (printed as Figure 2d) : slot A's dataset, by cancer type ----
## Cancer type comes from lineage_vec, which the loader aligns to
## rownames(expr_mat). Note the two datasets label it at different granularity:
##   CCLE -> OncotreeLineage from Models_DT.csv (broad: Lung, Breast, ...)
##   PDX  -> OncoTreeCode from the phenotype CSV (specific: LUAD, BRCA, ...)
## Codes are finer than lineages, so the same f2_min_n cut leaves fewer, smaller
## groups on PDX. The dropped list printed below tells you what was lost.
f2_panelB_ds <- f2_panelA_ds         # Figure 2d follows Figure 2e's dataset
FIG2$B <- tryCatch({
  pr <- PIPE[[f2_panelB_ds]]$pred_df
  if (is.null(pr)) stop("no pred_df for ", f2_panelB_ds, " — run 7.0")
  pr <- pr[pr$model == f2_model &
           as.numeric(pr$top_n_genes) == f2_top_n &
           pr$ranking_method == f2_ranking &
           pr$gene_direction == f2_direction, ]
  if (nrow(pr) == 0) stop("no predictions matched the filter")
  pr <- pr[!duplicated(pr$cell_line), ]

  lin_map <- setNames(as.character(PIPE[[f2_panelB_ds]]$lineage_vec),
                      rownames(PIPE[[f2_panelB_ds]]$expr_mat))
  d <- data.frame(x = as.numeric(pr$observed_DT),
                  y = as.numeric(pr$predicted_DT),
                  lineage = unname(lin_map[pr$cell_line]),
                  stringsAsFactors = FALSE)
  d <- d[is.finite(d$x) & is.finite(d$y) & !is.na(d$lineage) & nzchar(d$lineage), ]

  cnt  <- table(d$lineage)
  keep <- names(cnt)[cnt >= f2_min_n]
  if (length(keep) == 0)
    stop("no cancer type reaches f2_min_n = ", f2_min_n, " (largest is ", max(cnt), ")")
  if (length(setdiff(names(cnt), keep)) > 0)
    cat(sprintf("     B: dropped (n < %d): %s\n", f2_min_n,
                paste(sprintf("%s (n=%d)", setdiff(names(cnt), keep),
                              cnt[setdiff(names(cnt), keep)]), collapse = ", ")))
  d   <- d[d$lineage %in% keep, ]
  ord <- names(sort(cnt[keep], decreasing = TRUE))     # biggest type first
  d$lineage <- factor(d$lineage, levels = ord,
                      labels = sprintf("%s (n=%d)", ord, cnt[ord]))
  attr(d, "dataset") <- f2_panelB_ds
  cat(sprintf("     B: %d samples from %s across %d cancer types\n",
              nrow(d), f2_panelB_ds, nlevels(d$lineage)))
  d
}, error = function(e) { cat("     B FAILED:", conditionMessage(e), "\n"); NULL })

## ---- SLOTS C and D (printed as Figures 2a, 2b) : individual models -------
## PDX FIRST, then CCLE. The order is set HERE, at the data layer, so the letter
## and the dataset stay bound together everywhere downstream: panel block, the
## standalone panel PDF, and the supplementary table's cross-check column.
for (.p in list(list(id = "C", ds = "PDX_full_sig_3"),
                list(id = "D", ds = "CCLE_norm"))) {
  FIG2[[.p$id]] <- tryCatch({
    sm <- PIPE[[.p$ds]]$summary_metrics_flat
    if (is.null(sm)) stop("no summary_metrics_flat for ", .p$ds, " — run 7.0")
    sm <- sm[sm$ranking_method == f2_ranking & sm$gene_direction == f2_direction, ]
    sm$top_n_genes <- as.numeric(sm$top_n_genes)
    sm <- sm[is.finite(sm$Spearman_mean),
             c("top_n_genes", "model", "Spearman_mean", "Spearman_sd")]
    if (nrow(sm) == 0) stop("no rows after filtering")
    ## Fold counts for the SEM error bars. summary_metrics_flat was made with
    ## aggregate(cbind(RMSE, MAE, R2, Pearson, Spearman) ~ ...), which drops a
    ## fold if ANY of those five is NA, so count folds the same way.
    vm <- PIPE[[.p$ds]]$valid_metrics
    if (!is.null(vm)) {
      vm <- vm[vm$ranking_method == f2_ranking & vm$gene_direction == f2_direction &
               !is.na(vm$top_n_genes), ]
      vm <- vm[complete.cases(vm[, c("RMSE", "MAE", "R2", "Pearson", "Spearman")]), ]
      if (nrow(vm) > 0) {
        cnt <- aggregate(list(Spearman_n = vm$Spearman),
                         by = list(top_n_genes = as.numeric(vm$top_n_genes),
                                   model       = as.character(vm$model)),
                         FUN = length)
        sm$Spearman_n <- cnt$Spearman_n[match(
          paste(sm$model, sm$top_n_genes), paste(cnt$model, cnt$top_n_genes))]
      }
    }
    attr(sm, "dataset") <- .p$ds
    cat(sprintf("     %s: %d rows, %d models (%s)\n",
                .p$id, nrow(sm), length(unique(sm$model)), .p$ds))
    sm
  }, error = function(e) { cat("     ", .p$id, "FAILED:", conditionMessage(e), "\n"); NULL })
}

## ---- SLOT E (printed as Figure 2c) and the CCLE ensembles (Sup Fig 2c) ---
## The ensemble prediction is the unweighted mean of its members' PREDICTED
## values, computed per fold, then correlated against the observed doubling
## time. This is NOT the mean of the members' individual Spearman scores —
## averaging scores describes the typical member, not the ensemble.
##
## Slot E is the PDX ensemble curve (Figure 2c, Table 3.3); slot E_CCLE is the
## CCLE one (Supplementary Figure 2c, Table 3.4). Slot F is NOT filled here:
## panel 2f is the gene-level rho scatter, built in 7.4c from Supplemental
## Table 2.1.
for (.p in list(list(id = "E",      ds = "PDX_full_sig_3"),
                list(id = "E_CCLE", ds = "CCLE_norm"))) {
  FIG2[[.p$id]] <- tryCatch({
    pr <- PIPE[[.p$ds]]$pred_df
    if (is.null(pr)) stop("no pred_df for ", .p$ds, " — run 7.0")
    pr <- pr[pr$ranking_method == f2_ranking & pr$gene_direction == f2_direction, ]
    if (nrow(pr) == 0) stop("no predictions after filtering")

    sizes <- sort(unique(as.numeric(pr$top_n_genes)))
    folds <- sort(unique(pr$fold))

    out <- do.call(rbind, lapply(names(f2_combos), function(lab) {
      members <- f2_combos[[lab]]
      do.call(rbind, lapply(sizes, function(sz) {
        rho <- vapply(folds, function(f) {
          q <- pr[as.numeric(pr$top_n_genes) == sz & pr$fold == f, ]
          if (nrow(q) == 0) return(NA_real_)
          vecs <- lapply(members, function(m) {
            r <- q[q$model == m, ]
            r <- r[!duplicated(r$cell_line), ]
            setNames(r$predicted_DT, r$cell_line)
          })
          if (any(vapply(vecs, length, integer(1)) == 0)) return(NA_real_)
          common <- Reduce(intersect, lapply(vecs, names))
          if (length(common) < 3) return(NA_real_)
          ens <- rowMeans(do.call(cbind, lapply(vecs, function(v) v[common])))
          o   <- q[!duplicated(q$cell_line), ]
          obs <- setNames(o$observed_DT, o$cell_line)[common]
          ok  <- is.finite(ens) & is.finite(obs)
          if (sum(ok) < 3) return(NA_real_)
          suppressWarnings(cor(ens[ok], obs[ok], method = "spearman"))
        }, numeric(1))
        data.frame(top_n_genes = sz, ensemble = lab, n_members = length(members),
                   Spearman_mean = mean(rho, na.rm = TRUE),
                   Spearman_sd   = sd(rho,   na.rm = TRUE),
                   Spearman_n    = sum(is.finite(rho)),   # folds, for the SEM bars
                   stringsAsFactors = FALSE)
      }))
    }))
    out <- out[is.finite(out$Spearman_mean), ]
    if (nrow(out) == 0) stop("no ensemble rows — are all trio members present?")
    out$ensemble <- factor(out$ensemble, levels = names(f2_combos))
    attr(out, "dataset") <- .p$ds
    cat(sprintf("     %s: %d combinations x %d sizes (%s)\n",
                .p$id, nlevels(droplevels(out$ensemble)), length(sizes), .p$ds))
    out
  }, error = function(e) { cat("     ", .p$id, "FAILED:", conditionMessage(e), "\n"); NULL })
}

## ---- SWAP PANELS G and H : domain swap -------------------------------------
## G: CCLE measured vs the PDX-trained trio.  H: PDX measured vs CCLE-trained.
## These are Supplementary Figure 2 panels a and b, drawn in 7.12-7.14.
for (.p in list(list(id = "G", key = "PDX -> CCLE", target = "CCLE"),
                list(id = "H", key = "CCLE -> PDX", target = "PDX"))) {
  FIG2[[.p$id]] <- tryCatch({
    if (!exists("SWAP")) stop("SWAP not found — run 7.2 first")
    sw <- SWAP[[.p$key]]
    if (is.null(sw)) stop("the ", .p$key, " swap is unavailable")
    pg <- sw$predictions[sw$predictions$top_n_genes == f2_swap_n, ]
    if (nrow(pg) == 0) stop("no swap predictions at ", f2_swap_n, " genes")
    used <- unique(sw$summary$genes_used[sw$summary$top_n_genes == f2_swap_n])
    d <- data.frame(x = as.numeric(pg$observed), y = as.numeric(pg$Ensemble))
    d <- d[is.finite(d$x) & is.finite(d$y), ]
    attr(d, "genes_used") <- used
    attr(d, "dataset")    <- .p$key
    cat(sprintf("     %s: %d %s samples, %d genes used (of %d requested)\n",
                .p$id, nrow(d), .p$target, used, f2_swap_n))
    d
  }, error = function(e) { cat("     ", .p$id, "FAILED:", conditionMessage(e), "\n"); NULL })
}
rm(.p)

## ---- save ------------------------------------------------------------------
## Stamp the file with what produced it, so 7.4 can tell you what you are
## actually plotting. Stored as an attribute so it does not disturb the
## per-panel loops above.
attr(FIG2, "meta") <- list(
  built     = Sys.time(),
  panelA_ds = f2_panelA_ds,
  panelB_ds = f2_panelB_ds,
  model     = f2_model,
  top_n     = f2_top_n,
  ranking   = f2_ranking,
  direction = f2_direction,
  min_n     = f2_min_n,
  swap_n    = f2_swap_n)

saveRDS(FIG2, FIG2_RDS)
cat(sprintf("\n     Saved %s (%.0f KB) — panels below read only this.\n",
            basename(FIG2_RDS), file.size(FIG2_RDS) / 1024))
## What 7.3 is responsible for. F is absent on purpose — 7.4c fills it.
FIG2_SLOTS_73 <- c("A", "B", "C", "D", "E", "E_CCLE", "G", "H")
.built <- names(FIG2)[!vapply(FIG2, is.null, logical(1))]
cat(sprintf("     Slots with data: %s\n", paste(.built, collapse = ", ")))
if (!all(FIG2_SLOTS_73 %in% .built))
  cat(sprintf("     MISSING: %s\n",
              paste(setdiff(FIG2_SLOTS_73, .built), collapse = ", ")))
cat("     (F is filled by 7.4c, not here.)\n")
rm(.built)


################################################################################
################################################################################
##                                                                            ##
##   7.4  LOAD THE FIGURE 2 DATA                                              ##
##                                                                            ##
##   >>> START HERE to redraw Figure 2 without re-running anything else. <<<   ##
##                                                                            ##
##   IMPORTANT: the panels below read figure2_data.rds and NOTHING else. If    ##
##   you change WHICH DATA a panel shows (the f2_* settings in 7.3 — dataset,  ##
##   model, gene count, ranking), you must RE-RUN 7.3 to rebuild this file.    ##
##   Re-running only the panels replays whatever was cached last time.         ##
##                                                                            ##
##   The summary printed below tells you what is actually in the file.         ##
##                                                                            ##
##   From a fresh R session, run Sections 1-4 (libraries, parameters, paths,   ##
##   format_p_value), then this block, then any panel below. Nothing in        ##
##   Sections 5-7.3 is needed — no pipeline, no big .rds, no refitting.        ##
##                                                                            ##
##   Everything below this line reads FIG2 and nothing else.                   ##
##                                                                            ##
################################################################################
################################################################################

.need <- c("NEW_FIG_DIR", "format_p_value")
.miss <- .need[!vapply(.need, exists, logical(1))]
if (length(.miss) > 0)
  stop("Cannot start at 7.4 — missing: ", paste(.miss, collapse = ", "), "\n",
       "  Run Sections 1-5 first (top of the file down to just before\n",
       "  '# Section 6. PIPELINE RUN'), then jump back here.")
rm(.need, .miss)

FIG2_RDS <- file.path(NEW_FIG_DIR, "figure2_data.rds")

if (!file.exists(FIG2_RDS))
  stop("figure2_data.rds not found at ", FIG2_RDS,
       "\n  Run section 7.3 once to build it.")

FIG2 <- readRDS(FIG2_RDS)

cat(sprintf("\n7.4  Loaded %s\n", basename(FIG2_RDS)))

.meta <- attr(FIG2, "meta")
if (is.null(.meta)) {
  cat("     WARNING: this file has no build stamp, so it was not written by\n")
  cat("              7.3 of this script. Re-run 7.3 to rebuild it.\n")
} else {
  cat(sprintf("     built %s\n", format(.meta$built, "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("     panel a/b dataset: %s / %s | %s at %d genes | %s, %s\n",
              .meta$panelA_ds, .meta$panelB_ds, .meta$model, .meta$top_n,
              .meta$ranking, .meta$direction))
}
rm(.meta)
for (.k in names(FIG2)) {
  if (is.null(FIG2[[.k]])) { cat(sprintf("     %s : (empty)\n", .k)); next }
  cat(sprintf("     %s : %-28s %d rows\n", .k,
              if (is.null(attr(FIG2[[.k]], "dataset"))) "?" else attr(FIG2[[.k]], "dataset"),
              nrow(FIG2[[.k]])))
}
rm(.k)

################################################################################
# 7.4b  CHECK THE CACHED SLOTS
#
# 7.3 stamps every slot with the dataset it holds. The panels need C = PDX
# (Figure 2a), D = CCLE (Figure 2b) and E = the PDX ensembles (Figure 2c). A
# file that does not match was not written by 7.3 of this script — re-run 7.3.
################################################################################

FIG2_EXPECTED <- c(C = "PDX_full_sig_3", D = "CCLE_norm", E = "PDX_full_sig_3")

cat("\n7.4b Slot check\n")
for (.k in names(FIG2_EXPECTED)) {
  if (is.null(FIG2[[.k]])) next
  .ds <- attr(FIG2[[.k]], "dataset")
  if (!identical(.ds, unname(FIG2_EXPECTED[[.k]])))
    stop("FIG2$", .k, " holds '", if (is.null(.ds)) "(no stamp)" else .ds,
         "', expected '", FIG2_EXPECTED[[.k]], "'.\n",
         "  Re-run 7.3 to rebuild figure2_data.rds.")
}
cat("     OK — C = PDX, D = CCLE, E = PDX ensembles.\n")
rm(list = intersect(ls(), c(".k", ".ds")))

dir.create(NEW_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

################################################################################
# 7.4c  PANEL f DATA — gene-level Spearman rho, CCLE vs PDX
#
# WHY THIS IS NOT IN 7.3
# Every other panel's data comes out of the modelling pipeline, so 7.3 needs
# PIPE in memory. This one does not: it reads Table 2.1 of sup_tab2.xlsx,
# written by plot_figure_1.R, and nothing else. Keeping it here, in the load
# stage, means panel f can be rebuilt from a cold session in seconds and that
# an existing figure2_data.rds never has to be regenerated to get it.
#
# WHAT IT COMPUTES
#   1. read the pooled pan-cancer rho for each dataset, drop genes with no rho
#   2. keep the genes present in both
#   3. a gene is SAME-DIRECTION if its rho has the same sign in both datasets
#   4. rank the same-direction genes by |rho|, separately in each dataset, then
#      walk the depth n upward until the two top-n lists share at least
#      GENE_CORR$min_overlap genes. That intersection is the highlighted set.
#
# The sweep stops at the first depth that clears the threshold, so it is quick
# despite looking like a scan over every gene.
#
# Result -> FIG2$F, one row per gene common to both datasets:
#   gene, rho_ccle, rho_pdx, same_direction, top_set
# cached as GENE_CORR_RDS so repeated runs do not re-read the workbook. The
# cache is rebuilt automatically whenever sup_tab2.xlsx changes.
################################################################################

GENE_CORR_RDS <- file.path(NEW_FIG_DIR, "panel_f_gene_corr.rds")

## Set TRUE to ignore the cache and re-read the workbook.
if (!exists("GENE_CORR_REBUILD")) GENE_CORR_REBUILD <- FALSE

cat("\n7.4c Panel f data — gene-level rho, CCLE vs PDX\n")

FIG2$F <- tryCatch({
  
  ## The cache is used only if it was built from GENE_CORR$xlsx and that
  ## workbook has not been rewritten since.
  .gc_cache_ok <- !GENE_CORR_REBUILD && file.exists(GENE_CORR_RDS) &&
    identical(attr(readRDS(GENE_CORR_RDS), "source"),
              normalizePath(GENE_CORR$xlsx, mustWork = FALSE)) &&
    (!file.exists(GENE_CORR$xlsx) ||
       file.mtime(GENE_CORR$xlsx) <= file.mtime(GENE_CORR_RDS))
  if (.gc_cache_ok) {
    .gc <- readRDS(GENE_CORR_RDS)
    cat(sprintf("     Loaded %s — %d genes, %d highlighted (depth %d)\n",
                basename(GENE_CORR_RDS), nrow(.gc),
                sum(.gc$top_set), attr(.gc, "rank_depth")))
    .gc
    
  } else {
    if (!READXL_OK) stop("package 'readxl' is needed to read Supplemental Table 2")
    if (!file.exists(GENE_CORR$xlsx))
      stop("no Supplemental Table 2 at ", GENE_CORR$xlsx,
           " — run plot_figure_1.R first")
    .t21 <- as.data.frame(readxl::read_excel(GENE_CORR$xlsx,
                                             sheet = GENE_CORR$sheet))
    
    ## ---- one rho column of Table 2.1 -> named vector of rho ----------------
    ## Table 2.1 is merged on gene across analyses, so a gene tested in only
    ## one dataset has NA in the other's column; those are dropped here.
    .read_rho <- function(d, rho_col, what) {
      if (!all(c(GENE_CORR$gene_col, rho_col) %in% names(d)))
        stop("sheet '", GENE_CORR$sheet, "' has no '",
             GENE_CORR$gene_col, "' / '", rho_col, "' column")
      g <- as.character(d[[GENE_CORR$gene_col]])
      r <- as.numeric(d[[rho_col]])
      ok <- !is.na(g) & !duplicated(g) & is.finite(r)
      cat(sprintf("     %-4s %d genes with a usable rho (of %d rows)\n",
                  what, sum(ok), nrow(d)))
      setNames(r[ok], g[ok])
    }
    
    .rc <- .read_rho(.t21, GENE_CORR$ccle_rho_col, "CCLE")
    .rp <- .read_rho(.t21, GENE_CORR$pdx_rho_col,  "PDX")
    rm(.t21)
    
    ## ---- genes in both ------------------------------------------------------
    .common <- intersect(names(.rc), names(.rp))
    if (length(.common) == 0) stop("no genes shared by the two tables")
    .rc <- .rc[.common]; .rp <- .rp[.common]
    cat(sprintf("     %d genes in both tables\n", length(.common)))
    
    ## ---- same direction -----------------------------------------------------
    .same <- (.rc > 0 & .rp > 0) | (.rc < 0 & .rp < 0)
    .sd   <- .common[.same]
    cat(sprintf("     %d of them agree in sign\n", length(.sd)))
    
    ## ---- top-|rho| depth at which the two lists overlap enough --------------
    ## Ordered strongest-|rho| first, within the same-direction genes only.
    .ord_c <- .sd[order(abs(.rc[.sd]), decreasing = TRUE)]
    .ord_p <- .sd[order(abs(.rp[.sd]), decreasing = TRUE)]
    .depth_max <- min(GENE_CORR$max_depth, length(.sd))
    .top <- NULL; .depth <- NA_integer_
    for (.n in seq.int(10L, .depth_max)) {
      .hit <- intersect(.ord_c[seq_len(.n)], .ord_p[seq_len(.n)])
      if (length(.hit) >= GENE_CORR$min_overlap) {
        .top <- .hit; .depth <- .n; break
      }
    }
    if (is.null(.top))
      stop("the two |rho| rankings never share ", GENE_CORR$min_overlap,
           " genes within the first ", .depth_max,
           " — lower GENE_CORR$min_overlap or raise max_depth")
    cat(sprintf("     Highlighted set: %d genes, reached at depth %d\n",
                length(.top), .depth))
    
    ## ---- the table the panel plots -----------------------------------------
    .gc <- data.frame(
      gene           = .common,
      rho_ccle       = unname(.rc),
      rho_pdx        = unname(.rp),
      same_direction = .same,
      top_set        = .common %in% .top,
      stringsAsFactors = FALSE)
    .gc <- .gc[order(-abs(.gc$rho_ccle)), ]
    rownames(.gc) <- NULL
    
    attr(.gc, "dataset")     <- "CCLE_norm vs PDX_full_sig_3"
    attr(.gc, "rank_depth")  <- .depth
    attr(.gc, "min_overlap") <- GENE_CORR$min_overlap
    attr(.gc, "source")      <- normalizePath(GENE_CORR$xlsx, mustWork = FALSE)
    attr(.gc, "built")       <- Sys.time()
    
    saveRDS(.gc, GENE_CORR_RDS)
    cat(sprintf("     Saved %s (%.0f KB)\n", basename(GENE_CORR_RDS),
                file.size(GENE_CORR_RDS) / 1024))
    .gc
  }
  
}, error = function(e) {
  cat("     PANEL f DATA FAILED:", conditionMessage(e), "\n")
  NULL
})
rm(list = intersect(ls(), c(".gc", ".rc", ".rp", ".common", ".same", ".sd",
                            ".ord_c", ".ord_p", ".depth_max", ".top", ".depth",
                            ".n", ".hit", ".read_rho")))

################################################################################
# 7.4d  ENSEMBLE SCATTER — rebuild slots A and B (Figures 2e and 2d) from the
#       full trio
#
# Figure 2e and 2d show the full trio ensemble, RandomForest + Ridge +
# SVR_radial — the same ensemble whose curve is Figure 2c, so the scatter and
# the curve describe the same predictor. (7.3 fills slots A and B with a
# single-model scatter; this block replaces them with the ensemble.)
#
# HOW THE ENSEMBLE PREDICTION IS FORMED
# Each sample is held out in exactly one fold, and each of the three members
# predicts it in that fold. The ensemble value is the unweighted MEAN of those
# three predictions for that sample. This is the same definition 7.3 uses for
# the ensemble curves, so Figures 2e and 2c are consistent: it is NOT the
# mean of the three members' Spearman scores, which would describe the typical
# member rather than the ensemble.
#
# A sample is used only if all three members predicted it. Any sample missing a
# member is dropped, and the count is reported.
#
# WHY HERE AND NOT IN 7.3
# Overriding at load time means an existing figure2_data.rds does not have to
# be rebuilt: FIG2$A and FIG2$B are simply replaced after it is read. Every
# other panel, and every supplementary table, is untouched and still comes from
# the cached file. The result is cached in its own .rds, so the pipeline
# objects are read once and later runs are fast.
#
# WHAT IT NEEDS
# pipeline_objects.rds for the panel a/b dataset — the per-sample predictions
# are not in figure2_data.rds, which stores only the finished x/y pairs for the
# single model. If it is missing, this block says so and leaves FIG2$A and $B
# as they are, so the figure still builds with the SVR_radial scatter and the
# run does not fail.
################################################################################

ENS_AB_RDS <- file.path(NEW_FIG_DIR, "panel_ab_ensemble.rds")

## The members averaged for Figures 2e and 2d (slots A and B). Must match one
## of f2_combos in 7.3 for Figures 2e and 2c to describe the same predictor.
ENS_AB_MEMBERS <- c("RandomForest", "Ridge", "SVR_radial")

## Set TRUE to ignore the cache and re-read the pipeline objects.
if (!exists("ENS_AB_REBUILD")) ENS_AB_REBUILD <- FALSE

cat("\n7.4d Ensemble scatter for panels a and b\n")

.ens_ab <- tryCatch({

  if (!ENS_AB_REBUILD && file.exists(ENS_AB_RDS)) {
    .r <- readRDS(ENS_AB_RDS)
    cat(sprintf("     Loaded %s — %d samples, %d cancer types\n",
                basename(ENS_AB_RDS), nrow(.r$A), nlevels(.r$B$lineage)))
    .r

  } else {
    ## ---- settings, taken from the cached build stamp where possible --------
    .m    <- attr(FIG2, "meta")
    .ds   <- if (!is.null(.m$panelA_ds)) .m$panelA_ds else "PDX_full_sig_3"
    .topn <- if (!is.null(.m$top_n))     as.numeric(.m$top_n) else 100
    .rank <- if (!is.null(.m$ranking))   .m$ranking   else "cor_all"
    .dir  <- if (!is.null(.m$direction)) .m$direction else "both"
    .minn <- if (!is.null(.m$min_n))     as.numeric(.m$min_n) else 10
    cat(sprintf("     %s, %d genes, %s / %s, ensemble = %s\n",
                .ds, .topn, .rank, .dir,
                paste(ENS_AB_MEMBERS, collapse = " + ")))

    ## ---- the pipeline objects ----------------------------------------------
    .P <- if (exists("PIPE") && !is.null(PIPE[[.ds]])) PIPE[[.ds]] else {
      .hit <- vapply(DATASETS, function(d) identical(d$name, .ds), logical(1))
      if (!any(.hit)) stop(.ds, " is not in DATASETS")
      .f <- file.path(DATASETS[[which(.hit)[1]]]$out_dir, "pipeline_objects.rds")
      if (!file.exists(.f))
        stop("no pipeline_objects.rds at ", .f,
             "\n  The per-sample predictions are not in figure2_data.rds, so",
             "\n  this block cannot run without it. Point out_dir at the folder",
             "\n  the pipeline actually wrote to.")
      cat(sprintf("     reading %s\n", .f))
      readRDS(.f)
    }
    if (is.null(.P$pred_df)) stop("no pred_df for ", .ds)

    ## ---- the three members' predictions for this gene set -------------------
    .pr <- .P$pred_df
    .pr <- .pr[.pr$model %in% ENS_AB_MEMBERS &
               as.numeric(.pr$top_n_genes) == .topn &
               .pr$ranking_method == .rank &
               .pr$gene_direction == .dir, ]
    if (nrow(.pr) == 0)
      stop("no predictions matched: check the model names, ", .topn,
           " genes, ", .rank, " / ", .dir)

    ## one row per sample per member — each sample is validated in one fold
    .vecs <- lapply(ENS_AB_MEMBERS, function(m) {
      .q <- .pr[.pr$model == m, ]
      .q <- .q[!duplicated(.q$cell_line), ]
      setNames(as.numeric(.q$predicted_DT), as.character(.q$cell_line))
    })
    names(.vecs) <- ENS_AB_MEMBERS
    .n_each <- vapply(.vecs, length, integer(1))
    if (any(.n_each == 0))
      stop("no predictions for: ",
           paste(ENS_AB_MEMBERS[.n_each == 0], collapse = ", "))

    ## only samples all three members predicted
    .common <- Reduce(intersect, lapply(.vecs, names))
    if (length(.common) == 0) stop("the three members share no samples")
    cat(sprintf("     per member: %s -> %d samples common to all three\n",
                paste(sprintf("%s %d", names(.n_each), .n_each), collapse = ", "),
                length(.common)))
    if (length(.common) < max(.n_each))
      cat(sprintf("     dropped %d sample(s) missing at least one member\n",
                  max(.n_each) - length(.common)))

    ## THE ENSEMBLE: unweighted mean of the members' predictions, per sample
    .ens <- rowMeans(do.call(cbind, lapply(.vecs, function(v) v[.common])))

    ## measured DT is the same whichever member reported it
    .obs <- .pr[!duplicated(.pr$cell_line), ]
    .obs <- setNames(as.numeric(.obs$observed_DT), as.character(.obs$cell_line))

    ## ---- slot A (Figure 2e) ------------------------------------------------
    .A <- data.frame(x = unname(.obs[.common]), y = unname(.ens))
    .A <- .A[is.finite(.A$x) & is.finite(.A$y), ]
    attr(.A, "dataset") <- .ds
    cat(sprintf("     A: %d samples\n", nrow(.A)))

    ## ---- slot B (Figure 2d) : the same points, split by cancer type -------
    ## Same lineage source and same f2_min_n cut as 7.3.
    .lin <- setNames(as.character(.P$lineage_vec), rownames(.P$expr_mat))
    .B <- data.frame(x = unname(.obs[.common]), y = unname(.ens),
                     lineage = unname(.lin[.common]), stringsAsFactors = FALSE)
    .B <- .B[is.finite(.B$x) & is.finite(.B$y) &
             !is.na(.B$lineage) & nzchar(.B$lineage), ]
    .cnt  <- table(.B$lineage)
    .keep <- names(.cnt)[.cnt >= .minn]
    if (length(.keep) == 0)
      stop("no cancer type reaches min_n = ", .minn)
    if (length(setdiff(names(.cnt), .keep)) > 0)
      cat(sprintf("     B: dropped (n < %d): %s\n", .minn,
                  paste(sprintf("%s (n=%d)", setdiff(names(.cnt), .keep),
                                .cnt[setdiff(names(.cnt), .keep)]),
                        collapse = ", ")))
    .B   <- .B[.B$lineage %in% .keep, ]
    .ord <- names(sort(.cnt[.keep], decreasing = TRUE))
    .B$lineage <- factor(.B$lineage, levels = .ord,
                         labels = sprintf("%s (n=%d)", .ord, .cnt[.ord]))
    attr(.B, "dataset") <- .ds
    cat(sprintf("     B: %d samples across %d cancer types\n",
                nrow(.B), nlevels(.B$lineage)))

    .r <- list(A = .A, B = .B,
               members = ENS_AB_MEMBERS, dataset = .ds, n_genes = .topn,
               built = Sys.time())
    saveRDS(.r, ENS_AB_RDS)
    cat(sprintf("     Saved %s\n", basename(ENS_AB_RDS)))
    .r
  }

}, error = function(e) {
  cat("     ENSEMBLE SCATTER FAILED:", conditionMessage(e), "\n")
  cat("     Panels a and b will fall back to whatever figure2_data.rds holds\n")
  cat("     (the single-model scatter). The rest of the figure is unaffected.\n")
  NULL
})

## ---- the override -----------------------------------------------------------
## Only A and B are replaced. C through H come from the cached file untouched.
if (!is.null(.ens_ab)) {
  FIG2$A <- .ens_ab$A
  FIG2$B <- .ens_ab$B
  cat(sprintf("     FIG2$A and $B now hold the %s ensemble\n",
              paste(.ens_ab$members, collapse = " + ")))
}
rm(list = intersect(ls(), c(".ens_ab", ".r", ".m", ".ds", ".topn", ".rank",
                            ".dir", ".minn", ".P", ".hit", ".f", ".pr", ".vecs",
                            ".n_each", ".common", ".ens", ".obs", ".A", ".B",
                            ".lin", ".cnt", ".keep", ".ord", ".q")))

################################################################################
# ERROR BARS — standard error of the mean across the CV folds
#
# Used by Figures 2a, 2b, 2c and Supplementary Figure 2c. Each point is the
# mean validation rho over the folds; its bar is mean +/- SEM, where
# SEM = SD / sqrt(n) and n is the number of folds that gave a finite rho at
# that point (stored by 7.3 as Spearman_n). If the data has no Spearman_n
# column, SEM_N_FOLDS_DEFAULT is used for n and a message says so.
################################################################################
SEM_N_FOLDS_DEFAULT <- 5

add_sem <- function(df) {
  if ("Spearman_n" %in% names(df) && !all(is.na(df$Spearman_n))) {
    n <- df$Spearman_n
  } else {
    cat(sprintf("     (no fold counts in the data — SEM assumes %d folds; rebuild 7.3 to use the real counts)\n",
                SEM_N_FOLDS_DEFAULT))
    n <- rep(SEM_N_FOLDS_DEFAULT, nrow(df))
  }
  df$Spearman_sem <- ifelse(n > 1, df$Spearman_sd / sqrt(n), NA_real_)
  df
}

## One layer, shared by all four panels so the bars look the same everywhere.
## width is in log10(x) units because the x axis is log10.
sem_errorbar <- function()
  geom_errorbar(aes(ymin = Spearman_mean - Spearman_sem,
                    ymax = Spearman_mean + Spearman_sem),
                width = 0.04, linewidth = 0.5, na.rm = TRUE, show.legend = FALSE)


################################################################################
# 7.5  FIGURE 2e — PDX ENSEMBLE SCATTER (measured vs predicted DT)
#      Measured doubling time (x) vs ENSEMBLE predicted (y), all PDX samples,
#      100 genes. The prediction is the mean of RandomForest, Ridge and
#      SVR_radial (see 7.4d, which fills this data slot after the cache is read).
#      Printed as panel e (bottom left).  Data: FIG2$A
################################################################################

cat("\n7.5  Figure 2e — PDX ensemble scatter\n")

fig2e_PDX_ensemble_scatter <- NULL
tryCatch({

  h_df <- FIG2$A
  if (is.null(h_df)) stop("FIG2$A is empty — rebuild in 7.3")

  
  fig2e_PDX_ensemble_scatter <- ggplot(h_df, aes(x = log10(x), y = log10(y))) +
    geom_point(alpha = 0.4, colour = "#808080", size = 0.9) +
    geom_smooth(method = "lm", colour = "#005AB5", se = FALSE, na.rm = TRUE) +
    stat_cor(aes(label = paste(..r.label..,
                               after_stat(format_p_value(..p..)),
                               sep = "~`,`~")),
             method = "spearman", r.accuracy = 0.001, parse = TRUE, size = 4.5,
             label.y = 1.5) +
    xlab(TeX("$log_{10}(DT)\\ (days)$")) +
    ylab(TeX("$log_{10}(pDT)$")) +
    theme(legend.position  = "none",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line        = element_line(colour = "black"),
          axis.text        = element_text(size = 14,  face = "bold"),
          axis.title       = element_text(size = 16, face = "bold"))
  
  ggsave(file.path(NEW_FIG_DIR, "Figure_2e_PDX_ensemble_scatter.pdf"), fig2e_PDX_ensemble_scatter,
         width = 6, height = 6, device = cairo_pdf)
  cat("     Saved Figure_2e_PDX_ensemble_scatter.pdf\n")
  
}, error = function(e) cat("     FIGURE 2e FAILED:", conditionMessage(e), "\n"))



################################################################################
# 7.6  FIGURE 2d — PDX ENSEMBLE SCATTER BY CANCER TYPE
#      The same ensemble scatter as Figure 2e, one sub-panel per cancer type
#      (four types are chosen below). Uses Figure 2e's dataset.
#      Printed as panel d (middle right).  Data: FIG2$B
################################################################################

cat("\n7.6  Figure 2d — PDX ensemble scatter by cancer type\n")

fig2d_PDX_by_cancer_type <- NULL
tryCatch({
  
  h_df <- FIG2$B
  if (is.null(h_df)) stop("FIG2$B is empty — rebuild in 7.3")
  n_facets <- 4   # cancer types shown (lineage_plot below)
  h_df$lineage <- as.character(h_df$lineage)
  ## Shorten one facet title. The lineage labels carry their sample counts, so
  ## this line and lineage_plot below must match the counts in the data.
  h_df[which(h_df$lineage == "Head and Neck (n=311)"),"lineage"] <-"Head&Neck (n=311)"
  
  
  ## The four cancer types shown
  lineage_plot <- c("Bowel (n=438)","Pancreas (n=328)","Head&Neck (n=311)","Soft tissue (n=284)")
  
  ## scales = "free": each sub-panel has its own axis ranges, and the
  ## correlation labels are placed in npc units so they sit in the same corner
  ## of every sub-panel.
  fig2d_PDX_by_cancer_type <- ggplot(subset(h_df, lineage %in% lineage_plot), aes(x = log10(x), y = log10(y))) +
    geom_point(alpha = 0.4, colour = "#808080", size = 0.9) +
    geom_smooth(method = "lm", colour = "#005AB5", se = FALSE, na.rm = TRUE) +
    ## ---- correlation label, with a white halo ------------------------------
  ## The label sits on top of the point cloud here, so the glyphs need
  ## separating from what is behind them. bg.color draws a halo AROUND each
  ## glyph rather than a filled box behind the whole string the way
  ## geom = "label" would, which keeps the points visible through the gaps in
  ## the text. bg.r is the halo radius in units of font size: 0.1 is a
  ## hairline, 0.25 is heavy.
  ##
  ## PASS THE GEOM OBJECT, NOT THE STRING "text_repel". ggplot2 resolves a
  ## geom given as a string by searching the ATTACHED packages, and ggpubr
  ## only imports ggrepel rather than attaching it — so the string form dies
  ## with `Can't find geom called "text_repel"` unless library(ggrepel) has
  ## been run. ggrepel::GeomTextRepel works either way, and since ggpubr
  ## imports ggrepel it is always installed wherever stat_cor is.
  ##
  ## text_repel would normally MOVE the label to dodge collisions, which
  ## would quietly override the fixed label.y and put it somewhere different
  ## in each sub-panel. force = 0 turns the repulsion off, seed makes what
  ## little remains reproducible, and segment.color = NA drops the leader
  ## line it would otherwise draw back to the anchor point.
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
    facet_wrap(~ lineage, ncol = 2, scales = "free") +
    xlab(TeX("$log_{10}(DT)\\ (days)$")) +
    ylab(TeX("$log_{10}(pDT)$")) +
    theme(legend.position  = "none",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line        = element_line(colour = "black"),
          axis.text        = element_text(size = 14,  face = "bold"),
          axis.title       = element_text(size = 16, face = "bold"),
          strip.background = element_blank(),
          strip.text       = element_text(size = 13,  face = "bold"))
  
  ggsave(file.path(NEW_FIG_DIR, "Figure_2d_PDX_ensemble_scatter_by_cancer_type.pdf"), fig2d_PDX_by_cancer_type,
         width  = 6,
         height = 6,
         device = cairo_pdf)
  cat(sprintf("     Saved Figure_2d_PDX_ensemble_scatter_by_cancer_type.pdf (%d sub-panels)\n", n_facets))
  
}, error = function(e) cat("     FIGURE 2d FAILED:", conditionMessage(e), "\n"))



################################################################################
# 7.7  FIGURE 2a — PDX INDIVIDUAL MODELS (rho vs number of genes)
#      PDX: mean validation Spearman rho vs number of genes, every trained model.
#      Printed as panel a (top left).  Data: FIG2$C (PDX — set in 7.3).
################################################################################

cat("\n7.7  Figure 2a — PDX individual models\n")

fig2a_PDX_individual_models <- NULL
tryCatch({
  
  sm <- FIG2$C
  if (is.null(sm)) stop("FIG2$C is empty — rebuild in 7.3")
  sm <- add_sem(sm)
  
  ## Backstop (7.4b already checks this): plotting the wrong dataset under the
  ## right letter is the one failure that is invisible in the output.
  if (!identical(attr(sm, "dataset"), "PDX_full_sig_3"))
    stop("FIG2$C holds '", attr(sm, "dataset"), "', expected 'PDX_full_sig_3'.\n",
         "  Run 7.3 to rebuild figure2_data.rds, then 7.4.")
  
  ## one colour per model, same order in 2a and 2b so a colour always means the
  ## same model
  model_cols <- c(
    "ElasticNet_alpha_0.25" = "#33a02c", "ElasticNet_alpha_0.5"  = "#fb9a99",  # col 1
    "ElasticNet_alpha_0.75" = "#e31a1c", "RandomForest"          = "#6a3d9a",  # col 2
    "SVR_linear"            = "#b15928", "SVR_radial"            = "#000000",  # col 3
    "PLS_1_components"      = "#fdbf6f", "PLS_2_components"      = "#ff7f00",  # col 4
    "PLS_5_components"      = "#cab2d6", "XGBoost"               = "#8B4513",  # col 5
    "Linear"                = "#a6cee3", "Ridge"                 = "#1f78b4",  # col 6
    "Lasso"                 = "#b2df8a")                                       # col 7
  
  present <- names(model_cols)[names(model_cols) %in% as.character(sm$model)]
  sm$model <- factor(as.character(sm$model), levels = present)
  
  model_labs <- c(
    "ElasticNet_alpha_0.25" = "ElasticNet \u03b1 = 0.25",
    "ElasticNet_alpha_0.5"  = "ElasticNet \u03b1 = 0.5",
    "ElasticNet_alpha_0.75" = "ElasticNet \u03b1 = 0.75",
    "PLS_1_components"      = "PLS 1",
    "PLS_2_components"      = "PLS 2",
    "PLS_5_components"      = "PLS 5",
    "RandomForest"          = "RandomForest",
    "SVR_linear"            = "SVR linear",
    "SVR_radial"            = "SVR radial",
    "XGBoost"               = "XGBoost",
    "Linear"                = "Linear",
    "Ridge"                 = "Ridge",
    "Lasso"                 = "Lasso")
  
  fig2a_PDX_individual_models <- ggplot(sm, aes(x = top_n_genes, y = Spearman_mean,
                            colour = model, group = model)) +
    geom_vline(xintercept = 100, linetype = "dashed",
               colour = "grey30", linewidth = 0.5) +
    geom_line(linewidth = 0.8) +
    sem_errorbar() +                       # mean +/- SEM over folds
    geom_point(size = 2.2) +
    scale_x_continuous(breaks = c(10, 20, 50, 100, 200, 500), trans = "log10",
                       labels = c(10, 20, 50, 100, 200, 500)) +
    scale_colour_manual(values = model_cols[present],
                        labels = model_labs[present], name = NULL,
                        guide = guide_legend(nrow = 2, byrow = FALSE)) +
    labs(x = TeX("$Number\\ of\\ genes$"),
         y = TeX("$Mean\\ \\rho$")) +
    ## ---- custom text ---------------------------------------------------------
  ## x is a DATA value (10 = the leftmost break). y = Inf is the top of the
  ## panel, which is safe here because the y scale is linear — it is the LOG10
  ## x axis that turns -Inf into NaN and triggers "In log(x, base) : NaNs
  ## produced". vjust pushes the text back down inside the panel.
  annotate("text", x = 10, y = Inf,
           label = "PDX\n13 models",
           hjust = 0, vjust = 1.25,
           size = 4.5, fontface = "bold", colour = "grey20") +
    ## --------------------------------------------------------------------------
    theme_bw(base_size = 10) +
    theme(legend.position  = "bottom",
          legend.text = element_text(size = 12),
          panel.grid.minor = element_blank(),
          axis.text        = element_text(size = 14,  face = "bold"),
          axis.title       = element_text(size = 16, face = "bold"))
  
  ggsave(file.path(NEW_FIG_DIR, "Figure_2a_PDX_individual_models.pdf"), fig2a_PDX_individual_models,
         width = 12, height = 8.5, device = cairo_pdf)
  cat("     Saved Figure_2a_PDX_individual_models.pdf\n")
  
}, error = function(e) cat("     FIGURE 2a FAILED:", conditionMessage(e), "\n"))


################################################################################
# 7.8  FIGURE 2b — CCLE INDIVIDUAL MODELS (rho vs number of genes)
#      CCLE: same as Figure 2a.
#      Printed as panel b (top right).  Data: FIG2$D (CCLE — set in 7.3).
################################################################################

cat("\n7.8  Figure 2b — CCLE individual models\n")

fig2b_CCLE_individual_models <- NULL
tryCatch({
  
  sm <- FIG2$D
  if (is.null(sm)) stop("FIG2$D is empty — rebuild in 7.3")
  sm <- add_sem(sm)
  
  ## Same backstop as Figure 2a.
  if (!identical(attr(sm, "dataset"), "CCLE_norm"))
    stop("FIG2$D holds '", attr(sm, "dataset"), "', expected 'CCLE_norm'.\n",
         "  Run 7.3 to rebuild figure2_data.rds, then 7.4.")
  
  model_cols <- c(
    "ElasticNet_alpha_0.25" = "#33a02c", "ElasticNet_alpha_0.5"  = "#fb9a99",  # col 1
    "ElasticNet_alpha_0.75" = "#e31a1c", "RandomForest"          = "#6a3d9a",  # col 2
    "SVR_linear"            = "#b15928", "SVR_radial"            = "#000000",  # col 3
    "PLS_1_components"      = "#fdbf6f", "PLS_2_components"      = "#ff7f00",  # col 4
    "PLS_5_components"      = "#cab2d6", "XGBoost"               = "#8B4513",  # col 5
    "Linear"                = "#a6cee3", "Ridge"                 = "#1f78b4",  # col 6
    "Lasso"                 = "#b2df8a")                                       # col 7
  
  
  present <- names(model_cols)[names(model_cols) %in% as.character(sm$model)]
  sm$model <- factor(as.character(sm$model), levels = present)
  
  model_labs <- c(
    "ElasticNet_alpha_0.25" = "ElasticNet \u03b1 = 0.25",
    "ElasticNet_alpha_0.5"  = "ElasticNet \u03b1 = 0.5",
    "ElasticNet_alpha_0.75" = "ElasticNet \u03b1 = 0.75",
    "PLS_1_components"      = "PLS 1",
    "PLS_2_components"      = "PLS 2",
    "PLS_5_components"      = "PLS 5",
    "RandomForest"          = "RandomForest",
    "SVR_linear"            = "SVR linear",
    "SVR_radial"            = "SVR radial")
  
  fig2b_CCLE_individual_models <- ggplot(sm, aes(x = top_n_genes, y = Spearman_mean,
                            colour = model, group = model)) +
    geom_vline(xintercept = 100, linetype = "dashed",
               colour = "grey30", linewidth = 0.5) +
    geom_line(linewidth = 0.8) +
    sem_errorbar() +                       # mean +/- SEM over folds
    geom_point(size = 2.2) +
    scale_x_continuous(breaks = c(10, 20, 50, 100, 200, 500), trans = "log10",
                       labels = c(10, 20, 50, 100, 200, 500)) +
    scale_colour_manual(values = model_cols[present],
                        labels = model_labs[present], name = NULL,
                        guide = guide_legend(nrow = 2, byrow = FALSE)) +
    labs(x = TeX("$Number\\ of\\ genes$"),
         y = TeX("$Mean\\ \\rho$")) +
    ## ---- custom text ---------------------------------------------------------
  ## x is a DATA value (10 = the leftmost break); y = 0.3 is in rho units.
  annotate("text", x = 10, y = 0.3,
           label = "CCLE\n13 models",
           hjust = 0, vjust = 1.25,
           size = 4.5, fontface = "bold", colour = "grey20") +
    ## --------------------------------------------------------------------------
    theme_bw(base_size = 10) +
    theme(legend.position  = "bottom",
          legend.text = element_text(size = 12),
          panel.grid.minor = element_blank(),
          axis.text        = element_text(size = 14,  face = "bold"),
          axis.title       = element_text(size = 16, face = "bold"))
  
  ggsave(file.path(NEW_FIG_DIR, "Figure_2b_CCLE_individual_models.pdf"), fig2b_CCLE_individual_models,
         width = 12, height = 8.5, device = cairo_pdf)
  cat("     Saved Figure_2b_CCLE_individual_models.pdf\n")
  
}, error = function(e) cat("     FIGURE 2b FAILED:", conditionMessage(e), "\n"))

################################################################################
# 7.9  FIGURE 2c — PDX ENSEMBLES (rho vs number of genes)
#      PDX ensembles: the three pairs from the trio plus the full trio.
#      Printed as panel c (middle left).  Data: FIG2$E (per-fold prediction
#      averaging already done in 7.3).
#
# The matching CCLE curve is not a Figure 2 panel; it is panel c of the
# supplementary figure and Table 3.4, from FIG2$E_CCLE.
################################################################################

cat("\n7.9  Figure 2c — PDX ensembles\n")

fig2c_PDX_ensembles <- NULL
tryCatch({
  
  ens_df <- FIG2$E
  if (is.null(ens_df)) stop("FIG2$E is empty — rebuild in 7.3")
  ens_df <- add_sem(ens_df)
  
  ## Backstop (7.4b already checks this): the panel would otherwise plot CCLE
  ## where PDX belongs.
  if (!identical(attr(ens_df, "dataset"), "PDX_full_sig_3"))
    stop("FIG2$E holds '", attr(ens_df, "dataset"), "', expected 'PDX_full_sig_3'.\n",
         "  Run 7.3 to rebuild figure2_data.rds, then 7.4.")
  
  ens_cols <- c("RandomForest + Ridge"              = "#1f78b4",
                "RandomForest + SVR_radial"         = "#33a02c",
                "Ridge + SVR_radial"                = "#ff7f00",
                "RandomForest + Ridge + SVR_radial" = "#000000")
  present <- names(ens_cols)[names(ens_cols) %in% as.character(ens_df$ensemble)]
  ens_df$ensemble <- factor(as.character(ens_df$ensemble), levels = present)
  
  ens_labs <- c(
    "RandomForest + Ridge"              = "RF + Ridge",
    "RandomForest + SVR_radial"         = "RF + SVR radial",
    "Ridge + SVR_radial"                = "Ridge + SVR radial",
    "RandomForest + Ridge + SVR_radial" = "RF + Ridge + SVR radial")
  
  fig2c_PDX_ensembles <- ggplot(ens_df, aes(x = top_n_genes, y = Spearman_mean,
                                colour = ensemble, group = ensemble)) +
    geom_vline(xintercept = 100, linetype = "dashed",
               colour = "grey30", linewidth = 0.5) +
    geom_line(linewidth = 0.8) +
    ## full trio drawn thicker so it reads against the three pairs
    geom_line(data = ens_df[ens_df$n_members == 3, ], linewidth = 1.5) +
    sem_errorbar() +                       # mean +/- SEM over folds
    geom_point(size = 2.2) +
    scale_x_continuous(breaks = c(10, 20, 50, 100, 200, 500), trans = "log10",
                       labels = c(10, 20, 50, 100, 200, 500)) +
    scale_colour_manual(values = ens_cols[present],
                        labels = ens_labs[present], name = NULL,
                        guide  = guide_legend(nrow = 2, byrow = FALSE)) +
    labs(x = TeX("$Number\\ of\\ genes$"),
         y = TeX("$Mean\\ \\rho$")) +
    ## ---- custom text ---------------------------------------------------------
  ## x is a DATA value (10 = the leftmost break); y = 0.54 is in rho units.
  annotate("text", x = 10, y = 0.54,
           label = "PDX Ensemble models",
           hjust = 0, vjust = 1.25,
           size = 4.5, fontface = "bold", colour = "grey20") +
    ## --------------------------------------------------------------------------
    theme_bw(base_size = 10) +
    theme(legend.position  = "bottom",
          legend.text = element_text(size = 12),
          panel.grid.minor = element_blank(),
          axis.text        = element_text(size = 14,  face = "bold"),
          axis.title       = element_text(size = 16, face = "bold"))
  
  ggsave(file.path(NEW_FIG_DIR, "Figure_2c_PDX_ensembles.pdf"), fig2c_PDX_ensembles,
         width = 12, height = 8.5, device = cairo_pdf)
  cat("     Saved Figure_2c_PDX_ensembles.pdf\n")
  
}, error = function(e) cat("     FIGURE 2c FAILED:", conditionMessage(e), "\n"))


################################################################################
# 7.10  FIGURE 2f — GENE-LEVEL RHO, CCLE vs PDX
#      Gene-level agreement between the two datasets: pooled pan-cancer
#      Spearman rho against doubling time in CCLE (x) vs the same in PDX (y),
#      one point per gene.  Data: FIG2$F, built in 7.4c
#      Printed as panel f (bottom right).
#
# Not how well a model transfers, but how far the underlying gene-level signal
# itself is shared.
#
# GREY   every gene measured in both datasets
# BLUE   the highlighted set — same-sign genes that are top-|rho| in BOTH
#        datasets at the shallowest depth where the two lists overlap by at
#        least GENE_CORR$min_overlap genes (see 7.4c)
#
# The correlation statistic is computed on ALL the points, not just the
# highlighted ones, because stat_cor inherits the plot's data.
################################################################################

cat("\n7.10  Figure 2f — gene-level rho CCLE vs PDX\n")

fig2f_gene_rho_CCLE_vs_PDX <- NULL
tryCatch({
  
  gc_df <- FIG2$F
  if (is.null(gc_df)) stop("FIG2$F is empty — run 7.4c")
  if (!all(c("rho_ccle", "rho_pdx", "top_set") %in% names(gc_df)))
    stop("FIG2$F is not the gene-correlation table — run 7.4c")
  
  hl <- gc_df[gc_df$top_set, ]
  
  fig2f_gene_rho_CCLE_vs_PDX <- ggplot(gc_df, aes(x = rho_ccle, y = rho_pdx)) +
    geom_point(size = 0.9, alpha = 0.5, colour = "grey50") +
    ## highlighted genes drawn on top so they are not buried in the cloud
    geom_point(data = hl, size = 0.9, alpha = 1, colour = "blue") +
    geom_smooth(method = "lm", colour = "#005AB5", se = FALSE, na.rm = TRUE) +
    stat_cor(aes(label = paste(..r.label..,
                               after_stat(format_p_value(..p..)),
                               sep = "~`,`~")),
             method = "spearman", r.accuracy = 0.001, parse = TRUE, size = 4.5,
             geom = ggrepel::GeomTextRepel, bg.color = "white", bg.r = 0.25,
             force = 0, seed = 1, segment.color = NA) +
    xlab(TeX("$CCLE\\ \\rho$")) +
    ylab(TeX("$PDX\\ \\rho$")) +
    theme(legend.position  = "none",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line        = element_line(colour = "black"),
          axis.text        = element_text(size = 14,  face = "bold"),
          axis.title       = element_text(size = 16, face = "bold"))
  
  ggsave(file.path(NEW_FIG_DIR, "Figure_2f_gene_rho_CCLE_vs_PDX.pdf"), fig2f_gene_rho_CCLE_vs_PDX,
         width = 6, height = 6, device = cairo_pdf)
  cat(sprintf("     Saved Figure_2f_gene_rho_CCLE_vs_PDX.pdf (%d genes, %d highlighted)\n",
              nrow(gc_df), nrow(hl)))
  
}, error = function(e) cat("     FIGURE 2f FAILED:", conditionMessage(e), "\n"))
  
################################################################################
# 7.11  FIGURE 2 — ASSEMBLE
#
# Six panels in a 2-column grid, 3 rows, as printed:
#   a b   Figure 2a PDX individual models | Figure 2b CCLE individual models
#   c d   Figure 2c PDX ensembles         | Figure 2d PDX scatter by cancer type
#   e f   Figure 2e PDX ensemble scatter  | Figure 2f gene-level rho, CCLE vs PDX
#
# The plot objects and the standalone PDFs are named by their PRINTED letter
# (fig2a_..., Figure_2a_....pdf). If you change the layout or tag_levels below,
# rename them to match. The FIG2 slot letters are internal keys; the mapping is
#   FIG2$C -> 2a, $D -> 2b, $E -> 2c, $B -> 2d, $A -> 2e, $F -> 2f.
#
# The two domain-swap scatters are a supplementary figure of their own, built
# in 7.12-7.14 below, from FIG2$G and FIG2$H.
#
# plot_annotation(tag_levels = ...) draws the letters a-f. Tags are handed out
# in the order the plots are ADDED (2e, 2d, 2a, 2b, 2c, 2f), which is why the
# tag vector below is not alphabetical.
#
# The standalone panel PDFs carry no tag; their file name says which panel
# they are (Figure_2a_PDX_individual_models.pdf and so on).
################################################################################

cat("\n7.11  Assembling Figure_2.pdf\n")

## Every panel must exist — `+` chaining onto a NULL fails with an unhelpful
## error, so check first and say which one is missing.
.fig2_objs <- c("fig2a_PDX_individual_models", "fig2b_CCLE_individual_models",
                "fig2c_PDX_ensembles", "fig2d_PDX_by_cancer_type",
                "fig2e_PDX_ensemble_scatter", "fig2f_gene_rho_CCLE_vs_PDX")
.have <- vapply(.fig2_objs, function(k) !is.null(get0(k)), logical(1))
if (any(!.have))
  stop("missing panel(s): ", paste(.fig2_objs[!.have], collapse = ", "),
       " — rebuild them above before assembling")
rm(.have, .fig2_objs)

## Legend collection only merges IDENTICAL legends. Check the pair actually
## matches, so a silently duplicated legend does not come as a surprise.
## Only 2a and 2b share one — 2c carries its own and 2f has none.
for (.pair in list(list(a = "C", b = "D", la = "2a", lb = "2b", col = "model", what = "Model"))) {
  .x <- FIG2[[.pair$a]]; .y <- FIG2[[.pair$b]]
  if (is.null(.x) || is.null(.y)) next
  .lx <- sort(unique(as.character(.x[[.pair$col]])))
  .ly <- sort(unique(as.character(.y[[.pair$col]])))
  if (!identical(.lx, .ly)) {
    cat(sprintf("     NOTE: panels %s and %s have different %s entries, so the\n",
                .pair$la, .pair$lb, .pair$what))
    cat(sprintf("           legend will NOT merge. Only in %s: %s | only in %s: %s\n",
                .pair$la, paste(setdiff(.lx, .ly), collapse = ", "),
                .pair$lb, paste(setdiff(.ly, .lx), collapse = ", ")))
  } else {
    cat(sprintf("     %s legend shared by panels %s and %s (%d entries)\n",
                .pair$what, .pair$la, .pair$lb, length(.lx)))
  }
}
rm(.pair, .x, .y, .lx, .ly)

## ---- row 1: Figures 2a and 2b side by side, with ONE shared legend ----------
row_cd <- fig2a_PDX_individual_models + fig2b_CCLE_individual_models +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

## Figure 2c keeps its own ensemble legend and 2f has none, so they are added
## to the figure individually rather than as a legend-collecting row: a
## `guides = "collect"` with a trailing `&` would force an empty legend slot
## onto 2f and shrink it.

## ---- full figure -----------------------------------------------------------
## No guides argument at top level, so the 2a/2b legend stays under row 1 and
## 2c keeps its own.
fig2 <- fig2e_PDX_ensemble_scatter +
  wrap_elements(full = fig2d_PDX_by_cancer_type) +
  row_cd +
  fig2c_PDX_ensembles +
  wrap_elements(full = fig2f_gene_rho_CCLE_vs_PDX) +
  plot_layout(
    ## design letters are patchwork's, NOT panel names: A..E are the 1st..5th
    ## plots added above (2e, 2d, the 2a+2b row, 2c, 2f).
    design = paste(c("CC",
                     "DB",
                     "AE"),
                   collapse = "\n"),
    ## one entry per design row: row 1 = 2a/2b, row 2 = 2c/2d (2d holds a grid
    ## of sub-panels), row 3 = 2e/2f.
    heights = c(2, 3, 3)) +
  plot_annotation(tag_levels = list(c("e","d","a","b","c","f"))) &
  theme(plot.tag    = element_text(size = 16, face = "bold"))

## Figure size: 10 x 12 in (PDF) and 10 x 14 in (PNG). Change here if the rows
## look cramped.
## Rendering is where a panel that BUILT fine can still fail — a facet with no
## rows, a label outside the scale. Without this the error would abort the whole
## script and take the supplementary figure and table with it, so catch it: the
## six standalone panel PDFs are already on disk and are unaffected.
tryCatch({
  ggsave(file.path(NEW_FIG_DIR, "Figure_2.pdf"), fig2,
         width = 10, height = 12, limitsize = FALSE)
  ggsave(file.path(NEW_FIG_DIR, "Figure_2.png"), fig2,
         width = 10, height = 14, limitsize = FALSE)
  message(sprintf("\n>>> Figure_2 -> %s",
                  file.path(NEW_FIG_DIR, "Figure_2.pdf")))
}, error = function(e) {
  message("\n>>> ASSEMBLY FAILED: ", conditionMessage(e))
  message("    The individual panel PDFs were still written; check which panel")
  message("    reported a problem above. Carrying on.")
})


################################################################################
################################################################################
##                                                                            ##
##   SUPPLEMENTARY FIGURE 2  (7.12 - 7.14)                                    ##
##                                                                            ##
##     a   CCLE measured DT vs CCLE predicted by the PDX-trained trio         ##
##     b   PDX  measured DT vs PDX  predicted by the CCLE-trained trio        ##
##     c   CCLE ensembles vs number of genes (7.13b)                          ##
##                                                                            ##
##   Data: FIG2$G, FIG2$H and FIG2$E_CCLE, built in 7.3 and cached in         ##
##   figure2_data.rds.                                                        ##
##   The plot objects are called sup_swap_A / sup_swap_B so nothing collides   ##
##   with the fig2a_..fig2f_ objects of Figure 2.                             ##
##                                                                            ##
##   Rename the figure in ONE place, here. The panel file names follow this   ##
##   tag.                                                                     ##
##                                                                            ##
################################################################################
################################################################################

SUP_FIG2_TAG <- "Sup_Figure_2"   # <<< figure name lives here; the
                                 ## panel files follow it

################################################################################
# 7.12  SUPPLEMENTARY DOMAIN-SWAP FIGURE — PANEL a
#       CCLE measured doubling time (x) vs CCLE predicted by the PDX-trained
#       ensemble trio (y).  Data: FIG2$G
#
# Y AXIS: the trio was trained on PDX doubling times, so its predictions are on
# the PDX scale — not calibrated CCLE days, and the fitted slope has no
# absolute meaning. Read this as rank agreement, which is what Spearman reports.
################################################################################

cat("\n7.12  Domain-swap supplementary figure, panel a\n")

sup_swap_A <- NULL
tryCatch({
  
  h_df <- FIG2$G
  if (is.null(h_df)) stop("FIG2$G is empty — rebuild in 7.3")
  
  
  sup_swap_A <- ggplot(h_df, aes(x = log10(x), y = log10(y))) +
    geom_point(alpha = 0.4, colour = "#808080", size = 0.9) +
    geom_smooth(method = "lm", colour = "#005AB5", se = FALSE, na.rm = TRUE) +
    stat_cor(aes(label = paste(..r.label..,
                               after_stat(format_p_value(..p..)),
                               sep = "~`,`~")),
             method = "spearman", r.accuracy = 0.001, parse = TRUE, size = 2.8,
             label.y = 1.6) +
    xlab(TeX("$log_{10}(DT,\\ CCLE) (hours)$")) +
    ylab(TeX("$log_{10}(pDT,\\ PDX-trained)$")) +
    theme(legend.position  = "none",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line        = element_line(colour = "black"),
          axis.text        = element_text(size = 8,  face = "bold"),
          axis.title       = element_text(size = 10, face = "bold"))
  
  ggsave(file.path(NEW_FIG_DIR, sprintf("%s_panela.pdf", SUP_FIG2_TAG)),
         sup_swap_A, width = 6, height = 6, device = cairo_pdf)
  cat(sprintf("     Saved %s_panela.pdf (%d genes used)\n",
              SUP_FIG2_TAG, attr(h_df, "genes_used")))
  
}, error = function(e) cat("     SWAP PANEL a FAILED:", conditionMessage(e), "\n"))


################################################################################
# 7.13  SUPPLEMENTARY DOMAIN-SWAP FIGURE — PANEL b
#       PDX measured doubling time (x) vs PDX predicted by the CCLE-trained
#       ensemble trio (y).  Data: FIG2$H
#
# Same y-axis caveat as panel a, in the other direction.
################################################################################

cat("\n7.13  Domain-swap supplementary figure, panel b\n")

sup_swap_B <- NULL
tryCatch({
  
  h_df <- FIG2$H
  if (is.null(h_df)) stop("FIG2$H is empty — rebuild in 7.3")
  
  
  sup_swap_B <- ggplot(h_df, aes(x = log10(x), y = log10(y))) +
    geom_point(alpha = 0.4, colour = "#808080", size = 0.9) +
    geom_smooth(method = "lm", colour = "#005AB5", se = FALSE, na.rm = TRUE) +
    stat_cor(aes(label = paste(..r.label..,
                               after_stat(format_p_value(..p..)),
                               sep = "~`,`~")),
             method = "spearman", r.accuracy = 0.001, parse = TRUE, size = 2.8,
             label.y = 2.2) +
    xlab(TeX("$log_{10}(DT,\\ PDX) (days)$")) +
    ylab(TeX("$log_{10}(pDT,\\ CCLE-trained)$")) +
    theme(legend.position  = "none",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line        = element_line(colour = "black"),
          axis.text        = element_text(size = 8,  face = "bold"),
          axis.title       = element_text(size = 10, face = "bold"))
  
  ggsave(file.path(NEW_FIG_DIR, sprintf("%s_panelb.pdf", SUP_FIG2_TAG)),
         sup_swap_B, width = 6, height = 6, device = cairo_pdf)
  cat(sprintf("     Saved %s_panelb.pdf (%d genes used)\n",
              SUP_FIG2_TAG, attr(h_df, "genes_used")))
  
}, error = function(e) cat("     SWAP PANEL b FAILED:", conditionMessage(e), "\n"))


################################################################################
# 7.13b  SUPPLEMENTARY FIGURE 2 — PANEL c
#        CCLE: ensembles — the three pairs from the trio plus the full trio.
#        Data: FIG2$E_CCLE
#
# The CCLE companion to Figure 2c (PDX ensembles); its values are Table 3.4.
#
# Deliberately identical to 7.9 apart from the slot it reads and the note it
# draws, so the two can be laid side by side and compared without wondering
# whether a styling difference is doing any of the work.
################################################################################

cat("\n7.13b Supplementary Figure 2, panel c\n")

sup_ens_CCLE <- NULL
tryCatch({

  ens_df <- FIG2$E_CCLE
  if (is.null(ens_df))
    stop("FIG2$E_CCLE is empty — re-run 7.3 to build it")

  ## Backstop, as in the main panels: E_CCLE must hold CCLE, not PDX.
  if (!identical(attr(ens_df, "dataset"), "CCLE_norm"))
    stop("FIG2$E_CCLE holds '", attr(ens_df, "dataset"), "', expected 'CCLE_norm'.")
  ens_df <- add_sem(ens_df)

  ens_cols <- c("RandomForest + Ridge"              = "#1f78b4",   # col 1
                "RandomForest + SVR_radial"         = "#33a02c",
                "Ridge + SVR_radial"                = "#ff7f00",   # col 2
                "RandomForest + Ridge + SVR_radial" = "#000000")
  present <- names(ens_cols)[names(ens_cols) %in% as.character(ens_df$ensemble)]
  ens_df$ensemble <- factor(as.character(ens_df$ensemble), levels = present)

  ens_labs <- c(
    "RandomForest + Ridge"              = "RF + Ridge",
    "RandomForest + SVR_radial"         = "RF + SVR radial",
    "Ridge + SVR_radial"                = "Ridge + SVR radial",
    "RandomForest + Ridge + SVR_radial" = "RF + Ridge + SVR radial")

  sup_ens_CCLE <- ggplot(ens_df, aes(x = top_n_genes, y = Spearman_mean,
                                     colour = ensemble, group = ensemble)) +
    geom_vline(xintercept = 100, linetype = "dashed",
               colour = "grey30", linewidth = 0.5) +
    geom_line(linewidth = 0.8) +
    ## full trio drawn thicker so it reads against the three pairs
    geom_line(data = ens_df[ens_df$n_members == 3, ], linewidth = 1.5) +
    sem_errorbar() +                       # mean +/- SEM over folds
    geom_point(size = 2.2) +
    scale_x_continuous(breaks = c(10, 20, 50, 100, 200, 500), trans = "log10",
                       labels = c(10, 20, 50, 100, 200, 500)) +
    scale_colour_manual(values = ens_cols[present],
                        labels = ens_labs[present], name = NULL,
                        guide = guide_legend(nrow = 2, byrow = FALSE)) +
    labs(x = TeX("$Number\\ of\\ genes$"),
         y = TeX("$Mean\\ \\rho$")) +
    annotate("text", x = 10, y = Inf,
             label = "CCLE Ensemble models",
             hjust = 0, vjust = 1.25,
             size = 4.5, fontface = "bold", colour = "grey20") +
    theme_bw(base_size = 10) +
    theme(legend.position  = "bottom",
          legend.text      = element_text(size = 12),
          panel.grid.minor = element_blank(),
          axis.text        = element_text(size = 14, face = "bold"),
          axis.title       = element_text(size = 16, face = "bold"))

  ggsave(file.path(NEW_FIG_DIR, sprintf("%s_panelc.pdf", SUP_FIG2_TAG)),
         sup_ens_CCLE, width = 12, height = 8.5, device = cairo_pdf)
  cat(sprintf("     Saved %s_panelc.pdf\n", SUP_FIG2_TAG))

}, error = function(e) cat("     SUP PANEL c FAILED:", conditionMessage(e), "\n"))


################################################################################
# 7.14  SUPPLEMENTARY FIGURE 2 — ASSEMBLE
#
#   a b   the two domain-swap scatters, side by side
#   c     CCLE ensembles vs gene count, spanning the full width
#
# Panel c is a line plot with its own legend while a and b are plain scatters
# with none, so nothing is collected here — guides = "collect" would only force
# an empty legend slot onto the two scatters and shrink them.
#
# Panel c gets its own row rather than a third column: at a third of the width
# the six x-axis breaks and the four-entry legend would both be unreadable.
#
# The three standalone panel PDFs are written by 7.12, 7.13 and 7.13b
# regardless of whether this assembly succeeds.
################################################################################

cat("\n7.14  Assembling Supplementary Figure 2\n")

.sup_parts <- c(a = "sup_swap_A", b = "sup_swap_B", c = "sup_ens_CCLE")
.have <- vapply(.sup_parts, function(k) !is.null(get0(k)), logical(1))
if (any(!.have)) {
  cat(sprintf("     SKIPPED — missing panel(s): %s\n",
              paste(names(.sup_parts)[!.have], collapse = ", ")))
} else {

  fig_sup2 <- sup_swap_A + sup_swap_B + sup_ens_CCLE +
    plot_layout(design = paste(c("AB",
                                 "CC"),
                               collapse = "\n"),
                ## row 1 holds two square scatters, row 2 the wide curve panel
                ## plus its legend, so they get roughly equal height
                heights = c(1, 1)) +
    plot_annotation(
      tag_levels = "a",
      theme = theme(plot.tag = element_text(size = 15, face = "bold"))) &
    theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1))

  SUP_FIG2_FILE <- file.path(NEW_FIG_DIR, sprintf("%s.pdf", SUP_FIG2_TAG))
  tryCatch({
    ggsave(SUP_FIG2_FILE, fig_sup2,
           width = 6, height = 6, device = cairo_pdf, limitsize = FALSE)
    ggsave(file.path(NEW_FIG_DIR, sprintf("%s.png", SUP_FIG2_TAG)), fig_sup2,
           width = 6, height = 6, limitsize = FALSE)
    message(sprintf("\n>>> %s -> %s", SUP_FIG2_TAG, SUP_FIG2_FILE))
  }, error = function(e) {
    message("\n>>> SUPPLEMENTARY FIGURE 2 ASSEMBLY FAILED: ", conditionMessage(e))
    message("    Panels a, b and c were still written separately.")
  })
}
rm(.have, .sup_parts)


################################################################################
# Section 7.15  SUPPLEMENTARY TABLE 3  (sup_tab3.xlsx)
#
#   Contents     what each sheet holds
#   Table 3.1    Figure 2a        PDX,  per-model performance vs gene count
#   Table 3.2    Figure 2b        CCLE, the same
#   Table 3.3    Figure 2c        PDX,  validation Spearman by ensemble and gene count
#   Table 3.4    Sup. Figure 2c   CCLE, the same
#   Table 3.5    Figure 2f        the 100 genes used for prediction
#
# Figure 2d/2e (the ensemble scatters), Figure 2f's per-gene values and
# Supplementary Figure 2a/2b (domain swap) have no sheet; their data are cached
# in panel_ab_ensemble.rds, panel_f_gene_corr.rds and figure2_data.rds.
#
# Tables 3.1 and 3.2 are built from the out-of-fold predictions (pred_df in
# PIPE: in memory after 7.0, otherwise read from each dataset's
# pipeline_objects.rds). Tables 3.3 and 3.4 read FIG2 only. Any sheet whose
# source is missing is skipped and Contents says so. Only the cor_all ranking
# is used, as in every Figure 2 panel.
#
# ---------------------------------------------------------------------------
# NO P-VALUES ARE REPORTED. Tables 3.1 and 3.2 give the five per-fold rho
# values and their mean and SD, which is what Figures 2a and 2b are built from.
#
# There is no valid p-value to attach to that mean: averaging the folds'
# p-values is not a legitimate combination, and the five folds share 4/5 of
# their training data, so they are not independent tests either. calc_metrics()
# uses cor() rather than cor.test() for the same reason, and the p_value /
# p_adj columns that do exist in the pipeline belong to the GENE-level
# rankings — a different quantity.
# ---------------------------------------------------------------------------
################################################################################

cat("\n7.15  Supplementary table\n")

SUP_TAB_FILE <- file.path(NEW_FIG_DIR, "sup_tab3.xlsx")
SUP_XL_MAX_ROWS <- 1048575L   # 1,048,576 minus the header row

## Gene-ranking method used in the curve tables (as in every Figure 2 panel).
SUP_RANKINGS <- "cor_all"

## Gene direction for Tables 3.1 and 3.2 — must match what Figures 2a and 2b
## plot. Taken from the FIG2 build stamp when present.
SUP_PANEL_DIR <- {
  .m <- attr(FIG2, "meta")
  if (!is.null(.m) && !is.null(.m$direction)) .m$direction else "both"
}

## Short dataset labels. The DATASETS names carry preprocessing detail that
## does not belong in a supplement.
SUP_DS_LABEL <- c("CCLE_norm" = "CCLE", "PDX_full_sig_3" = "PDX")

.sup_label <- function(x) {
  x <- as.character(x)
  ifelse(is.na(x), NA_character_,
         ifelse(x %in% names(SUP_DS_LABEL), unname(SUP_DS_LABEL[x]), x))
}

################################################################################
# pred_CCLE and pred_PDX — every out-of-fold prediction, per dataset
#
# Inputs, not sheets: Tables 3.1 and 3.2 are derived from them.
################################################################################

## Load one dataset's pipeline objects, from PIPE if present, else off disk.
.sup_get_pipe <- function(ds_name) {
  if (exists("PIPE") && !is.null(PIPE[[ds_name]])) return(PIPE[[ds_name]])
  .hit <- vapply(DATASETS, function(d) identical(d$name, ds_name), logical(1))
  if (!any(.hit)) return(NULL)
  .f <- file.path(DATASETS[[which(.hit)[1]]]$out_dir, "pipeline_objects.rds")
  if (!file.exists(.f)) return(NULL)
  cat(sprintf("     reading %s\n", .f))
  readRDS(.f)
}

.sup_predictions <- function(ds_name) {
  P <- .sup_get_pipe(ds_name)
  if (is.null(P))         return(NULL)
  if (is.null(P$pred_df)) return(NULL)

  pr <- P$pred_df
  pr <- pr[pr$gene_direction %in% c("negative", "positive", "both"), ]
  pr <- pr[pr$ranking_method %in% SUP_RANKINGS, ]
  if (nrow(pr) == 0) return(NULL)

  ## Cancer type is not in pred_df; lineage_vec is aligned to rownames(expr_mat).
  lin <- if (!is.null(P$lineage_vec) && !is.null(P$expr_mat))
    setNames(as.character(P$lineage_vec), rownames(P$expr_mat)) else character(0)

  out <- data.frame(
    dataset         = .sup_label(ds_name),
    sample          = as.character(pr$cell_line),
    cancer_type     = if (length(lin) > 0) unname(lin[as.character(pr$cell_line)])
                      else NA_character_,
    validation_fold = as.integer(pr$fold),
    gene_selection  = as.character(pr$ranking_method),
    gene_direction  = as.character(pr$gene_direction),
    n_genes         = as.integer(pr$top_n_genes),
    model           = as.character(pr$model),
    measured_DT     = as.numeric(pr$observed_DT),
    predicted_DT    = as.numeric(pr$predicted_DT),
    stringsAsFactors = FALSE)

  out <- out[order(out$gene_selection, out$gene_direction, out$n_genes,
                   out$model, out$validation_fold, out$sample), ]
  rownames(out) <- NULL
  out
}

pred_CCLE <- .sup_predictions("CCLE_norm")
pred_PDX <- .sup_predictions("PDX_full_sig_3")

for (.nm in c("pred_CCLE", "pred_PDX")) {
  .t <- get(.nm)
  if (is.null(.t)) {
    cat(sprintf("     %s: no predictions available\n", .nm))
  } else {
    cat(sprintf("     %s: %s — %d rows, %d samples, %d models, sizes %s\n",
                .nm, .t$dataset[1], nrow(.t), length(unique(.t$sample)),
                length(unique(.t$model)),
                paste(sort(unique(.t$n_genes)), collapse = "/")))
  }
}
rm(.nm, .t)

################################################################################
# Tables 3.1 and 3.2 — Figures 2a and 2b, per-model performance vs gene count
#
# One row per model x gene-set size. Built from the prediction tables above so
# the per-fold correlations come from the same numbers the figure does.
# The same means are also read from FIG2, by a different route, and compared:
# the console prints max |recomputed - figure| for each table. That comparison
# is NOT a column of the sheet — it rides along as an attribute.
################################################################################

.sup_curve <- function(pred_tab, fig_df) {
  if (is.null(pred_tab)) return(NULL)

  d <- pred_tab[pred_tab$gene_direction == SUP_PANEL_DIR, ]
  if (nrow(d) == 0) return(NULL)

  folds <- sort(unique(d$validation_fold))
  grps  <- split(seq_len(nrow(d)), list(d$n_genes, d$model), drop = TRUE)

  out <- do.call(rbind, lapply(grps, function(idx) {
    g <- d[idx, ]
    g <- g[!duplicated(g$sample), ]     # each sample is held out exactly once

    ## per-fold rho — the five points the plotted mean and SD come from
    rho_f <- vapply(folds, function(f) {
      q  <- g[g$validation_fold == f, ]
      ok <- is.finite(q$measured_DT) & is.finite(q$predicted_DT)
      if (sum(ok) < 3) return(NA_real_)
      suppressWarnings(cor(q$measured_DT[ok], q$predicted_DT[ok],
                           method = "spearman"))
    }, numeric(1))

    ## gene_selection, gene_direction and dataset are not output columns: each
    ## is a single value per table (SUP_RANKINGS fixes the ranking to cor_all,
    ## the filter above fixes the direction to SUP_PANEL_DIR, and there is one
    ## table per dataset). The dataset is named in the Contents sheet; state
    ## the ranking and direction in the table caption.
    row <- data.frame(
      n_genes         = g$n_genes[1],
      model           = g$model[1],
      n_folds         = sum(!is.na(rho_f)),
      n_samples       = sum(is.finite(g$measured_DT) & is.finite(g$predicted_DT)),
      Spearman_mean   = mean(rho_f, na.rm = TRUE),
      Spearman_sd     = sd(rho_f,   na.rm = TRUE),
      stringsAsFactors = FALSE)

    per_fold <- as.data.frame(as.list(rho_f))
    names(per_fold) <- sprintf("rho_fold%d", folds)

    cbind(row, per_fold)
  }))

  ## CROSS-CHECK AGAINST FIG2, matched on model and gene count.
  ##
  ## Not a published column. This check is the one thing that catches a
  ## figure2_data.rds and a pipeline_objects.rds from different runs: the two
  ## numbers reach this point by completely different routes (recomputed from
  ## pred_df here, read from summary_metrics_flat there) and must agree to
  ## rounding. It rides along as an attribute, which saveRDS and the console
  ## report can see but write.xlsx will not put in a cell.
  .fig_mean <- rep(NA_real_, nrow(out))
  if (!is.null(fig_df)) {
    .key <- paste(out$model, out$n_genes, sep = "|")
    .fk  <- paste(as.character(fig_df$model),
                  as.numeric(fig_df$top_n_genes), sep = "|")
    .fig_mean <- fig_df$Spearman_mean[match(.key, .fk)]
  }
  attr(out, "figure_mean") <- .fig_mean

  ## Subsetting drops custom attributes, so reorder the check alongside the
  ## rows and re-attach it afterwards.
  .ord <- order(out$model, out$n_genes)
  .fm  <- attr(out, "figure_mean")[.ord]
  out <- out[.ord, ]
  rownames(out) <- NULL
  attr(out, "figure_mean") <- .fm
  out
}

## PDX first, matching the figure: Figure 2a is PDX, Figure 2b is CCLE.
sup_tab_3_1 <- .sup_curve(pred_PDX,  FIG2$C)   # PDX  -> Figure 2a
sup_tab_3_2 <- .sup_curve(pred_CCLE, FIG2$D)   # CCLE -> Figure 2b

## Report the cross-check. A mismatch means the curve table and the plotted
## panel disagree, which is worth knowing before the supplement goes out.
## The dataset is not a column of these tables, so it is named here.
.sup_curve_tabs <- c("sup_tab_3_1" = "PDX", "sup_tab_3_2" = "CCLE")
for (.i in seq_along(.sup_curve_tabs)) {
  .nm <- names(.sup_curve_tabs)[.i]
  .t <- get(.nm)
  if (is.null(.t)) { cat(sprintf("     %s: not built\n", .nm)); next }
  .dv <- abs(.t$Spearman_mean - attr(.t, "figure_mean"))
  .dv <- .dv[is.finite(.dv)]
  cat(sprintf("     %s: %s — %d rows (%d models x %d sizes), %s\n",
              .nm, unname(.sup_curve_tabs[.i]), nrow(.t),
              length(unique(.t$model)),
              length(unique(.t$n_genes)),
              if (length(.dv) == 0) "no figure values to compare"
              else sprintf("max |recomputed - figure| = %.2e", max(.dv))))
}
rm(list = intersect(ls(), c(".nm", ".t", ".dv", ".i", ".sup_curve_tabs")))

################################################################################
# Tables 3.3 and 3.4 — the two ensemble curves
################################################################################

## ---- Table 3.3  Figure 2c : PDX ensembles ----------------------------------
## From slot E.
## As in 3.1 / 3.2: no `dataset` column, one table per dataset.
sup_tab_3_3 <- if (is.null(FIG2$E)) NULL else data.frame(
  n_genes          = as.numeric(FIG2$E$top_n_genes),
  ensemble         = as.character(FIG2$E$ensemble),
  n_members        = FIG2$E$n_members,
  Spearman_mean    = FIG2$E$Spearman_mean,
  Spearman_sd      = FIG2$E$Spearman_sd,
  stringsAsFactors = FALSE)
if (!is.null(sup_tab_3_3))
  sup_tab_3_3 <- sup_tab_3_3[order(sup_tab_3_3$ensemble, sup_tab_3_3$n_genes), ]

## ---- Table 3.4  Sup Figure 2 panel c : CCLE ensembles ----------------------
## From slot E_CCLE: Supplementary Figure 2c, the companion to Figure 2c.
sup_tab_3_4 <- if (is.null(FIG2$E_CCLE)) NULL else data.frame(
  n_genes          = as.numeric(FIG2$E_CCLE$top_n_genes),
  ensemble         = as.character(FIG2$E_CCLE$ensemble),
  n_members        = FIG2$E_CCLE$n_members,
  Spearman_mean    = FIG2$E_CCLE$Spearman_mean,
  Spearman_sd      = FIG2$E_CCLE$Spearman_sd,
  stringsAsFactors = FALSE)
if (!is.null(sup_tab_3_4))
  sup_tab_3_4 <- sup_tab_3_4[order(sup_tab_3_4$ensemble, sup_tab_3_4$n_genes), ]

## Note on the domain swap (Supplementary Figure 2a/2b): its predictions are on
## the TRAINING dataset's doubling-time scale, so only their rank agreement with
## the measured DT is meaningful, which is what the figure reports.

## ---- Table 3.5  Figure 2f : the 100-gene prediction signature -------------
## The gene set the prediction models are trained on, written out here so the
## workbook is reproducible from this script alone.
##
## Order is as supplied: ranked by |rho|, strongest first. It is NOT sorted
## alphabetically, because the rank is information.
##
## Panel f plots every gene measured in both datasets and highlights the
## top-|rho| overlap set in blue; this list is that signature. If you ever
## change GENE_CORR$min_overlap or the rankings behind it, check the two still
## agree — nothing in the code forces them to.
SIG_100_GENES <- c(
  "CYP2J2", "POLR3G", "HNRNPC", "ENPP5", "TMC4",
  "CCT4", "CCNB1", "ADSL", "IDUA", "YPEL3",
  "NPM1", "MGAT4A", "SAAL1", "SLC40A1", "EPB41L4A",
  "VMAC", "KIF18A", "PGAP3", "PBK", "ZBTB7C",
  "C1QBP", "SELENOP", "CKS2", "VTCN1", "MYH14",
  "MUC20", "PRMT3", "FAAH2", "CRACDL", "MAB21L4",
  "DEPDC1", "ILDR1", "PCLAF", "CDC20", "B3GALT4",
  "INHBB", "PLK1", "RAI2", "SH3YL1", "AURKA",
  "CENPW", "NME1", "SOWAHB", "FAM216A", "OSBPL2",
  "THNSL2", "SPINT1", "CCNA2", "ORMDL3", "PRSS8",
  "PTTG1", "CTNNAL1", "DLGAP5", "PAICS", "EPS8L1",
  "NUF2", "TOMM5", "MYO6", "ARHGEF10L", "RAN",
  "CYP4B1", "DPH5", "NCAPG", "EXOSC8", "KNSTRN",
  "MRTO4", "PIK3IP1", "HAUS1", "AMY2B", "ZBTB7B",
  "KRT23", "SPSB3", "RSAD2", "GSDMB", "PLXNB1",
  "RASEF", "BAG2", "DUSP8", "CCDC159", "PATJ",
  "AVIL", "TNFSF13", "PARPBP", "CENPN", "SUSD6",
  "RAB11FIP4", "DENND2D", "DDR1", "C17orf107", "NEURL1B",
  "MRPL11", "YPEL2", "NEURL3", "MRPL50", "ABCA1",
  "RNF103", "ERP27", "SRSF7", "SNRPB", "DTYMK")

sup_tab_3_5 <- data.frame(gene_symbol = SIG_100_GENES, stringsAsFactors = FALSE)

################################################################################
# Contents
################################################################################

.sup_sheets <- list(
  list(sheet = "Table 3.1", obj = "sup_tab_3_1",
       desc = "PDX: per-model performance vs gene count. The five per-fold Spearman rho values with their mean and SD; Spearman_mean is the value Figures 2a and 2b plot. PDX is Figure 2a."),
  list(sheet = "Table 3.2", obj = "sup_tab_3_2",
       desc = "CCLE: per-model performance vs gene count, same columns as Table 3.1. CCLE is Figure 2b."),
  list(sheet = "Table 3.3", obj = "sup_tab_3_3",
       desc = "PDX: validation Spearman by ensemble and number of genes (mean, SD over 5 folds). This is Figure 2c."),
  list(sheet = "Table 3.4", obj = "sup_tab_3_4",
       desc = "CCLE: validation Spearman by ensemble and number of genes, same columns as Table 3.3. This is panel c of Supplementary Figure 2 — the companion to Figure 2c (Table 3.3)."),
  list(sheet = "Table 3.5", obj = "sup_tab_3_5",
       desc = "The 100 genes used for prediction. Figure 2f."))

## Contents lists ONLY the tables that actually made it into the workbook.
##
## A table can be absent for two reasons: its source was missing, or it was too
## large for a sheet and went out as a gzipped CSV instead. Either way, listing
## it here would point a reader at a tab that is not in the file. The console
## says what happened to each one as the workbook is written, and the CSV
## fallback names itself on disk.
sup_tab_3_contents <- do.call(rbind, lapply(.sup_sheets, function(s) {
  x <- get(s$obj)
  if (is.null(x) || nrow(x) > SUP_XL_MAX_ROWS) return(NULL)
  data.frame(sheet = s$sheet, description = s$desc,
             stringsAsFactors = FALSE)
}))

################################################################################
# Build the workbook
################################################################################

sup_tab_3 <- createWorkbook()

.sup_hdr <- createStyle(textDecoration = "bold", halign = "left",
                        border = "bottom", borderStyle = "thin")

## Contents is added FIRST, so it is the sheet the workbook opens on. The
## table sheets follow in the loop below, in .sup_sheets order.
addWorksheet(sup_tab_3, "Contents")
writeData(sup_tab_3, sheet = "Contents", x = sup_tab_3_contents,
          headerStyle = .sup_hdr)
setColWidths(sup_tab_3, "Contents", cols = 1:2, widths = c(11, 110))
freezePane(sup_tab_3, "Contents", firstRow = TRUE)

for (.s in .sup_sheets) {
  .x <- get(.s$obj)

  if (is.null(.x)) {
    cat(sprintf("     %-10s: SKIPPED — no data\n", .s$sheet))
    next
  }

  ## Oversized tables go to a gzipped CSV instead of a sheet.
  if (nrow(.x) > SUP_XL_MAX_ROWS) {
    .csv <- file.path(NEW_FIG_DIR,
                      sprintf("sup_tab3_%s.csv.gz",
                              gsub("[^A-Za-z0-9]", "", .s$sheet)))
    .con <- gzfile(.csv, "wt")
    write.csv(.x, .con, row.names = FALSE)
    close(.con)
    cat(sprintf("     %-10s: %d rows exceeds the Excel limit -> %s\n",
                .s$sheet, nrow(.x), basename(.csv)))
    next
  }

  addWorksheet(sup_tab_3, .s$sheet)
  writeData(sup_tab_3, sheet = .s$sheet, x = .x, headerStyle = .sup_hdr)
  ## "auto" widths measure every cell, so only do it on short sheets.
  if (nrow(.x) <= 50000)
    setColWidths(sup_tab_3, .s$sheet, cols = seq_len(ncol(.x)), widths = "auto")
  freezePane(sup_tab_3, .s$sheet, firstRow = TRUE)
  cat(sprintf("     %-10s: %d rows x %d columns\n",
              .s$sheet, nrow(.x), ncol(.x)))
}

dir.create(dirname(SUP_TAB_FILE), recursive = TRUE, showWarnings = FALSE)
saveWorkbook(sup_tab_3, SUP_TAB_FILE, overwrite = TRUE)
cat(sprintf("\n     Saved %s (%.0f KB)\n",
            basename(SUP_TAB_FILE), file.size(SUP_TAB_FILE) / 1024))

rm(list = intersect(ls(), c(".sup_label", ".sup_sheets",
                            ".sup_hdr", ".sup_get_pipe", ".sup_predictions",
                            ".sup_curve", ".sup_ens_predictions", ".m", ".s",
                            ".x", ".csv", ".con", ".hit", ".f", ".key", ".fk")))


cat("\n================================================================\n")
cat(" FIGURES COMPLETE\n")
cat(sprintf(" New Figure 2 (a-f) + panels: %s\n", NEW_FIG_DIR))
if (exists("SUP_FIG2_TAG"))
  cat(sprintf(" Supplementary Figure 2: %s.pdf (+ panels a, b, c) in the same folder\n",
              SUP_FIG2_TAG))
if (exists("GENE_CORR_RDS") && file.exists(GENE_CORR_RDS))
  cat(sprintf(" Panel f gene-correlation cache: %s\n", basename(GENE_CORR_RDS)))
if (exists("SUP_TAB_FILE") && file.exists(SUP_TAB_FILE))
  cat(sprintf(" Supplementary table:   %s\n", SUP_TAB_FILE))
cat("================================================================\n")

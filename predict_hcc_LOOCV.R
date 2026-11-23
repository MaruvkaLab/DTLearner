################################################################################
#  predict_hcc_LOOCV.R
#  Leave-one-out cross-validation (LOOCV) of the pDT models within the HCC
#  cohort (GSE54236). Output is used by Supplemental Figure 3 (panel a) in
#  plot_figure_3.R.
#
#  ── INPUTS (inside WORK_DIR) ──────────────────────────────────────────────
#    input_files/hcc/GSE54236_qnorm.tsv                    expression
#    input_files/hcc/GSE54236_clinical_metadata_compact.tsv measured DT
#
#  ── WHAT IT DOES ──────────────────────────────────────────────────────────
#  Uses ONLY the HCC samples (no PDX training data). For each HCC tumour with
#  a measured doubling time, the three models (Ridge, RandomForest,
#  SVR_radial) are trained on the other N-1 tumours using the 100-gene
#  signature (Section 1D) and predict the held-out tumour. Ensemble_Fixed is
#  the mean of the three.
#
#  Normalisation (NORM_MODES = "zscore"): each fold's N-1 training samples are
#  z-scored on their own mean/sd, and the held-out sample is standardised with
#  those same parameters.
#
#  OUTPUTS (in OUT_DIR):
#    pDT_loocv_testonly_predictions.csv — per (norm mode, model, sample)
#    pDT_loocv_testonly_metrics.csv     — Pearson r / Spearman rho / RMSE /
#                                          MAE / R² per (norm mode, model)
#    pDT_loocv_testonly_genes_used.csv  — which of your requested genes were
#                                          matched in the expression matrix
#
# ── REQUIRED PACKAGES ────────────────────────────────────────────────────────
#  glmnet  ranger  e1071  future  future.apply  (plus tidyr, used via ::)
#  Missing packages are installed automatically (Section 2).
################################################################################


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 1 — USER CONFIG    ← EDIT ONLY THIS SECTION                   ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

# ── WORKING DIRECTORY — defaults to the folder this script is in ──────────
#
#  Every input is read from, and every output written to, a folder inside
#  WORK_DIR — the same WORK_DIR as all the other scripts. Only the two HCC
#  files below (WORK_DIR/input_files/hcc/) are read.
#
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
  stop(sprintf("WORK_DIR does not exist: %s — set WORK_DIR in Section 1.", WORK_DIR))

# ── 1A  HCC samples: RNA expression TSV ──────────────────────────────────────
#  Genes × samples (genes as rows); the same file predict_hcc.R uses.
NEW_RNA_FILE <- file.path(WORK_DIR, "input_files", "hcc", "GSE54236_qnorm.tsv")

# ── 1B  HCC samples: measured doubling time ──────────────────────────────────
#  Only tumour samples with a usable doubling time are used — both for
#  training each fold AND for being predicted/scored.
NEW_META_FILE <- file.path(WORK_DIR, "input_files", "hcc", "GSE54236_clinical_metadata_compact.tsv")

# ── 1C  Output directory ──────────────────────────────────────────────────────
OUT_DIR <- file.path(WORK_DIR, "loocv_hcc")

# ── 1D  Gene panel ───────────────────────────────────────────────────────────
#
#  Gene names are matched against the columns of NEW_RNA_FILE: exact match
#  first, then a whitespace-trimmed / case-insensitive fallback. Anything
#  that still can't be found is reported and skipped. (GENE_LIST_FILE can
#  instead point to a text/CSV file of gene names; if both are set, the two
#  are combined.)
#
#  The 100-gene prediction signature: the same list, in the same order
#  (ranked by |rho|, strongest first), as SIG_100_GENES in the Figure 2
#  scripts and Supplemental Table 3.5. Keep the three in step if it changes.
GENE_LIST      <- c(
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
GENE_LIST_FILE <- NULL

#  Minimum usable genes required to proceed. Ridge needs >= 2 columns, so
#  don't set this below 2.
MIN_GENES   <- 2L

#  Per-fold normalisation: "zscore" (fold-wise z-scoring, see the header) or
#  "none" (raw values). Several modes can be listed; each is run separately.
NORM_MODES  <- c("zscore")

# ── 1E  Fixed ensemble composition ───────────────────────────────────────────
FIXED_ENSEMBLE_MODELS <- c("Ridge", "RandomForest", "SVR_radial")

# ── 1F  Minimum labeled samples required to proceed ──────────────────────────
MIN_SAMPLES <- 10L   # hard stop below this; a strong warning is issued below 20

# ── 1G  Parallel workers (NULL = all logical cores - 1) ──────────────────────
N_WORKERS <- NULL


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2 — PACKAGES                                                   ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

# ── Install missing packages ───────────────────────────────────────────────────
# Installs only the packages that are not already installed; packages that are
# present are left as they are (no updates). All are on CRAN.
CRAN_PKGS <- c("glmnet", "ranger", "e1071", "future", "future.apply", "tidyr")

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

suppressPackageStartupMessages({
  library(glmnet); library(ranger); library(e1071)
  library(future);  library(future.apply)
})


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 3 — HELPER FUNCTIONS                                           ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

# ── Read a user gene list from a plain text / CSV / TSV file ──────────────────
#  Accepts one gene per line, or a single-column file, or a multi-column file
#  in which case the first column is used. A header row whose first field
#  looks like a column label ("gene", "genes", "symbol", ...) is skipped.
read_gene_list_file <- function(path) {
  if (is.null(path) || !nzchar(path)) return(character(0))
  if (!file.exists(path))
    stop("GENE_LIST_FILE not found: ", path)

  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (!length(lines)) stop("GENE_LIST_FILE is empty: ", path)

  # take the first field of each line (handles CSV / TSV with extra columns)
  first_field <- trimws(sub("[,\t;].*$", "", lines))
  first_field <- gsub('^"|"$', "", first_field)

  header_words <- c("gene", "genes", "gene_symbol", "genesymbol", "symbol",
                    "gene_id", "geneid", "feature", "id", "x")
  if (tolower(first_field[1L]) %in% header_words) first_field <- first_field[-1L]

  unique(first_field[nzchar(first_field)])
}

# ── Match requested gene names onto the columns of the expression matrix ─────
#  Exact match first; then a trimmed / case-insensitive fallback for the
#  leftovers. Returns the matched column names (in the order requested) plus
#  a per-gene report of how each one resolved.
match_genes_to_matrix <- function(requested, available) {
  requested <- unique(trimws(as.character(requested)))
  requested <- requested[nzchar(requested)]

  key_avail <- toupper(trimws(available))
  lookup    <- available[!duplicated(key_avail)]
  names(lookup) <- key_avail[!duplicated(key_avail)]

  matched <- rep(NA_character_, length(requested))
  how     <- rep("not_found", length(requested))

  exact <- requested %in% available
  matched[exact] <- requested[exact]
  how[exact]     <- "exact"

  fuzzy <- !exact & (toupper(requested) %in% names(lookup))
  matched[fuzzy] <- unname(lookup[toupper(requested[fuzzy])])
  how[fuzzy]     <- "case_insensitive"

  report <- data.frame(requested_gene = requested,
                       matched_column = matched,
                       match_type     = how,
                       stringsAsFactors = FALSE)

  list(genes  = unique(matched[!is.na(matched)]),
       report = report)
}

# ── z-score helpers ────────────────────────────────────────────────────────────
#  store_zparams: per-gene mean and SD (SD floored at 1e-8).
#  apply_zparams: z-score a matrix with stored parameters.
store_zparams <- function(mat) {
  mat <- as.matrix(mat); storage.mode(mat) <- "numeric"
  list(means = colMeans(mat, na.rm = TRUE),
       sds   = pmax(apply(mat, 2, sd, na.rm = TRUE), 1e-8))
}

apply_zparams <- function(mat, zp) {
  mat <- as.matrix(mat[, names(zp$means), drop = FALSE])
  storage.mode(mat) <- "numeric"
  sweep(sweep(mat, 2, zp$means, "-"), 2, zp$sds, "/")
}

# ── Replace non-finite values with the per-gene median in `medians` ────────────
impute_medians <- function(mat, medians) {
  bad <- !is.finite(mat)
  if (any(bad)) {
    idx <- which(bad, arr.ind = TRUE)
    mat[bad] <- medians[colnames(mat)[idx[, 2L]]]
  }
  mat
}

# ── Fit the three ensemble members on ALREADY-PREPARED data ────────────────────
#  Unlike predict_*.R, the caller (run_one_fold) normalises the data; this only
#  imputes and fits. Hyperparameters match predict_*.R: Ridge (cv.glmnet,
#  alpha = 0, lambda.min), RandomForest (1000 trees, mtry = sqrt(p),
#  min.node.size = 5), SVR_radial (tune.svm over cost x gamma, 3-fold CV).
#  A model that errors is returned as NULL.
fit_all_models <- function(x_train, y_train, genes, seed_base = 0L) {
  nf   <- max(3L, min(10L, floor(length(y_train) / 5L)))
  xtr  <- as.matrix(x_train[, genes, drop = FALSE])
  storage.mode(xtr) <- "numeric"

  med  <- setNames(
    apply(xtr, 2, function(v) { v <- v[is.finite(v)]; if (length(v)) median(v) else 0 }),
    genes)
  xtr  <- impute_medians(xtr, med)

  safe <- function(expr) tryCatch(expr, error = function(e) NULL)

  fits <- list()

  # The 3 models that make up FIXED_ENSEMBLE_MODELS.
  fits$Ridge <- safe({
    list(cv = cv.glmnet(xtr, y_train, alpha = 0, nfolds = nf, standardize = FALSE))
  })

  set.seed(1000L + seed_base)
  fits$RandomForest <- safe({
    df <- data.frame(y = y_train, xtr, check.names = FALSE)
    list(model = ranger(y ~ ., data = df,
                        num.trees     = 1000L,
                        mtry          = max(1L, floor(sqrt(ncol(xtr)))),
                        min.node.size = 5L,
                        importance    = "impurity",
                        seed          = 1000L + seed_base))
  })

  set.seed(3000L + seed_base)
  fits$SVR_radial <- safe({
    tuned <- tune.svm(xtr, y_train, kernel = "radial",
                      cost  = c(0.1, 1, 10),
                      gamma = c(0.01, 0.1, 1 / ncol(xtr)),
                      epsilon = 0.1,
                      tunecontrol = tune.control(cross = 3))
    list(model = tuned$best.model)
  })

  list(fits = fits, genes = genes, medians = med)
}

# ── Predict prepared new data with each fitted model ───────────────────────────
predict_bundle <- function(bundle, x_new_prepared) {
  x   <- as.matrix(x_new_prepared[, bundle$genes, drop = FALSE])
  storage.mode(x) <- "numeric"
  x   <- impute_medians(x, bundle$medians)
  df  <- data.frame(x, check.names = FALSE)
  safe <- function(expr) tryCatch(expr, error = function(e) NULL)
  preds <- list()
  for (nm in names(bundle$fits)) {
    obj <- bundle$fits[[nm]]
    if (is.null(obj)) { preds[[nm]] <- NULL; next }
    preds[[nm]] <- safe({
      if (nm == "Ridge") {
        as.numeric(predict(obj$cv, x, s = "lambda.min"))
      } else if (nm == "RandomForest") {
        predict(obj$model, df)$predictions
      } else if (nm == "SVR_radial") {
        as.numeric(predict(obj$model, df))
      }
    })
    if (!is.null(preds[[nm]])) names(preds[[nm]]) <- rownames(x_new_prepared)
  }
  preds
}

# ── Pearson r, Spearman rho, RMSE, MAE and R^2 of predicted vs actual ──────────
#  (NA for fewer than 3 finite pairs)
compute_metrics <- function(predicted, actual) {
  ok <- is.finite(predicted) & is.finite(actual)
  predicted <- predicted[ok]; actual <- actual[ok]
  n <- length(predicted)
  if (n < 3L) {
    return(data.frame(n = n, pearson_r = NA_real_, spearman_rho = NA_real_,
                      rmse = NA_real_, mae = NA_real_, r_squared = NA_real_))
  }
  resid <- predicted - actual
  ss_res <- sum(resid^2)
  ss_tot <- sum((actual - mean(actual))^2)
  data.frame(
    n            = n,
    pearson_r    = suppressWarnings(cor(predicted, actual, method = "pearson")),
    spearman_rho = suppressWarnings(cor(predicted, actual, method = "spearman")),
    rmse         = sqrt(mean(resid^2)),
    mae          = mean(abs(resid)),
    r_squared    = if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  )
}


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 4 — LOAD DATA                                                  ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

message("\n>>> Validating inputs...")
for (.f in c(NEW_RNA_FILE, NEW_META_FILE)) {
  if (!is.character(.f) || length(.f) != 1L || nchar(.f) == 0L || !file.exists(.f))
    stop("File not found or invalid path: ", .f,
         "\n  Update the paths in Section 1 (USER CONFIG) and re-source.")
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("\n>>> Loading data...")

new_rna_raw <- as.data.frame(t(read.table(NEW_RNA_FILE, sep = "\t",
                                           header = TRUE, row.names = 1,
                                           check.names = FALSE)))
new_mat <- as.matrix(new_rna_raw); storage.mode(new_mat) <- "numeric"
rm(new_rna_raw)
message(sprintf("    Test samples: %d  |  genes provided: %d",
                nrow(new_mat), ncol(new_mat)))

# ── Resolve YOUR gene panel against the expression matrix ────────────────────
requested_genes <- unique(c(
  if (!is.null(GENE_LIST))      trimws(as.character(GENE_LIST)) else character(0),
  if (!is.null(GENE_LIST_FILE)) read_gene_list_file(GENE_LIST_FILE) else character(0)
))
requested_genes <- requested_genes[nzchar(requested_genes)]

if (!length(requested_genes))
  stop("No gene panel supplied. Set GENE_LIST and/or GENE_LIST_FILE in Section 1D ",
       "(USER CONFIG) — this script does not select genes on its own.")

.m          <- match_genes_to_matrix(requested_genes, colnames(new_mat))
GENE_SET    <- .m$genes
GENE_REPORT <- .m$report
.missing    <- .m$report$requested_gene[.m$report$match_type == "not_found"]
.fuzzy      <- sum(.m$report$match_type == "case_insensitive")

message(sprintf("    Gene panel: %d requested, %d matched%s",
                nrow(.m$report), length(GENE_SET),
                if (.fuzzy > 0L) sprintf(" (%d by case-insensitive match)", .fuzzy) else ""))
if (length(.missing)) {
  message(sprintf("    %d requested gene(s) not present in NEW_RNA_FILE and will be skipped:",
                  length(.missing)))
  message("      ", paste(utils::head(.missing, 25L), collapse = ", "),
          if (length(.missing) > 25L) sprintf(" ... (+%d more)", length(.missing) - 25L) else "")
}
if (length(GENE_SET) < MIN_GENES)
  stop(sprintf("Only %d of your %d requested gene(s) were found in NEW_RNA_FILE — need at least %d. ",
               length(GENE_SET), nrow(.m$report), MIN_GENES),
       "Check gene naming (symbols vs Ensembl IDs vs probe IDs) against the row names of NEW_RNA_FILE.")
rm(.m, .missing, .fuzzy)

# ── Measured doubling time for the HCC tumours ────────────────────────────────
#  Same parsing and exclusions as plot_figure_1.R / plot_figure_3.R: split
#  Sample_title into type / patient / biopsy, keep tumours ("T"), and drop
#  biopsies 047, 127 and 107.
new_meta <- subset(read.table(NEW_META_FILE,
                              sep = "\t",
                              header = TRUE,
                              row = 1)[,c("Sample_title","doubling.time..days.")],
                   complete.cases(doubling.time..days.))

new_meta <- rbind(tidyr::separate(new_meta[c(1:81),],col = (Sample_title), sep = " ", into = c("type", "sample"), remove = FALSE),
                  tidyr::separate(new_meta[c(82:161),],col = (Sample_title), sep = " ", into = c("type", "patient","sample"), remove = FALSE)[,-c(3)])[,-c(2)]

new_meta <- tidyr::separate(new_meta,col = (sample), sep = "_", into = c("type", "patient","biopsy"), remove = FALSE)

new_meta <- subset(new_meta, (type == "T") & !(biopsy %in% c("047","127","107")))[,-c(2,3,4,5)]

# Keep the title and DT, renamed to the column names used below
colnames(new_meta) <- c("Sample_title","dt_days")


if (!"dt_days" %in% colnames(new_meta))
  stop("NEW_META_FILE must contain a 'dt_days' column.")

common_new <- intersect(rownames(new_mat), rownames(new_meta))
if (length(common_new) == 0L)
  stop("No row names shared between NEW_RNA_FILE and NEW_META_FILE.")

# DT for the samples present in both files; these are the LOOCV samples
.new_dt_raw <- setNames(as.numeric(new_meta[common_new, "dt_days"]), common_new)
NEW_DT <- .new_dt_raw[is.finite(.new_dt_raw)]
message(sprintf("    Test samples with known dt_days: %d / %d",
                length(NEW_DT), nrow(new_mat)))

if (length(NEW_DT) < MIN_SAMPLES)
  stop(sprintf(
    "Only %d test sample(s) have a known dt_days — need at least %d (MIN_SAMPLES) ",
    length(NEW_DT), MIN_SAMPLES), "to run LOOCV meaningfully.")
if (length(NEW_DT) < 20L)
  warning(sprintf(
    "Only %d labeled test samples for a %d-gene panel. Your gene set is fixed ",
    length(NEW_DT), length(GENE_SET)),
    "(so there's no selection instability), but the models themselves are ",
    "being fit on very few samples — treat results as directional, not definitive.")

rm(new_meta, common_new, .new_dt_raw)
eval_samples <- names(NEW_DT)


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 5 — PARALLELISATION SETUP                                      ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

message("\n>>> Setting up parallel workers...")

if (is.null(N_WORKERS))
  N_WORKERS <- max(1L, parallel::detectCores(logical = TRUE) - 1L)

.mem <- max(512 * 1024^2, min(2 * 1024^3, floor(12 * 1024^3 / N_WORKERS)))
options(future.globals.maxSize = .mem)
message(sprintf("    Workers     : %d  (detected %d logical cores)",
                N_WORKERS, parallel::detectCores(logical = TRUE)))
rm(.mem)

plan(multisession, workers = N_WORKERS)


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 6 — LOOCV: normalization + training + prediction redone per     ║
# ║              held-out sample, per normalization mode (fixed gene panel)  ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

message(sprintf("\n>>> Running LOOCV: %d sample(s) × %d norm mode(s) = %d fold(s)...",
                length(eval_samples), length(NORM_MODES),
                length(eval_samples) * length(NORM_MODES)))
message(sprintf("    (fixed %d-gene panel; 3-model training [Ridge, RandomForest, SVR_radial] repeated for every fold)",
                length(GENE_SET)))

# One job per held-out sample x normalisation mode
fold_jobs <- list()
for (s in eval_samples) {
  for (nm_mode in NORM_MODES) {
    fold_jobs[[paste0(s, "__", nm_mode)]] <- list(held_out = s, mode = nm_mode)
  }
}

# One LOOCV fold: train on every labelled sample except `held_out`, predict it.
# Everything fitted from data (variance filter, medians, z-score parameters,
# models) uses the N-1 training samples only.
run_one_fold <- function(job, new_mat, NEW_DT, eval_samples,
                         GENE_SET, MIN_GENES = 2L) {
  s        <- job$held_out
  mode     <- job$mode
  fold_ids <- setdiff(eval_samples, s)

  fold_expr <- new_mat[fold_ids, , drop = FALSE]
  fold_dt   <- NEW_DT[fold_ids]

  # Your panel is used as-is in every fold — no selection happens here. The
  # only thing dropped is genes that are flat within this fold's N-1 samples
  # (they carry no information and can destabilise the fits). That filter is
  # unsupervised and excludes the held-out sample, so it leaks nothing.
  present  <- intersect(GENE_SET, colnames(fold_expr))
  if (!length(present)) return(NULL)
  var_fold <- apply(fold_expr[, present, drop = FALSE], 2, sd, na.rm = TRUE)
  selected <- present[!is.na(var_fold) & var_fold > 1e-5]
  if (length(selected) < MIN_GENES) return(NULL)

  # Impute missing values in both the training rows and the held-out row with
  # the training rows' medians
  xtr_raw     <- fold_expr[, selected, drop = FALSE]
  fold_med    <- setNames(
    apply(xtr_raw, 2, function(v) { v <- v[is.finite(v)]; if (length(v)) median(v) else 0 }),
    selected)
  xtr_raw     <- impute_medians(xtr_raw, fold_med)

  one_row_raw <- new_mat[s, selected, drop = FALSE]
  one_row_raw <- impute_medians(one_row_raw, fold_med)

  if (mode == "none") {
    xtr <- xtr_raw; one_row <- one_row_raw
  } else if (mode == "zscore") {
    zp      <- store_zparams(xtr_raw)
    xtr     <- apply_zparams(xtr_raw, zp)
    one_row <- apply_zparams(one_row_raw, zp)
  }

  bnd   <- fit_all_models(xtr, fold_dt, selected, seed_base = 0L)
  preds <- predict_bundle(bnd, one_row)

  list(sample = s, mode = mode, n_genes = length(selected), preds = preds)
}

# Run all folds in parallel. future.seed = TRUE gives each fold its own
# reproducible random stream.
fold_results <- future_lapply(
  fold_jobs,
  run_one_fold,
  new_mat = new_mat, NEW_DT = NEW_DT, eval_samples = eval_samples,
  GENE_SET = GENE_SET, MIN_GENES = MIN_GENES,
  future.seed     = TRUE,
  future.packages = c("glmnet","ranger","e1071")
)
fold_results <- Filter(Negate(is.null), fold_results)
message(sprintf("    Completed %d / %d fold(s).", length(fold_results), length(fold_jobs)))


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 7 — AGGREGATE PER-FOLD PREDICTIONS                            ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

message("\n>>> Aggregating fold predictions...")

# Per normalisation mode: one row per held-out sample, one column per model,
# plus Ensemble_Fixed = mean of the FIXED_ENSEMBLE_MODELS predictions.
by_mode <- split(fold_results, sapply(fold_results, `[[`, "mode"))

all_loocv_preds <- list()
for (nm_mode in names(by_mode)) {
  mode_folds <- by_mode[[nm_mode]]
  all_models <- unique(unlist(lapply(mode_folds, function(f) names(f$preds))))

  pred_mat <- matrix(NA_real_, nrow = length(mode_folds), ncol = length(all_models),
                     dimnames = list(sapply(mode_folds, `[[`, "sample"), all_models))
  for (f in mode_folds) {
    for (mn in names(f$preds)) {
      if (!is.null(f$preds[[mn]])) pred_mat[f$sample, mn] <- f$preds[[mn]][1]
    }
  }

  fixed_avail <- intersect(FIXED_ENSEMBLE_MODELS, colnames(pred_mat))
  if (length(fixed_avail) >= 2L) {
    pred_mat <- cbind(pred_mat,
                      Ensemble_Fixed = rowMeans(pred_mat[, fixed_avail, drop = FALSE],
                                                 na.rm = TRUE))
  }

  all_loocv_preds[[nm_mode]] <- do.call(rbind, lapply(colnames(pred_mat), function(mn) {
    data.frame(
      norm_mode = nm_mode, model = mn, sample_id = rownames(pred_mat),
      predicted = pred_mat[, mn], actual = NEW_DT[rownames(pred_mat)],
      row.names = NULL, stringsAsFactors = FALSE)
  }))
}


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 8 — SAVE OUTPUTS                                              ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

message("\n>>> Saving outputs to: ", OUT_DIR)

loocv_preds_all <- do.call(rbind, all_loocv_preds)
rownames(loocv_preds_all) <- NULL
write.csv(loocv_preds_all,
          file.path(OUT_DIR, "pDT_loocv_testonly_predictions.csv"),
          row.names = FALSE, quote = FALSE)

# Metrics per mode x model over all held-out predictions, best Spearman first
loocv_metrics_all <- do.call(rbind, lapply(
  split(loocv_preds_all, list(loocv_preds_all$norm_mode, loocv_preds_all$model), drop = TRUE),
  function(d) {
    m <- compute_metrics(d$predicted, d$actual)
    cbind(norm_mode = d$norm_mode[1L], model = d$model[1L], m)
  }))
rownames(loocv_metrics_all) <- NULL
loocv_metrics_all <- loocv_metrics_all[order(loocv_metrics_all$norm_mode,
                                             -loocv_metrics_all$spearman_rho), ]
write.csv(loocv_metrics_all,
          file.path(OUT_DIR, "pDT_loocv_testonly_metrics.csv"),
          row.names = FALSE, quote = FALSE)

message("    pDT_loocv_testonly_predictions.csv")
message("    pDT_loocv_testonly_metrics.csv")

write.csv(GENE_REPORT,
          file.path(OUT_DIR, "pDT_loocv_testonly_genes_used.csv"),
          row.names = FALSE, quote = FALSE)
message("    pDT_loocv_testonly_genes_used.csv")


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 9 — SUMMARY + RELEASE WORKERS                                 ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

plan(sequential)

message("\n═══════════════════════════════════════════════════════════")
message(">>> Done.")
message(sprintf("    Labeled test samples : %d", length(eval_samples)))
.ng <- unique(vapply(fold_results, function(f) as.integer(f$n_genes), integer(1L)))
message(sprintf("    Genes                : your panel, %d matched (%s used per fold)",
                length(GENE_SET),
                if (length(.ng) == 0L) "no folds completed"
                else if (length(.ng) == 1L) as.character(.ng)
                else sprintf("%d-%d after per-fold variance filter", min(.ng), max(.ng))))
rm(.ng)
message(sprintf("    Norm modes           : %s", paste(names(by_mode), collapse = ", ")))
message(sprintf("    Output               : %s/pDT_loocv_testonly_metrics.csv", OUT_DIR))
message("═══════════════════════════════════════════════════════════")

################################################################################
# predict_ccle.R
# Predict doubling time (pDT) for CCLE cell lines without a measured DT.
#
# ── INPUTS (inside WORK_DIR) ─────────────────────────────────────────────────
#    input_files/ccle/DT_Expression_CCLE_Data.csv  expression + Avg_2DT
#                                                   (cell lines in rows)
#    input_files/ccle/Models_DT.csv                 OncotreeLineage
#
# ── QUICK START ──────────────────────────────────────────────────────────────
#  1. Put this script in WORK_DIR (next to input_files/), or set WORK_DIR
#     in SECTION 1 (USER CONFIG).
#  2. Source the file:  Ctrl / Cmd + Shift + S
#  3. A data frame called  pDT_results  appears in your environment.
#     The same table is also saved to OUT_DIR as CSV and Excel files.
#
# ── WHAT IT DOES ─────────────────────────────────────────────────────────────
#  Trains on the CCLE cell lines that have a measured doubling time (Avg_2DT)
#  and predicts pDT for the CCLE cell lines that do not. Output is used by
#  plot_figure_5.R.
#
#  Three models (Ridge, RandomForest, SVR_radial) are trained on the 100-gene
#  signature (SECTION 1E); pDT_ensemble is the mean of the three.
#
# ── REQUIRED PACKAGES ────────────────────────────────────────────────────────
#  future  future.apply  glmnet  ranger  e1071  openxlsx
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
#  WORK_DIR — the same WORK_DIR as all the other scripts. The input files are
#  listed at the top of this script.
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

# ── 1A  New (prediction) samples: CCLE cell lines without a DT ────────────────
#
#  The CCLE file holds expression and Avg_2DT for all cell lines (cell lines
#  in rows). Cell lines WITHOUT a measured Avg_2DT are predicted; those with
#  one form the training set (1B/1C). Values are log2-normalised expression.
#

NEW_RNA_FILE <- file.path(WORK_DIR, "input_files", "ccle", "DT_Expression_CCLE_Data.csv")

# ── 1B  Training data: CCLE cell lines with a measured DT ─────────────────────
#
#  Same file as NEW_RNA_FILE: the cell lines WITH a measured Avg_2DT.
#

TRAIN_RNA_FILE <- file.path(WORK_DIR, "input_files", "ccle", "DT_Expression_CCLE_Data.csv")

# ── 1C  Training data: metadata CSV ──────────────────────────────────────────
#
#  Rows = cell lines, matched to the expression file by cell-line ID.
#  Required columns:
#    Avg_2DT      — numeric doubling time in days (used as the response variable)
#    OncotreeLineage — tumour-type / lineage label
#  Column 1 is the row/sample ID and will be used to align with the RNA CSV.
#

TRAIN_META_FILE <- file.path(WORK_DIR, "input_files", "ccle", "Models_DT.csv")

# ── 1D  Output directory ──────────────────────────────────────────────────────
OUT_DIR <- file.path(WORK_DIR, "ccle_pred")

# ── 1E  Gene list ─────────────────────────────────────────────────────────────
#
#  The gene list the models are trained on (HGNC symbols), and the model name
#  used as the prefix of the output columns (pDT_Ridge, ..., pDT_ensemble).
#  plot_figure_3/4/5.R read the pDT_ensemble column, so keep MODEL_NAME "pDT".
#
#  Genes that are absent from the training/new-data overlap (or fail the
#  variance filter in Section 7) are dropped automatically, with a message
#  telling you which ones were dropped — you do not need to pre-filter.
#
#  The 100-gene prediction signature: the same list, in the same order
#  (ranked by |rho|, strongest first), as SIG_100_GENES in the Figure 2
#  scripts and Supplemental Table 3.5. Keep the three in step if it changes.
MODEL_NAME <- "pDT"
GENE_LIST  <- c(
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


# Internal: gene-set run definition built from the list above
# (do not edit this block — edit MODEL_NAME / GENE_LIST instead)
GENE_SET_RUNS <- setNames(
  list(list(genes = GENE_LIST)),
  MODEL_NAME
)

if (!is.character(MODEL_NAME) || length(MODEL_NAME) != 1L || is.na(MODEL_NAME) ||
    nchar(MODEL_NAME) == 0L)
  stop("MODEL_NAME must be a single non-empty character string.")

if (length(GENE_LIST) == 0L || all(is.na(GENE_LIST)))
  stop("GENE_LIST is empty — check the gene list in Section 1E.")

# Ensemble composition — these are also the only models that get fitted
ENSEMBLE_MODELS <- c("Ridge", "RandomForest", "SVR_radial")

# Training strategy. "train_set_Transfer" (fit on the training set, predict
# the new samples) is the only one; the strategy and gene-set loops below
# therefore run once.
STRATEGIES <- c("train_set_Transfer")


# ── 1F  Parallelisation ──────────────────────────────────────────────────────
#  NULL → auto: all logical cores − 1, capped at 11
#  1    → sequential — use this first to confirm the script runs without errors
N_WORKERS <- NULL


# ── 1G  Global seed ──────────────────────────────────────────────────────────
SEED <- 42L


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2 — LIBRARIES                                                  ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

# ── Install missing packages ───────────────────────────────────────────────────
# Installs only the packages that are not already installed; packages that are
# present are left as they are (no updates). All are on CRAN.
CRAN_PKGS <- c("future", "future.apply", "glmnet", "ranger", "e1071",
               "openxlsx")

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

library(future)
library(future.apply)
library(glmnet)
library(ranger)
library(e1071)
library(openxlsx)

message(">>> predict_ccle.R")
message(sprintf("    R version    : %s", R.version.string))
message(sprintf("    Architecture : %s", Sys.info()["machine"]))
message(sprintf("    OS           : %s %s",
                Sys.info()["sysname"], Sys.info()["release"]))
for (.p in c("future","future.apply","glmnet","ranger",
             "e1071","openxlsx"))
  message(sprintf("    %-14s : %s", .p, as.character(packageVersion(.p))))
rm(.p)


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 3 — HELPER FUNCTIONS                                           ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

# ── 3A  z-score helpers ───────────────────────────────────────────────────────
#  store_zparams: per-gene mean and SD (SD floored at 1e-8 so a flat gene does
#  not divide by zero). apply_zparams: z-score a matrix with stored parameters,
#  reordering its columns to match them.
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

# ── 3B  Impute non-finite values with column medians ─────────────────────────
#  Replaces NA / NaN / Inf in each column with that gene's entry in `medians`
#  (0 if the gene has no median).
impute_medians <- function(mat, medians) {
  mat <- as.matrix(mat); storage.mode(mat) <- "numeric"
  for (j in seq_len(ncol(mat))) {
    bad <- !is.finite(mat[, j])
    if (any(bad)) {
      fill <- medians[colnames(mat)[j]]
      mat[bad, j] <- if (!is.na(fill)) fill else 0
    }
  }
  mat
}

# ── 3C  Align new data to training gene set ───────────────────────────────────
#  Returns the new samples x training genes, in training-gene order. Genes the
#  new data lacks are filled with the training median (and reported).
align_to_genes <- function(new_mat, genes, medians) {
  new_mat <- as.matrix(new_mat); storage.mode(new_mat) <- "numeric"
  out <- matrix(NA_real_, nrow(new_mat), length(genes),
                dimnames = list(rownames(new_mat), genes))
  present <- intersect(genes, colnames(new_mat))
  absent  <- setdiff(genes, colnames(new_mat))
  if (length(present)) out[, present] <- new_mat[, present, drop = FALSE]
  if (length(absent)) {
    message(sprintf(
      "    [align] %d / %d training genes absent from new data — imputed with training median.",
      length(absent), length(genes)))
    for (g in absent) out[, g] <- medians[g]
  }
  out
}

# ── 3D  Fit the 3 ensemble models ─────────────────────────────────────────────
#  Only Ridge, RandomForest and SVR_radial are fitted — these are the three
#  members of the ensemble (see ENSEMBLE_MODELS in Section 1E).
#
#  Input is raw-scale expression. It is z-scored on the training samples, then
#  non-finite values are set to the per-gene median of the z-scored data.
#  The bundle returned holds the fits plus the z-score parameters and medians,
#  so predict_bundle() can prepare new data the same way.
#
#  Hyperparameters:
#    Ridge         cv.glmnet, alpha = 0, lambda.min; nf CV folds
#                  (n / 5, clamped to 3-10)
#    RandomForest  ranger, 1000 trees, mtry = sqrt(p), min.node.size = 5
#    SVR_radial    tune.svm over cost {0.1, 1, 10} x gamma {0.01, 0.1, 1/p},
#                  epsilon 0.1, 3-fold CV
#  A model that errors is skipped (NULL) with a message.
fit_all_models <- function(x_train_raw, y_train, genes, seed_base = 0L) {
  nf   <- max(3L, min(10L, floor(length(y_train) / 5L)))
  xtr  <- as.matrix(x_train_raw[, genes, drop = FALSE])
  storage.mode(xtr) <- "numeric"
  zp   <- store_zparams(xtr)
  xtr  <- apply_zparams(xtr, zp)
  med  <- setNames(
    apply(xtr, 2, function(v) { v <- v[is.finite(v)]; if (length(v)) median(v) else 0 }),
    genes)
  xtr  <- impute_medians(xtr, med)

  safe <- function(expr) tryCatch(expr, error = function(e) {
    message("    [fit] SKIPPED: ", conditionMessage(e)); NULL })

  fits <- list()

  # Ridge
  fits$Ridge <- safe({
    list(cv = cv.glmnet(xtr, y_train, alpha = 0,
                        nfolds = nf, standardize = FALSE))
  })

  # Random Forest
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

  # SVR radial
  set.seed(3000L + seed_base)
  fits$SVR_radial <- safe({
    tuned <- tune.svm(xtr, y_train, kernel = "radial",
                      cost  = c(0.1, 1, 10),
                      gamma = c(0.01, 0.1, 1 / ncol(xtr)),
                      epsilon = 0.1,
                      tunecontrol = tune.control(cross = 3))
    list(model = tuned$best.model)
  })

  list(fits = fits, genes = genes, zp = zp, medians = med)
}

# ── 3E  Apply a fitted bundle to new raw-scale expression data ────────────────
#  Aligns the new data to the bundle's genes, imputes, z-scores with the
#  TRAINING parameters, and returns one named prediction vector per model.
predict_bundle <- function(bundle, new_mat_raw) {
  x   <- align_to_genes(new_mat_raw, bundle$genes, bundle$medians)
  x   <- impute_medians(x, bundle$medians)
  x   <- apply_zparams(x, bundle$zp)
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
    if (!is.null(preds[[nm]])) names(preds[[nm]]) <- rownames(new_mat_raw)
  }
  preds
}

################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 4 — INPUT VALIDATION                                           ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

message("\n>>> Validating inputs...")

for (.f in c(NEW_RNA_FILE, TRAIN_RNA_FILE, TRAIN_META_FILE)) {
  if (!is.character(.f) || length(.f) != 1L || nchar(.f) == 0L)
    stop("A file path in Section 1 is not a valid string.\n",
         "  Check NEW_RNA_FILE, TRAIN_RNA_FILE, and TRAIN_META_FILE.")
  if (!file.exists(.f))
    stop("File not found: ", .f,
         "\n  Update the path in Section 1 (USER CONFIG) and re-source.")
}

STRATEGIES <- match.arg(STRATEGIES, c("train_set_Transfer"), several.ok = TRUE)

message(sprintf("    New RNA file    : %s", NEW_RNA_FILE))
message(sprintf("    Train RNA file  : %s", TRAIN_RNA_FILE))
message(sprintf("    Train meta file : %s", TRAIN_META_FILE))
message(sprintf("    Strategies      : %s", paste(STRATEGIES, collapse = ", ")))
message(sprintf("    Gene set        : %s", paste(names(GENE_SET_RUNS), collapse = ", ")))
message(sprintf("    Ensemble        : %s", paste(ENSEMBLE_MODELS, collapse = " + ")))
message(sprintf("    Output dir      : %s", OUT_DIR))


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 5 — PARALLELISATION SETUP                                      ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

message("\n>>> Setting up parallel workers...")

if (is.null(N_WORKERS))
  N_WORKERS <- max(1L, min(parallel::detectCores(logical = TRUE) - 1L, 11L))

.mem <- max(512 * 1024^2, min(2 * 1024^3, floor(12 * 1024^3 / N_WORKERS)))
options(future.globals.maxSize = .mem)
message(sprintf("    Workers     : %d  (detected %d logical cores)",
                N_WORKERS, parallel::detectCores(logical = TRUE)))
message(sprintf("    globals.max : %.1f GB per worker", .mem / 1024^3))
rm(.mem)

plan(multisession, workers = N_WORKERS)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 6 — LOAD DATA                                                  ║
# ║                                                                          ║
# ║  Training set:                                                           ║
# ║    · TRAIN_RNA_FILE  — expression matrix (column 1 = sample IDs)        ║
# ║    · TRAIN_META_FILE — metadata with Avg_2DT + OncotreeLineage cols
# ║      Rows are matched by the sample ID in column 1 of each file.        ║
# ║                                                                          ║
# ║  New samples:                                                            ║
# ║    · NEW_RNA_FILE    — expression matrix (column 1 = sample IDs)        ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

message("\n>>> Loading data...")

# ── New (prediction) samples ──────────────────────────────────────────────────
#  One CCLE file holds both sets: cell lines without Avg_2DT are predicted,
#  cell lines with Avg_2DT are used for training (below).

# Cell-line IDs are in column 2; the two remaining non-gene columns (1 and 3)
# are dropped. Avg_2DT stays in the matrix as a column; it is not in GENE_LIST,
# so it is never used as a model feature.
rna_raw <- read.csv(TRAIN_RNA_FILE, row = 2)[, -c(1, 3)]
new_rna_raw <- subset(rna_raw, is.na(Avg_2DT))
new_mat     <- as.matrix(new_rna_raw); storage.mode(new_mat) <- "numeric"
if (is.null(rownames(new_mat))) stop("NEW_RNA_FILE must have row names (sample IDs).")
if (is.null(colnames(new_mat))) stop("NEW_RNA_FILE must have column names (gene symbols).")
message(sprintf("    New samples : %d  |  genes provided : %d",
                nrow(new_mat), ncol(new_mat)))
rm(new_rna_raw)

# ── Training expression ───────────────────────────────────────────────────────
train_rna_raw <- subset(rna_raw, complete.cases(Avg_2DT))
TRAIN_expr_raw <- as.matrix(train_rna_raw); storage.mode(TRAIN_expr_raw) <- "numeric"
if (is.null(rownames(TRAIN_expr_raw))) stop("TRAIN_RNA_FILE must have row names (sample IDs).")
if (is.null(colnames(TRAIN_expr_raw))) stop("TRAIN_RNA_FILE must have column names (gene symbols).")
message(sprintf("    Training RNA: %d samples  |  %d genes",
                nrow(TRAIN_expr_raw), ncol(TRAIN_expr_raw)))
rm(train_rna_raw)

# ── Training metadata ─────────────────────────────────────────────────────────
#  Column 1 = sample ID → row names.
#  Required: Avg_2DT (numeric), OncotreeLineage (character).
train_meta <- read.csv(TRAIN_META_FILE, header = TRUE, row.names = 1)
train_meta <- subset(train_meta, complete.cases(Avg_2DT))   # cell lines with a measured DT

if (!"Avg_2DT" %in% colnames(train_meta))
  stop("TRAIN_META_FILE must contain a 'Avg_2DT' column.")
if (!"OncotreeLineage" %in% colnames(train_meta))
  stop("TRAIN_META_FILE must contain an 'OncotreeLineage' column.")

# Align training expression and metadata by shared sample IDs
common_train <- intersect(rownames(TRAIN_expr_raw), rownames(train_meta))
if (length(common_train) == 0L)
  stop("No row names shared between TRAIN_RNA_FILE and TRAIN_META_FILE.\n",
       "  Check that column 1 of both CSVs contains the same sample IDs.")

TRAIN_expr_raw <- TRAIN_expr_raw[common_train, , drop = FALSE]
train_meta     <- train_meta[common_train, , drop = FALSE]

TRAIN_DT      <- setNames(as.numeric(train_meta[["Avg_2DT"]]),   common_train)
TRAIN_lineage <- setNames(as.character(train_meta[["OncotreeLineage"]]), common_train)

# Remove samples with missing DT or lineage
ok_train      <- is.finite(TRAIN_DT) & !is.na(TRAIN_lineage) & TRAIN_lineage != ""
TRAIN_expr_raw <- TRAIN_expr_raw[ok_train, , drop = FALSE]
TRAIN_DT      <- TRAIN_DT[ok_train]
TRAIN_lineage <- TRAIN_lineage[ok_train]

message(sprintf("    Training set after filtering: %d samples  |  %d genes",
                nrow(TRAIN_expr_raw), ncol(TRAIN_expr_raw)))
rm(train_meta, common_train, ok_train)
gc()


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 7 — COMMON GENE UNIVERSE + VARIANCE FILTER                    ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

message("\n>>> Common gene universe + variance filter...")

# Genes present in both the training and the new data, keeping those with
# training variance > 1e-5. GENE_LIST is intersected with this set in
# Section 8.
common_genes     <- intersect(colnames(TRAIN_expr_raw), colnames(new_mat))
var_train_all    <- apply(TRAIN_expr_raw[, common_genes, drop = FALSE], 2,
                          var, na.rm = TRUE)
keep_genes_train <- common_genes[!is.na(var_train_all) & var_train_all > 1e-5]

message(sprintf(
  "    Train ∩ new : %d common genes  →  %d after variance filter",
  length(common_genes), length(keep_genes_train)))

if (length(keep_genes_train) < 10L)
  stop("Fewer than 10 genes pass the variance filter. ",
       "Check that column names are HGNC gene symbols with no version suffixes.")

keep_genes <- keep_genes_train

# Training-set z-score parameters over the filtered genes. These are only
# written to pDT_gene_info.csv for reference; the models use the parameters
# fit_all_models() computes on the selected genes (identical values for those
# genes).
train_means  <- colMeans(TRAIN_expr_raw[, keep_genes, drop = FALSE], na.rm = TRUE)
train_sds    <- pmax(apply(TRAIN_expr_raw[, keep_genes, drop = FALSE],
                           2, sd, na.rm = TRUE), 1e-8)
train_zp_sel <- list(means = train_means, sds = train_sds)

rm(train_means, train_sds)
gc()


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTIONS 8–11 — TRAIN MODELS + PREDICT                                   ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

set.seed(SEED)

# Accumulators across gene-set runs (there is one run, named MODEL_NAME)
all_pred_cols       <- list()
all_ensemble_cols   <- list()
all_ensemble_members <- list()
all_gene_infos      <- list()
all_trained         <- list()
all_cor_matrices    <- list()

for (.gs_name in names(GENE_SET_RUNS)) {
  .gs_def <- GENE_SET_RUNS[[.gs_name]]

  message(sprintf("\n╔══════════════════════════════════════════════════╗"))
  message(sprintf("║  Gene-set run: %-33s║", .gs_name))
  message(sprintf("╚══════════════════════════════════════════════════╝"))

  # ── 8. Resolve selected genes (user-supplied list) ──────────────────────────
  .genes_requested <- unique(.gs_def$genes)
  selected_train   <- intersect(.genes_requested, keep_genes_train)
  .missing_genes   <- setdiff(.genes_requested, keep_genes_train)

  message(sprintf(
    "    '%s': %d / %d requested genes usable (present in train ∩ new-data ∩ variance filter)",
    .gs_name, length(selected_train), length(.genes_requested)))
  if (length(.missing_genes) > 0L) {
    message(sprintf("    [warning] %d gene(s) dropped: %s",
                    length(.missing_genes), paste(.missing_genes, collapse = ", ")))
  }

  .gene_info <- data.frame(
    gene_set  = .gs_name, strategy = "train_set_Transfer",
    gene      = .genes_requested,
    selected  = .genes_requested %in% selected_train,
    stringsAsFactors = FALSE)
  all_gene_infos[[.gs_name]] <- .gene_info

  if (length(selected_train) < 2L) {
    warning(sprintf("Gene-set '%s': fewer than 2 usable genes — skipped.", .gs_name))
    next
  }

  # ── 9. Build parallel job table ────────────────────────────────────────────
  #  One job per strategy: the training expression for the selected genes,
  #  the training DT, and a seed offset.
  message(sprintf("    Building job table for '%s'...", .gs_name))
  .jobs <- list()

  if ("train_set_Transfer" %in% STRATEGIES) {
    .x_train <- TRAIN_expr_raw[, selected_train, drop = FALSE]


    .jobs[["train_set_Transfer"]] <- list(
      x_raw     = .x_train,
      y         = TRAIN_DT,
      seed_base = 0L,
      norm      = "train_set"
    )
  }

  message(sprintf("    %d job(s): %s", length(.jobs), paste(names(.jobs), collapse=", ")))

  # ── 10. Parallel model training ─────────────────────────────────────────────
  #  future.seed = TRUE gives each job its own reproducible random stream.
  message(sprintf("    Dispatching %d job(s) across %d worker(s)...",
                  length(.jobs), N_WORKERS))

  .trained_raw <- future_lapply(
    .jobs,
    function(job) {
      tryCatch({
        bnd <- fit_all_models(
          x_train_raw = job$x_raw,
          y_train     = job$y,
          genes       = colnames(job$x_raw),
          seed_base   = job$seed_base
        )
        if (is.null(bnd)) return(NULL)
        bnd$norm     <- job$norm
        bnd
      }, error = function(e) {
        message(sprintf("    [ERROR] %s", conditionMessage(e))); NULL
      })
    },
    future.seed     = TRUE,
    future.globals  = c(".jobs",
                        "fit_all_models", "store_zparams", "apply_zparams",
                        "impute_medians"),
    future.packages = c("glmnet","ranger","e1071")
  )
  names(.trained_raw) <- names(.jobs)
  .trained <- Filter(Negate(is.null), .trained_raw)
  if (length(.trained) == 0L) {
    warning(sprintf("Gene-set '%s': all training jobs failed — skipped.", .gs_name))
    next
  }
  message(sprintf("    Trained strategies: %s", paste(names(.trained), collapse = ", ")))

  for (.sn in names(.trained))
    all_trained[[if (length(.trained) == 1L) .gs_name
                 else paste0(.gs_name, ".", .sn)]] <- .trained[[.sn]]

  # ── 11. Predict + ensembles ─────────────────────────────────────────────────
  #  Output column names are  <MODEL_NAME>_<model>  (e.g. pDT_Ridge,
  #  pDT_ensemble).
  message(sprintf("    Predicting %d sample(s)...", nrow(new_mat)))
  .gs_pred_cols    <- list()
  .strat_pred_list <- list()   # per-strategy predictions, keyed by model name

  for (.strat in names(.trained)) {
    .bnd        <- .trained[[.strat]]
    .col_prefix <- paste0(.gs_name, "_")
    .strat_preds <- predict_bundle(.bnd, new_mat)
    .strat_pred_list[[.strat]] <- .strat_preds
    for (.mn in names(.strat_preds))
      .gs_pred_cols[[paste0(.col_prefix, .mn)]] <- .strat_preds[[.mn]]
  }

  all_pred_cols <- c(all_pred_cols, .gs_pred_cols)

  # Ensembles: per sample, the mean of the ENSEMBLE_MODELS predictions that
  # are available (at least 2 are required). The Spearman correlations between
  # the member models are kept for pDT_ensemble_correlations.csv.
  message(sprintf("    Building ensembles for '%s'...", .gs_name))

  for (.strat in names(.trained)) {
    .preds_s <- Filter(Negate(is.null), .strat_pred_list[[.strat]])
    if (length(.preds_s) < 2L) next
    .pm   <- do.call(cbind, .preds_s)
    colnames(.pm) <- names(.preds_s)
    .pm   <- .pm[, colSums(is.finite(.pm)) > 0, drop = FALSE]
    if (ncol(.pm) < 2L) next

    .ok_rows <- rowSums(!is.finite(.pm)) == 0
    .cor_mat <- cor(.pm[.ok_rows, , drop = FALSE], method = "spearman",
                    use = "complete.obs")
    all_cor_matrices[[paste0(.gs_name, "__", .strat)]] <- .cor_mat

    .ens_avail <- intersect(ENSEMBLE_MODELS, colnames(.pm))
    if (length(.ens_avail) >= 2L) {
      .ens_key <- paste0(.gs_name, "_ensemble")
      all_ensemble_members[[.ens_key]] <- .ens_avail
      all_ensemble_cols[[.ens_key]] <-
        setNames(rowMeans(.pm[, .ens_avail, drop = FALSE], na.rm = TRUE),
                 rownames(new_mat))
      message(sprintf("    %-30s ensemble → %s",
                      paste0(.gs_name, " / ", .strat),
                      paste(.ens_avail, collapse = " | ")))
    } else {
      message(sprintf("    %-30s ensemble — too few models available (%d), skipped.",
                      paste0(.gs_name, " / ", .strat), length(.ens_avail)))
    }
  }

}   # end gene-set loop


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 12 — ASSEMBLE OUTPUT DATA FRAME                                ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

message("\n>>> Assembling output data frame...")

# One row per new sample: sample_id, then pDT_<model> and pDT_ensemble.
all_cols <- c(all_pred_cols, all_ensemble_cols)

pDT_results <- data.frame(
  sample_id = rownames(new_mat),
  do.call(cbind, lapply(all_cols, function(p)
    if (is.null(p)) rep(NA_real_, nrow(new_mat)) else p)),
  row.names = NULL,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

message(sprintf("    pDT_results : %d rows × %d columns  (%d individual + %d ensemble)",
                nrow(pDT_results),
                ncol(pDT_results) - 1L,
                length(all_pred_cols),
                length(all_ensemble_cols)))


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 13 — SAVE OUTPUTS                                              ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

message("\n>>> Saving outputs to: ", OUT_DIR)

# pDT_results.csv is the file the plot_figure_*.R scripts read.
write.csv(pDT_results,
          file.path(OUT_DIR, "pDT_results.csv"),
          row.names = FALSE, quote = FALSE)

# Which requested genes were used, with their training mean / SD
gene_info_all <- do.call(rbind, all_gene_infos)
gene_info_all$zscore_mean_train <- train_zp_sel$means[gene_info_all$gene]
gene_info_all$zscore_sd_train   <- train_zp_sel$sds[gene_info_all$gene]
write.csv(gene_info_all,
          file.path(OUT_DIR, "pDT_gene_info.csv"),
          row.names = FALSE, quote = FALSE)

# Member-model correlation matrices in long format
cor_long <- do.call(rbind, lapply(names(all_cor_matrices), function(run_key) {
  mat   <- all_cor_matrices[[run_key]]
  parts <- strsplit(run_key, "__", fixed = TRUE)[[1L]]
  gene_set <- parts[1L]
  strategy <- paste(parts[-1L], collapse = "__")
  idx <- which(!is.na(mat), arr.ind = TRUE)
  data.frame(
    gene_set   = gene_set,
    strategy   = strategy,
    model_row  = rownames(mat)[idx[, 1L]],
    model_col  = colnames(mat)[idx[, 2L]],
    spearman_r = mat[idx],
    stringsAsFactors = FALSE
  )
}))
write.csv(cor_long,
          file.path(OUT_DIR, "pDT_ensemble_correlations.csv"),
          row.names = FALSE, quote = FALSE)

# Fitted models (with their z-score parameters), and an Excel summary
saveRDS(all_trained, file.path(OUT_DIR, "pDT_model_objects.rds"))

wb <- createWorkbook()
addWorksheet(wb, "pDT_results");  writeData(wb, "pDT_results",  pDT_results)
addWorksheet(wb, "Gene_Info");    writeData(wb, "Gene_Info",    gene_info_all)
ens_tbl <- do.call(rbind, lapply(names(all_ensemble_members), function(k)
  data.frame(ensemble  = k,
             model     = all_ensemble_members[[k]],
             stringsAsFactors = FALSE)))
addWorksheet(wb, "Ensemble_Members"); writeData(wb, "Ensemble_Members", ens_tbl)
addWorksheet(wb, "Ensemble_Correlations"); writeData(wb, "Ensemble_Correlations", cor_long)
saveWorkbook(wb, file.path(OUT_DIR, "pDT_predictions_summary.xlsx"),
             overwrite = TRUE)

message("    pDT_results.csv")
message("    pDT_gene_info.csv")
message("    pDT_ensemble_correlations.csv")
message("    pDT_model_objects.rds")
message("    pDT_predictions_summary.xlsx")


################################################################################
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 14 — SUMMARY + RELEASE WORKERS                                 ║
# ╚══════════════════════════════════════════════════════════════════════════╝
################################################################################

message("\n═══════════════════════════════════════════════════════════")
message(">>> Done.")
message(sprintf("    Samples      : %d", nrow(pDT_results)))
message(sprintf("    Training src : %s  (%d samples)", TRAIN_RNA_FILE, nrow(TRAIN_expr_raw)))
message(sprintf("    Gene set     : %s", paste(names(GENE_SET_RUNS), collapse=", ")))
message(sprintf("    Strategies   : %s", paste(STRATEGIES, collapse=", ")))
message(sprintf("    Ensemble     : %s", paste(ENSEMBLE_MODELS, collapse=" + ")))
message(sprintf("    Columns      : %d individual  +  %d ensemble",
                length(all_pred_cols), length(all_ensemble_cols)))
message(sprintf("    Output       : %s/pDT_results.csv", OUT_DIR))
message("    R object     : pDT_results")
message("═══════════════════════════════════════════════════════════")

plan(sequential)   # shut down the parallel workers

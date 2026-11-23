################################################################################
# plot_figure_3.R
# Figure 3, Supplemental Figure 3 and Supplemental Table 4: predicted doubling
# time (pDT) in two patient cohorts with an independent growth measurement.
#
#   HCC (GSE54236)  measured tumour DT, MKI67 and overall survival
#   LGG             subtype (1p/19q codel vs non-codel) and tumour-volume
#                   doubling time (TVDT)
#
# Run after predict_hcc.R, predict_lgg.R and predict_hcc_LOOCV.R, which write
# the pDT tables read here.
#
# WHAT THIS SCRIPT WRITES (all in WORK_DIR/figure_3/)
#   Figure_3.pdf / .png          panels a-f
#   sup_Figure_3.pdf / .png      panels a-c
#   sup_tab4.xlsx                Contents + Tables 4.1 (HCC) and 4.2 (LGG)
#
# Figure 3 panels
#   a  HCC: measured DT vs pDT
#   b  HCC: measured DT vs MKI67 expression
#   c  HCC: Kaplan-Meier, short vs long MEASURED DT (lowest quartile vs rest)
#   d  HCC: Kaplan-Meier, short vs long PREDICTED DT (lowest quartile vs rest)
#   e  LGG: pDT by subtype
#   f  LGG: MKI67 by subtype
# Supplemental Figure 3 panels
#   a  HCC leave-one-out CV: measured vs predicted DT (HCC-only training)
#   b  LGG codel: TVDT vs pDT
#   c  LGG non-codel: TVDT vs pDT
################################################################################

################################################################################
# WORKING DIRECTORY — defaults to the folder this script is in
#
# Same WORK_DIR as the other figure and predict_*.R scripts. Inputs are read
# from WORK_DIR/input_files/ and from the predict_*.R output folders; outputs
# go to WORK_DIR/figure_3/.
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
OUT_DIR <- file.path(WORK_DIR, "figure_3")

# ── Install missing packages ───────────────────────────────────────────────────
# Installs only the packages that are not already installed; packages that are
# present are left as they are (no updates). All are on CRAN.
CRAN_PKGS <- c("ggplot2", "ggpubr", "latex2exp", "survival", "survminer",
               "patchwork", "openxlsx", "readxl", "tibble", "tidyr")

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

# readxl, tibble, tidyr: used via :: only, so installed but not attached.
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(survival)
library(survminer)
library(patchwork)
library(openxlsx)


# ── format_p_value() ──────────────────────────────────────────────────────────
# Formats p-values as plotmath strings, e.g. "italic(P) == 1.23 %*% 10^-04".
# Use with parse = TRUE (stat_cor / stat_compare_means labels below).
# Values below `threshold` are shown with no decimal places on the mantissa.
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

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

################################################################################
# HCC COHORT (GSE54236) — expression, clinical data and pDT
################################################################################

# Expression: genes x samples on disk, transposed so samples are rows.
rna_hcc <- as.data.frame(t(read.table(file.path(WORK_DIR, "input_files", "hcc", "GSE54236_qnorm.tsv"), sep = "\t", header = TRUE,row = 1)))
# Clinical metadata: sample title, measured doubling time (days), overall
# survival (months) and sex.
meta_hcc <- read.table(file.path(WORK_DIR, "input_files", "hcc", "GSE54236_clinical_metadata_compact.tsv"),
                       sep = "\t",
                       header = TRUE,
                       row = 1)[,c("Sample_title","doubling.time..days.","survival.time.months.","gender")]

# Parse Sample_title into sample type / patient / biopsy. The titles come in two
# formats: rows 1-81 are "<type> <sample>" and rows 82-161 are
# "<type> <patient> <sample>". The row ranges assume the file's row order and
# row count are as supplied (161 rows); check them if the file changes.
meta_hcc <- rbind(tidyr::separate(meta_hcc[c(1:81),],col = (Sample_title), sep = " ", into = c("type", "sample"), remove = FALSE),
                  tidyr::separate(meta_hcc[c(82:161),],col = (Sample_title), sep = " ", into = c("type", "patient","sample"), remove = FALSE)[,-c(3)])[,-c(2)]

# The sample field is "<type>_<patient>_<biopsy>"; split it into its parts
# (type "T" = tumour).
meta_hcc <- tidyr::separate(meta_hcc,col = (sample), sep = "_", into = c("type", "patient","biopsy"), remove = FALSE)

# pDT from predict_hcc.R (PDX-trained ensemble), keyed by sample ID.
pDT_results_hcc <- read.csv(file.path(WORK_DIR, "hcc_pred", "pDT_results.csv"))

pDT_results_hcc <- tibble::column_to_rownames(pDT_results_hcc, var = "sample_id")
# Add MKI67 expression and the clinical columns, matched on sample ID.
pDT_results_hcc$mki67 <- rna_hcc[row.names(pDT_results_hcc),"MKI67"]

pDT_results_hcc <- cbind(pDT_results_hcc,meta_hcc[rownames(pDT_results_hcc),])

# Tumour samples only. Biopsies 047, 127 and 107 are excluded here and in
# plot_figure_1.R and predict_hcc_LOOCV.R, so all three use the same tumour set.
pDT_results_hcc_tumor <- subset(pDT_results_hcc, (type == "T") & !(biopsy %in% c("047","127","107")))

# Vital status in 2018 from the extended clinical table (header on row 2, sample
# IDs in the first, unnamed, column). Column 9 is dropped.
surv_stat <- tibble::column_to_rownames(as.data.frame(readxl::read_xlsx(file.path(WORK_DIR, "input_files", "hcc", "GSE54236_clinical_meta_compact_EV.xlsx"),
                                                                        skip = 1)),
                                        var = "...1")[,-c(9)]
pDT_results_hcc_tumor <- cbind(pDT_results_hcc_tumor,surv_stat[rownames(pDT_results_hcc_tumor),c("Status in 2018"), drop = FALSE])
# Recode vital status to the "1:DECEASED" / "0:LIVING" labels used by the
# survival models below. In the source sheet 0 = deceased.
pDT_results_hcc_tumor$`Status in 2018` <- ifelse(test = pDT_results_hcc_tumor$`Status in 2018` == 0, yes = "1:DECEASED", no = "0:LIVING")

# Short vs long DT groups for the Kaplan-Meier panels: samples at or below the
# 25th percentile are "Short_DT", the rest "Long_DT". Done separately for the
# measured DT (TDT) and the predicted DT (pDT).
pDT_results_hcc_tumor$TDT_group_quantile <- factor(ifelse(pDT_results_hcc_tumor$doubling.time..days. <=  quantile(pDT_results_hcc_tumor$doubling.time..days.)["25%"], "Short_DT", "Long_DT"),
                                                   levels = c("Short_DT", "Long_DT"))

pDT_results_hcc_tumor$pDT_group_quantile <- factor(ifelse(pDT_results_hcc_tumor$pDT_ensemble <=  quantile(pDT_results_hcc_tumor$pDT_ensemble)["25%"], "Short_DT", "Long_DT"),
                                                   levels = c("Short_DT", "Long_DT"))

################################################################################
# LGG COHORT — pDT, MKI67 and tumour-volume growth
################################################################################

# pDT from predict_lgg.R. Sample IDs use "." where the volume table uses "_",
# so both are normalised to "_".
pDT_results_lgg <- read.csv(file.path(WORK_DIR, "lgg_pred", "pDT_results.csv"), row = 1)
rownames(pDT_results_lgg) <- gsub("\\.", "_", rownames(pDT_results_lgg))

# LGG expression (genes x samples, log2(TPM + 1)); used here for MKI67 only.
log_2_tpm_short <- read.csv(file.path(WORK_DIR, "input_files", "lgg", "rna_input_lgg.csv"), row = 1)
colnames(log_2_tpm_short) <- gsub("\\.", "_", colnames(log_2_tpm_short))

# Tumour-volume growth table (Table S1, sheet TableS1_TVGR_estimation): one row
# per sample with subtype and tumour-volume doubling time (TVDTyears).
TVGR_tab <- tibble::column_to_rownames(as.data.frame(readxl::read_excel(file.path(WORK_DIR, "input_files", "lgg", "TableS1_GRvolume.xlsx"), sheet = "TableS1_TVGR_estimation")),var = "ID")

# Add pDT and MKI67 to the volume table, matched on sample ID.
TVGR_tab$pDT <- pDT_results_lgg[rownames(TVGR_tab),"pDT_ensemble"]
TVGR_tab$mki67 <- as.data.frame(t(log_2_tpm_short))[rownames(TVGR_tab),"MKI67"]

################################################################################
# FIGURE 3 PANELS
################################################################################

# ── Panel a — HCC: measured DT vs pDT (log10 both axes, Spearman) ─────────────
panel_A <- ggplot(data=pDT_results_hcc_tumor,
                  aes(x=log10(`doubling.time..days.`),y=log10(pDT_ensemble)))+
  geom_point(size=1, alpha=0.5) +
  geom_smooth(method="lm", se=FALSE, na.rm = TRUE) +
  stat_cor(aes(label = paste(..r.label.., after_stat(format_p_value(..p..)),sep = "~`,`~")), r.accuracy = 0.001,method = "spearman",parse = TRUE,size=5, label.y = 1.48, label.x = 1.8)+
  xlab(TeX('$log_{10}(Measured\\ DT\\ (days))$'))+
  ylab(TeX('$log_{10}(pDT)$'))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"))

# ── Panel b — HCC: measured DT vs MKI67 expression (Spearman) ─────────────────
panel_B <- ggplot(data=pDT_results_hcc_tumor,
                  aes(x=log10(`doubling.time..days.`),y=mki67))+
  geom_point(size=1, alpha=0.5) +
  geom_smooth(method="lm", se=FALSE, na.rm = TRUE) +
  stat_cor(aes(label = paste(..r.label.., after_stat(format_p_value(..p..)),sep = "~`,`~")), r.accuracy = 0.001,method = "spearman",parse = TRUE,size=5, label.x = 1.8)+
  xlab(TeX('$log_{10}(Measured\\ DT\\ (days))$'))+
  ylab(TeX("$log_{2}(mki67_{TPM}+1)$")) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"))

# ── Panel c — HCC overall survival by MEASURED DT group ───────────────────────
# Kaplan-Meier curves with log-rank p and a risk table. Red = short DT.
panel_C <- ggsurvplot(
  survfit(Surv(time  = survival.time.months.,
               event = `Status in 2018` == "1:DECEASED") ~ TDT_group_quantile,
          data = pDT_results_hcc_tumor),
  data = pDT_results_hcc_tumor,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("#DC3220", "#005AB5"),
  ylab = "Survival Probability",
  legend.title = "Measured DT",
  legend.labs = c("Short DT", "Long DT"),
  ggtheme = theme_classic(base_size = 12) +
    theme(legend.position = "none",
          axis.text  = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"))
)

# Label the curve panel "Measured" (x is in months)
panel_C$plot <- panel_C$plot +
  annotate("text",
           x = 30,                      # months; must fall inside the time range
           y = 0.9,
           label = "Measured",
           size = 5,
           fontface = "bold",
           hjust = 0)

# ── Panel d — HCC overall survival by PREDICTED DT group ──────────────────────
# Same layout as panel c, grouped by pDT instead of measured DT.
panel_D <- ggsurvplot(
  survfit(Surv(time  = pDT_results_hcc_tumor$survival.time.months. ,event = pDT_results_hcc_tumor$`Status in 2018` == "1:DECEASED") ~ pDT_group_quantile,
          data = pDT_results_hcc_tumor),
  data = pDT_results_hcc_tumor,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("Short DT" = "#DC3220","Long DT"  = "#005AB5"),
  ylab = "Survival Probability",
  legend.title = "Predicted DT",
  legend.labs = c("Short DT", "Long DT"),
  ggtheme = theme_classic(base_size = 12)+
    theme(legend.position = "none",
          axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"))
)


# Label the curve panel "Predicted" (x is in months)
panel_D$plot <- panel_D$plot +
  annotate("text",
           x = 30,                      # months; must fall inside the time range
           y = 0.9,
           label = "Predicted",
           size = 5,
           fontface = "bold",
           hjust = 0)

# ── Panel e — LGG: pDT by subtype (violin + box + points, Wilcoxon p) ─────────
panel_E <- ggplot(data = TVGR_tab,
                  aes(x = subtype,y = log10(pDT),colour = subtype)) +
  geom_violin(width = 0.8,trim = TRUE,scale = "width") +
  geom_jitter()+
  geom_boxplot(width = 0.2,outlier.shape = NA,colour = "black") +
  xlab(" ") +
  stat_compare_means(na.rm = TRUE, size=5, aes(label = after_stat(format_p_value(..p..))), parse = TRUE, vjust = 0.4, hjust = 0.3, label.y = 1.46)+
  ylab(TeX('$log_{10}(pDT)$'))+
  scale_colour_manual(labels = c("codel" = "Codel", "non-codel" = "Non Codel"),
                      values = c("codel" = "#2CA02C","non-codel" = "#8E44AD"))+
  scale_x_discrete(labels = c("Codel", "Non Codel")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold")
  )

# ── Panel f — LGG: MKI67 by subtype (violin + box + points, Wilcoxon p) ───────
panel_F <- ggplot(data = TVGR_tab,
                  aes(x = subtype,y = mki67,colour = subtype)) +
  geom_violin(width = 0.8,trim = TRUE,scale = "width") +
  geom_jitter()+
  geom_boxplot(width = 0.2,outlier.shape = NA,colour = "black") +
  xlab(" ") +
  stat_compare_means(na.rm = TRUE, size=5, aes(label = after_stat(format_p_value(..p..))), parse = TRUE, vjust = 0.4, hjust = 0.2, label.y = 5)+
  ylab(TeX("$log_{2}(mki67_{TPM}+1)$")) +
  scale_colour_manual(labels = c("codel" = "Codel", "non-codel" = "Non Codel"),
                      values = c("codel" = "#2CA02C","non-codel" = "#8E44AD"))+
  scale_x_discrete(labels = c("Codel", "Non Codel")) +
  scale_y_continuous(limits = c(1,5.5),breaks = seq(1,5,1)) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold")
  )


################################################################################
# ASSEMBLE FIGURE 3
################################################################################

# ggsurvplot returns the curves ($plot) and the risk table ($table) as separate
# ggplots. Strip the x label from the curves and the title from the tables,
# then stack each pair (curves 2/3 of the height, table 1/3).
km_C   <- panel_C$plot  + labs(x = NULL, y = "Survival \n Probability") + theme(legend.position = "none")
km_D   <- panel_D$plot  + labs(x = NULL, y = "Survival \n Probability") + theme(legend.position = "none")
tab_C  <- panel_C$table + labs(x = "Time (months)", y = NULL, title = NULL)
tab_D  <- panel_D$table + labs(x = "Time (months)", y = NULL, title = NULL)

pan_c <- km_C + tab_C +
  plot_layout(
    design = paste(c("A",
                     "B"),
                   collapse = "\n"),
    heights = c(2, 1))

pan_d <- km_D + tab_D +
  plot_layout(
    design = paste(c("A",
                     "B"),
                   collapse = "\n"),
    heights = c(2, 1))

# 3 rows x 2 columns: a | b, c | d, e | f. Rows 1 and 3 are 2 units high, the
# survival row (curves + risk tables) 3 units.
fig3 <- wrap_elements(full = panel_A) + wrap_elements(full = panel_B) +
  wrap_elements(full = pan_c) + wrap_elements(full = pan_d) +
  wrap_elements(full = panel_E) + wrap_elements(full = panel_F) + 
  plot_layout(
    design = paste(c("AAABBB",
                     "CCCDDD",
                     "EEEFFF"),
                   collapse = "\n"),
    heights = c(2, 3, 2)) + 
  plot_annotation(tag_levels = list(c("a", "b", "c", "d","e","f"))) & 
  theme(plot.tag = element_text(size = 15, face  = "bold"), plot.tag.position = c(0, 1.00))

ggsave(file.path(OUT_DIR, "Figure_3.pdf"), fig3, width = 9, height = 10)
ggsave(file.path(OUT_DIR, "Figure_3.png"), fig3, width = 9, height = 10)

################################################################################
# SUPPLEMENTAL TABLE 4 (sup_tab4.xlsx)
################################################################################

## ---- Contents sheet, written by hand ---------------------------------------
## One row per sheet. Edit the text freely; keep the two vectors the same
## length, or data.frame() will fail or recycle values.
sup_tab_4_contents <- data.frame(
  sheet = c(
    "Table 4.1",
    "Table 4.2"),
  
  description = c(
    "Clinical data for the HCC cohort.",
    "Clinical data for the LGG cohort."),
  stringsAsFactors = FALSE)

sup_tab_4_contents <- sup_tab_4_contents[, c("sheet", "description")]

# Table 4.1 — HCC tumours: sample, MKI67, pDT, measured DT, survival and the
# DT groups used in panels c and d.
pDT_results_hcc_tumor$sample_name <- rownames(pDT_results_hcc_tumor)
sup_table_4_1 <- pDT_results_hcc_tumor[,c("sample_name","Sample_title","mki67","pDT_ensemble","doubling.time..days.","survival.time.months.","Status in 2018","pDT_group_quantile","TDT_group_quantile")]
colnames(sup_table_4_1) <- c("sample_name","sample_title","mki67_expression","pDT","doubling_time_days","survival_time_months","survival_status_in_2018","pDT_group","TDT_group")

# Table 4.2 — LGG: sample, MKI67, pDT, TVDT and subtype (samples with a pDT).
TVGR_tab$sample_name <- rownames(TVGR_tab)
sup_table_4_2 <- subset(TVGR_tab[,c("sample_name","mki67", "pDT", "TVDTyears","subtype")], complete.cases(pDT))
colnames(sup_table_4_2) <- c("sample_name","mki67_expression","pDT","TVDTyears","subtype")

# Write the workbook: Contents first, then one sheet per table.
sup_tab_4 <- createWorkbook()

addWorksheet(sup_tab_4, "Contents")
addWorksheet(sup_tab_4, "sup_tab_4.1")
addWorksheet(sup_tab_4, "sup_tab_4.2")

writeData(sup_tab_4, sheet = "Contents", x = sup_tab_4_contents)
writeData(sup_tab_4, sheet = "sup_tab_4.1", x = sup_table_4_1)
writeData(sup_tab_4, sheet = "sup_tab_4.2", x = sup_table_4_2)
saveWorkbook(sup_tab_4, file.path(OUT_DIR, "sup_tab4.xlsx"), overwrite = TRUE)

################################################################################
# SUPPLEMENTAL FIGURE 3 PANELS
################################################################################

# ── Supp. panel a — HCC leave-one-out CV ──────────────────────────────────────
# Predictions from predict_hcc_LOOCV.R: models trained on the other HCC tumours
# only (no PDX data), z-scored per fold. Ensemble_Fixed = mean of Ridge,
# RandomForest and SVR_radial.
sup_panel_A <- ggplot(data=subset(read.csv(file.path(WORK_DIR, "loocv_hcc", "pDT_loocv_testonly_predictions.csv")), (model == "Ensemble_Fixed") & (norm_mode == "zscore")),
                      aes(x=log10(actual),y=log10(predicted)))+
  geom_point(size=1, alpha=0.5) +
  geom_smooth(method="lm", se=FALSE, na.rm = TRUE) +
  stat_cor(aes(label = paste(..r.label.., after_stat(format_p_value(..p..)),sep = "~`,`~")), r.accuracy = 0.001,method = "spearman",parse = TRUE,size=5)+
  xlab(TeX('$log_{10}(TDT)(days)$'))+
  ylab(TeX('$log_{10}(pDT)$'))+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"))

# ── Supp. panel b — LGG codel: tumour-volume DT vs pDT ────────────────────────
sup_panel_B <- ggplot(data=subset(TVGR_tab, subtype == "codel"),
                      aes(x=log10(TVDTyears),y=log10(pDT)))+
  geom_point(size=1, alpha=0.5) +
  geom_smooth(method="lm", se=FALSE, na.rm = TRUE) +
  stat_cor(aes(label = paste(..r.label.., after_stat(format_p_value(..p..)),sep = "~`,`~")), r.accuracy = 0.001,method = "spearman",parse = TRUE,size=5)+
  xlab(TeX('$log_{10}(TVDT)(Years)$'))+
  ylab(TeX('$log_{10}(pDT)$'))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"))

# ── Supp. panel c — LGG non-codel: tumour-volume DT vs pDT ────────────────────
sup_panel_C <- ggplot(data=subset(TVGR_tab, subtype == "non-codel"),
                      aes(x=log10(TVDTyears),y=log10(pDT)))+
  geom_point(size=1, alpha=0.5) +
  geom_smooth(method="lm", se=FALSE, na.rm = TRUE) +
  stat_cor(aes(label = paste(..r.label.., after_stat(format_p_value(..p..)),sep = "~`,`~")), r.accuracy = 0.001,method = "spearman",parse = TRUE,size=5, label.y = 1.46 )+
  xlab(TeX('$log_{10}(TVDT)(Years)$'))+
  ylab(TeX('$log_{10}(pDT)$'))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"))

# ── Assemble Supplemental Figure 3: three panels in one row ───────────────────
sup_fig3 <- sup_panel_A + sup_panel_B + sup_panel_C +
  plot_layout(
    design = paste(c("ABC"),
                   collapse = "\n"),
    heights = c(1)) + 
  plot_annotation(tag_levels = list(c("a", "b", "c"))) & 
  theme(plot.tag = element_text(size = 15, face  = "bold"), plot.tag.position = c(0, 1.00))

ggsave(file.path(OUT_DIR, "sup_Figure_3.pdf"), sup_fig3, width = 10, height = 3.3)
ggsave(file.path(OUT_DIR, "sup_Figure_3.png"), sup_fig3, width = 10, height = 3.3)

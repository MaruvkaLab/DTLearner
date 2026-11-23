################################################################################
# plot_figure_4.R
# Figure 4, Supplemental Figures 4 and 5 and Supplemental Table 5: predicted doubling
# time (pDT) across TCGA, and its association with mutation count, stage and
# survival.
#
# Run after predict_tcga.R, which writes tcga_pred/pDT_results.csv.
#
# WHAT THIS SCRIPT WRITES (all in WORK_DIR/figure_4/)
#   Figure_4.pdf / .png          panels a-e
#   sup_Figure_4.pdf / .png      panels a-g
#   sup_Figure_5.pdf / .png      panels a-e
#   sup_tab5.xlsx                Contents + Tables 5.1-5.8
#   DSS/, OS/, PFS/, DFS/        one Kaplan-Meier PDF per cancer type and
#                                survival endpoint
#
# Figure 4 panels (letters as printed; set by tag_levels in the assembly)
#   a  pDT by cancer type (types with > 300 samples)            panel_A
#   b  per cancer type: BH q-value of mutation count vs pDT     panel_C
#   c  BRCA: mutation count vs pDT                              panel_B
#   d  per cancer type: BH q-value of the DSS log-rank test     panel_E
#   e  MESO: disease-specific survival, short vs long pDT       panel_D
# Supplemental Figure 4 panels
#   a     pDT by cancer type, all types
#   b-d   pDT by AJCC stage / lymph-node stage / T stage in one cancer type
#   e-g   per cancer type q-values for those three stage correlations
# Supplemental Figure 5 panels
#   a     DSS forest plot: Cox HR (95% CI), long vs short pDT, per cancer type
#         and Pan_Cancer (the same HRs as Figure 4d / Table 5.6)
#   b     DSS forest plot, Cox model stratified by cancer type and adjusted
#         for age and T stage: unadjusted vs
#         adjusted HR on the same complete-case patients
#   c-e   per cancer type q-values for OS, PFS and DFS (median pDT split)
################################################################################

################################################################################
# WORKING DIRECTORY — defaults to the folder this script is in
#
# Same WORK_DIR as the other figure and predict_*.R scripts. Inputs are read
# from WORK_DIR/input_files/ and from the predict_*.R output folders; outputs
# go to WORK_DIR/figure_4/.
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
OUT_DIR <- file.path(WORK_DIR, "figure_4")

# ── Install missing packages ───────────────────────────────────────────────────
# Installs only the packages that are not already installed; packages that are
# present are left as they are (no updates). All are on CRAN.
CRAN_PKGS <- c("dplyr", "ggplot2", "ggpubr", "latex2exp", "survival",
               "survminer", "patchwork", "openxlsx")

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

library(dplyr)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(survival)
library(survminer)
library(patchwork)
library(openxlsx)


# ── format_p_value() ──────────────────────────────────────────────────────────
# Formats p-values as plotmath strings, e.g. "italic(P) == 1.23 %*% 10^-04".
# Use with parse = TRUE (stat_cor labels below). Values below `threshold` are
# shown with no decimal places on the mantissa.
format_p_value <- function(p, threshold = 0.001, scientific = TRUE) {
  
  to_sci <- function(x, digits = 2) {
    sci <- formatC(x, format = "e", digits = digits)
    gsub("e", " %*% 10^", sci)
  }
  
  vapply(p, function(pi) {
    
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
# DATA — pDT joined to the TCGA clinical table
################################################################################

# pDT from predict_tcga.R, keyed by TCGA sample ID.
pDT_results <- read.csv(file.path(WORK_DIR, "tcga_pred", "pDT_results.csv"), row = 1)

# TCGA PanCancer Atlas clinical data (one row per sample). Only samples with a
# pDT are kept.
Meta_TCGA <- read.table(file = file.path(WORK_DIR, "input_files", "tcga", "TCGA_clinical_data.tsv"),
                        header=T,sep="\t",comment.char = "",quote = "")
Meta_TCGA$DT_pred <- pDT_results[Meta_TCGA$Sample.ID,"pDT_ensemble"]
Meta_TCGA <- subset(Meta_TCGA, complete.cases(DT_pred))

# Figure 4a shows only cancer types with more than 300 samples.
panel_a_ttypes <- names(which(table(Meta_TCGA$TCGA.PanCanAtlas.Cancer.Type.Acronym) > 300))

# x tick labels for panel a: "BRCA\n(n=123)", one per cancer type, counted
# after the > 300 filter
lin_n   <- table(subset(Meta_TCGA, TCGA.PanCanAtlas.Cancer.Type.Acronym %in% panel_a_ttypes)$TCGA.PanCanAtlas.Cancer.Type.Acronym)
x_labs  <- setNames(paste0(sort(subset(Meta_TCGA, TCGA.PanCanAtlas.Cancer.Type.Acronym %in% panel_a_ttypes)$TCGA.PanCanAtlas.Cancer.Type.Acronym), "\n(n=",
                           lin_n[sort(subset(Meta_TCGA, TCGA.PanCanAtlas.Cancer.Type.Acronym %in% panel_a_ttypes)$TCGA.PanCanAtlas.Cancer.Type.Acronym)], ")"),
                    sort(subset(Meta_TCGA, TCGA.PanCanAtlas.Cancer.Type.Acronym %in% panel_a_ttypes)$TCGA.PanCanAtlas.Cancer.Type.Acronym))


# Same for supplemental panel a, which shows every cancer type
lin_n_sup   <- table(Meta_TCGA$TCGA.PanCanAtlas.Cancer.Type.Acronym)
x_labs_sup  <- setNames(paste0(sort(Meta_TCGA$TCGA.PanCanAtlas.Cancer.Type.Acronym),"(n=",
                               lin_n_sup[sort(Meta_TCGA$TCGA.PanCanAtlas.Cancer.Type.Acronym)], ")"),
                        sort(Meta_TCGA$TCGA.PanCanAtlas.Cancer.Type.Acronym))

# ── Figure 4a — pDT by cancer type (> 300 samples), ordered by median pDT ─────
panel_A <- ggplot(subset(Meta_TCGA, TCGA.PanCanAtlas.Cancer.Type.Acronym %in% panel_a_ttypes),
                  aes(x = reorder(x=TCGA.PanCanAtlas.Cancer.Type.Acronym,X=DT_pred,FUN=median),
                      y = log10(DT_pred),
                      fill = TCGA.PanCanAtlas.Cancer.Type.Acronym)) +
  geom_jitter(aes(colour = TCGA.PanCanAtlas.Cancer.Type.Acronym), width = 0.22, alpha = 1.0,size = 1, shape = 16, na.rm = TRUE) +
  geom_violin(alpha = 0.55, linewidth = 0.25, scale = "width",bw = "nrd0", na.rm = TRUE) +
  geom_boxplot(width = 0.18, fill = "white", outlier.shape = NA,linewidth = 0.70, colour = "black", na.rm = TRUE) +
  scale_x_discrete(labels = x_labs)+
  xlab(" ") +
  ylab(TeX('$log_{10}(pDT)$'))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"))

# ── Supp. 4a — pDT by cancer type, all types, ordered by median pDT ───────────
sup_panel_A <- ggplot(Meta_TCGA,
                      aes(x = reorder(x=TCGA.PanCanAtlas.Cancer.Type.Acronym,X=DT_pred,FUN=median),
                          y = log10(DT_pred),
                          fill = TCGA.PanCanAtlas.Cancer.Type.Acronym)) +
  geom_jitter(aes(colour = TCGA.PanCanAtlas.Cancer.Type.Acronym), width = 0.22, alpha = 1.0,size = 1, shape = 16, na.rm = TRUE) +
  geom_violin(alpha = 0.55, linewidth = 0.25, scale = "width",bw = "nrd0", na.rm = TRUE) +
  geom_boxplot(width = 0.18, fill = "white", outlier.shape = NA,linewidth = 0.70, colour = "black", na.rm = TRUE) +
  xlab(" ") +
  scale_x_discrete(labels = x_labs_sup)+
  ylab(TeX('$log_{10}(pDT)$'))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"))

# Add a "Pan_Cancer" copy of every sample, so each per-cancer-type analysis
# below (correlations and survival) also gets a pooled row. From here on
# Meta_TCGA holds every sample twice: once under its cancer type, once under
# Pan_Cancer. Pan_Cancer is removed again before Table 5.1 is written.
Meta_TCGA_all <- Meta_TCGA
Meta_TCGA_all$TCGA.PanCanAtlas.Cancer.Type.Acronym <- "Pan-Cancer"
Meta_TCGA <- rbind(Meta_TCGA,Meta_TCGA_all)

# Spearman correlation of mutation count with pDT, per cancer type (plus
# Pan_Cancer), with BH-adjusted q-values across cancer types.
cor_results <- Meta_TCGA %>%
  group_by(TCGA.PanCanAtlas.Cancer.Type.Acronym) %>%
  summarise(
    # samples with both a mutation count and a pDT (those used in the test)
    N = sum(!is.na(Mutation.Count) & !is.na(DT_pred)),
    Spearman_r = cor(Mutation.Count, DT_pred, method = "spearman", use = "pairwise"),
    p_value = cor.test(Mutation.Count, DT_pred, method = "spearman")$p.value
  ) %>%
  arrange(p_value)

# Add BH-adjusted p-values
cor_results <- cor_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# ── Figure 4c — BRCA: mutation count vs pDT (log10 both axes) ────────────────
panel_B <- ggplot(data=subset(Meta_TCGA, TCGA.PanCanAtlas.Cancer.Type.Acronym == "BRCA"),
                  aes(x=log10(Mutation.Count),y=log10(DT_pred)))+
  geom_point(alpha = 0.4, colour = "#808080") +
  geom_smooth(method = "lm", colour = "#005AB5", se = FALSE, na.rm = TRUE) +
  stat_cor(aes(label = paste(..r.label.., after_stat(format_p_value(..p..)),sep = "~`,`~")),
           r.accuracy = 0.001,
           method = "spearman",
           parse = TRUE,
           size=5,
           label.y = 1.65)+
  ylab(TeX('$log_{10}(pDT)$'))+
  xlab(TeX('$log_{10}(Mutation\\ Count)$'))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        strip.text = element_text(size = 12),
        ggside.panel.scale = 0.2,
        ggside.axis.text = element_text(size = 10))

# ── Figure 4b — q-value of mutation count vs pDT, per cancer type ─────────────
# Bars are -log10(q); fill = log10(sample count). Dashed line at q = 0.1.
panel_C <- ggplot(cor_results, aes(x= -log10(p_adj), y=reorder(x=TCGA.PanCanAtlas.Cancer.Type.Acronym, X= - log10(p_adj)), fill=log10(N))) +
  geom_col() +
  theme_classic(base_size=14) +
  ylab(" ")+
  xlab(TeX('$-log_{10}(Q_{value})$'))+
  labs(fill = TeX('$log_{10}(n_{Sample})$'))+
  geom_vline(xintercept = 1, linetype = "dashed", colour = "red")+
  theme(legend.position = c(0.6,0.2),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold")
  )

################################################################################
# STAGE — collapse AJCC sub-stages to main stages
#
# Each staging column gets two versions: *_simple (a label, e.g. "STAGE II")
# and *_num (an ordinal number used for the Spearman correlations). Values not
# listed (e.g. "STAGE X", "NX", "TX") stay NA.
################################################################################

# Overall AJCC stage: I-IV -> 1-4
Meta_TCGA$disease_stage_simple <- NA
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code %in% c("STAGE IA","STAGE IB","STAGE IS","STAGE I")),"disease_stage_simple"] <- "STAGE I"
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code %in% c("STAGE IIA","STAGE IIB","STAGE IIC","STAGE II")),"disease_stage_simple"] <- "STAGE II"
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code %in% c("STAGE IIIA","STAGE IIIB","STAGE IIIC","STAGE III")),"disease_stage_simple"] <- "STAGE III"
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code %in% c("STAGE IVA","STAGE IVB","STAGE IVC","STAGE IV")),"disease_stage_simple"] <- "STAGE IV"

Meta_TCGA$disease_stage_num <- NA
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code %in% c("STAGE IA","STAGE IB","STAGE IS","STAGE I")),"disease_stage_num"] <- 1
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code %in% c("STAGE IIA","STAGE IIB","STAGE IIC","STAGE II")),"disease_stage_num"] <- 2
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code %in% c("STAGE IIIA","STAGE IIIB","STAGE IIIC","STAGE III")),"disease_stage_num"] <- 3
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code %in% c("STAGE IVA","STAGE IVB","STAGE IVC","STAGE IV")),"disease_stage_num"] <- 4

# Lymph-node stage: N0-N3 -> 1-4
Meta_TCGA$lymph_node_stage_simple <- NA
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code %in% c("N0","N0 (I+)","N0 (I-)","N0 (MOL+)")),"lymph_node_stage_simple"] <- "N0"
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code %in% c("N1","N1A","N1B","N1C","N1MI")),"lymph_node_stage_simple"] <- "N1"
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code %in% c("N2","N2A","N2B","N2C")),"lymph_node_stage_simple"] <- "N2"
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code %in% c("N3","N3A","N3B","N3C")),"lymph_node_stage_simple"] <- "N3"

Meta_TCGA$lymph_node_stage_num <- NA
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code %in% c("N0","N0 (I+)","N0 (I-)","N0 (MOL+)")),"lymph_node_stage_num"] <- 1
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code %in% c("N1","N1A","N1B","N1C","N1MI")),"lymph_node_stage_num"] <- 2
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code %in% c("N2","N2A","N2B","N2C")),"lymph_node_stage_num"] <- 3
Meta_TCGA[which(Meta_TCGA$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code %in% c("N3","N3A","N3B","N3C")),"lymph_node_stage_num"] <- 4

# Primary-tumour stage: T0-T4 -> 0-4
Meta_TCGA$tumor_stage_simple <- NA
Meta_TCGA[which(Meta_TCGA$American.Joint.Committee.on.Cancer.Tumor.Stage.Code %in% c("T0")),"tumor_stage_simple"] <- "T0"
Meta_TCGA[which(Meta_TCGA$American.Joint.Committee.on.Cancer.Tumor.Stage.Code %in% c("T1","T1A","T1A1","T1B","T1B1","T1B2","T1C")),"tumor_stage_simple"] <- "T1"
Meta_TCGA[which(Meta_TCGA$American.Joint.Committee.on.Cancer.Tumor.Stage.Code %in% c("T2","T2A","T2A1","T2A2","T2B","T2C")),"tumor_stage_simple"] <- "T2"
Meta_TCGA[which(Meta_TCGA$American.Joint.Committee.on.Cancer.Tumor.Stage.Code %in% c("T3","T3A","T3B","T3C")),"tumor_stage_simple"] <- "T3"
Meta_TCGA[which(Meta_TCGA$American.Joint.Committee.on.Cancer.Tumor.Stage.Code %in% c("T4","T4A","T4B","T4D")),"tumor_stage_simple"] <- "T4"

Meta_TCGA$tumor_stage_num <- NA
Meta_TCGA[which(Meta_TCGA$American.Joint.Committee.on.Cancer.Tumor.Stage.Code %in% c("T0")),"tumor_stage_num"] <- 0
Meta_TCGA[which(Meta_TCGA$American.Joint.Committee.on.Cancer.Tumor.Stage.Code %in% c("T1","T1A","T1A1","T1B","T1B1","T1B2","T1C")),"tumor_stage_num"] <- 1
Meta_TCGA[which(Meta_TCGA$American.Joint.Committee.on.Cancer.Tumor.Stage.Code %in% c("T2","T2A","T2A1","T2A2","T2B","T2C")),"tumor_stage_num"] <- 2
Meta_TCGA[which(Meta_TCGA$American.Joint.Committee.on.Cancer.Tumor.Stage.Code %in% c("T3","T3A","T3B","T3C")),"tumor_stage_num"] <- 3
Meta_TCGA[which(Meta_TCGA$American.Joint.Committee.on.Cancer.Tumor.Stage.Code %in% c("T4","T4A","T4B","T4D")),"tumor_stage_num"] <- 4

# Spearman correlation of AJCC stage with pDT, per cancer type. Cancer types
# with no staged samples are left out; q-values are BH across cancer types.
cor_results_disease <- subset(Meta_TCGA, !(`TCGA.PanCanAtlas.Cancer.Type.Acronym` %in%
                                             names(which(rowSums(table(Meta_TCGA[,c("TCGA.PanCanAtlas.Cancer.Type.Acronym","disease_stage_num")])) == 0)))) %>%
  group_by(TCGA.PanCanAtlas.Cancer.Type.Acronym) %>%
  summarise(
    # samples with both a stage and a pDT (those used in the test)
    N = sum(!is.na(disease_stage_num) & !is.na(DT_pred)),
    Spearman_r = cor.test(disease_stage_num, DT_pred, method = "spearman")$estimate,
    p_value = cor.test(disease_stage_num, DT_pred, method = "spearman")$p.value
  ) %>%
  arrange(p_value)

# Add BH-adjusted p-values
cor_results_disease <- cor_results_disease %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))


# Same for lymph-node stage
cor_results_lymph_node <- subset(Meta_TCGA, !(`TCGA.PanCanAtlas.Cancer.Type.Acronym` %in%
                                                names(which(rowSums(table(Meta_TCGA[,c("TCGA.PanCanAtlas.Cancer.Type.Acronym","lymph_node_stage_simple")])) == 0)))) %>%
  group_by(TCGA.PanCanAtlas.Cancer.Type.Acronym) %>%
  summarise(
    # samples with both a stage and a pDT (those used in the test)
    N = sum(!is.na(lymph_node_stage_num) & !is.na(DT_pred)),
    Spearman_r = cor.test(lymph_node_stage_num, DT_pred, method = "spearman")$estimate,
    p_value = cor.test(lymph_node_stage_num, DT_pred, method = "spearman")$p.value
  ) %>%
  arrange(p_value)

# Add BH-adjusted p-values
cor_results_lymph_node <- cor_results_lymph_node %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Same for primary-tumour (T) stage
cor_results_tumor <- subset(Meta_TCGA, !(`TCGA.PanCanAtlas.Cancer.Type.Acronym` %in%
                                           names(which(rowSums(table(Meta_TCGA[,c("TCGA.PanCanAtlas.Cancer.Type.Acronym","tumor_stage_num")])) == 0)))) %>%
  group_by(TCGA.PanCanAtlas.Cancer.Type.Acronym) %>%
  summarise(
    # samples with both a stage and a pDT (those used in the test)
    N = sum(!is.na(tumor_stage_num) & !is.na(DT_pred)),
    Spearman_r = cor.test(tumor_stage_num, DT_pred, method = "spearman")$estimate,
    p_value = cor.test(tumor_stage_num, DT_pred, method = "spearman")$p.value
  ) %>%
  arrange(p_value)

# Add BH-adjusted p-values
cor_results_tumor <- cor_results_tumor %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# ── Supp. 4b — KIRC: pDT by AJCC stage ────────────────────────────────────────
sup_panel_B <- ggplot(data = subset(Meta_TCGA, complete.cases(disease_stage_simple) & TCGA.PanCanAtlas.Cancer.Type.Acronym == "KIRC"), aes(x = disease_stage_num ,y = log10(DT_pred),colour = disease_stage_simple)) +
  geom_jitter()+
  geom_violin(width = 0.8,trim = TRUE,scale = "width",alpha = 0.5) +
  geom_boxplot(width = 0.2,outlier.shape = NA) +
  xlab(" ") +
  ylab(TeX('$log_{10}(pDT)$'))+
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels =c("STAGE I","STAGE II","STAGE III","STAGE IV"))+
  stat_cor(aes(label = paste(..r.label.., after_stat(format_p_value(..p..)),sep = "~`,`~")),
           r.accuracy = 0.001,
           method = "spearman",
           parse = TRUE,
           size=5,
           label.y = 1.7,
           color = "black")+
  geom_smooth(method="lm", se=FALSE, na.rm = TRUE, colour = "black", linetype = "dashed") +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 16),
    strip.placement = "outside",
    strip.background = element_blank()
  )

# ── Supp. 4e — q-value of AJCC stage vs pDT, per cancer type (fill = rho) ─────
# (object names follow the order the panels were written, not their printed
# letters; see the assembly below)
sup_panel_C <- ggplot(cor_results_disease, aes(x= -log10(p_adj), y=reorder(x=TCGA.PanCanAtlas.Cancer.Type.Acronym, X= - log10(p_adj)), fill=Spearman_r)) +
  geom_col() +
  theme_classic(base_size=14) +
  ylab(" ")+
  xlab(TeX('$-log_{10}(Q_{value})$'))+
  labs(fill = TeX('Rho'))+
  geom_vline(xintercept = 1, linetype = "dashed", colour = "red")+
  theme(legend.position = c(0.8,0.2),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        strip.text = element_text(size = 12),
        ggside.panel.scale = 0.2,
        ggside.axis.text = element_text(size = 10))

# ── Supp. 4c — KIRP: pDT by lymph-node stage ──────────────────────────────────
sup_panel_D <- ggplot(data = subset(Meta_TCGA, complete.cases(lymph_node_stage_simple) & TCGA.PanCanAtlas.Cancer.Type.Acronym == "KIRP"), aes(x = lymph_node_stage_num ,y = log10(DT_pred),colour = lymph_node_stage_simple)) +
  geom_jitter()+
  geom_violin(width = 0.8,trim = TRUE,scale = "width",alpha = 0.5) +
  geom_boxplot(width = 0.2,outlier.shape = NA) +
  xlab(" ") +
  ylab(TeX('$log_{10}(pDT)$'))+
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("N0", "N1", "N2", "N3"))+
  stat_cor(aes(label = paste(..r.label.., after_stat(format_p_value(..p..)),sep = "~`,`~")),
           r.accuracy = 0.001,
           method = "spearman",
           parse = TRUE,
           size=5,
           label.y = 1.7,
           color = "black")+
  geom_smooth(method="lm", se=FALSE, na.rm = TRUE, colour = "black", linetype = "dashed") +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 16),
    strip.placement = "outside",
    strip.background = element_blank()
  )

# ── Supp. 4f — q-value of lymph-node stage vs pDT, per cancer type ────────────
sup_panel_E <- ggplot(cor_results_lymph_node, aes(x= -log10(p_adj), y=reorder(x=TCGA.PanCanAtlas.Cancer.Type.Acronym, X= - log10(p_adj)), fill=Spearman_r)) +
  geom_col() +
  theme_classic(base_size=14) +
  ylab(" ")+
  xlab(TeX('$-log_{10}(Q_{value})$'))+
  labs(fill = TeX('Rho'))+
  geom_vline(xintercept = 1, linetype = "dashed", colour = "red")+
  theme(legend.position = c(0.8,0.2),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        strip.text = element_text(size = 12),
        ggside.panel.scale = 0.2,
        ggside.axis.text = element_text(size = 10))

# ── Supp. 4d — KIRC: pDT by T stage ───────────────────────────────────────────
sup_panel_F <- ggplot(data = subset(Meta_TCGA, complete.cases(tumor_stage_simple) & TCGA.PanCanAtlas.Cancer.Type.Acronym == "KIRC"), aes(x = tumor_stage_num ,y = log10(DT_pred),colour = tumor_stage_simple)) +
  geom_jitter()+
  geom_violin(width = 0.8,trim = TRUE,scale = "width",alpha = 0.5) +
  geom_boxplot(width = 0.2,outlier.shape = NA) +
  xlab(" ") +
  ylab(TeX('$log_{10}(pDT)$'))+
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4), labels =c("T0","T1","T2","T3","T4"))+
  stat_cor(aes(label = paste(..r.label.., after_stat(format_p_value(..p..)),sep = "~`,`~")),
           r.accuracy = 0.001,
           method = "spearman",
           parse = TRUE,
           size=5,
           label.y = 1.7,
           color = "black")+
  geom_smooth(method="lm", se=FALSE, na.rm = TRUE, colour = "black", linetype = "dashed") +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 16),
    strip.placement = "outside",
    strip.background = element_blank()
  )

# ── Supp. 4g — q-value of T stage vs pDT, per cancer type ─────────────────────
sup_panel_G <- ggplot(cor_results_tumor, aes(x= -log10(p_adj), y=reorder(x=TCGA.PanCanAtlas.Cancer.Type.Acronym, X= - log10(p_adj)), fill=Spearman_r)) +
  geom_col() +
  theme_classic(base_size=14) +
  ylab(" ")+
  xlab(TeX('$-log_{10}(Q_{value})$'))+
  labs(fill = TeX('Rho'))+
  geom_vline(xintercept = 1, linetype = "dashed", colour = "red")+
  theme(legend.position = c(0.8,0.2),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        strip.text = element_text(size = 12),
        ggside.panel.scale = 0.2,
        ggside.axis.text = element_text(size = 10))

################################################################################
# SURVIVAL — Kaplan-Meier per cancer type and endpoint
#
# For every cancer type (and Pan_Cancer) with >= 30 samples, samples are split
# at the median pDT into Short_DT / Long_DT, and each of four endpoints is
# tested (log-rank p, Cox hazard ratio of Long_DT vs Short_DT). An endpoint is
# skipped when fewer than 10 samples have it or when it has no events.
################################################################################

# Survival endpoints: time column, status column, and the status value(s)
# that count as an event.
surv_list <- list(
  
  DSS = list(
    time  = "Months.of.disease.specific.survival",
    event = "Disease.specific.Survival.status",
    event_codes = c("1:DEAD WITH TUMOR")
  ),
  
  OS = list(
    time  = "Overall.Survival..Months.",
    event = "Overall.Survival.Status",
    event_codes = c("1:DECEASED")
  ),
  
  PFS = list(
    time  = "Progress.Free.Survival..Months.",
    event = "Progression.Free.Status",
    event_codes = c("1:PROGRESSION")
  ),
  
  DFS = list(
    time  = "Disease.Free..Months.",
    event = "Disease.Free.Status",
    event_codes = c("1:Recurred/Progressed")
  )
)


# One row per cancer type x endpoint
KM_results <- data.frame()
# Same rows, with the 95% CI of the HR and event counts (for Supp. Figure 5a)
KM_ci <- data.frame()

for (tumor in unique(Meta_TCGA$TCGA.PanCanAtlas.Cancer.Type.Acronym)) {
  
  df_tumor <- Meta_TCGA %>% 
    filter(TCGA.PanCanAtlas.Cancer.Type.Acronym == tumor)
  if (nrow(df_tumor) < 30) next
  
  # Median split of pDT within this cancer type
  median_DT <- median(df_tumor$DT_pred, na.rm = TRUE)
  df_tumor$DT_group <- factor(ifelse(df_tumor$DT_pred <= median_DT, "Short_DT", "Long_DT"),levels = c("Short_DT", "Long_DT"))
  
  for (stype in names(surv_list)) {
    
    time_col     <- surv_list[[stype]]$time
    event_col    <- surv_list[[stype]]$event
    event_codes  <- surv_list[[stype]]$event_codes
    
    # Samples with both a time and a status for this endpoint
    df_sub <- df_tumor %>%
      filter(!is.na(.data[[time_col]]), !is.na(.data[[event_col]]))
    
    if (nrow(df_sub) < 10) next
    
    df_sub$event_binary <- df_sub[[event_col]] %in% event_codes
    if (length(unique(df_sub$event_binary)) < 2) next
    
    surv_obj <- Surv(
      time  = df_sub[[time_col]],
      event = df_sub$event_binary
    )
    
    # KM fit, log-rank test, and Cox HR (Long_DT relative to Short_DT)
    fit  <- survfit(surv_obj ~ DT_group, data = df_sub)
    test <- survdiff(surv_obj ~ DT_group, data = df_sub)
    pval <- pchisq(test$chisq, df = length(test$n) - 1, lower.tail = FALSE)
    hr   <- summary(coxph(surv_obj ~ DT_group, data = df_sub))$coefficients[, "exp(coef)"]
    
    # Kaplan-Meier plots are saved in one subfolder per survival type
    # (DSS / OS / PFS / DFS)
    surv_dir <- file.path(OUT_DIR, stype)
    dir.create(surv_dir, showWarnings = FALSE, recursive = TRUE)
    
    g <- ggsurvplot(
      fit, data = df_sub,
      pval = TRUE,
      risk.table = TRUE,
      palette = c("Short DT" = "#DC3220","Long DT"  = "#005AB5"),
      xlab = paste0(stype, " (Months)"),
      ylab = "Survival Probability",
      legend.title = "Predicted DT",
      legend.labs = c("Short DT", "Long DT"),
      ggtheme = theme_classic(base_size = 12)+
        theme(axis.text = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 14, face = "bold"))
    )
    
    ggsave(filename =file.path(surv_dir, paste0("panel_H_KM_", tumor, "_", stype, ".pdf")),
           plot = g$plot, width = 6, height = 6)
    
    # Keep the MESO disease-specific survival plot for Figure 4e
    if(stype == "DSS" & tumor == "MESO"){
      panel_D <- g
    }
    
    KM_results <- dplyr::bind_rows(
      KM_results,
      data.frame(
        Study.ID_short = tumor,
        Survival_Type  = stype,
        N              = nrow(df_sub),
        Median_DT_pred = median_DT,
        p_value        = pval,
        hazard_ratio   = hr
      )
    )
    
    # 95% CI of the same Cox HR, and events per group, for Supp. Figure 5a
    hr_ci <- summary(coxph(surv_obj ~ DT_group, data = df_sub))$conf.int
    KM_ci <- dplyr::bind_rows(
      KM_ci,
      data.frame(
        Study.ID_short = tumor,
        Survival_Type  = stype,
        N              = nrow(df_sub),
        Events         = sum(df_sub$event_binary),
        Events_short   = sum(df_sub$event_binary[df_sub$DT_group == "Short_DT"]),
        Events_long    = sum(df_sub$event_binary[df_sub$DT_group == "Long_DT"]),
        hazard_ratio   = hr_ci[1, "exp(coef)"],
        HR_lower       = hr_ci[1, "lower .95"],
        HR_upper       = hr_ci[1, "upper .95"]
      )
    )
  }
}

# BH-adjust the log-rank p-values across cancer types, separately per endpoint
surv_type_vector <- unique(KM_results$Survival_Type)
for (surv_type in surv_type_vector) {
  
  surv_index <- which(KM_results$Survival_Type == surv_type)
  KM_results[surv_index,"p_adj"] <- p.adjust(KM_results[surv_index,"p_value"], method = "BH")
}

# ── Broken x-axis for the survival q-value bar plots (Fig. 4d, Supp. 5c-e) ───
# Pan_Cancer's -log10(q) is much larger than that of any single cancer type,
# so its bar is shortened: the axis runs normally up to just above the largest
# cancer type, then a break (//, and a white gap in the Pan_Cancer bar), then a
# single tick at the end of the Pan_Cancer bar labelled with its actual value.
# The part of the Pan_Cancer bar beyond the break is not to scale. The break is
# only used when Pan_Cancer exceeds twice the normal axis range; otherwise the
# axis is drawn as before.
#   df:    KM_results rows of one survival type
#   gap:   space between the normal axis and the break, as a fraction of it
#   extra: length of the Pan_Cancer bar beyond the break, as a fraction of it
q_axis <- function(df, gap = 0.10, extra = 0.25) {
  df$neg_log_q <- -log10(df$p_adj)
  is_pan  <- df$Study.ID_short == "Pan-Cancer"
  others  <- df$neg_log_q[!is_pan]
  ticks   <- pretty(c(0, max(others, 1.5)))      # normal axis ticks
  cut     <- max(ticks)                          # end of the normal axis
  pan     <- df$neg_log_q[is_pan]
  df$x_plot <- df$neg_log_q
  if (length(pan) != 1 || !is.finite(pan) || pan <= 2 * cut) {
    return(list(data = df, layers = list()))      # no break needed
  }
  break_at <- cut * (1 + gap)                    # position of the break
  pan_end  <- cut * (1 + gap + extra)            # where the Pan_Cancer bar ends
  df$x_plot[is_pan] <- pan_end
  n_bars   <- nrow(df)                           # Pan_Cancer is the top bar
  half_gap <- 0.012 * pan_end                    # half-width of the white gap
  layers <- list(
    scale_x_continuous(breaks = c(ticks, pan_end),
                       labels = c(ticks, formatC(pan, format = "f", digits = 0)),
                       limits = c(0, pan_end * 1.03),
                       expand = expansion(mult = c(0, 0))),
    # white gap through the Pan_Cancer bar at the break
    annotate("rect", xmin = break_at - half_gap, xmax = break_at + half_gap,
             ymin = n_bars - 0.5, ymax = n_bars + 0.5, fill = "white", colour = NA),
    # break symbol on the x-axis line
    annotate("text", x = break_at, y = -Inf, label = "//", size = 5,
             fontface = "bold", vjust = 0.5),
    coord_cartesian(clip = "off"))
  list(data = df, layers = layers)
}

q_axis_DSS <- q_axis(subset(KM_results, Survival_Type == "DSS"))
# ── Figure 4d — DSS q-value per cancer type, fill = hazard ratio ──────────────
panel_E <- ggplot(q_axis_DSS$data,
                  aes(x= x_plot,
                      y=reorder(x=Study.ID_short, X= - log10(p_adj)),
                      fill=hazard_ratio)) +
  geom_col() +
  theme_classic(base_size=14) +
  ylab(" ")+
  xlab(TeX('$-log_{10}(Q_{value})$'))+
  labs(fill = TeX('HR'))+
  geom_vline(xintercept = 1, linetype = "dashed", colour = "red")+
  theme(legend.position = c(0.8,0.2),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        strip.text = element_text(size = 12),
        ggside.panel.scale = 0.2,
        ggside.axis.text = element_text(size = 10)) +
  q_axis_DSS$layers

# ── Supp. 5a — DSS forest plot: Cox HR (95% CI) per cancer type ──────────────
# The HRs are those of Figure 4d / Table 5.6 (unadjusted Cox model, long vs
# short pDT, Pan_Cancer = pooled median split). Red = log-rank q < 0.1 (the
# q-values of Figure 4d); the right-hand column gives HR (95% CI) and q.
# CIs are drawn within 0.1-10; an x marks an HR that cannot be estimated
# because one pDT group has no events.
forest_x_lim <- c(0.1, 10)
forest_DSS <- merge(KM_ci,
                    KM_results[, c("Study.ID_short", "Survival_Type", "p_adj")],
                    by = c("Study.ID_short", "Survival_Type")) %>%
  filter(Survival_Type == "DSS") %>%
  arrange(hazard_ratio) %>%
  mutate(label      = sprintf("%s  (%d/%d)", Study.ID_short, Events, N),
         label      = factor(label, levels = label),
         estimable  = Events_short > 0 & Events_long > 0,
         hr_txt     = ifelse(estimable,
                             sprintf("%.2f (%.2f-%.2f)   q=%s", hazard_ratio, HR_lower, HR_upper,
                                     formatC(p_adj, format = "g", digits = 2)),
                             sprintf("n.e.   q=%s", formatC(p_adj, format = "g", digits = 2))),
         # no CI bar for an HR that cannot be estimated
         lo_plot    = ifelse(estimable, pmax(HR_lower, forest_x_lim[1]), NA),
         hi_plot    = ifelse(estimable, pmin(HR_upper, forest_x_lim[2]), NA),
         hr_plot    = pmin(pmax(hazard_ratio, forest_x_lim[1]), forest_x_lim[2]),
         sig        = factor(ifelse(p_adj < 0.1, "q < 0.1", "n.s."), levels = c("q < 0.1", "n.s.")))

sup5_panel_A <- ggplot(forest_DSS, aes(y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lo_plot, xmax = hi_plot), width = 0.3, orientation = "y",
                na.rm = TRUE) +
  geom_point(aes(x = hr_plot, colour = sig, shape = estimable), size = 3) +
  geom_text(aes(x = forest_x_lim[2] * 1.5, label = hr_txt, colour = sig),
            hjust = 0, size = 4, show.legend = FALSE) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10)) +
  scale_colour_manual(values = c("q < 0.1" = "#DC3220", "n.s." = "black"),
                      drop = FALSE, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 15, `FALSE` = 4),
                     labels = c(`TRUE` = "estimable", `FALSE` = "not estimable (no events in one group)"),
                     name = NULL) +
  coord_cartesian(xlim = forest_x_lim, clip = "off") +
  labs(x = "Hazard ratio, long vs short pDT (DSS)", y = "Cancer type (events / N)") +
  theme_classic(base_size = 14) +
  theme(plot.margin = margin(5.5, 230, 5.5, 5.5),
        legend.position = "bottom",
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12))

# ── Supp. 5b — DSS forest plot, adjusted for age and T stage ──────────────────
# DSS: Cox model of long vs short pDT
# (median split within each cancer type), adjusted for age (per 10 years) and
# T stage (tumor_stage_num, T0-T4 = 0-4). Pan_Cancer: short/long split at the
# median of each sample's own cancer type, and stratified by cancer type.
# Only patients with age and T stage (complete cases); cancer types without
# usable T stage (non-missing in < 50% or a single value) have no row.
# Grey = unadjusted HR on the same patients, blue = adjusted HR; a square is
# red when its own q < 0.1 (its CI keeps its grey or blue colour). The right-hand
# column gives, for each of the two, HR (95% CI) and q (Wald p of the pDT term,
# BH across cancer types), in the colour of its CI; a red * marks q < 0.1.
# CIs are drawn within 0.1-10.

# cancer type of every sample, to stratify the Pan_Cancer copy
sample_type <- with(subset(Meta_TCGA, TCGA.PanCanAtlas.Cancer.Type.Acronym != "Pan-Cancer"),
                    setNames(TCGA.PanCanAtlas.Cancer.Type.Acronym, Sample.ID))

# Fit a Cox model and return the long-vs-short HR, 95% CI, p and any warning
fit_cox_adj <- function(formula, data) {
  warn <- NA_character_
  fit <- withCallingHandlers(
    coxph(formula, data = data),
    warning = function(w) { warn <<- conditionMessage(w); invokeRestart("muffleWarning") })
  s <- summary(fit)
  list(hr = s$conf.int["DT_groupLong_DT", "exp(coef)"],
       lo = s$conf.int["DT_groupLong_DT", "lower .95"],
       hi = s$conf.int["DT_groupLong_DT", "upper .95"],
       p  = s$coefficients["DT_groupLong_DT", "Pr(>|z|)"],
       n_par = length(coef(fit)),
       warn = warn)
}

adj_DSS <- data.frame()
for (tumor in unique(Meta_TCGA$TCGA.PanCanAtlas.Cancer.Type.Acronym)) {

  df <- Meta_TCGA %>% filter(TCGA.PanCanAtlas.Cancer.Type.Acronym == tumor)
  if (nrow(df) < 30) next
  is_pan <- tumor == "Pan-Cancer"
  df$Cancer_Type <- sample_type[df$Sample.ID]

  # median split: within each cancer type for Pan_Cancer
  if (is_pan) {
    df <- df %>% group_by(Cancer_Type) %>% mutate(median_DT = median(DT_pred)) %>%
      ungroup() %>% as.data.frame()
  } else {
    df$median_DT <- median(df$DT_pred)
  }
  df$DT_group <- factor(ifelse(df$DT_pred <= df$median_DT, "Short_DT", "Long_DT"),
                        levels = c("Short_DT", "Long_DT"))

  # patients with a DSS time and status
  df <- df %>% filter(!is.na(Months.of.disease.specific.survival),
                      !is.na(Disease.specific.Survival.status),
                      Months.of.disease.specific.survival >= 0)
  if (nrow(df) < 10) next
  df$time_        <- df$Months.of.disease.specific.survival
  df$event_binary <- df$Disease.specific.Survival.status %in% surv_list$DSS$event_codes
  if (length(unique(df$event_binary)) < 2 || n_distinct(df$DT_group) < 2) next

  # covariates: T stage must be usable; age is used if available
  df$age_10y <- suppressWarnings(as.numeric(df$Diagnosis.Age)) / 10
  if (mean(!is.na(df$tumor_stage_num)) < 0.5 || n_distinct(na.omit(df$tumor_stage_num)) < 2) next
  covars <- c(if (mean(!is.na(df$age_10y)) >= 0.5) "age_10y", "tumor_stage_num")
  df_cc  <- df %>% filter(if_all(all_of(covars), ~ !is.na(.x)))
  if (nrow(df_cc) < 10 || n_distinct(df_cc$tumor_stage_num) < 2 ||
      n_distinct(df_cc$DT_group) < 2 || length(unique(df_cc$event_binary)) < 2) next

  strat <- if (is_pan) "+ strata(Cancer_Type)" else ""
  unadj <- fit_cox_adj(as.formula(paste("Surv(time_, event_binary) ~ DT_group", strat)), df_cc)
  adj   <- fit_cox_adj(as.formula(paste("Surv(time_, event_binary) ~ DT_group +",
                                   paste(covars, collapse = " + "), strat)), df_cc)

  adj_DSS <- bind_rows(adj_DSS, data.frame(
    Study.ID_short = tumor, N_cc = nrow(df_cc), Events = sum(df_cc$event_binary),
    N_excluded = nrow(df) - nrow(df_cc),
    Covariates = paste(c(if ("age_10y" %in% covars) "age", "T stage"), collapse = " + "),
    HR_unadj = unadj$hr, HR_unadj_lower = unadj$lo, HR_unadj_upper = unadj$hi,
    p_unadj = unadj$p,
    HR_adj = adj$hr, HR_adj_lower = adj$lo, HR_adj_upper = adj$hi, p_adj_model = adj$p))
}
adj_DSS$q_adj_model <- p.adjust(adj_DSS$p_adj_model, method = "BH")
adj_DSS$q_unadj     <- p.adjust(adj_DSS$p_unadj,     method = "BH")

d_adj <- adj_DSS %>%
  arrange(HR_adj) %>%
  mutate(row     = row_number(),
         label   = sprintf("%s  (%d/%d)", Study.ID_short, Events, N_cc),
         txt_unadj = sprintf("%.2f (%.2f-%.2f)   q=%s", HR_unadj, HR_unadj_lower, HR_unadj_upper,
                             formatC(q_unadj, format = "g", digits = 2)),
         txt_adj   = sprintf("%.2f (%.2f-%.2f)   q=%s", HR_adj, HR_adj_lower, HR_adj_upper,
                             formatC(q_adj_model, format = "g", digits = 2)))
long_adj <- bind_rows(
  d_adj %>% transmute(row, model = "Unadjusted", hr = HR_unadj,
                     lo = HR_unadj_lower, hi = HR_unadj_upper, offset = 0.18,
                     square = ifelse(q_unadj < 0.1, "q < 0.1", "Unadjusted")),
  d_adj %>% transmute(row, model = "Adjusted (age + T stage)", hr = HR_adj,
                     lo = HR_adj_lower, hi = HR_adj_upper, offset = -0.18,
                     square = ifelse(q_adj_model < 0.1, "q < 0.1",
                                     "Adjusted (age + T stage)"))) %>%
  mutate(y  = row + offset,
         lo = pmax(lo, forest_x_lim[1]), hi = pmin(hi, forest_x_lim[2]),
         hr = pmin(pmax(hr, forest_x_lim[1]), forest_x_lim[2]),
         model  = factor(model, levels = c("Unadjusted", "Adjusted (age + T stage)")),
         square = factor(square, levels = c("Unadjusted", "Adjusted (age + T stage)",
                                            "q < 0.1")))

sup5_panel_B <- ggplot(long_adj, aes(x = hr, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  # CI bars: grey (unadjusted) or blue (adjusted)
  geom_errorbar(aes(xmin = lo, xmax = hi, colour = model), width = 0.15, orientation = "y") +
  # squares: grey (unadjusted), blue (adjusted), or red for either with q < 0.1
  geom_point(aes(fill = square), size = 3.4, shape = 22, stroke = 0) +
  # right-hand text: unadjusted (grey, upper line) and adjusted (blue, lower
  # line), each with a red * in front when q < 0.1
  annotate("text", x = forest_x_lim[2] * 1.6, y = d_adj$row + 0.18, label = d_adj$txt_unadj,
           colour = "grey45", hjust = 0, size = 3.6) +
  annotate("text", x = forest_x_lim[2] * 1.6, y = d_adj$row - 0.18, label = d_adj$txt_adj,
           colour = "#005AB5", hjust = 0, size = 3.6) +
  annotate("text", x = forest_x_lim[2] * 1.3, y = d_adj$row[d_adj$q_unadj < 0.1] + 0.18,
           label = "*", colour = "#DC3220", hjust = 0, size = 5, fontface = "bold") +
  annotate("text", x = forest_x_lim[2] * 1.3, y = d_adj$row[d_adj$q_adj_model < 0.1] - 0.18,
           label = "*", colour = "#DC3220", hjust = 0, size = 5, fontface = "bold") +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10)) +
  scale_y_continuous(breaks = d_adj$row, labels = d_adj$label, expand = expansion(add = 0.6)) +
  scale_colour_manual(values = c("Unadjusted" = "grey55", "Adjusted (age + T stage)" = "#005AB5"),
                      guide = "none") +
  scale_fill_manual(values = c("Unadjusted" = "grey55", "Adjusted (age + T stage)" = "#005AB5",
                               "q < 0.1" = "#DC3220"),
                    drop = FALSE, name = NULL) +
  coord_cartesian(xlim = forest_x_lim, clip = "off") +
  labs(x = "Hazard ratio, long vs short pDT (DSS, adjusted)",
       y = "Cancer type (events / N, complete cases)") +
  theme_classic(base_size = 14) +
  theme(plot.margin = margin(5.5, 230, 5.5, 5.5),
        legend.position = "bottom",
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12))

# ── Supp. 5c-e — the same for OS, PFS and DFS ─────────────────────────────────
q_axis_OS  <- q_axis(subset(KM_results, Survival_Type == "OS"))
q_axis_PFS <- q_axis(subset(KM_results, Survival_Type == "PFS"))
q_axis_DFS <- q_axis(subset(KM_results, Survival_Type == "DFS"))
sup5_panel_C <- ggplot(q_axis_OS$data,
                      aes(x= x_plot,
                          y=reorder(x=Study.ID_short, X= - log10(p_adj)),
                          fill=hazard_ratio)) +
  geom_col() +
  theme_classic(base_size=14) +
  ylab(" ")+
  xlab(TeX('$-log_{10}(Q_{value})$'))+
  labs(fill = TeX('Hazard Ratio'))+
  geom_vline(xintercept = 1, linetype = "dashed", colour = "red")+
  theme(legend.position = c(0.8,0.2),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        strip.text = element_text(size = 12),
        ggside.panel.scale = 0.2,
        ggside.axis.text = element_text(size = 10)) +
  q_axis_OS$layers

sup5_panel_D <- ggplot(q_axis_PFS$data,
                      aes(x= x_plot,
                          y=reorder(x=Study.ID_short, X= - log10(p_adj)),
                          fill=hazard_ratio)) +
  geom_col() +
  theme_classic(base_size=14) +
  ylab(" ")+
  xlab(TeX('$-log_{10}(Q_{value})$'))+
  labs(fill = TeX('Hazard Ratio'))+
  geom_vline(xintercept = 1, linetype = "dashed", colour = "red")+
  theme(legend.position = c(0.8,0.2),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        strip.text = element_text(size = 12),
        ggside.panel.scale = 0.2,
        ggside.axis.text = element_text(size = 10)) +
  q_axis_PFS$layers

sup5_panel_E <- ggplot(q_axis_DFS$data,
                      aes(x= x_plot,
                          y=reorder(x=Study.ID_short, X= - log10(p_adj)),
                          fill=hazard_ratio)) +
  geom_col() +
  theme_classic(base_size=14) +
  ylab(" ")+
  xlab(TeX('$-log_{10}(Q_{value})$'))+
  labs(fill = TeX('Hazard Ratio'))+
  geom_vline(xintercept = 1, linetype = "dashed", colour = "red")+
  theme(legend.position = c(0.8,0.2),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        strip.text = element_text(size = 12),
        ggside.panel.scale = 0.2,
        ggside.axis.text = element_text(size = 10)) +
  q_axis_DFS$layers

################################################################################
# ASSEMBLE FIGURE 4
################################################################################

# Panel e: stack the MESO KM curves over their risk table (4:1).
km_D   <- panel_D$plot  + labs(x = NULL, y = "Survival \n Probability")+ theme(legend.position = "none")
tab_D  <- panel_D$table + labs(x = "Time (months)", y = NULL, title = NULL)

pan_d <- wrap_elements(full = km_D) + tab_D +
  plot_layout(
    design = paste(c("A",
                     "B"),
                   collapse = "\n"),
    heights = c(4, 1))

# Design areas are filled in alphabetical order by the plots in the order they
# are added: panel_A -> A, panel_C -> B, panel_B -> C, pan_d -> D,
# panel_E -> E. Tags go in the same order, so c("a","b","c","e","d") prints:
#   row 1     a (full width)
#   rows 2-3  b (left, both rows) | c over e (middle) | d (right, both rows)
fig4 <- wrap_elements(full = panel_A) +
  wrap_elements(full = panel_C) +
  wrap_elements(full = panel_B) +
  wrap_elements(full = pan_d) +
  wrap_elements(full = panel_E) +
  plot_layout(
    design = paste(c("AAAAAA",
                     "BBCCEE",
                     "BBDDEE"),
                   collapse = "\n"),
    heights = c(0.8, 1, 1)) + 
  plot_annotation(tag_levels = list(c("a","b","c","e","d"))) & 
  theme(plot.tag = element_text(size = 15, face  = "bold"))
ggsave(file.path(OUT_DIR, "Figure_4.pdf"), fig4, width = 10, height = 10)
ggsave(file.path(OUT_DIR, "Figure_4.png"), fig4, width = 10, height = 10)

################################################################################
# ASSEMBLE SUPPLEMENTAL FIGURE 4
#
# Row 1: a (pDT by cancer type). Row 2: stage examples b-d. Row 3: stage
# q-values e-g. (The OS / PFS / DFS q-values, formerly h-j, are now
# Supplemental Figure 5c-e.) The height is reduced from 20 to 16.5 so that the
# remaining rows keep the same size as before (20 x 4.8 / 5.8).
################################################################################
sup_fig4 <-  wrap_elements(full = sup_panel_A) +
  wrap_elements(full = sup_panel_B) +  wrap_elements(full = sup_panel_D) +  wrap_elements(full = sup_panel_F) +
  wrap_elements(full = sup_panel_C) +  wrap_elements(full = sup_panel_E) +  wrap_elements(full = sup_panel_G) +
  plot_layout(
    design = paste(c("AAA",
                     "BCD",
                     "EFG"),
                   collapse = "\n"),
    heights = c(1, 1, 1.8)) + 
  plot_annotation(tag_levels = "a")
sup_fig4
ggsave(file.path(OUT_DIR, "sup_Figure_4.pdf"), sup_fig4, width = 12, height = 16.5)
ggsave(file.path(OUT_DIR, "sup_Figure_4.png"), sup_fig4, width = 12, height = 16.5)

################################################################################
# ASSEMBLE SUPPLEMENTAL FIGURE 5
#
# Panels c-e: per cancer type q-values of the OS, PFS and DFS log-rank tests
# (formerly Supplemental Figure 4h-j), same size as before (row height 20 x
# 2 / 5.8, i.e. about 7).
#
# Left column: a (DSS forest plot of Figure 4d / Table 5.6) above b (DSS
# forest plot adjusted for age and T stage), 12 inches wide.
# Right column: c, d, e (OS, PFS, DFS q-values) stacked, 7 inches tall each
# and 4.06 inches wide, i.e. the same proportions as in the old Supp. Figure 4
# (4 x 6.9 inches).
################################################################################
sup5_left  <- wrap_elements(full = sup5_panel_A) /
  wrap_elements(full = sup5_panel_B) +
  plot_layout(heights = c(10, 9))
sup5_right <- wrap_elements(full = sup5_panel_C) /
  wrap_elements(full = sup5_panel_D) /
  wrap_elements(full = sup5_panel_E)
sup_fig5 <- (sup5_left | sup5_right) +
  plot_layout(widths = c(12, 4.06)) +
  plot_annotation(tag_levels = "a")
sup_fig5
ggsave(file.path(OUT_DIR, "sup_Figure_5.pdf"), sup_fig5, width = 16.06, height = 21)
ggsave(file.path(OUT_DIR, "sup_Figure_5.png"), sup_fig5, width = 16.06, height = 21)

################################################################################
# ANALYSES FOR SUPPLEMENTAL TABLE 5.8
#
# Table 5.8: pDT as a continuous variable (untransformed; HR per one-unit
#   increase in pDT), per cancer type and endpoint, unadjusted and adjusted
#   for age (per 10 years) and T stage (tumor_stage_num, T0-T4 = 0-4).
#   Pan_Cancer models are stratified by cancer type.
# Sample selection, endpoints and covariate rules are the same as for
# Supplemental Figure 5b (adj_DSS).
################################################################################

# Fit a Cox model and return HR, 95% CI and Wald p of one term
cox_term <- function(formula, data, term) {
  warn <- NA_character_
  fit <- withCallingHandlers(
    coxph(formula, data = data),
    warning = function(w) { warn <<- conditionMessage(w); invokeRestart("muffleWarning") })
  s <- summary(fit)
  c(HR = s$conf.int[term, "exp(coef)"], lower = s$conf.int[term, "lower .95"],
    upper = s$conf.int[term, "upper .95"], p = s$coefficients[term, "Pr(>|z|)"])
}

# Patients of one cohort with a time and status for one endpoint, with the
# Short/Long pDT split (within cancer type for Pan_Cancer if within_type),
# or NULL if the endpoint cannot be analysed
endpoint_data <- function(tumor, stype, within_type = TRUE) {
  df <- Meta_TCGA %>% filter(TCGA.PanCanAtlas.Cancer.Type.Acronym == tumor)
  if (nrow(df) < 30) return(NULL)
  df$Cancer_Type <- sample_type[df$Sample.ID]
  if (tumor == "Pan-Cancer" && within_type) {
    df <- df %>% group_by(Cancer_Type) %>% mutate(median_DT = median(DT_pred)) %>%
      ungroup() %>% as.data.frame()
  } else {
    df$median_DT <- median(df$DT_pred)
  }
  df$DT_group <- factor(ifelse(df$DT_pred <= df$median_DT, "Short_DT", "Long_DT"),
                        levels = c("Short_DT", "Long_DT"))
  time_col  <- surv_list[[stype]]$time
  event_col <- surv_list[[stype]]$event
  df <- df %>% filter(!is.na(.data[[time_col]]), !is.na(.data[[event_col]]),
                      .data[[time_col]] >= 0)
  if (nrow(df) < 10) return(NULL)
  df$time_        <- df[[time_col]]
  df$event_binary <- df[[event_col]] %in% surv_list[[stype]]$event_codes
  if (length(unique(df$event_binary)) < 2 || n_distinct(df$DT_group) < 2) return(NULL)
  df$age_10y <- suppressWarnings(as.numeric(df$Diagnosis.Age)) / 10
  df
}

# Complete cases for the age + T-stage adjustment, or NULL if T stage is not
# usable (non-missing in < 50% or a single value) or too few patients remain
complete_cases <- function(df) {
  if (mean(!is.na(df$tumor_stage_num)) < 0.5 || n_distinct(na.omit(df$tumor_stage_num)) < 2)
    return(NULL)
  covars <- c(if (mean(!is.na(df$age_10y)) >= 0.5) "age_10y", "tumor_stage_num")
  cc <- df %>% filter(if_all(all_of(covars), ~ !is.na(.x)))
  if (nrow(cc) < 10 || n_distinct(cc$tumor_stage_num) < 2 || n_distinct(cc$DT_group) < 2 ||
      length(unique(cc$event_binary)) < 2) return(NULL)
  attr(cc, "covars") <- covars
  cc
}

surv_formula <- function(exposure, covars = character(0), strat = FALSE)
  as.formula(paste("Surv(time_, event_binary) ~", paste(c(exposure, covars), collapse = " + "),
                   if (strat) "+ strata(Cancer_Type)" else ""))

## ---- Table 5.8: continuous pDT ----------------------------------------------
cont_results <- data.frame()
for (tumor in unique(Meta_TCGA$TCGA.PanCanAtlas.Cancer.Type.Acronym)) {
  for (stype in names(surv_list)) {
    df <- endpoint_data(tumor, stype)
    if (is.null(df)) next
    strat <- tumor == "Pan-Cancer"
    u  <- cox_term(surv_formula("DT_pred", strat = strat), df, "DT_pred")
    cc <- complete_cases(df)
    a  <- if (!is.null(cc)) cox_term(surv_formula("DT_pred", attr(cc, "covars"), strat), cc, "DT_pred")
          else c(HR = NA, lower = NA, upper = NA, p = NA)
    q  <- quantile(df$DT_pred, c(0.25, 0.5, 0.75), names = FALSE)
    cont_results <- bind_rows(cont_results, data.frame(
      Study.ID_short = tumor, Survival_Type = stype,
      N = nrow(df), Events = sum(df$event_binary),
      pDT_median = q[2], pDT_Q1 = q[1], pDT_Q3 = q[3],
      HR = u[["HR"]], HR_lower = u[["lower"]], HR_upper = u[["upper"]], p = u[["p"]],
      N_cc = if (!is.null(cc)) nrow(cc) else NA, Events_cc = if (!is.null(cc)) sum(cc$event_binary) else NA,
      HR_adj = a[["HR"]], HR_adj_lower = a[["lower"]], HR_adj_upper = a[["upper"]], p_adj = a[["p"]]))
  }
}
# BH across cancer types within each endpoint (adjusted: rows with a model)
cont_results <- cont_results %>%
  group_by(Survival_Type) %>%
  mutate(q     = p.adjust(p, method = "BH"),
         q_adj = { qq <- rep(NA_real_, n()); ok <- !is.na(p_adj)
                   qq[ok] <- p.adjust(p_adj[ok], method = "BH"); qq }) %>%
  ungroup() %>% as.data.frame()

################################################################################
# SUPPLEMENTAL TABLE 5 (sup_tab5.xlsx)
################################################################################

## ---- Contents sheet, written by hand ---------------------------------------
## One row per sheet. Edit the text freely; keep the two vectors the same
## length, or data.frame() will fail or recycle values.
sup_tab_5_contents <- data.frame(
  sheet = c(
    "Table 5.1",
    "Table 5.2",
    "Table 5.3",
    "Table 5.4",
    "Table 5.5",
    "Table 5.6",
    "Table 5.7",
    "Table 5.8"),
  
  description = c(
    "Clinical data for the TCGA cohort.",
    "Spearman Correlation results of Mutation Count vs predicted doubling time.",
    "Spearman Correlation results of Cancer Tumor Stage vs predicted doubling time.",
    "Spearman Correlation results of Neoplasm Disease Stage vs predicted doubling time.",
    "Spearman Correlation results of Lymph Node Stage vs predicted doubling time.",
    paste("Survival analysis of samples with long vs short predicted doubling time (median split within each",
          "tumor type), per tumor type and endpoint (DSS, OS, PFS, DFS): hazard ratio (long vs short) with 95% CI",
          "from a Cox model; p- and q-values from log-rank tests, Benjamini-Hochberg-corrected across tumor types",
          "within each endpoint."),
    paste("Disease-specific survival (DSS) of samples with long vs short predicted doubling time, Cox models",
          "adjusted for age (per 10 years) and T stage (T0-T4 as a numerical covariate), on patients with both",
          "covariates available (Supplemental Figure 5b). Pan_Cancer: median split within each tumor type and",
          "model stratified by tumor type. Unadjusted = same patients without covariates. q-values: Wald test,",
          "Benjamini-Hochberg-corrected across tumor types.",
          "Tumor types without usable T stage (non-missing in < 50% of patients or a single value) are not",
          "included."),
    paste("Survival analysis with pDT as a continuous variable (untransformed; hazard ratio per one-unit increase",
          "in pDT), per tumor type and endpoint. Unadjusted: all patients with the endpoint. Adjusted: age (per 10",
          "years) and T stage (T0-T4 as a numerical covariate), patients with both covariates available; empty for",
          "tumor types without usable T stage. Pan_Cancer models are stratified by tumor type. p-values: Wald test;",
          "q-values: Benjamini-Hochberg-corrected across tumor types within each endpoint. pDT median and",
          "interquartile range describe the pDT values in the unadjusted model.")),
  
  
  stringsAsFactors = FALSE)

sup_tab_5_contents <- sup_tab_5_contents[, c("sheet", "description")]

# Table 5.1 — clinical data with pDT, one row per sample (Pan_Cancer copies
# removed)
sup_tab_5_1 <- subset(Meta_TCGA, TCGA.PanCanAtlas.Cancer.Type.Acronym != "Pan-Cancer")[,c("TCGA.PanCanAtlas.Cancer.Type.Acronym","Patient.ID","Sample.ID","Diagnosis.Age",
                                                                                          "DT_pred",
                                                                                          "Mutation.Count",
                                                                                          "Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code",
                                                                                          "Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code",
                                                                                          "American.Joint.Committee.on.Cancer.Tumor.Stage.Code",
                                                                                          "Months.of.disease.specific.survival",
                                                                                          "Disease.specific.Survival.status",
                                                                                          "Overall.Survival..Months.",
                                                                                          "Overall.Survival.Status",
                                                                                          "Progress.Free.Survival..Months.",
                                                                                          "Progression.Free.Status",
                                                                                          "Disease.Free..Months.",
                                                                                          "Disease.Free.Status")]

colnames(sup_tab_5_1) <- c("TCGA cohort", "Patient ID","Sample ID","Diagnosis Age",
                           "pDT",
                           "Mutation Count",
                           "Neoplasm Disease Stage American Joint Committee on Cancer Code",
                           "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code",
                           "American Joint Committee on Cancer Tumor Stage Code",
                           "Months of disease specific survival",
                           "Disease specific Survival status",
                           "Overall Survival Months",
                           "Overall Survival Status",
                           "Progress Free Survival Months",
                           "Progression Free Status",
                           "Disease Free Months",
                           "Disease Free Status")

# Tables 5.2-5.5 — the per-cancer-type correlation results (incl. Pan_Cancer)
sup_tab_5_2 <- cor_results
colnames(sup_tab_5_2) <- c("TCGA cohort","sample count","correlation","p_value","q_value")

sup_tab_5_3 <- cor_results_tumor
colnames(sup_tab_5_3) <- c("TCGA cohort","sample count","correlation","p_value","q_value")

sup_tab_5_4 <- cor_results_disease
colnames(sup_tab_5_4) <- c("TCGA cohort","sample count","correlation","p_value","q_value")

sup_tab_5_5 <- cor_results_lymph_node
colnames(sup_tab_5_5) <- c("TCGA cohort","sample count","correlation","p_value","q_value")

# Table 5.6 — survival results, one row per cancer type x endpoint, with the
# 95% CI of the hazard ratio (from KM_ci) as the last two columns
sup_tab_5_6 <- (KM_results)
colnames(sup_tab_5_6) <- c("TCGA cohort","Survival Type","sample count","Median_pDT","p_value","hazard_ratio","q_value")
rownames(sup_tab_5_6) <- NULL
ci_row <- match(paste(KM_results$Study.ID_short, KM_results$Survival_Type),
                paste(KM_ci$Study.ID_short, KM_ci$Survival_Type))
sup_tab_5_6$"HR 95% CI lower" <- KM_ci$HR_lower[ci_row]
sup_tab_5_6$"HR 95% CI upper" <- KM_ci$HR_upper[ci_row]
# an HR that cannot be estimated (no events in one group, e.g. PCPG DSS) has an
# infinite upper limit; leave that cell empty
sup_tab_5_6$"HR 95% CI upper"[!is.finite(sup_tab_5_6$"HR 95% CI upper")] <- NA

# Table 5.7 — DSS, Cox models adjusted for age and T stage (Supp. Figure 5b)
sup_tab_5_7 <- adj_DSS %>%
  arrange(p_adj_model) %>%
  transmute(`TCGA cohort`                = Study.ID_short,
            `sample count`               = N_cc,
            `events`                     = Events,
            `samples without age or T stage` = N_excluded,
            `covariates`                 = Covariates,
            `hazard_ratio unadjusted`    = HR_unadj,
            `HR unadjusted 95% CI lower` = HR_unadj_lower,
            `HR unadjusted 95% CI upper` = HR_unadj_upper,
            `p_value unadjusted`         = p_unadj,
            `q_value unadjusted`         = q_unadj,
            `hazard_ratio adjusted`      = HR_adj,
            `HR adjusted 95% CI lower`   = HR_adj_lower,
            `HR adjusted 95% CI upper`   = HR_adj_upper,
            `p_value adjusted`           = p_adj_model,
            `q_value adjusted`           = q_adj_model)

# Table 5.8 — continuous pDT, per cancer type and endpoint
sup_tab_5_8 <- cont_results %>%
  arrange(Survival_Type, p) %>%
  transmute(`TCGA cohort`                = Study.ID_short,
            `Survival Type`              = Survival_Type,
            `sample count`               = N,
            `events`                     = Events,
            `pDT median`                 = pDT_median,
            `pDT Q1`                     = pDT_Q1,
            `pDT Q3`                     = pDT_Q3,
            `hazard_ratio per unit pDT`  = HR,
            `HR 95% CI lower`            = HR_lower,
            `HR 95% CI upper`            = HR_upper,
            `p_value`                    = p,
            `q_value`                    = q,
            `sample count adjusted`      = N_cc,
            `events adjusted`            = Events_cc,
            `hazard_ratio per unit pDT adjusted` = HR_adj,
            `HR adjusted 95% CI lower`   = HR_adj_lower,
            `HR adjusted 95% CI upper`   = HR_adj_upper,
            `p_value adjusted`           = p_adj,
            `q_value adjusted`           = q_adj)

# Write the workbook: Contents first, then one sheet per table.
sup_tab_5 <- createWorkbook()

addWorksheet(sup_tab_5, "Contents")
addWorksheet(sup_tab_5, "sup_tab_5.1")
addWorksheet(sup_tab_5, "sup_tab_5.2")
addWorksheet(sup_tab_5, "sup_tab_5.3")
addWorksheet(sup_tab_5, "sup_tab_5.4")
addWorksheet(sup_tab_5, "sup_tab_5.5")
addWorksheet(sup_tab_5, "sup_tab_5.6")
addWorksheet(sup_tab_5, "sup_tab_5.7")
addWorksheet(sup_tab_5, "sup_tab_5.8")

writeData(sup_tab_5, sheet = "Contents", x = sup_tab_5_contents)
writeData(sup_tab_5, sheet = "sup_tab_5.1", x = sup_tab_5_1)
writeData(sup_tab_5, sheet = "sup_tab_5.2", x = sup_tab_5_2)
writeData(sup_tab_5, sheet = "sup_tab_5.3", x = sup_tab_5_3)
writeData(sup_tab_5, sheet = "sup_tab_5.4", x = sup_tab_5_4)
writeData(sup_tab_5, sheet = "sup_tab_5.5", x = sup_tab_5_5)
writeData(sup_tab_5, sheet = "sup_tab_5.6", x = sup_tab_5_6)
writeData(sup_tab_5, sheet = "sup_tab_5.7", x = sup_tab_5_7)
writeData(sup_tab_5, sheet = "sup_tab_5.8", x = sup_tab_5_8)

saveWorkbook(sup_tab_5, file.path(OUT_DIR, "sup_tab5.xlsx"), overwrite = TRUE)

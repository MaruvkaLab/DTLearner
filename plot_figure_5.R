################################################################################
# plot_figure_5.R
# Figure 5, Supplemental Figure 6 and Supplemental Tables 6 and 7: predicted
# doubling time (pDT) and drug response.
#
#   Patient cohort  pre-treatment pDT in responders vs non-responders, per
#                   cancer type x regimen (Wilcoxon rank-sum)
#   DepMap / CCLE   Spearman correlation of doubling time with drug effect,
#                   per compound (PRISM Repurposing 24Q2)
#
# Run after predict_drug_cohort.R and predict_ccle.R, which write the pDT
# tables read here.
#
# WHAT THIS SCRIPT WRITES (all in WORK_DIR/figure_5/)
#   Figure_5.pdf / .png          panels a-d
#   sup_Figure_6.pdf / .png      panels a-c
#   sup_tab6.xlsx                patient cohort: clinical data and Wilcoxon
#                                results
#   sup_tab7.xlsx                DepMap: drug correlations and per-mechanism
#                                summaries, with and without lung lines
#
# Figure 5 panels
#   a  pDT, responders vs non-responders, paclitaxel + doxorubicin
#   b  volcano of the chemotherapy cohorts: log2(mean pDT R / NR) vs p
#   c  CCLE: DT vs effect of one example compound (BRD-K18518344)
#   d  volcano of DT vs drug effect, all compounds
# Supplemental Figure 5 panels
#   a-c  as panel b, for immune checkpoint blockade, targeted therapy and
#        chemo + targeted cohorts
################################################################################

################################################################################
# WORKING DIRECTORY — defaults to the folder this script is in
#
# Same WORK_DIR as the other figure and predict_*.R scripts. Inputs are read
# from WORK_DIR/input_files/ and from the predict_*.R output folders; outputs
# go to WORK_DIR/figure_5/.
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
OUT_DIR <- file.path(WORK_DIR, "figure_5")

# ── Install missing packages ───────────────────────────────────────────────────
# Installs only the packages that are not already installed; packages that are
# present are left as they are (no updates). All are on CRAN.
CRAN_PKGS <- c("dplyr", "ggplot2", "ggpubr", "latex2exp", "ggrepel",
               "patchwork", "openxlsx", "data.table", "stringr", "tibble")

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

# data.table, stringr, tibble: used via :: only, so installed but not attached.
library(dplyr)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(ggrepel)
library(patchwork)
library(openxlsx)


# ── format_p_value() ──────────────────────────────────────────────────────────
# Formats p-values as plotmath strings, e.g. "italic(P) == 1.23 %*% 10^-04".
# Use with parse = TRUE (stat_cor / stat_compare_means labels below). Values
# below `threshold` are shown with no decimal places on the mantissa.
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
# PATIENT DRUG-RESPONSE COHORT — clinical data
################################################################################

# One row per sample: dataset, sample, timepoint, regimen, therapy class and
# binary response.
clinical_data <- read.delim(file.path(WORK_DIR, "input_files", "patient_drug_cohort", "clinical.tsv"))

# Expression file for each dataset, and a sample ID that is unique across
# datasets (<resource>_<dataset_id>_<sample_id>, as in predict_drug_cohort.R).
clinical_data$file_name <- paste0(clinical_data$resource, "__" ,clinical_data$dataset_id,".tsv.gz")
clinical_data$unique_sample_id <- paste0(clinical_data$resource, "_" ,clinical_data$dataset_id,"_",clinical_data$sample_id)

# Expression data type, parsed from the platform description. Later rules
# overwrite earlier ones when a description matches more than one.
clinical_data$data_type <- NA
clinical_data[which(grepl(pattern = "raw counts", x = clinical_data$platform, ignore.case = TRUE)),"data_type"] <- "raw counts"
clinical_data[which(grepl(pattern = "RPKM", x = clinical_data$platform, ignore.case = TRUE, fixed = TRUE)),"data_type"] <- "RPKM"
clinical_data[which(grepl(pattern = "DESeq", x = clinical_data$platform, ignore.case = TRUE)),"data_type"] <- "DESeq"
clinical_data[which(grepl(pattern = "TPM", x = clinical_data$platform, ignore.case = TRUE)),"data_type"] <- "TPM"
clinical_data[which(grepl(pattern = "TIGER", x = clinical_data$platform, ignore.case = TRUE)),"data_type"] <- "TIGER_TPM"
clinical_data[which(grepl(pattern = "CPM", x = clinical_data$platform, ignore.case = TRUE)),"data_type"] <- "CPM"
clinical_data[which(grepl(pattern = "microarray", x = clinical_data$platform, ignore.case = TRUE)),"data_type"] <- "microarray"
clinical_data[which(grepl(pattern = "Ion\\ ", x = clinical_data$platform, ignore.case = FALSE)),"data_type"] <- "Ion"

# Targeted therapies are reclassified with regimen_classification.csv (other
# therapy classes are kept as given), and regimen names are mapped to one
# canonical name per regimen with regimens_deduplicated.csv.
regimen_classification <- read.csv(file.path(WORK_DIR, "input_files", "patient_drug_cohort", "regimen_classification.csv"))
idx <- match(
  clinical_data$drug_or_regimen,
  regimen_classification$regimen
)

regimen_classification <- regimen_classification[idx, ]
clinical_data$therapy_class_updated <- ifelse(test = clinical_data$therapy_class == "targeted", yes = regimen_classification$classification, no = clinical_data$therapy_class)
clinical_data$drug_or_regimen_updated <- read.csv(file.path(WORK_DIR, "input_files", "patient_drug_cohort", "regimens_deduplicated.csv"), row = 2)[clinical_data$drug_or_regimen,"Canonical_name"]

## The 100-gene prediction signature: the same list, in the same order, as
## SIG_100_GENES in the Figure 2 scripts, Supplemental Table 3.5 and the
## predict_*.R scripts.
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

# Read the signature genes + MKI67 from each dataset's harmonised expression
# file and combine them into one genes x samples table. Files that fail to read
# (or lack a requested gene) are printed and skipped.
file_list_index <- unique(clinical_data[,c("file_name","data_type","resource","dataset_id")])
rownames(file_list_index) <- file_list_index$file_name

rna_tab_pat_drag_comb <- data.frame(row.names = c(GENE_LIST,"MKI67"))

for (file in rownames(file_list_index)) {
  
  tryCatch(
    {
      rna_tab_pat_drag <- tibble::column_to_rownames(data.table::fread(file.path(WORK_DIR, "input_files", "patient_drug_cohort", "harmonized", file),
                                                                       header = TRUE, check.names = FALSE),
                                                     var = "gene_id")[c(GENE_LIST,"MKI67"),,drop = F]
      colnames(rna_tab_pat_drag) <- paste0(file_list_index[file,"resource"], "_" ,file_list_index[file,"dataset_id"],"_",colnames(rna_tab_pat_drag))
      rna_tab_pat_drag_comb <- cbind(rna_tab_pat_drag_comb, rna_tab_pat_drag)
    },
    error = function(e) {
      print(file)
    }
  )
}

# pDT from predict_drug_cohort.R
pDT_results_pat <- read.csv(file.path(WORK_DIR, "drug_cohort_pred", "pDT_results.csv"), row = 1)

# Number of the 100 signature genes (the first 100 rows) that are measured in
# each sample, and the sample's pDT.
clinical_data$valid_gene_num <- colSums(as.data.frame(sapply(rna_tab_pat_drag_comb[c(1:100),],function(x) ifelse(test =  is.na(x) , yes = 0, no = 1))))[clinical_data$unique_sample_id]
clinical_data$DT_pred <- pDT_results_pat[clinical_data$unique_sample_id,"pDT_ensemble"]

# Keep pre-treatment samples with more than 60 of the 100 signature genes.
clinical_data <- subset(clinical_data, timepoint == "pre")
clinical_data <- subset(clinical_data, valid_gene_num > 60)

# Keep regimens with more than n_respon responders AND more than n_respon
# non-responders.
n_respon <- 15
coh_treat_index <- rownames(subset(as.data.frame.matrix(table(clinical_data[,c("drug_or_regimen_updated","response_binary")])),
                                   (`non-responder` > n_respon) & (responder > n_respon)))

clinical_data <-subset(clinical_data, drug_or_regimen_updated %in% coh_treat_index)

# MKI67 expression, a cancer type x regimen key, and the short regimen label
# used on the plots.
clinical_data$MKI67 <- as.vector(unlist(rna_tab_pat_drag_comb["MKI67",clinical_data$unique_sample_id]))
clinical_data$type_drag <- paste0(clinical_data$cancer_type_harmonized, "_", clinical_data$drug_or_regimen_updated)
clinical_data$short_name_no_cancer <- read.csv(file.path(WORK_DIR, "input_files", "patient_drug_cohort", "regimen_short_names.csv"), row =1)[clinical_data$type_drag,c("Short_name_no_cancer")]

################################################################################
# WILCOXON TEST PER CANCER TYPE x REGIMEN
#
# For each cohort with more than n_respon responders and non-responders (with
# MKI67 measured): Wilcoxon rank-sum test of pDT by response, and
# Mean_diff = log2(mean pDT responders / mean pDT non-responders).
################################################################################
wilcox_results_dr <- data.frame()
drug_list_2 <- unique(clinical_data$type_drag)

for (coh2 in drug_list_2) {
  
  df <- subset(clinical_data, type_drag == coh2)
  
  if((dim(subset(df, complete.cases(MKI67) & response_binary == "non-responder"))[1] > n_respon) & (dim(subset(df, complete.cases(MKI67) & response_binary == "responder"))[1] > n_respon)){
    test <- wilcox.test(formula = DT_pred ~ response_binary, data = df)
    
    wilcox_results_dr[coh2,"statistic"] <- test$statistic
    wilcox_results_dr[coh2,"p_value"] <- test$p.value
    wilcox_results_dr[coh2,"Mean_diff"] <- log2(mean(subset(df, response_binary == "responder")$DT_pred, na.rm = TRUE)/
                                                  mean(subset(df, response_binary == "non-responder")$DT_pred, na.rm = TRUE))
    wilcox_results_dr[coh2,"num_responder"] <- dim(subset(df, response_binary == "responder"))[1]
    wilcox_results_dr[coh2,"num_non-responder"] <- dim(subset(df, response_binary == "non-responder"))[1]
  }
}

# BH q-values across cohorts. q < 0.1 cohorts are labelled by direction:
# high_dt_responder (responders have the longer pDT) or high_dt_non-responder.
wilcox_results_dr$q_value<-p.adjust(wilcox_results_dr$p_value,method = "BH")
wilcox_results_dr$significance <- "not_significant"
wilcox_results_dr[which((wilcox_results_dr$q_value < 0.1) & ((wilcox_results_dr$Mean_diff) >= 0)),"significance"] <- "high_dt_responder"
wilcox_results_dr[which((wilcox_results_dr$q_value < 0.1) & ((wilcox_results_dr$Mean_diff) <= 0)),"significance"] <- "high_dt_non-responder"
# Add cancer type, regimen, therapy class and short name for each cohort.
wilcox_results_dr <- cbind(wilcox_results_dr,
                           tibble::column_to_rownames(dplyr::distinct(clinical_data[,c("type_drag",
                                                                                       "cancer_type_harmonized",
                                                                                       "drug_or_regimen_updated",
                                                                                       "therapy_class_updated",
                                                                                       "short_name_no_cancer")]),
                                                      var = "type_drag")[rownames(wilcox_results_dr),,drop = F])
# Plot label: the short name for significant cohorts, NA otherwise.
wilcox_results_dr$sig_column <- NA
wilcox_results_dr[which((wilcox_results_dr$q_value < 0.1)),"sig_column"] <- (wilcox_results_dr[which((wilcox_results_dr$q_value < 0.1)),"short_name_no_cancer"])

# Chemotherapy cohorts for panel b. Cohorts with p = 0 (infinite -log10 p) are
# dropped. The labels are reset to the 3 lowest and 3 highest Mean_diff
# cohorts, whatever their significance.
wilcox_results_chemo <- subset(x = wilcox_results_dr, (!(is.infinite(-log10(p_value))) & therapy_class_updated == "chemo"))
index_low <- (base::sort(wilcox_results_chemo$Mean_diff, decreasing = FALSE, index.return = TRUE)$ix[c(1:3)])
index_high <- (base::sort(wilcox_results_chemo$Mean_diff, decreasing = TRUE, index.return = TRUE)$ix[c(1:3)])
wilcox_results_chemo$sig_column <- NA
wilcox_results_chemo[c(index_low,index_high),"sig_column"] <-wilcox_results_chemo[c(index_low,index_high),"short_name_no_cancer"]

lbl_chemo <- wilcox_results_chemo[c(index_low, index_high), "short_name_no_cancer"]
lbl_chemo <- gsub("([/,;+])", "\\1 ", lbl_chemo)          # allow breaks after / , ; +
lbl_chemo <- gsub("\\s+", " ", trimws(lbl_chemo))

wilcox_results_chemo[c(index_low, index_high), "sig_column"] <-
  stringr::str_wrap(
    lbl_chemo,
    width = 40   # characters per label line
  )

################################################################################
# FIGURE 5 PANELS
################################################################################

# ── Figure 5a — pDT by response, paclitaxel + doxorubicin cohort ──────────────
panel_A <- ggplot(data = subset(clinical_data, drug_or_regimen_updated == "paclitaxel,doxorubicin"),
                  aes(x = response_binary,y = log10(DT_pred),colour = response_binary)) +
  geom_violin(width = 0.8,trim = TRUE,scale = "width") +
  geom_jitter()+
  geom_boxplot(width = 0.2,outlier.shape = NA,colour = "black") +
  xlab(" ") +
  stat_compare_means(na.rm = TRUE, size=6, aes(label = after_stat(format_p_value(..p..))), parse = TRUE, vjust = 0.4, hjust = 0.5, label.y = 1.6)+
  ylab(TeX('$log_{10}(pDT)$'))+
  scale_colour_manual(labels = c("responder" = "Responder", "non-responder" = "Non Responder"),
                      values = c("responder" = "#DC3220","non-responder" = "#005AB5"))+
  scale_x_discrete(labels = c("Non \n Responder", "Responder")) +
  scale_y_continuous(limits = c(1.35,1.61),breaks = seq(1.35,1.6,0.05)) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20, face = "bold"),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title = element_text(size=20, face = "bold")
  )

# ── Figure 5b — chemotherapy cohorts: log2(R/NR) vs -log10(p) ─────────────────
panel_B <- ggplot(data = wilcox_results_chemo,
                  aes(x = (Mean_diff), y = -log10(p_value), colour = significance))+
  geom_point(aes(alpha = significance))+
  xlab(TeX('$log_{2}(R/NR)$'))+
  ylab(TeX('$-log_{10}(P_{value})$'))+
  geom_text_repel(aes(label = wilcox_results_chemo$sig_column,
                      colour = wilcox_results_chemo$significance),
                  bg.color = "white", bg.r = 0.2,
                  size = 6.5,box.padding = 1,
                  max.overlaps = Inf,point.size = 2, segment.alpha = 0.5)+
  scale_colour_manual(values = c("not_significant" = "#808080", "high_dt_responder" = "#DC3220","high_dt_non-responder" = "#005AB5"))+
  scale_alpha_manual(values = c("not_significant" = 0.5, "high_dt_responder" = 1,"high_dt_non-responder" = 1), guide="none")+
  scale_x_continuous(limits = c(-0.15,0.1),breaks = c(-0.15,0,0.1)) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=20, face = "bold"),
        axis.title = element_text(size=20, face = "bold")
  )

################################################################################
# DEPMAP / CCLE — doubling time vs drug effect
#
# Each cell line's doubling time is its measured Avg_2DT where available,
# otherwise its pDT from predict_ccle.R (Final_DT). The analysis is run on all
# cell lines and again without lung lines.
################################################################################

Models <- read.csv(file.path(WORK_DIR, "input_files", "ccle", "Models_DT.csv"), header = TRUE, row.names = 1)
rna_raw <- read.csv(file.path(WORK_DIR, "input_files", "ccle", "DT_Expression_CCLE_Data.csv"), row = 2)[, -c(1, 3)]
Models$mki67 <- rna_raw[rownames(Models),"MKI67"]

pDT_results_ccle <- read.csv(file.path(WORK_DIR, "ccle_pred", "pDT_results.csv"), row = 1)

Models_with_pred <- cbind(Models,pDT_results_ccle[rownames(Models),])

Models_with_pred$Final_DT <- ifelse(
  is.na(Models_with_pred$Avg_2DT),
  Models_with_pred$pDT_ensemble,
  Models_with_pred$Avg_2DT
)

Models_with_pred <- subset(Models_with_pred, complete.cases(Final_DT))
Models_with_pred_nolung <- subset(Models_with_pred, OncotreeLineage != "Lung")

# PRISM drug matrix: compounds in rows, cell lines in columns, so it is
# transposed. Cell-line IDs have "." for "-" after read.csv and are restored.
drug_df_CCLE <- read.csv(file = file.path(WORK_DIR, "input_files", "ccle", "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv"),
                         sep = ",", header = TRUE, row.names = 1)
drug_df_CCLE <- as.data.frame(t(drug_df_CCLE))
row.names(drug_df_CCLE) <- gsub(pattern = "[.]", "-", row.names(drug_df_CCLE))
Models_drug <- Models_with_pred[row.names(drug_df_CCLE),]
Models_drug_nolung <- Models_with_pred_nolung[row.names(drug_df_CCLE),]

# Spearman correlation of one compound's effect with Final_DT, all cell lines
# and without lung lines. Applied to every column of drug_df_CCLE.
correlation_function_dt_drug <- function(col1) {
  test <- cor.test((col1), Models_drug$Final_DT, method = "spearman")
  list(correlation = test$estimate, p_value = test$p.value)
}
correlation_function_dt_drug_nolung <- function(col1) {
  test <- cor.test((col1), Models_drug_nolung$Final_DT, method = "spearman")
  list(correlation = test$estimate, p_value = test$p.value)
}
drug_dt_corr_ress <- mapply(correlation_function_dt_drug,drug_df_CCLE, SIMPLIFY = FALSE)
drug_dt_corr_ress_nolung <- mapply(correlation_function_dt_drug_nolung,drug_df_CCLE, SIMPLIFY = FALSE)

# Collect rho and p per compound; BH q-values across compounds.
correlation_results_drug_dtime <- data.frame(
  column = names(drug_dt_corr_ress),
  correlation = sapply(drug_dt_corr_ress, function(x) x$correlation),
  p_value = sapply(drug_dt_corr_ress, function(x) x$p_value)
)
rownames(correlation_results_drug_dtime) <- correlation_results_drug_dtime$column
correlation_results_drug_dtime$q_value<-p.adjust(correlation_results_drug_dtime$p_value,method = "BH")

correlation_results_drug_dtime_nolung <- data.frame(
  column = names(drug_dt_corr_ress_nolung),
  correlation = sapply(drug_dt_corr_ress_nolung, function(x) x$correlation),
  p_value = sapply(drug_dt_corr_ress_nolung, function(x) x$p_value)
)
rownames(correlation_results_drug_dtime_nolung) <- correlation_results_drug_dtime_nolung$column
correlation_results_drug_dtime_nolung$q_value<-p.adjust(correlation_results_drug_dtime_nolung$p_value,method = "BH")

# Add compound names and synonyms (compound list), then targets and mechanism
# (PortalCompounds.csv), matched on compound ID.
drugs_description <- read.csv(file.path(WORK_DIR, "input_files", "ccle", "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv"),
                              sep = ",", header = TRUE, row.names = 5)
correlation_results_drug_dtime$Drug_Name<-drugs_description[row.names(correlation_results_drug_dtime),"Drug.Name"]
correlation_results_drug_dtime$Synonyms<-drugs_description[row.names(correlation_results_drug_dtime),"Synonyms"]
correlation_results_drug_dtime_nolung$Drug_Name<-drugs_description[row.names(correlation_results_drug_dtime_nolung),"Drug.Name"]
correlation_results_drug_dtime_nolung$Synonyms<-drugs_description[row.names(correlation_results_drug_dtime_nolung),"Synonyms"]

drugs_description2 <- read.csv(file.path(WORK_DIR, "input_files", "ccle", "PortalCompounds.csv"),
                               sep = ",", header = TRUE, row.names = 10)
correlation_results_drug_dtime$GeneSymbolOfTargets<-drugs_description2[row.names(correlation_results_drug_dtime),"GeneSymbolOfTargets"]
correlation_results_drug_dtime$TargetOrMechanism<-drugs_description2[row.names(correlation_results_drug_dtime),"TargetOrMechanism"]
correlation_results_drug_dtime_nolung$GeneSymbolOfTargets<-drugs_description2[row.names(correlation_results_drug_dtime_nolung),"GeneSymbolOfTargets"]
correlation_results_drug_dtime_nolung$TargetOrMechanism<-drugs_description2[row.names(correlation_results_drug_dtime_nolung),"TargetOrMechanism"]

# ── Figure 5c — DT vs effect of BRD-K18518344, all cell lines ─────────────────
panel_C <- ggplot(data.frame(x = log10(Models_drug$Final_DT), y = drug_df_CCLE$`BRD:BRD-K18518344-001-15-9`), aes(x = x, y = y)) +
  geom_point(alpha = 0.4, colour = "#808080") +
  geom_smooth(method = "lm", colour = "#005AB5", se = FALSE, na.rm = TRUE) +
  xlab(TeX('$log_{10}(DT)$'))+
  ylab(TeX('Drug Effect')) +
  stat_cor(aes(label = paste(..r.label..)), r.accuracy = 0.001,method = "spearman",parse = TRUE,label.y = 2.4,size=6.5)+
  stat_cor(aes(label = paste(after_stat(format_p_value(..p..)))), r.accuracy = 0.001,method = "spearman",parse = TRUE,label.y = 1.4,size=6.5)+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"))

# Volcano data, all cell lines: compounds with a defined rho. The 3 most
# negative and 3 most positive rho get their mechanism as a label. q < 0.1
# compounds are high_corr (rho >= 0) or low_corr (rho < 0).
correlation_results_drug_dtime_plot <- correlation_results_drug_dtime[which(complete.cases(correlation_results_drug_dtime$correlation)),]
index_low <- (base::sort(correlation_results_drug_dtime_plot$correlation, decreasing = FALSE, index.return = TRUE)$ix[c(1:3)])
index_high <- (base::sort(correlation_results_drug_dtime_plot$correlation, decreasing = TRUE, index.return = TRUE)$ix[c(1:3)])

correlation_results_drug_dtime_plot$sig_column <- NA

correlation_results_drug_dtime_plot[c(index_low, index_high), "sig_column"] <-
  stringr::str_wrap(
    correlation_results_drug_dtime_plot[c(index_low, index_high), "TargetOrMechanism"],
    width = 20   # characters per label line
  )

correlation_results_drug_dtime_plot$significance <- "not_significant"
correlation_results_drug_dtime_plot[which((correlation_results_drug_dtime_plot$q_value < 0.1) & ((correlation_results_drug_dtime_plot$correlation) >= 0)),"significance"] <- "high_corr"
correlation_results_drug_dtime_plot[which((correlation_results_drug_dtime_plot$q_value < 0.1) & ((correlation_results_drug_dtime_plot$correlation) < 0)),"significance"] <- "low_corr"


# Same, lung lines removed (used for Tables 7.3 and 7.4)
correlation_results_drug_dtime_nolung_plot <- correlation_results_drug_dtime_nolung[which(complete.cases(correlation_results_drug_dtime_nolung$correlation)),]
index_low <- (base::sort(correlation_results_drug_dtime_nolung_plot$correlation, decreasing = FALSE, index.return = TRUE)$ix[c(1:3)])
index_high <- (base::sort(correlation_results_drug_dtime_nolung_plot$correlation, decreasing = TRUE, index.return = TRUE)$ix[c(1:3)])

correlation_results_drug_dtime_nolung_plot$sig_column <- NA

correlation_results_drug_dtime_nolung_plot[c(index_low, index_high), "sig_column"] <-
  stringr::str_wrap(
    correlation_results_drug_dtime_nolung_plot[c(index_low, index_high), "TargetOrMechanism"],
    width = 20   # characters per label line
  )

correlation_results_drug_dtime_nolung_plot$significance <- "not_significant"
correlation_results_drug_dtime_nolung_plot[which((correlation_results_drug_dtime_nolung_plot$q_value < 0.1) & ((correlation_results_drug_dtime_nolung_plot$correlation) >= 0)),"significance"] <- "high_corr"
correlation_results_drug_dtime_nolung_plot[which((correlation_results_drug_dtime_nolung_plot$q_value < 0.1) & ((correlation_results_drug_dtime_nolung_plot$correlation) < 0)),"significance"] <- "low_corr"


# ── Figure 5d — rho vs -log10(p), all compounds (p = 0 compounds dropped) ─────
panel_D <- ggplot(data = subset(x = correlation_results_drug_dtime_plot, !(is.infinite(-log10(p_value)))),
                  aes(x = (correlation), y = -log10(p_value), colour = significance))+
  geom_point(aes(alpha = significance))+
  xlab(TeX("$\\rho$")) +
  ylab(TeX('$-log_{10}(P_{value})$'))+
  geom_text_repel(aes(label = subset(x = correlation_results_drug_dtime_plot,!(is.infinite(-log10(p_value))))$sig_column,
                      colour = subset(x = correlation_results_drug_dtime_plot,!(is.infinite(-log10(p_value))))$significance),
                  bg.color = "white", bg.r = 0.2,
                  size = 6.5,box.padding = 1,
                  max.overlaps = Inf,point.size = 2, segment.alpha = 0.5,
                  show.legend = FALSE)+
  scale_colour_manual(values = c("not_significant" = "#808080", "high_corr" = "#DC3220","low_corr" = "#005AB5"))+
  scale_alpha_manual(values = c("not_significant" = 0.1, "high_corr" = 1,"low_corr" = 1), guide="none")+
  scale_x_continuous(limits=c(-0.3,0.3), breaks = round(seq(-0.3,0.3,0.3), digits = 1)) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=20, face = "bold"),
        axis.title = element_text(size=20, face = "bold"))

# ── Assemble Figure 5: a | b over c | d ────────────────────────────────────────
fig5 <- wrap_elements(full =panel_A) + wrap_elements(full =panel_B) +
  wrap_elements(full =panel_C) + wrap_elements(full =panel_D) +
  plot_layout(
    design = paste(c("AB",
                     "CD"),
                   collapse = "\n"),
    heights = c(1, 1)) + 
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 20, face  = "bold"))

ggsave(file.path(OUT_DIR, "Figure_5.pdf"), fig5, width = 8, height = 8)
ggsave(file.path(OUT_DIR, "Figure_5.png"), fig5, width = 8, height = 8)

################################################################################
# SUPPLEMENTAL FIGURE 6 — the same volcano for the other therapy classes
#
# For ICB, targeted and chemo + targeted cohorts, label only significant
# cohorts that are also among the 10 lowest / 10 highest Mean_diff.
################################################################################
wilcox_results_icb <- subset(x = wilcox_results_dr, (!(is.infinite(-log10(p_value))) & therapy_class_updated == "ICB"))
index_low <- intersect((base::sort(wilcox_results_icb$Mean_diff, decreasing = FALSE, index.return = TRUE)$ix[c(1:10)]),
                       which(wilcox_results_icb$significance == "high_dt_non-responder"))
index_high <- intersect((base::sort(wilcox_results_icb$Mean_diff, decreasing = TRUE, index.return = TRUE)$ix[c(1:10)]),
                        which(wilcox_results_icb$significance == "high_dt_responder"))
wilcox_results_icb$sig_column <- NA
wilcox_results_icb[c(index_low,index_high),"sig_column"] <-wilcox_results_icb[c(index_low,index_high),"short_name_no_cancer"]

lbl_icb <- wilcox_results_icb[c(index_low, index_high), "short_name_no_cancer"]
lbl_icb <- gsub("([/,;+])", "\\1 ", lbl_icb)          # allow breaks after / , ; +
lbl_icb <- gsub("\\s+", " ", trimws(lbl_icb))

wilcox_results_icb[c(index_low, index_high), "sig_column"] <-
  stringr::str_wrap(
    lbl_icb,
    width = 40   # characters per label line
  )

wilcox_results_targeted <- subset(x = wilcox_results_dr, (!(is.infinite(-log10(p_value))) & therapy_class_updated == "targeted"))
index_low <- intersect((base::sort(wilcox_results_targeted$Mean_diff, decreasing = FALSE, index.return = TRUE)$ix[c(1:10)]),
                       which(wilcox_results_targeted$significance == "high_dt_non-responder"))
index_high <- intersect((base::sort(wilcox_results_targeted$Mean_diff, decreasing = TRUE, index.return = TRUE)$ix[c(1:10)]),
                        which(wilcox_results_targeted$significance == "high_dt_responder"))
wilcox_results_targeted$sig_column <- NA
wilcox_results_targeted[c(index_low,index_high),"sig_column"] <-wilcox_results_targeted[c(index_low,index_high),"short_name_no_cancer"]

lbl_targeted <- wilcox_results_targeted[c(index_low, index_high), "short_name_no_cancer"]
lbl_targeted <- gsub("([/,;+])", "\\1 ", lbl_targeted)          # allow breaks after / , ; +
lbl_targeted <- gsub("\\s+", " ", trimws(lbl_targeted))

wilcox_results_targeted[c(index_low, index_high), "sig_column"] <-
  stringr::str_wrap(
    lbl_targeted,
    width = 40   # characters per label line
  )

wilcox_results_chemo_targeted <- subset(x = wilcox_results_dr, (!(is.infinite(-log10(p_value))) & therapy_class_updated == "chemo+targeted"))
index_low <- intersect((base::sort(wilcox_results_chemo_targeted$Mean_diff, decreasing = FALSE, index.return = TRUE)$ix[c(1:10)]),
                       which(wilcox_results_chemo_targeted$significance == "high_dt_non-responder"))
index_high <- intersect((base::sort(wilcox_results_chemo_targeted$Mean_diff, decreasing = TRUE, index.return = TRUE)$ix[c(1:10)]),
                        which(wilcox_results_chemo_targeted$significance == "high_dt_responder"))
wilcox_results_chemo_targeted$sig_column <- NA
wilcox_results_chemo_targeted[c(index_low,index_high),"sig_column"] <-wilcox_results_chemo_targeted[c(index_low,index_high),"short_name_no_cancer"]

lbl_chemo_targeted <- wilcox_results_chemo_targeted[c(index_low, index_high), "short_name_no_cancer"]
lbl_chemo_targeted <- gsub("([/,;+])", "\\1 ", lbl_chemo_targeted)          # allow breaks after / , ; +
lbl_chemo_targeted <- gsub("\\s+", " ", trimws(lbl_chemo_targeted))

wilcox_results_chemo_targeted[c(index_low, index_high), "sig_column"] <-
  stringr::str_wrap(
    lbl_chemo_targeted,
    width = 40   # characters per label line
  )
# ── Supp. 6a — immune checkpoint blockade cohorts ─────────────────────────────
sup_panel_icb_pdt <- ggplot(data = wilcox_results_icb,
                            aes(x = (Mean_diff), y = -log10(p_value), colour = significance))+
  geom_point(aes(alpha = significance))+
  xlab(TeX('$log_{2}(R/NR)$'))+
  ylab(TeX('$-log_{10}(P_{value})$'))+
  geom_text_repel(aes(label = wilcox_results_icb$sig_column,
                      colour = wilcox_results_icb$significance),
                  bg.color = "white", bg.r = 0.2,
                  size = 6.5,box.padding = 1,
                  max.overlaps = Inf,point.size = 2, segment.alpha = 0.5)+
  scale_colour_manual(values = c("not_significant" = "#808080", "high_dt_responder" = "#DC3220","high_dt_non-responder" = "#005AB5"))+
  scale_alpha_manual(values = c("not_significant" = 0.5, "high_dt_responder" = 1,"high_dt_non-responder" = 1), guide="none")+
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=20, face = "bold"),
        axis.title = element_text(size=20, face = "bold"))

# ── Supp. 6b — targeted-therapy cohorts ───────────────────────────────────────
sup_panel_targeted_pdt <- ggplot(data = wilcox_results_targeted,
                                 aes(x = (Mean_diff), y = -log10(p_value), colour = significance))+
  geom_point(aes(alpha = significance))+
  xlab(TeX('$log_{2}(R/NR)$'))+
  ylab(TeX('$-log_{10}(P_{value})$'))+
  geom_text_repel(aes(label = wilcox_results_targeted$sig_column,
                      colour = wilcox_results_targeted$significance),
                  bg.color = "white", bg.r = 0.2,
                  size = 6.5,box.padding = 1,
                  max.overlaps = Inf,point.size = 2, segment.alpha = 0.5)+
  scale_colour_manual(values = c("not_significant" = "#808080", "high_dt_responder" = "#DC3220","high_dt_non-responder" = "#005AB5"))+
  scale_alpha_manual(values = c("not_significant" = 0.5, "high_dt_responder" = 1,"high_dt_non-responder" = 1), guide="none")+
  scale_x_continuous(limits = c(-0.1,0.1),breaks = c(-0.1,0,0.1)) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=20, face = "bold"),
        axis.title = element_text(size=20, face = "bold"))

# ── Supp. 6c — chemo + targeted cohorts ───────────────────────────────────────
sup_panel_chemo_targeted_pdt <- ggplot(data = wilcox_results_chemo_targeted,
                                       aes(x = (Mean_diff), y = -log10(p_value), colour = significance))+
  geom_point(aes(alpha = significance))+
  xlab(TeX('$log_{2}(R/NR)$'))+
  ylab(TeX('$-log_{10}(P_{value})$'))+
  geom_text_repel(aes(label = wilcox_results_chemo_targeted$sig_column,
                      colour = wilcox_results_chemo_targeted$significance),
                  bg.color = "white", bg.r = 0.2,
                  size = 6.5,box.padding = 1,
                  max.overlaps = Inf,point.size = 2, segment.alpha = 0.5)+
  scale_colour_manual(values = c("not_significant" = "#808080", "high_dt_responder" = "#DC3220","high_dt_non-responder" = "#005AB5"))+
  scale_alpha_manual(values = c("not_significant" = 0.5, "high_dt_responder" = 1,"high_dt_non-responder" = 1), guide="none")+
  scale_x_continuous(limits = c(-0.1,0.05),breaks = c(-0.1,-0.05,0,0.05)) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=20, face = "bold"),
        axis.title = element_text(size=20, face = "bold"))


################################################################################
# SUPPLEMENTAL TABLE 6 (sup_tab6.xlsx) — patient cohort
################################################################################

## ---- Contents sheet, written by hand ---------------------------------------
## One row per sheet. Edit the text freely; keep the two vectors the same
## length, or data.frame() will fail or recycle values.
sup_tab_6_contents <- data.frame(
  sheet = c(
    "Table 6.1",
    "Table 6.2"),
  
  description = c(
    "Clinical data for the patient drug response cohort.",
    "Wilcoxon rank-sum test results for binary drug response vs measured doubling time."),
  
  stringsAsFactors = FALSE)

sup_tab_6_contents <- sup_tab_6_contents[, c("sheet", "description")]

# Table 6.1 — clinical data of the samples used. Internal helper columns are
# dropped by position, and the reclassified therapy class / canonical regimen
# columns are renamed to the plain names. The positions depend on the column
# order of clinical.tsv; check them if that file changes.
sup_tab_6_1 <- clinical_data[,-c(8,9,10,12,13,14,15,16,17,18,21,24,25)]
rownames(sup_tab_6_1) <- NULL
colnames(sup_tab_6_1)[c(9,10,11,12)] <- c("therapy_class","drug_or_regimen","pDT","MKI67 expression")

# Table 6.2 — Wilcoxon results per cohort (helper columns dropped by position)
sup_tab_6_2 <- wilcox_results_dr[,-c(1,7,12)]
rownames(sup_tab_6_2) <- NULL
colnames(sup_tab_6_2)[c(2,7,8,9)] <- c("log2_fc","drug_or_regimen","therapy_class","plot_abbreviation")

# Write the workbook: Contents first, then one sheet per table.
sup_tab_6 <- createWorkbook()

addWorksheet(sup_tab_6, "Contents")
addWorksheet(sup_tab_6, "Table 6.1")
addWorksheet(sup_tab_6, "Table 6.2")

writeData(sup_tab_6, sheet = "Contents", x = sup_tab_6_contents)
writeData(sup_tab_6, sheet = "Table 6.1", x = sup_tab_6_1)
writeData(sup_tab_6, sheet = "Table 6.2", x = sup_tab_6_2)

saveWorkbook(sup_tab_6, file.path(OUT_DIR, "sup_tab6.xlsx"), overwrite = TRUE)

################################################################################
# SUPPLEMENTAL TABLE 7 (sup_tab7.xlsx) — DepMap drug correlations
################################################################################

# Per mechanism of action: number of compounds, how many correlate significantly
# with DT in each direction, and those counts as fractions. The first row of the
# table (the first mechanism in sort order) is dropped.
correlation_results_drug_dtime_summary <- (as.data.frame.matrix(table(correlation_results_drug_dtime_plot[,c("TargetOrMechanism","significance")])))[-c(1),]
correlation_results_drug_dtime_summary$drug_num <- rowSums(correlation_results_drug_dtime_summary)
correlation_results_drug_dtime_summary$sig <- correlation_results_drug_dtime_summary$high_corr + correlation_results_drug_dtime_summary$low_corr
correlation_results_drug_dtime_summary <- cbind(correlation_results_drug_dtime_summary, (correlation_results_drug_dtime_summary/correlation_results_drug_dtime_summary$drug_num)[,c(1,2,5)])
colnames(correlation_results_drug_dtime_summary)[c(6,7,8)] <- paste0(colnames(correlation_results_drug_dtime_summary)[c(6,7,8)], "_frac")

# Same, lung lines removed
correlation_results_drug_dtime_nolung_summary <- (as.data.frame.matrix(table(correlation_results_drug_dtime_nolung_plot[,c("TargetOrMechanism","significance")])))[-c(1),]
correlation_results_drug_dtime_nolung_summary$drug_num <- rowSums(correlation_results_drug_dtime_nolung_summary)
correlation_results_drug_dtime_nolung_summary$sig <- correlation_results_drug_dtime_nolung_summary$high_corr + correlation_results_drug_dtime_nolung_summary$low_corr
correlation_results_drug_dtime_nolung_summary <- cbind(correlation_results_drug_dtime_nolung_summary, (correlation_results_drug_dtime_nolung_summary/correlation_results_drug_dtime_nolung_summary$drug_num)[,c(1,2,5)])
colnames(correlation_results_drug_dtime_nolung_summary)[c(6,7,8)] <- paste0(colnames(correlation_results_drug_dtime_nolung_summary)[c(6,7,8)], "_frac")

## ---- Contents sheet, written by hand ---------------------------------------
## One row per sheet. Edit the text freely; keep the two vectors the same
## length, or data.frame() will fail or recycle values.
sup_tab_7_contents <- data.frame(
  sheet = c(
    "Table 7.1",
    "Table 7.2",
    "Table 7.3",
    "Table 7.4"),
  
  description = c(
    "Spearman Correlation results of DepMap. drug response vs combined (measured and predicted) doubling time.",
    "Summary of drugs enriched for association with doubling time in DepMap.",
    "Spearman Correlation results of DepMap. drug response vs combined (measured and predicted) doubling time, Lung samples removed.",
    "Summary of drugs enriched for association with doubling time in DepMap, Lung samples removed."
    ),
  
  stringsAsFactors = FALSE)


sup_tab_7_contents <- sup_tab_7_contents[, c("sheet", "description")]


# Tables 7.1 / 7.3 — per-compound results (label columns dropped);
# Tables 7.2 / 7.4 — per-mechanism summaries
sup_tab_7_1 <- correlation_results_drug_dtime_plot[,-c(9,10)]
colnames(sup_tab_7_1)[1] <- "Drug_ID"
rownames(sup_tab_7_1) <- NULL

sup_tab_7_2 <- tibble::rownames_to_column(correlation_results_drug_dtime_summary[,-c(3,5,8)])
colnames(sup_tab_7_2) <- c("TargetOrMechanism","Significant Positive Correlation","Significant Negative Correlation","Total Drugs","Fraction of Significant Positive ","Fraction of Significant Negative")

sup_tab_7_3 <- correlation_results_drug_dtime_nolung_plot[,-c(9,10)]
colnames(sup_tab_7_3)[1] <- "Drug_ID"
rownames(sup_tab_7_3) <- NULL

sup_tab_7_4 <- tibble::rownames_to_column(correlation_results_drug_dtime_nolung_summary[,-c(3,5,8)])
colnames(sup_tab_7_4) <- c("TargetOrMechanism","Significant Positive Correlation","Significant Negative Correlation","Total Drugs","Fraction of Significant Positive ","Fraction of Significant Negative")


# Write the workbook: Contents first, then one sheet per table.
sup_tab_7 <- createWorkbook()

addWorksheet(sup_tab_7, "Contents")
addWorksheet(sup_tab_7, "Table 7.1")
addWorksheet(sup_tab_7, "Table 7.2")
addWorksheet(sup_tab_7, "Table 7.3")
addWorksheet(sup_tab_7, "Table 7.4")

writeData(sup_tab_7, sheet = "Contents", x = sup_tab_7_contents)
writeData(sup_tab_7, sheet = "Table 7.1", x = sup_tab_7_1)
writeData(sup_tab_7, sheet = "Table 7.2", x = sup_tab_7_2)
writeData(sup_tab_7, sheet = "Table 7.3", x = sup_tab_7_3)
writeData(sup_tab_7, sheet = "Table 7.4", x = sup_tab_7_4)

saveWorkbook(sup_tab_7, file.path(OUT_DIR, "sup_tab7.xlsx"), overwrite = TRUE)

# ── Assemble Supplemental Figure 5: a | b | c ──────────────────────────────────
sup_fig6 <- wrap_elements(full = sup_panel_icb_pdt) +
  wrap_elements(full = sup_panel_targeted_pdt) +
  wrap_elements(full = sup_panel_chemo_targeted_pdt) +
  plot_layout(
    design = paste(c("ABC"),
                   collapse = "\n"),
    heights = c(1)) + 
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 20, face  = "bold"))

ggsave(file.path(OUT_DIR, "sup_Figure_6.pdf"), sup_fig6, width = 12, height = 4)
ggsave(file.path(OUT_DIR, "sup_Figure_6.png"), sup_fig6, width = 12, height = 4)

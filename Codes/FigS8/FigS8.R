# =============================================================================================
# Loss of Y, chromosome-arm loss and gain, and gene essentiality in the Sanger tumour-organoid biobank
# Supplementary figure (panels B-F) and supplementary tables
#
# Data: Herranz-Ors et al., Nature 2026, doi:10.1038/s41586-026-10830-y
#   Supplementary Tables 2, 3, 4 and 6, and CNV-segments_organoids.tsv (Figshare,
#   doi:10.6084/m9.figshare.28339340). Reference files (hg38 gaps, cytobands, ENCODE blacklist,
#   Ensembl GRCh38 release-98 genes) are downloaded into reference_dir on the first run.
#
# Sections
#   0.  Paths and settings
#   1.  Packages and helper functions
#   2.  CRISPR dependency scores (one screen per organoid)
#   3.  Organoid annotations (clinical, WGS summary, screen quality)
#   4.  Copy-number segments and DNA quality control
#   5.  Loss-of-Y calls (DNA + RNA)
#   6.  Chromosome-arm calls (loss/neutral/gain) for every organoid
#   7.  Comparison cohorts: altered vs neutral organoids for every arm (and Y)
#   8.  Differential dependency: one regression per gene for every tested arm
#   9.  Expression of genes on the altered arm (panel C)
#   10. Figures (panels B-G)
#   11. Supplementary tables
#   12. Provenance (input checksums, settings, session info)
#
# Panels
#   B  LOY and Y-WT organoids in the comparison cohort
#   C  Expression of genes on the altered arm, altered vs neutral organoids
#   D  Volcano plot of every gene x arm-loss and gene x arm-gain test
#   E  DDX3X dependency in LOY vs Y-WT organoids
#   F  DDX3X dependency vs DDX3Y expression
#   G  Biomarker-target relationships
#
# Conventions that I used throughout:
#   -Dependency score = the normalized log fold change in Supplementary Table 6 (lower = more essential).
#   -Relative dependency = neutral - altered, i.e. the regression coefficient times -1. POSITIVE means
#     more essential in organoids with the arm altered. 
#   -A hit is p < 0.01 and relative dependency >= 0.2.
# =============================================================================================


## ==== 0. Paths and settings ===================================================================
project_dir <- "/Users/shaunak/Documents/Crispr-screens_Organoid"
supp_dir <- file.path(project_dir, "41586_2026_10830_MOESM3_ESM")

# ---- Input files ----
crispr_file <- file.path(supp_dir, "supplementary_table_6_revision.csv") # space-separated
clinical_xlsx <- file.path(supp_dir, "supplementary_table_2_revision.xlsx")
molecular_xlsx <- file.path(supp_dir, "supplementary_table_3_revision.xlsx") # WGS summary + RNA TPM
screen_qc_xlsx <- file.path(supp_dir, "supplementary_table_4_revision.xlsx") # CRISPR screen QC
cnv_segments_file <- file.path(project_dir, "CNV-segments_organoids.tsv") # comma-separated
reference_dir <- file.path(project_dir, "reference") # downloaded references

clinical_sheet <- "clinical_data"
wgs_sheet <- "data_obtained_from_WGS"
rna_sheet <- "RNAseq_TPM"
screen_qc_sheet <- "summary_QC_metrics"

gap_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz"
cytoband_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz"
blacklist_url <- "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
gtf_url <- "https://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.chr.gtf.gz"

# ---- Output folders (everything is written here) ----
output_dir <- file.path(project_dir, "results_supplementary_figure")
figure_dir <- file.path(output_dir, "figures") # panels + their source data
table_dir <- file.path(output_dir, "tables") # intermediate and full results
supp_table_dir <- file.path(output_dir, "supplementary_tables")  # the supplementary table workbook

# ---- Loss-of-Y calls (section 5) ----
# DNA: absolute copy number over the masked male-specific region of chrY
y_dna_min_coverage   <- 0.8    # supported fraction of the region needed for a call
y_dna_loss_median    <- 0.25   # loss: weighted median copy number below this ...
y_dna_loss_fraction  <- 0.9    # ... and at least this fraction of supported bases below 0.5
y_dna_keep_median    <- 0.75   # retained: weighted median at or above this ...
y_dna_keep_fraction  <- 0.8    # ... and at least this fraction of supported bases at or above 0.5
# RNA: seven broadly expressed Y genes; threshold per gene = max(1 TPM, female 99th percentile)
y_rna_floor_tpm      <- 1      # low Y: median marker TPM at or below this ...
y_rna_min_low        <- 6      # ... and at least this many of the 7 markers below their threshold
min_purity           <- 0.8    # tumour purity needed for DNA calls
# Stability: a high-confidence LOY call must hold under every one of these settings
dna_stability_grid <- data.table::CJ(coverage = c(0.7, 0.8, 0.9), median_cutoff = c(0.1, 0.25, 0.5),
                                     fraction = c(0.8, 0.9, 0.95))
rna_stability_grid <- data.table::CJ(floor = c(0.5, 1, 2), number = c(5, 6, 7))

# ---- Chromosome-arm calls (section 6), on copy number relative to ploidy ----
arm_min_coverage <- 0.90           # >= 90% of the masked arm covered by supported segments
loss_ratio       <- 0.75           # loss: median < 0.75 and >= 80% of the arm below 0.75
loss_fraction    <- 0.80
neutral_range    <- c(0.85, 1.15)  # neutral: median in [0.85, 1.15] and >= 80% of the arm in that range
neutral_fraction <- 0.80
gain_ratio       <- 1.25           # gain: median > 1.25 and >= 80% of the arm above 1.25
gain_fraction    <- 0.80
partial_fraction <- 0.20           # partial loss: 20-80% of the arm below 0.75

# ---- Which comparisons are tested (sections 7-8) ----
min_per_lineage <- 3    # a cancer type is included only with >= 3 altered AND >= 3 neutral organoids
min_altered     <- 25   # altered organoids needed for an arm to be tested
min_neutral     <- 20   # neutral organoids needed
min_donors      <- 40   # distinct donors needed
min_residual_df <- 20   # degrees of freedom left after fitting the model

# ---- Hit rule (section 8) ----
p_cut      <- 0.01
effect_cut <- 0.2       # relative dependency (neutral - altered) needed

# ---- Panel C ----
min_tpm_neutral  <- 1      # a gene is used for an arm if its median TPM in the neutral group is >= this
include_Y_panelC <- FALSE  # the LOY call itself requires low Y expression, so a Y box would be circular
show_p_panelC    <- TRUE   # Wilcoxon p value (genes on the arm vs genes elsewhere) over each pair
p_digits_panelC  <- 2      # significant digits for those p values (Yi uses 4, which is crowded here)

# ---- Panels E and F ----
focus_x_gene      <- "DDX3X"      # X gene (dependency)
focus_y_gene      <- "DDX3Y"      # its Y paralog (expression); must be one of the seven Y markers
panelE_p_method   <- "wilcoxon"   # p over the boxes: "wilcoxon" (rank-sum test) or "model" (regression)

# P values for plot labels in Nature style: italic P, "P = 2.9 × 10^−8" (a plotmath string, so use parse = TRUE)
fmt_P <- function(p) {
  if (is.na(p)) return('italic(P)~"="~"NA"')
  if (p < 2.2e-16) return('italic(P)~"<"~"2.2"~"\u00d7"~10^"\u221216"')
  if (p >= 1e-3) return(sprintf('italic(P)~"="~"%s"', formatC(p, format = "g", digits = 2)))
  e <- floor(log10(p)); m <- round(p / 10^e, 1)
  if (m >= 10) { m <- m / 10; e <- e + 1 }
  sprintf('italic(P)~"="~"%.1f"~"\u00d7"~10^"\u2212%d"', m, -e)
}

# ---- Figure style (manuscript style: Arial, 7 pt, axis lines that stop at the outer ticks) ----
fig_font      <- "Arial"
fig_base_pt   <- 7
x_label_angle <- 90               # angle of the group labels under the boxes and bars (panels B, E)
label_ywt     <- "Y-WT"           # Y-retained male organoids
label_loy     <- "LOY"
col_ywt <- "#8E73C2"; col_ywt_dark <- "#5E4396"   # Y-WT purple (panels B, E)
col_loy <- "#F5AC3D"; col_loy_dark <- "#D9881A"   # LOY orange (panels B, E)
col_yi_orange <- "#EF8536"        # Yi's orange: hits (D) and lost arms (C)
col_yi_blue   <- "#3A76AF"        # Yi's blue: gained arms (C)
col_ns        <- "grey80"         # genes that are not hits (D)
col_elsewhere <- "grey55"         # genes on other arms (C)

# Volcano (panel D)
volcano_hit_size  <- 1.0          # point size of hits
volcano_focus_size <- 2           # point size for the labelled genes
volcano_ns_size   <- 0.25         # point size of everything else
volcano_label_n   <- 0            # most significant hits labelled
volcano_always    <- c("Y: DDX3X", "Y: RPS4X")   # "event: gene", labelled even if not a hit

# Panel sizes in mm (width, height)
size_B <- c(62, 58); size_C <- c(114, 58); size_D <- c(56, 50); size_E <- c(32, 50); size_F <- c(44, 50)

# ---- Fixed inputs (copied from the earlier scripts; you shouldn't need to change these) ----
# 16 organoids were screened with both libraries; the authors keep one screen each. Exclusion list from
# the authors' code (github.com/Garnett-Lab/SangerOrganoidBiobank, commit 0acbbd4, CRISPR_biomarker_analysis.R)
duplicate_screens_to_drop <- data.table::data.table(
  sample_ID = c("HCM-SANG-0267-D12", "HCM-SANG-0268-C18", "HCM-SANG-0270-C20", "HCM-SANG-0272-C20",
                "HCM-SANG-0273-C18", "HCM-SANG-0276-C18", "HCM-SANG-0282-C18", "HCM-SANG-0529-C18",
                "HCM-SANG-0284-C18", "HCM-SANG-0288-C18", "WTSI-COLO_278",     "WTSI-COLO_376",
                "HCM-SANG-0291-C15", "HCM-SANG-0299-C15-A", "HCM-SANG-0299-C15-B", "HCM-SANG-0306-C15"),
  library   = c("minLib", "v1_1", "minLib", "v1_1", "minLib", "minLib", "minLib", "v1_1",
                "minLib", "v1_1", "minLib", "minLib", "v1_1", "minLib", "minLib", "minLib"))
y_markers <- c(RPS4Y1 = "ENSG00000129824", DDX3Y = "ENSG00000067048", KDM5D = "ENSG00000012817",
               USP9Y  = "ENSG00000114374", EIF1AY = "ENSG00000198692", UTY  = "ENSG00000183878",
               ZFY    = "ENSG00000067646")
mixed_known  <- c("HCM-SANG-0270-C20", "HCM-SANG-0276-C18", "HCM-SANG-0308-C15", "WTSI-OVAR_018")
acrocentric  <- c("13p", "14p", "15p", "21p", "22p")   # < 1 Mb of usable sequence: never called
chrY_length  <- 57227415                               # hg38
par_regions  <- data.table::data.table(start0 = c(10000, 56887902), end0 = c(2781479, 57217415))  # PAR1, PAR2


## ==== 1. Packages and helper functions ========================================================
required <- c("data.table", "readxl", "R.utils", "ggplot2", "ggrepel", "svglite", "scales", "writexl")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Install these packages first: install.packages(c(\"", paste(missing_pkgs, collapse = "\", \""), "\"))",
       call. = FALSE)
}
if (packageVersion("ggplot2") < "3.5.0") {
  stop("ggplot2 3.5.0 or later is needed (capped axes, inset legend): install.packages(\"ggplot2\")", call. = FALSE)
}
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
})
inputs <- c(crispr_file, clinical_xlsx, molecular_xlsx, screen_qc_xlsx, cnv_segments_file)
if (any(!file.exists(inputs))) {
  stop("Missing input files:\n  ", paste(inputs[!file.exists(inputs)], collapse = "\n  "), call. = FALSE)
}
for (folder in c(reference_dir, figure_dir, table_dir, supp_table_dir)) dir.create(folder, showWarnings = FALSE, recursive = TRUE)

# ---- Reading files ----

# Return the first column of `dt` whose name matches one of `candidates` (case-insensitive)
find_col <- function(dt, candidates, what) {
  hit <- names(dt)[tolower(names(dt)) %in% tolower(candidates)]
  if (length(hit) == 0) {
    stop("Could not find the ", what, " column. Columns are: ", paste(names(dt), collapse = ", "), call. = FALSE)
  }
  hit[1]
}

# Supplementary sheets can have notes above the table: use the row containing `header_word` as the
# header. Everything is read as text and converted to numbers later, where needed.
read_sheet_table <- function(path, sheet, header_word) {
  raw <- as.data.table(read_excel(path, sheet = sheet, col_names = FALSE, col_types = "text",
                                  .name_repair = "minimal"))
  first_rows <- as.matrix(raw[seq_len(min(50, nrow(raw)))])
  is_header  <- apply(first_rows, 1, function(r) any(tolower(trimws(r)) == tolower(header_word), na.rm = TRUE))
  header_row <- which(is_header)[1]
  if (is.na(header_row)) stop("No '", header_word, "' header found in sheet ", sheet, call. = FALSE)
  header <- trimws(as.character(unlist(raw[header_row])))
  keep   <- !is.na(header) & header != ""
  out    <- raw[(header_row + 1):nrow(raw), which(keep), with = FALSE]
  setnames(out, make.unique(header[keep]))
  out[rowSums(!is.na(out)) > 0]   # drop empty rows
}

# Download a reference file once
get_reference <- function(url, dest) {
  if (!file.exists(dest)) {
    message("Downloading ", url)
    download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  dest
}

# ---- Genomic intervals (all half-open, 0-based: [start0, end0)) ----

# Merge overlapping intervals
union_intervals <- function(start0, end0) {
  if (length(start0) == 0) return(data.table(start0 = numeric(0), end0 = numeric(0)))
  o <- order(start0, end0)
  s <- start0[o]
  e <- end0[o]
  out_s <- numeric(0); out_e <- numeric(0)
  cur_s <- s[1]; cur_e <- e[1]
  for (i in seq_along(s)[-1]) {
    if (s[i] <= cur_e) {
      cur_e <- max(cur_e, e[i])
    } else {
      out_s <- c(out_s, cur_s); out_e <- c(out_e, cur_e)
      cur_s <- s[i]; cur_e <- e[i]
    }
  }
  data.table(start0 = c(out_s, cur_s), end0 = c(out_e, cur_e))
}

# The part of [a, b) that is NOT covered by the masked intervals
complement_within <- function(a, b, mask) {
  m <- mask[end0 > a & start0 < b]
  if (nrow(m) == 0) return(data.table(start0 = a, end0 = b))
  m <- union_intervals(pmax(m$start0, a), pmin(m$end0, b))
  gaps <- data.table(start0 = c(a, m$end0), end0 = c(m$start0, b))
  gaps[end0 > start0]
}

# Number of bases of [s, e) that fall inside the territory intervals
overlap_bp <- function(s, e, territory) {
  sum(pmax(0, pmin(e, territory$end0) - pmax(s, territory$start0)))
}

# Weighted median: the first value at which the cumulative weight reaches half the total weight
weighted_median <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  x <- x[keep]; w <- w[keep]
  o <- order(x)
  x[o][which(cumsum(w[o]) >= sum(w) / 2)[1]]
}

# ---- Regression ----

# Design matrix for one comparison:
#   dependency ~ altered + cancer type + sex/Y group + library + MSI + ploidy + screen quality (AUC_PR)
# The sex/Y group is left out for the Y comparison itself (all organoids are male and Y status is the
# comparison). Continuous covariates are z-scored within the cohort; adjustment terms that are constant
# in the cohort are dropped. Returns X = NULL with a reason if the design is unstable.
make_design <- function(frame, use_sex_y) {
  X <- cbind(intercept = 1, altered = as.numeric(frame$altered))

  add_dummies <- function(X, values, prefix) {
    levels_present <- sort(unique(values))
    for (lev in levels_present[-1]) {             # first level is the reference
      X <- cbind(X, as.numeric(values == lev))
      colnames(X)[ncol(X)] <- paste0(prefix, lev)
    }
    X
  }
  X <- add_dummies(X, frame$primary_tumour_type, "type_")
  if (use_sex_y) X <- add_dummies(X, frame$sex_Y_group, "sexY_")
  X <- add_dummies(X, frame$library, "library_")
  X <- cbind(X, MSI = as.numeric(frame$msStatus == "MSI"))

  for (v in c("ploidy", "AUC_PR")) {
    x <- frame[[v]]
    z <- if (sd(x) > 0) (x - mean(x)) / sd(x) else rep(0, length(x))
    X <- cbind(X, z)
    colnames(X)[ncol(X)] <- paste0(v, "_z")
  }

  # Drop adjustment columns that don't vary in this cohort (e.g. cancer type in a one-lineage cohort)
  nuisance <- setdiff(colnames(X), c("intercept", "altered"))
  is_constant <- vapply(nuisance, function(col) length(unique(X[, col])) == 1, logical(1))
  X <- X[, setdiff(colnames(X), nuisance[is_constant]), drop = FALSE]

  # Stability checks
  if (qr(X)$rank < ncol(X)) return(list(X = NULL, reason = "design_not_full_rank"))
  singular_values <- svd(X)$d
  condition_number <- max(singular_values) / min(singular_values)
  others <- X[, colnames(X) != "altered", drop = FALSE]
  altered_resid <- X[, "altered"] - others %*% qr.solve(others, X[, "altered"])
  vif <- sum((X[, "altered"] - mean(X[, "altered"]))^2) / sum(altered_resid^2)
  residual_df <- nrow(X) - ncol(X)
  reasons <- c(if (condition_number > 100) "condition_number_above_100",
               if (vif > 10) "altered_VIF_above_10",
               if (residual_df < min_residual_df) paste0("residual_df_below_", min_residual_df))
  if (length(reasons) > 0) return(list(X = NULL, reason = paste(reasons, collapse = ";")))
  list(X = X, reason = "")
}

# OLS for all genes at once (identical to running lm() gene by gene), with a donor-clustered (CR1)
# standard error for the "altered" coefficient. t tests and 95% CIs use (number of donors - 1) df.
fit_genes <- function(X, Y, donors) {
  n <- nrow(X); k <- ncol(X); g <- length(unique(donors))
  XtX_inv   <- solve(crossprod(X))
  B         <- XtX_inv %*% crossprod(X, Y)          # coefficients, one column per gene
  residuals <- Y - X %*% B
  a         <- (XtX_inv %*% t(X))["altered", ]       # row of the hat-matrix for "altered"
  by_donor  <- rowsum(residuals * a, group = donors)
  se   <- sqrt((g / (g - 1)) * ((n - 1) / (n - k)) * colSums(by_donor^2))
  beta <- B["altered", ]
  ok   <- apply(Y, 2, function(y) diff(range(y)) > 1e-12) & colSums(residuals^2) > 1e-20 & se > 0
  t_value  <- ifelse(ok, beta / se, NA_real_)
  t_crit   <- qt(0.975, df = g - 1)
  data.table(ensembl_id = colnames(Y), beta = beta, SE = ifelse(ok, se, NA_real_),
             CI95_low = ifelse(ok, beta - t_crit * se, NA_real_),
             CI95_high = ifelse(ok, beta + t_crit * se, NA_real_),
             p = 2 * pt(-abs(t_value), df = g - 1), n_models = n, n_donors = g)
}

# ---- Figures ----

pt_to_mm <- 1 / ggplot2::.pt
txt <- function(pt) pt * pt_to_mm            # text size in points -> ggplot size units

pub_theme <- theme_classic(base_size = fig_base_pt, base_family = fig_font) +
  theme(axis.line  = element_line(linewidth = 0.2, colour = "grey40"),
        axis.ticks = element_line(linewidth = 0.2, colour = "grey40"),
        axis.ticks.length = unit(1, "mm"),
        axis.text  = element_text(colour = "black", size = fig_base_pt),
        axis.title = element_text(colour = "black", size = fig_base_pt),
        strip.background = element_blank(), strip.text = element_text(size = fig_base_pt + 1),
        legend.text = element_text(size = fig_base_pt), legend.key.size = unit(2.5, "mm"),
        legend.background = element_blank(),
        plot.margin = margin(2, 3, 2, 2, "mm"))
capped_axes <- guides(y = guide_axis(cap = "both"))  # like Yi's range frame

# p values in the style of Yi's figures: "p = 4.91e-08", "p = 0.0123", "p < 2.2e-16"
format_p <- function(p, digits = 3) {
  out <- character(length(p))
  for (i in seq_along(p)) {
    if (is.na(p[i])) {
      out[i] <- "p = NA"
    } else if (p[i] < 2.2e-16) {
      out[i] <- "p < 2.2e-16"
    } else if (p[i] < 1e-3) {
      out[i] <- sprintf(paste0("p = %.", digits - 1, "e"), p[i])
    } else {
      out[i] <- paste0("p = ", formatC(p[i], format = "g", digits = digits))
    }
  }
  out
}

# Save a figure as SVG (editable text), PDF and PNG
save_figure <- function(plot, stem, size_mm) {
  w <- size_mm[1] / 25.4
  h <- size_mm[2] / 25.4
  out <- file.path(figure_dir, stem)
  ggsave(paste0(out, ".svg"), plot, width = w, height = h,
         device = function(filename, ...) svglite::svglite(filename, ..., fix_text_size = FALSE))
  if (capabilities("cairo")) {
    ggsave(paste0(out, ".pdf"), plot, width = w, height = h, device = cairo_pdf)
  } else {   # base pdf(): Helvetica instead of Arial
    ggsave(paste0(out, ".pdf"), plot + theme(text = element_text(family = "Helvetica")), width = w, height = h,
           device = function(filename, ...) grDevices::pdf(file = filename, ..., useDingbats = FALSE))
  }
  ggsave(paste0(out, ".png"), plot, width = w, height = h, dpi = 600)
  cat("Saved", paste0(stem, c(".svg", ".pdf", ".png"), collapse = ", "), "\n")
}


## ==== 2. CRISPR dependency scores (one screen per organoid) ===================================
lfc <- fread(crispr_file, sep = " ")   # space-separated despite the .csv name
cat(sprintf("Read %d organoids and %d screens   (expected 162 organoids, 178 screens)\n",
            uniqueN(lfc$sample_ID), nrow(unique(lfc[, .(sample_ID, library)]))))

# Keep one screen per organoid
lfc <- lfc[!duplicate_screens_to_drop, on = .(sample_ID, library)]
screens <- unique(lfc[, .(sample_ID, library)])
if (anyDuplicated(screens$sample_ID) > 0) stop("Some organoids still have more than one screen")
if (anyDuplicated(lfc[, .(sample_ID, ensembl_id)]) > 0) stop("Repeated organoid/gene rows in Table 6")

# Organoids x genes matrix of dependency scores
fitness_wide <- dcast(lfc, sample_ID ~ ensembl_id, value.var = "LFC")
fitness <- as.matrix(fitness_wide[, -1])
rownames(fitness) <- fitness_wide$sample_ID
gene_names <- unique(lfc[, .(ensembl_id, gene)])
symbol <- setNames(gene_names$gene, gene_names$ensembl_id)   # Ensembl ID -> gene symbol
rm(lfc, fitness_wide); invisible(gc())

# Genes with a score in every organoid
complete_genes <- colSums(!is.finite(fitness)) == 0
fitness <- fitness[, complete_genes, drop = FALSE]
cat(sprintf("Dependency matrix: %d organoids x %s genes   (expected 162 x 16,562)\n",
            nrow(fitness), format(ncol(fitness), big.mark = ",")))


## ==== 3. Organoid annotations ==================================================================
cat("\n==== 3. Organoid annotations ====\n")

# Clinical annotations (Supplementary Table 2)
clin <- read_sheet_table(clinical_xlsx, clinical_sheet, "sample_ID")
meta <- data.table(
  sample_ID           = clin[[find_col(clin, "sample_ID", "sample ID")]],
  individual_ID       = clin[[find_col(clin, "individual_ID", "donor ID")]],
  primary_tumour_type = clin[[find_col(clin, "primary_tumour_type", "cancer type")]],
  clinical_sex        = trimws(clin[[find_col(clin, "sex", "clinical sex")]]))
if (anyDuplicated(meta$sample_ID) > 0) stop("Duplicate sample_IDs in the clinical table")

# WGS summary (Supplementary Table 3): organoid rows only
wgs_raw <- read_sheet_table(molecular_xlsx, wgs_sheet, "sample_ID")
wgs_raw <- wgs_raw[tolower(get(find_col(wgs_raw, "model", "model"))) == "organoid"]
wgs <- data.table(
  sample_ID  = wgs_raw[[find_col(wgs_raw, "sample_ID", "sample ID")]],
  wgs_gender = wgs_raw[[find_col(wgs_raw, "gender", "WGS gender")]],
  purity     = suppressWarnings(as.numeric(wgs_raw[[find_col(wgs_raw, "purity", "purity")]])),
  ploidy     = suppressWarnings(as.numeric(wgs_raw[[find_col(wgs_raw, "ploidy", "ploidy")]])),
  msStatus   = toupper(trimws(wgs_raw[[find_col(wgs_raw, "msStatus", "MSI status")]])))
if (anyDuplicated(wgs$sample_ID) > 0) stop("Duplicate organoid rows in the WGS summary")
meta <- merge(meta, wgs, by = "sample_ID", all.x = TRUE)

# Clinical sex and WGS-inferred sex must agree; organoids where they don't are never used in comparisons
meta[, sex_annotation_conflict := toupper(trimws(fcoalesce(clinical_sex, ""))) !=
                                  toupper(trimws(fcoalesce(wgs_gender, "")))]
meta[, crispr_screened := sample_ID %in% screens$sample_ID]
cat(sprintf("Organoids: %d | male %d | female %d | sex conflicts %d | screened %d   (expected 256 | 171 | 85 | 6 | 162)\n",
            nrow(meta), sum(meta$clinical_sex == "Male", na.rm = TRUE), sum(meta$clinical_sex == "Female", na.rm = TRUE),
            sum(meta$sex_annotation_conflict), sum(meta$crispr_screened)))

# Screen quality (Supplementary Table 4): AUC_PR, used as a covariate
qc_raw <- read_sheet_table(screen_qc_xlsx, screen_qc_sheet, "sample_ID")
screen_qc <- data.table(
  sample_ID = qc_raw[[find_col(qc_raw, "sample_ID", "sample ID")]],
  library   = qc_raw[[find_col(qc_raw, "library", "library")]],
  AUC_PR    = as.numeric(qc_raw[[find_col(qc_raw, c("AUC_PR", "PR_AUC", "PRRC_AUC", "AUC_PRC"), "AUC_PR")]]))


## ==== 4. Copy-number segments and DNA quality control ==========================================
cat("\n==== 4. Copy-number segments ====\n")
seg <- fread(cnv_segments_file)   # comma-separated despite the .tsv name
standard_names <- c("sample_ID", "model", "chromosome", "start", "end", "copyNumber", "depthWindowCount", "QCStatus")
for (nm in standard_names) setnames(seg, find_col(seg, nm, nm), nm)
seg <- unique(seg[model == "organoid"])   # baseline organoid profiles; identical rows counted once
seg[, chromosome := as.character(chromosome)]
seg[!startsWith(chromosome, "chr"), chromosome := paste0("chr", chromosome)]

# Mixed profiles: four known, plus any organoid with overlapping segments on a chromosome
setorder(seg, sample_ID, chromosome, start, end)
overlap_check <- seg[, .(has_overlap = .N > 1 && any(start[-1] <= cummax(end)[-.N])), by = .(sample_ID, chromosome)]
mixed_ids <- union(mixed_known, overlap_check[has_overlap == TRUE, unique(sample_ID)])
meta[, ambiguous_baseline_DNA_profile := sample_ID %in% mixed_ids]
cat(sprintf("Mixed baseline profiles: %d   (expected 4)\n", length(mixed_ids)))

# Source QC status, summarised per organoid ("PASS" only if every segment is PASS)
source_qc <- seg[, .(dna_source_qc = paste(sort(unique(QCStatus)), collapse = ";")), by = sample_ID]
meta <- merge(meta, source_qc, by = "sample_ID", all.x = TRUE)

# DNA QC for arm calls: not mixed, source QC PASS, purity >= 0.8 (missing purity fails)
meta[, dna_base_qc := !ambiguous_baseline_DNA_profile &
                      fcoalesce(dna_source_qc == "PASS", FALSE) &
                      fcoalesce(purity >= min_purity, FALSE)]

# Regions masked from every copy-number measurement: assembly gaps and blacklisted regions
gap <- fread(get_reference(gap_url, file.path(reference_dir, "hg38_gap.txt.gz")), header = FALSE)
blacklist <- fread(get_reference(blacklist_url, file.path(reference_dir, "hg38-blacklist.v2.bed.gz")), header = FALSE)
mask_all <- rbind(gap[, .(chromosome = V2, start0 = V3, end0 = V4)],
                  blacklist[, .(chromosome = V1, start0 = V2, end0 = V3)])


## ==== 5. Loss-of-Y calls ========================================================================
# DNA state (absolute copy number over the masked male-specific region of chrY):
#   loss_candidate  weighted median CN < 0.25 and >= 90% of supported bases with CN < 0.5
#   retained        weighted median CN >= 0.75 and >= 80% of supported bases with CN >= 0.5
#   unassessed      supported coverage < 0.8 or no usable chrY data; intermediate = anything else
# RNA state (seven Y genes): low_Y if >= 6 of 7 below max(1 TPM, female 99th percentile) and median <= 1 TPM
# High-confidence LOY: DNA loss + low Y RNA, both assays pass QC, and the call holds under all 27 DNA
# and 9 RNA threshold settings. Y-WT (Y-retained): DNA retained + RNA not low.

dna_state <- function(coverage, median_cn, frac_lt05, frac_ge05,
                      min_coverage = y_dna_min_coverage, median_cutoff = y_dna_loss_median,
                      fraction = y_dna_loss_fraction) {
  fcase(is.na(coverage) | coverage < min_coverage, "unassessed",
        median_cn < median_cutoff & frac_lt05 >= fraction, "loss_candidate",
        median_cn >= y_dna_keep_median & frac_ge05 >= y_dna_keep_fraction, "retained",
        default = "intermediate")
}
rna_low <- function(values, thresholds, floor = y_rna_floor_tpm, number = y_rna_min_low) {
  if (!all(is.finite(values))) return(FALSE)
  n_low <- sum(values < pmax(floor, thresholds))
  n_low >= number && median(values) <= floor
}

# ---- 5a. DNA ----
y_mask <- rbind(par_regions, mask_all[chromosome == "chrY", .(start0, end0)])
y_territory <- complement_within(0, chrY_length, y_mask)
y_territory_bp <- sum(y_territory$end0 - y_territory$start0)
cat(sprintf("Masked male-specific chrY: %s bp   (expected 22,344,554)\n", format(y_territory_bp, big.mark = ",")))

segY <- seg[chromosome == "chrY"]
y_row_counts <- segY[, .(dna_Y_segment_rows = .N), by = sample_ID]
meta <- merge(meta, y_row_counts, by = "sample_ID", all.x = TRUE)
meta[is.na(dna_Y_segment_rows), dna_Y_segment_rows := 0L]

# Supported segments (depthWindowCount > 0) of non-mixed profiles; 1-based coordinates -> BED
y_supported <- segY[depthWindowCount > 0 & !(sample_ID %in% mixed_ids)]
y_supported[, bp := mapply(overlap_bp, start - 1, end, MoreArgs = list(territory = y_territory))]
y_supported <- y_supported[bp > 0]
dna_y <- y_supported[, .(dna_Y_supported_bases = sum(bp),
                         dna_Y_median          = weighted_median(copyNumber, bp),
                         dna_Y_fraction_lt05   = sum(bp[copyNumber < 0.5]) / sum(bp),
                         dna_Y_fraction_ge05   = sum(bp[copyNumber >= 0.5]) / sum(bp)), by = sample_ID]
meta <- merge(meta, dna_y, by = "sample_ID", all.x = TRUE)
meta[is.na(dna_Y_supported_bases), dna_Y_supported_bases := 0]
meta[, dna_Y_coverage := dna_Y_supported_bases / y_territory_bp]
meta[, dna_measurement_state := dna_state(dna_Y_coverage, dna_Y_median, dna_Y_fraction_lt05, dna_Y_fraction_ge05)]

# Stability: a loss candidate under every one of the 27 DNA settings
meta[, dna_loss_stable := TRUE]
for (i in seq_len(nrow(dna_stability_grid))) {
  setting <- dna_stability_grid[i]
  state_i <- dna_state(meta$dna_Y_coverage, meta$dna_Y_median, meta$dna_Y_fraction_lt05, meta$dna_Y_fraction_ge05,
                       min_coverage = setting$coverage, median_cutoff = setting$median_cutoff,
                       fraction = setting$fraction)
  meta[, dna_loss_stable := dna_loss_stable & (state_i == "loss_candidate")]
}

# DNA QC for the Y call (missing purity is allowed here, as in the original LOY analysis)
meta[, dna_qc_mixed_profile := ambiguous_baseline_DNA_profile]
meta[, dna_qc_source_fail   := is.na(dna_source_qc) | dna_source_qc != "PASS"]
meta[, dna_qc_low_purity    := !is.na(purity) & purity < min_purity]
meta[, dna_qc_low_Y_cover   := dna_Y_segment_rows == 0 | dna_Y_coverage < y_dna_min_coverage]
meta[, dna_QC_eligible := !dna_qc_mixed_profile & !dna_qc_source_fail & !dna_qc_low_purity & !dna_qc_low_Y_cover]

# ---- 5b. RNA (the full TPM table is kept for panels C and F) ----
rna_raw <- read_sheet_table(molecular_xlsx, rna_sheet, "ENSEMBL_ID")
rna_ensembl <- sub("\\.[0-9]+$", "", rna_raw[[find_col(rna_raw, "ENSEMBL_ID", "Ensembl ID")]])
rna_samples <- intersect(names(rna_raw), meta$sample_ID)
tpm <- sapply(rna_raw[, ..rna_samples], function(v) suppressWarnings(as.numeric(v)))
rownames(tpm) <- rna_ensembl
rm(rna_raw); invisible(gc())
cat(sprintf("RNA: %s genes x %d organoids   (expected 19,962 x 255)\n", format(nrow(tpm), big.mark = ","), ncol(tpm)))

# RNA QC: finite, positive total TPM and at least 1,000 detected genes
rna_qc <- data.table(sample_ID = rna_samples,
                     rna_total_TPM = colSums(tpm),
                     rna_detected_genes = colSums(tpm > 0, na.rm = TRUE))
rna_qc[, rna_QC_pass := is.finite(rna_total_TPM) & rna_total_TPM > 0 & rna_detected_genes >= 1000]

# The seven Y markers
marker_rows <- match(y_markers, rna_ensembl)
if (anyNA(marker_rows)) stop("Y markers missing from the RNA table: ", paste(names(y_markers)[is.na(marker_rows)], collapse = ", "))
marker_tpm <- t(tpm[marker_rows, , drop = FALSE])       # organoids x 7 markers
colnames(marker_tpm) <- names(y_markers)

# Per-marker thresholds from sex-consistent female organoids (99th percentile)
female_ids <- intersect(meta[clinical_sex == "Female" & !sex_annotation_conflict, sample_ID], rna_samples)
female_ids <- female_ids[rowSums(!is.finite(marker_tpm[female_ids, , drop = FALSE])) == 0]
female_q99 <- apply(marker_tpm[female_ids, , drop = FALSE], 2, quantile, probs = 0.99, names = FALSE)
cat(sprintf("Female controls for the RNA thresholds: %d   (expected 82)\n", length(female_ids)))

rna_call_list <- list()
for (id in rna_samples) {
  values <- marker_tpm[id, ]
  stable <- TRUE
  for (i in seq_len(nrow(rna_stability_grid))) {
    stable <- stable && rna_low(values, female_q99, floor = rna_stability_grid$floor[i],
                                number = rna_stability_grid$number[i])
  }
  rna_call_list[[id]] <- data.table(
    sample_ID             = id,
    rna_Y_median_TPM      = median(values),
    rna_low_marker_count  = sum(values < pmax(y_rna_floor_tpm, female_q99)),
    rna_measurement_state = if (rna_low(values, female_q99)) "low_Y" else "not_low",
    rna_loss_stable       = stable)
}
rna_calls <- rbindlist(rna_call_list)
marker_cols <- as.data.table(round(marker_tpm[rna_samples, , drop = FALSE], 4))
setnames(marker_cols, paste0("rna_", names(y_markers), "_TPM"))
rna_calls <- cbind(rna_calls, marker_cols)
rna_calls <- merge(rna_calls, rna_qc, by = "sample_ID")
meta <- merge(meta, rna_calls, by = "sample_ID", all.x = TRUE)
meta[, rna_available := sample_ID %in% rna_samples]
meta[rna_available == FALSE, `:=`(rna_measurement_state = "unassessed", rna_loss_stable = FALSE)]
meta[, rna_QC_eligible := rna_available & fcoalesce(rna_QC_pass, FALSE)]

# ---- 5c. Evidence tiers ----
meta[, is_male := fcoalesce(clinical_sex, "") == "Male"]
meta[, evidence_tier := fcase(
  !is_male & sex_annotation_conflict, "female_control_annotation_review",
  !is_male,                           "female_control",
  sex_annotation_conflict,            "male_sex_annotation_review",
  dna_measurement_state == "loss_candidate" & rna_measurement_state == "low_Y" &
    dna_QC_eligible & rna_QC_eligible & dna_loss_stable & rna_loss_stable,
                                      "concordant_stable_high_confidence_candidate",
  dna_measurement_state == "loss_candidate" & rna_measurement_state == "low_Y",      "concordant_lower_confidence_candidate",
  dna_measurement_state == "loss_candidate" & rna_measurement_state == "unassessed", "DNA_supported_RNA_unavailable",
  dna_measurement_state == "loss_candidate" & rna_measurement_state == "not_low",    "discordant_DNA_loss_RNA_not_low",
  dna_measurement_state == "retained"       & rna_measurement_state == "low_Y",      "discordant_DNA_retained_RNA_low",
  dna_measurement_state == "unassessed"     & rna_measurement_state == "low_Y",      "RNA_only_low_Y_candidate",
  dna_measurement_state == "intermediate"   & rna_measurement_state == "low_Y",      "RNA_low_DNA_intermediate",
  dna_measurement_state == "retained"       & rna_measurement_state == "not_low",    "Y_retained",
  dna_measurement_state == "retained"       & rna_measurement_state == "unassessed", "DNA_retained_RNA_unavailable",
  default = "indeterminate")]
meta[, high_confidence_LOY := evidence_tier == "concordant_stable_high_confidence_candidate"]
meta[, Y_retained          := evidence_tier == "Y_retained"]

# Short Y label for tables: LOY, Y-WT, female, or "unresolved" for every other male tier
meta[, Y_status := fcase(!is_male, "female",
                         high_confidence_LOY, label_loy,
                         Y_retained, label_ywt,
                         default = "unresolved_male")]

# Sex/Y group: adjustment covariate in the autosomal comparisons, so that LOY's own effects
# (e.g. on DDX3X) can't show up as arm effects
meta[, sex_Y_group := fcase(clinical_sex == "Female", "female_derived",
                            high_confidence_LOY,      "male_high_confidence_LOY",
                            Y_retained,               "male_Y_retaining",
                            default = "male_other_Y_status")]
cat(sprintf("High-confidence LOY: %d | Y-retained: %d   (expected 62 | 69)\n",
            sum(meta$high_confidence_LOY), sum(meta$Y_retained)))


## ==== 6. Chromosome-arm calls ===================================================================
# Arms are split at the centromere (cytoband "acen" bands); centromeres, assembly gaps and blacklisted
# regions are masked. Relative copy number = total copy number / organoid ploidy. For each organoid and
# arm, the base-pair-weighted median and the fraction of the arm below / within / above the thresholds
# decide the call. X arms are called only in sex-consistent female organoids. Y uses the LOY calls.
cat("\n==== 6. Chromosome-arm calls ====\n")

# ---- 6a. Arm boundaries and the usable (masked) territory of each arm ----
cyto <- fread(get_reference(cytoband_url, file.path(reference_dir, "hg38_cytoBandIdeo.txt.gz")),
              header = FALSE, col.names = c("chr", "start", "end", "band", "stain"))
arm_list <- list()
for (this_chr in c(paste0("chr", 1:22), "chrX")) {
  bands      <- cyto[chr == this_chr]
  centromere <- bands[stain == "acen"]
  arm_list[[this_chr]] <- data.table(event      = paste0(sub("chr", "", this_chr), c("p", "q")),
                                     chromosome = this_chr,
                                     start0     = c(0, max(centromere$end)),
                                     end0       = c(min(centromere$start), max(bands$end)))
}
arm_ref <- rbindlist(arm_list)

territories <- list()
for (i in seq_len(nrow(arm_ref))) {
  this_chrom <- arm_ref$chromosome[i]
  chrom_mask <- mask_all[chromosome == this_chrom]
  territories[[arm_ref$event[i]]] <- complement_within(arm_ref$start0[i], arm_ref$end0[i], chrom_mask)
}
arm_ref[, masked_bases := vapply(event, function(ev) sum(territories[[ev]]$end0 - territories[[ev]]$start0), numeric(1))]
odd_arms <- arm_ref[(masked_bases < 1e6) != (event %in% acrocentric), event]
if (length(odd_arms) > 0) warning("Arms with unexpected usable length (check the cytoband file): ", paste(odd_arms, collapse = ", "))
all_events <- c(arm_ref$event, "Y")   # genome order; used to order every table and figure

# ---- 6b. Measure every organoid x arm ----
supported_seg <- seg[depthWindowCount > 0, .(key = paste(sample_ID, chromosome, sep = "|"), start, end, copyNumber)]
segs_by_key   <- split(supported_seg, by = "key", keep.by = FALSE)   # one table per organoid x chromosome

measure_arm <- function(rows, territory, ploidy) {
  result <- list(coverage = 0, relative_median = NA_real_, fraction_low = NA_real_,
                 fraction_neutral = NA_real_, fraction_gain = NA_real_)
  if (is.null(rows) || nrow(rows) == 0 || !is.finite(ploidy)) return(result)
  bp <- mapply(overlap_bp, rows$start - 1, rows$end, MoreArgs = list(territory = territory))  # 1-based -> BED
  keep <- bp > 0
  if (!any(keep)) return(result)
  relative <- rows$copyNumber[keep] / ploidy
  w <- bp[keep]
  total <- sum(w)
  result$coverage         <- total / sum(territory$end0 - territory$start0)
  result$relative_median  <- weighted_median(relative, w)
  result$fraction_low     <- sum(w[relative < loss_ratio]) / total
  result$fraction_neutral <- sum(w[relative >= neutral_range[1] & relative <= neutral_range[2]]) / total
  result$fraction_gain    <- sum(w[relative > gain_ratio]) / total
  result
}

call_arm <- function(coverage, med, f_low, f_neutral, f_gain) {
  fcase(is.na(coverage) | coverage < arm_min_coverage, "unassessable",
        med < loss_ratio & f_low >= loss_fraction, "loss",
        med >= neutral_range[1] & med <= neutral_range[2] & f_neutral >= neutral_fraction, "neutral",
        med > gain_ratio & f_gain >= gain_fraction, "gain",
        f_low >= partial_fraction & f_low < loss_fraction, "partial_loss",
        default = "intermediate_or_mixed")
}

not_measured <- list(coverage = NA_real_, relative_median = NA_real_, fraction_low = NA_real_,
                     fraction_neutral = NA_real_, fraction_gain = NA_real_)
arm_rows <- vector("list", nrow(arm_ref) * nrow(meta))
k <- 0
for (a in seq_len(nrow(arm_ref))) {
  ev       <- arm_ref$event[a]
  chrom    <- arm_ref$chromosome[a]
  is_x_arm <- startsWith(ev, "X")
  for (i in seq_len(nrow(meta))) {
    id <- meta$sample_ID[i]
    if (is_x_arm) {
      in_population <- fcoalesce(meta$clinical_sex[i], "") == "Female" && !meta$sex_annotation_conflict[i]
    } else {
      in_population <- TRUE
    }
    attempted <- in_population && !(ev %in% acrocentric) && !meta$ambiguous_baseline_DNA_profile[i]
    if (attempted) {
      m <- measure_arm(segs_by_key[[paste(id, chrom, sep = "|")]], territories[[ev]], meta$ploidy[i])
    } else {
      m <- not_measured
    }
    k <- k + 1
    arm_rows[[k]] <- c(list(sample_ID = id, event = ev, in_population = in_population), m)
  }
}
arm_status <- rbindlist(arm_rows)

# ---- 6c. Final status. Later lines override earlier ones ----
arm_status <- merge(arm_status, meta[, .(sample_ID, ambiguous_baseline_DNA_profile, dna_base_qc)], by = "sample_ID")
arm_status[, measured_state := call_arm(coverage, relative_median, fraction_low, fraction_neutral, fraction_gain)]
arm_status[, status := measured_state]
arm_status[dna_base_qc == FALSE, status := "QC_limited"]
arm_status[ambiguous_baseline_DNA_profile == TRUE | event %in% acrocentric, status := "unassessable"]
arm_status[in_population == FALSE, status := "out_of_population"]
arm_status[, c("ambiguous_baseline_DNA_profile", "dna_base_qc") := NULL]

# ---- 6d. Y as one more event: loss = high-confidence LOY, neutral = Y-retained with both QCs passing ----
y_rows <- meta[, .(sample_ID, event = "Y",
                   in_population    = is_male & !sex_annotation_conflict,
                   coverage         = dna_Y_coverage,
                   relative_median  = NA_real_, fraction_low = NA_real_,
                   fraction_neutral = NA_real_, fraction_gain = NA_real_,
                   measured_state   = evidence_tier,
                   status = fcase(!(is_male & !sex_annotation_conflict), "out_of_population",
                                  high_confidence_LOY, "loss",
                                  Y_retained & dna_QC_eligible & rna_QC_eligible, "neutral",
                                  Y_retained, "QC_limited",
                                  default = "intermediate_or_mixed"))]
arm_status <- rbind(arm_status, y_rows, use.names = TRUE)
arm_status[, event := factor(event, levels = all_events)]
setorder(arm_status, sample_ID, event)
cat(sprintf("8p across all organoids: %d loss, %d neutral   (expected 83, 75)\n",
            arm_status[event == "8p" & status == "loss", .N], arm_status[event == "8p" & status == "neutral", .N]))


## ==== 7. Comparison cohorts =====================================================================
# For each arm (and Y): screened organoids with the arm altered vs copy-number neutral, with complete
# covariates and no sex conflict, restricted to cancer types with >= 3 altered AND >= 3 neutral organoids.
# A comparison is tested if it has >= 25 altered, >= 20 neutral organoids and >= 40 donors, and its
# design matrix is stable.
cat("\n==== 7. Comparison cohorts ====\n")
screened <- merge(screens, screen_qc, by = c("sample_ID", "library"), all.x = TRUE)
screened <- merge(screened, meta[, .(sample_ID, individual_ID, primary_tumour_type, clinical_sex, sex_Y_group,
                                     sex_annotation_conflict, msStatus, ploidy)], by = "sample_ID")

build_cohort <- function(ev, altered_state) {
  calls <- arm_status[event == ev, .(sample_ID, in_population, status)]
  d <- merge(screened, calls, by = "sample_ID")
  d[, altered := fcase(status == altered_state, 1L, status == "neutral", 0L, default = NA_integer_)]
  d[, eligible := in_population & !is.na(altered) & !sex_annotation_conflict & !is.na(AUC_PR) &
                  is.finite(ploidy) & msStatus %in% c("MSI", "MSS") & !is.na(individual_ID) &
                  !is.na(primary_tumour_type) & !is.na(library) & !is.na(sex_Y_group)]
  d <- d[eligible == TRUE]

  # Keep cancer types that have enough organoids in both groups
  per_type   <- d[, .(n_altered = sum(altered == 1), n_neutral = sum(altered == 0)), by = primary_tumour_type]
  keep_types <- per_type[n_altered >= min_per_lineage & n_neutral >= min_per_lineage, primary_tumour_type]
  d <- d[primary_tumour_type %in% keep_types]
  setorder(d, sample_ID)
  d
}

comparisons <- rbind(data.table(type = "Arm loss", event = all_events,    altered_state = "loss"),
                     data.table(type = "Arm gain", event = arm_ref$event, altered_state = "gain"))
comparisons[, key := paste(type, event)]

cohorts <- list()   # tested comparisons only, by key (e.g. "Arm loss 8p")
designs <- list()
eligibility_rows <- list()
for (k in seq_len(nrow(comparisons))) {
  this_key   <- comparisons$key[k]
  this_event <- comparisons$event[k]
  this_state <- comparisons$altered_state[k]

  frame <- build_cohort(this_event, this_state)
  n_alt <- sum(frame$altered == 1)
  n_neu <- sum(frame$altered == 0)
  n_don <- uniqueN(frame$individual_ID)

  reasons <- c(if (n_alt < min_altered) paste0("altered_below_", min_altered),
               if (n_neu < min_neutral) paste0("neutral_below_", min_neutral),
               if (n_don < min_donors)  paste0("donors_below_", min_donors))
  if (length(reasons) == 0) {
    design <- make_design(frame, use_sex_y = this_event != "Y")
    if (is.null(design$X)) {
      reasons <- design$reason
    } else {
      cohorts[[this_key]] <- frame
      designs[[this_key]] <- design$X
    }
  }

  calls_all <- arm_status[event == this_event]
  eligibility_rows[[k]] <- data.table(
    type = comparisons$type[k], event = this_event,
    biobank_altered  = calls_all[status == this_state, .N],
    biobank_neutral  = calls_all[status == "neutral", .N],
    screened_altered = calls_all[status == this_state & sample_ID %in% screens$sample_ID, .N],
    n_altered = n_alt, n_neutral = n_neu, n_donors = n_don,
    cancer_types = paste(sort(unique(frame$primary_tumour_type)), collapse = ";"),
    tested = length(reasons) == 0,
    reasons_not_tested = paste(reasons, collapse = ";"))
}
eligibility <- rbindlist(eligibility_rows)

tested_loss <- eligibility[type == "Arm loss" & tested == TRUE, event]
tested_gain <- eligibility[type == "Arm gain" & tested == TRUE, event]
cat("Tested arm losses:", paste(tested_loss, collapse = ", "), "\n   (expected 4p, 4q, 8p, 17p, 18p, 18q, 21q, 22q, Y)\n")
cat("Tested arm gains: ", paste(tested_gain, collapse = ", "), "\n   (expected 7p, 7q, 8q, 13q, 20p, 20q)\n")
print(eligibility[tested == TRUE, .(type, event, n_altered, n_neutral, n_donors, cancer_types)])
y_cohort <- cohorts[["Arm loss Y"]]
if (is.null(y_cohort)) stop("The LOY comparison (Arm loss Y) was not tested; panels B, E and F need it")
cat(sprintf("LOY comparison: %d LOY vs %d Y-WT organoids from %d donors   (expected 40 vs 42, 81 donors)\n",
            sum(y_cohort$altered == 1), sum(y_cohort$altered == 0), uniqueN(y_cohort$individual_ID)))


## ==== 8. Differential dependency ================================================================
# For each tested comparison and every gene:
#   dependency ~ altered + cancer type + sex/Y group + library + MSI + ploidy + screen quality
# OLS with donor-clustered standard errors. Relative dependency = -coefficient (positive = more
# essential in the altered organoids). Hit = p < 0.01 and relative dependency >= 0.2.
cat("\n==== 8. Differential dependency ====\n")
screen_results <- list()
for (this_key in names(designs)) {
  info  <- comparisons[match(this_key, comparisons$key)]
  frame <- cohorts[[this_key]]
  dep   <- fitness[frame$sample_ID, , drop = FALSE]
  fit   <- fit_genes(designs[[this_key]], dep, frame$individual_ID)
  fit[, `:=`(type = info$type, event = info$event,
             n_altered = sum(frame$altered == 1), n_neutral = sum(frame$altered == 0))]
  screen_results[[this_key]] <- fit
}
res <- rbindlist(screen_results)
res <- res[!is.na(p)]

# ---- Gene positions (Ensembl GRCh38 release 98), to note genes that sit on the altered arm ----
gene_pos_file <- file.path(reference_dir, "ensembl98_gene_positions.csv")
if (!file.exists(gene_pos_file)) {
  gtf_path <- get_reference(gtf_url, file.path(reference_dir, "Homo_sapiens.GRCh38.98.chr.gtf.gz"))
  con <- gzfile(gtf_path, "r")
  gene_lines <- character(0)
  repeat {   # read in chunks and keep only "gene" lines, so memory stays small
    lines <- readLines(con, n = 250000)
    if (length(lines) == 0) break
    gene_lines <- c(gene_lines, lines[grepl("\tgene\t", lines, fixed = TRUE)])
  }
  close(con)
  parts <- tstrsplit(gene_lines, "\t", fixed = TRUE)
  gene_pos <- data.table(ensembl_id = sub('.*gene_id "([^"]+)".*', "\\1", parts[[9]]),
                         chromosome = paste0("chr", parts[[1]]),
                         start = as.numeric(parts[[4]]), end = as.numeric(parts[[5]]),
                         feature = parts[[3]])
  gene_pos <- gene_pos[feature == "gene", .(ensembl_id, chromosome, start, end)]
  fwrite(gene_pos, gene_pos_file)
}
gene_pos <- fread(gene_pos_file)

# A gene is "on" an arm if it overlaps it (genes across a centromere count for both arms); Y = chrY genes
gene_arm_rows <- list()
for (i in seq_len(nrow(arm_ref))) {
  arm_i <- arm_ref[i]
  gene_arm_rows[[i]] <- gene_pos[chromosome == arm_i$chromosome & start - 1 < arm_i$end0 & end > arm_i$start0,
                                 .(ensembl_id, arm = arm_i$event)]
}
gene_arm_rows[[length(gene_arm_rows) + 1]] <- gene_pos[chromosome == "chrY", .(ensembl_id, arm = "Y")]
gene_arms <- rbindlist(gene_arm_rows)

gene_annot <- data.table(ensembl_id = colnames(fitness), gene = unname(symbol[colnames(fitness)]))
gene_chrom <- gene_pos[, .(chromosome = paste(unique(chromosome), collapse = ";")), by = ensembl_id]
gene_arm_text <- gene_arms[, .(gene_arm = paste(unique(arm), collapse = ";")), by = ensembl_id]
gene_annot <- merge(gene_annot, gene_chrom, by = "ensembl_id", all.x = TRUE)
gene_annot <- merge(gene_annot, gene_arm_text, by = "ensembl_id", all.x = TRUE)
gene_annot[is.na(chromosome), `:=`(chromosome = "unknown", gene_arm = "unknown")]
gene_annot[is.na(gene_arm), gene_arm := "centromere_or_unplaced"]

res <- merge(res, gene_annot, by = "ensembl_id", all.x = TRUE)
res[, on_altered_arm := mapply(function(arms, ev) ev %in% strsplit(arms, ";")[[1]], gene_arm, event)]

# ---- Direction, hits, and BH q values (reference only; not used for hits) ----
res[, rel_dep         := -beta]        # neutral - altered: positive = more essential when altered
res[, rel_dep_CI_low  := -CI95_high]
res[, rel_dep_CI_high := -CI95_low]
res[, hit := p < p_cut & rel_dep >= effect_cut]
res[, q_BH_reference_only := p.adjust(p, method = "BH"), by = type]
res[, event := factor(event, levels = all_events)]
setorder(res, type, event, p)
res[, event := as.character(event)]

hits_by_comparison <- res[, .(tests = .N, nominal_p01 = sum(p < p_cut), expected_by_chance = round(.N * p_cut),
                              hits = sum(hit)), by = .(type, event)]
cat("nominal_p01 counts both directions; about 1% of tests pass p < 0.01 by chance alone.\n")
print(hits_by_comparison)

# ---- Checks against the earlier LOY analysis ----
x_id <- names(symbol)[symbol == focus_x_gene][1]
if (is.na(x_id) || !(x_id %in% colnames(fitness))) stop(focus_x_gene, " is not among the screened genes")
print(res[type == "Arm loss" & event == "Y" & gene %in% c("DDX3X", "RPS4X"),
          .(gene, beta, CI95_low, CI95_high, p, rel_dep)], digits = 3)
cat("   (expected DDX3X beta -0.603 [-0.793, -0.413]; RPS4X beta -0.668 [-0.901, -0.436])\n")

# The fast fit must agree with lm() (shown for the focus gene in the LOY comparison)
X_y <- designs[["Arm loss Y"]]
lm_fit <- lm(fitness[y_cohort$sample_ID, x_id] ~ X_y - 1)   # "- 1": X_y already has an intercept column
lm_beta <- unname(coef(lm_fit)["X_yaltered"])
fast_beta <- res[type == "Arm loss" & event == "Y" & ensembl_id == x_id, beta]
if (abs(lm_beta - fast_beta) > 1e-8) stop("Fast fit disagrees with lm() for ", focus_x_gene)


## ==== 9. Expression of genes on the altered arm (panel C data) =================================
# Same comparison cohorts as the regressions. For every tested comparison and gene:
#   relative expression = median log2(TPM + 1) in altered organoids - median in neutral organoids
# (Yi's "relative expression", manuscript Figs 2L, 3G, 4J). Genes are split into those on the altered
# arm and those on all other autosomal arms. Only genes expressed in the neutral group are used.
cat("\n==== 9. Expression of genes on the altered arm ====\n")
panelC_keys <- names(cohorts)
if (!include_Y_panelC) panelC_keys <- setdiff(panelC_keys, "Arm loss Y")

panelC_samples <- character(0)
for (this_key in panelC_keys) panelC_samples <- union(panelC_samples, cohorts[[this_key]]$sample_ID)
n_without_rna  <- length(setdiff(panelC_samples, colnames(tpm)))
panelC_samples <- intersect(panelC_samples, colnames(tpm))
cat(sprintf("%d cohort organoids have RNA; %d without RNA are left out of panel C only\n",
            length(panelC_samples), n_without_rna))

expr <- log2(tpm[, panelC_samples, drop = FALSE] + 1)
keep_genes <- !duplicated(rownames(expr)) & rowSums(!is.finite(expr)) == 0
expr <- expr[keep_genes, , drop = FALSE]

# Genes wholly inside one arm (genes across a centromere are left out of this panel)
gene_pos_expr <- unique(gene_pos[ensembl_id %in% rownames(expr)], by = "ensembl_id")
contained_rows <- list()
for (i in seq_len(nrow(arm_ref))) {
  arm_i <- arm_ref[i]
  contained_rows[[i]] <- gene_pos_expr[chromosome == arm_i$chromosome & start - 1 >= arm_i$start0 & end <= arm_i$end0,
                                       .(ensembl_id, arm = arm_i$event, chromosome)]
}
genes_contained <- rbindlist(contained_rows)
autosomes <- paste0("chr", 1:22)

expr_rows <- list()
for (this_key in panelC_keys) {
  info  <- comparisons[match(this_key, comparisons$key)]
  frame <- cohorts[[this_key]][sample_ID %in% colnames(expr)]
  altered_ids <- frame[altered == 1, sample_ID]
  neutral_ids <- frame[altered == 0, sample_ID]

  if (info$event == "Y") {
    on_ids <- gene_pos_expr[chromosome == "chrY", ensembl_id]
  } else {
    on_ids <- genes_contained[arm == info$event, ensembl_id]
  }
  other_ids <- genes_contained[arm != info$event & chromosome %in% autosomes, ensembl_id]
  other_ids <- setdiff(other_ids, on_ids)
  use_ids   <- c(on_ids, other_ids)

  altered_median <- apply(expr[use_ids, altered_ids, drop = FALSE], 1, median)
  neutral_median <- apply(expr[use_ids, neutral_ids, drop = FALSE], 1, median)
  expr_rows[[this_key]] <- data.table(
    type = info$type, event = info$event, ensembl_id = use_ids,
    location = c(rep("On altered arm", length(on_ids)), rep("Elsewhere", length(other_ids))),
    altered_median_log2TPM1 = altered_median, neutral_median_log2TPM1 = neutral_median,
    rel_expr = altered_median - neutral_median,
    n_altered_rna = length(altered_ids), n_neutral_rna = length(neutral_ids))
}
expr_res <- rbindlist(expr_rows)
expr_res <- expr_res[neutral_median_log2TPM1 >= log2(min_tpm_neutral + 1)]   # expressed genes only
expr_res[, gene := unname(symbol[ensembl_id])]                                 # NA if not in the CRISPR library

# One row per comparison: medians and a two-sided Wilcoxon test, on-arm genes vs genes elsewhere
expr_summary_rows <- list()
for (this_key in panelC_keys) {
  info <- comparisons[match(this_key, comparisons$key)]
  d <- expr_res[type == info$type & event == info$event]
  on_arm    <- d[location == "On altered arm", rel_expr]
  elsewhere <- d[location == "Elsewhere", rel_expr]
  p_value <- NA_real_
  if (length(on_arm) > 0 && length(elsewhere) > 0) p_value <- wilcox.test(on_arm, elsewhere, exact = FALSE)$p.value
  expr_summary_rows[[this_key]] <- data.table(
    type = info$type, event = info$event,
    n_altered_rna = d$n_altered_rna[1], n_neutral_rna = d$n_neutral_rna[1],
    genes_on_arm = length(on_arm), genes_elsewhere = length(elsewhere),
    median_on_arm = median(on_arm), median_elsewhere = median(elsewhere), p_wilcoxon = p_value)
}
expr_summary <- rbindlist(expr_summary_rows)
expr_summary[, concordant := fifelse(type == "Arm loss", median_on_arm < median_elsewhere,
                                     median_on_arm > median_elsewhere)]
print(expr_summary, digits = 3)
cat(sprintf("Concordant direction in %d of %d comparisons\n", sum(expr_summary$concordant), nrow(expr_summary)))


## ==== 10. Figures ===============================================================================
cat("\n==== 10. Figures ====\n")
group_fill <- setNames(c(col_ywt, col_loy), c(label_ywt, label_loy))
group_dark <- setNames(c(col_ywt_dark, col_loy_dark), c(label_ywt, label_loy))
x_label_theme <- theme(axis.text.x = element_text(angle = x_label_angle,
                                                  hjust = if (x_label_angle == 0) 0.5 else 1,
                                                  vjust = if (x_label_angle == 90) 0.5 else 1))
dependency_title <- bquote(italic(.(focus_x_gene)) ~ "dependency (CRISPR)")

# The LOY comparison cohort, with group labels and the focus gene's scores (panels B, E, F)
y_data <- y_cohort[, .(sample_ID, individual_ID, primary_tumour_type, library, altered)]
y_data[, group := factor(fifelse(altered == 1, label_loy, label_ywt), levels = c(label_ywt, label_loy))]
y_data[, dependency := fitness[sample_ID, x_id]]

# ---- Panel B: LOY and Y-WT organoids in the comparison cohort ----
# Bars in Yi's style: fill and outline in the same colour, fill at 70% (geom_col(..., alpha = 0.7))
b_data <- y_data[, .(n = .N), by = group]
b_data[, fraction := n / sum(n)]
b_data[, label := sprintf("%d (%.0f%%)", n, 100 * fraction)]
setorder(b_data, group)
cat(sprintf("Panel B: %s\n", paste(b_data$group, b_data$label, collapse = " | ")))

p_B <- ggplot(b_data, aes(group, fraction)) +
  geom_col(aes(fill = group, colour = group), width = 0.5, alpha = 0.7, linewidth = 0.4) +
  geom_text(aes(label = label), vjust = -0.5, size = txt(fig_base_pt), family = fig_font) +
  scale_fill_manual(values = group_fill, guide = "none") +
  scale_colour_manual(values = group_fill, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, max(b_data$fraction) * 1.2),
                     breaks = scales::pretty_breaks(4), expand = expansion(mult = c(0, 0))) +
  labs(x = sprintf("Male organoids in the\ncomparison cohort (n = %d)", sum(b_data$n)), y = "Organoids") +
  capped_axes + pub_theme + x_label_theme
save_figure(p_B, "panel_B_loy_fraction", size_B)
fwrite(b_data, file.path(figure_dir, "panel_B_source_data.csv"))

# ---- Panel B (alternative): Y expression of every organoid in the comparison cohort, ranked ----
# One dot per organoid. x = median expression of the seven Y marker genes, log2(TPM + 1), which is the
# quantity behind the RNA part of the LOY call
b2_data <- y_data[sample_ID %in% rownames(marker_tpm), .(sample_ID, primary_tumour_type, group)]
b2_log_markers <- log2(marker_tpm[b2_data$sample_ID, , drop = FALSE] + 1)   # organoids x 7 markers
b2_data[, y_expression := apply(b2_log_markers, 1, median)]
b2_data[, rank := frank(y_expression, ties.method = "first")]               # 1 = lowest Y expression
print(b2_data[, .N, by = .(primary_tumour_type, group)][order(primary_tumour_type, group)])

# Colours by cancer type (any cancer type not listed here is drawn grey)
b2_cancer_cols <- c(Colorectal = "#2166AC", Oesophageal = "#B2182B")   
b2_types <- sort(unique(b2_data$primary_tumour_type))
b2_palette <- setNames(rep("grey55", length(b2_types)), b2_types)
b2_known <- intersect(b2_types, names(b2_cancer_cols))
b2_palette[b2_known] <- b2_cancer_cols[b2_known]

# x-axis ticks that bracket the data, so the capped axis line covers every dot
b2_x_breaks <- pretty(c(0, b2_data$y_expression))

b2_alpha <- 0.5

p_B2 <- ggplot(b2_data, aes(y_expression, rank)) +
  geom_vline(xintercept = log2(0.5 + 1), linetype = "dotted", linewidth = 0.4, colour = "grey40") +   # 0.5 TPM
  geom_point(aes(colour = primary_tumour_type), size = 1.2, alpha = b2_alpha, stroke = 0) +
  scale_colour_manual(values = b2_palette, name = NULL) +
  scale_x_continuous(breaks = b2_x_breaks) +
  scale_y_reverse(breaks = NULL, expand = expansion(mult = 0.02)) +      # rank 1 (lowest Y) at the top
  coord_cartesian(xlim = range(b2_x_breaks)) +
  guides(x = guide_axis(cap = "both"),
         colour = guide_legend(override.aes = list(alpha = b2_alpha, size = 2))) +
  labs(x = "Y gene expression\n(median of 7 genes, log2[TPM + 1])",
       y = sprintf("Organoids, ranked (n = %d)", nrow(b2_data))) +
  pub_theme +
  theme(legend.position = "inside", legend.position.inside = c(0.97, 0.70),
        legend.justification.inside = c(1, 0))
save_figure(p_B2, "panel_B_y_expression_ranked", size_B)
fwrite(b2_data[, .(sample_ID, primary_tumour_type, group, y_expression_log2TPM1 = y_expression, rank)],
       file.path(figure_dir, "panel_B_ranked_source_data.csv"))

# ---- Panel C: expression of genes on the altered arm ----
expr_res[, event_f    := factor(event, levels = all_events)]
expr_res[, type_f     := factor(type, levels = c("Arm loss", "Arm gain"))]
expr_res[, fill_group := factor(fifelse(location == "Elsewhere", "Elsewhere", type),
                                levels = c("Elsewhere", "Arm loss", "Arm gain"))]   
box_cols   <- c("Arm loss" = col_yi_orange, "Arm gain" = col_yi_blue, Elsewhere = col_elsewhere)
box_labels <- c("Arm loss" = "Genes on the lost arm", "Arm gain" = "Genes on the gained arm",
                Elsewhere = "Genes on other autosomal arms")

whiskers <- expr_res[, {
  q <- quantile(rel_expr, c(0.25, 0.75))
  iqr <- q[2] - q[1]
  .(lo = max(min(rel_expr), q[1] - 1.5 * iqr), hi = min(max(rel_expr), q[2] + 1.5 * iqr))
}, by = .(type, event, location)]
c_ylim <- range(c(whiskers$lo, whiskers$hi, 0))
c_span <- diff(c_ylim)

p_C <- ggplot(expr_res, aes(event_f, rel_expr)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_boxplot(aes(fill = fill_group, colour = fill_group), position = position_dodge2(padding = 0.15),
               width = 0.75, outlier.shape = NA, linewidth = 0.4, alpha = 0.4) +
  facet_grid(. ~ type_f, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = box_cols, labels = box_labels, breaks = c("Arm loss", "Arm gain", "Elsewhere"), name = NULL) +
  scale_colour_manual(values = box_cols, labels = box_labels, breaks = c("Arm loss", "Arm gain", "Elsewhere"), name = NULL) +
  labs(x = NULL, y = "Rel. gene expression\n(log2 fold change, altered vs. neutral)") +
  guides(x = guide_axis(cap = "both"), y = guide_axis(cap = "both")) +
  pub_theme + theme(legend.position = "none", panel.spacing = unit(4, "mm"))

p_stars <- function(p) fcase(p < 1e-4, "****", p < 1e-3, "***", p < 0.01, "**", p < 0.05, "*", default = "ns")
if (show_p_panelC) {
  # One bracket per comparison, over its two boxes
  brackets <- whiskers[, .(y = max(hi) + 0.05 * c_span), by = .(type, event)]
  brackets[, y := max(whiskers$hi) + 0.05 * c_span]
  brackets[, y := max(whiskers$hi) + 0.05 * c_span]   # one height for every bracket, in both panels
  brackets <- merge(brackets, expr_summary[, .(type, event, p_wilcoxon)], by = c("type", "event"))
  brackets[, type_f := factor(type, levels = levels(expr_res$type_f))]
  brackets[, label := p_stars(p_wilcoxon)]
  # x position of each arm within its panel (each panel numbers its arms 1, 2, ... in genome order)
  brackets[, x := as.integer(factor(event, levels = intersect(all_events, event))), by = type]
  c_ylim[2] <- max(brackets$y) + 0.12 * c_span
  p_C <- p_C +
    geom_segment(data = brackets, aes(x = x - 0.19, xend = x + 0.19, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.3) +
    geom_text(data = brackets, aes(x = x, y = y, label = label), inherit.aes = FALSE, vjust = -0.4,
              size = txt(fig_base_pt - 2), family = fig_font)
}

c_key <- data.table(type_f = factor(c("Arm loss", "Arm loss", "Arm gain", "Arm gain"), levels = levels(expr_res$type_f)),
                    label  = c("Lost arm genes", "Other arm genes", "Gained arm genes", "Other arm genes"),
                    col    = c(col_yi_orange, col_elsewhere, col_yi_blue, col_elsewhere),
                    vj     = c(-1.6, -0.3, -1.6, -0.3))
p_C <- p_C + geom_text(data = c_key, aes(x = Inf, y = -Inf, label = label, colour = I(col), vjust = vj),
                       inherit.aes = FALSE, hjust = 1, size = txt(fig_base_pt), family = fig_font)

c_ylim[1] <- c_ylim[1] - 0.05 * c_span
c_breaks <- pretty(c_ylim)
p_C <- p_C + scale_y_continuous(breaks = c_breaks, expand = expansion(0)) +
  coord_cartesian(ylim = range(c_breaks), clip = "off")

save_figure(p_C, "panel_C_expression_concordance", size_C)
fwrite(expr_summary, file.path(figure_dir, "panel_C_source_data_summary.csv"))
fwrite(expr_res[, .(type, event, ensembl_id, gene, location, rel_expr, altered_median_log2TPM1,
                    neutral_median_log2TPM1, n_altered_rna, n_neutral_rna)],
       file.path(figure_dir, "panel_C_source_data_genes.csv.gz"))

# ---- Panels B and C side by side, with aligned plot areas ----
p_row2 <- cowplot::plot_grid(p_B2, p_C, nrow = 1, align = "h", axis = "tb",
                             rel_widths = c(size_B[1], size_C[1]))
save_figure(p_row2, "panels_BC_combined", c(size_B[1] + size_C[1], size_C[2]))

# ---- Panel D: volcano plot of every gene x arm-loss and gene x arm-gain test ----
d_data <- copy(res)
d_data[, neg_log10_p := -log10(p)]
d_data[, lab_key := paste0(event, ": ", gene)]
d_ns   <- d_data[hit == FALSE]
d_hits <- d_data[hit == TRUE]

# Non-hits as a rasterised layer (keeps the SVG/PDF small); without ggrastr, thin the p > 0.05 cloud
has_ragg   <- requireNamespace("ragg", quietly = TRUE)
use_raster <- requireNamespace("ggrastr", quietly = TRUE) && (has_ragg || requireNamespace("Cairo", quietly = TRUE))
if (!use_raster) {
  set.seed(1)
  dense <- d_ns[p > 0.05]
  d_ns  <- rbind(d_ns[p <= 0.05], dense[sample(nrow(dense), ceiling(0.2 * nrow(dense)))])
  message("ggrastr not available: drawing all non-hits with p <= 0.05 and 20% of the rest. ",
          "install.packages(c(\"ggrastr\", \"ragg\")) to rasterise instead.")
}
ns_layer <- geom_point(data = d_ns, colour = col_ns, size = volcano_ns_size, stroke = 0)
if (use_raster) ns_layer <- ggrastr::rasterise(ns_layer, dpi = 600, dev = if (has_ragg) "ragg" else "cairo")

d_labels <- rbind(head(d_hits[order(p)], volcano_label_n), d_data[lab_key %in% volcano_always])
d_labels <- unique(d_labels, by = c("type", "event", "ensembl_id"))
d_xlim   <- c(-1, 1) * max(abs(d_data$rel_dep))   # symmetric about zero, as in Yi's figure

d_breaks <- pretty(c(-1, 1) * max(abs(d_data$rel_dep)))
d_ybreaks <- pretty(c(0, max(d_data$neg_log10_p)))
p_D <- ggplot(d_data, aes(rel_dep, neg_log10_p)) +
  ns_layer +
  geom_point(data = d_hits, colour = col_yi_orange, shape = 16, size = volcano_hit_size,
             alpha = 0.75, stroke = 0) +
  geom_point(data = d_labels, colour = col_yi_orange, shape = 16, size = volcano_focus_size, stroke = 0) +  # new: DDX3X and RPS4X, larger and opaque
  ggrepel::geom_text_repel(data = d_labels, aes(label = gene), size = txt(fig_base_pt + 1), family = fig_font,
                           segment.size = 0.3, segment.colour = "black", box.padding = 0.6,
                           point.size = volcano_focus_size,   # new: keeps labels clear of the larger dots
                           max.overlaps = Inf, seed = 1) +
  scale_x_continuous(limits = range(d_breaks), breaks = d_breaks) +
  scale_y_continuous(breaks = d_ybreaks, expand = expansion(0)) +
  coord_cartesian(ylim = range(d_ybreaks), clip = "off") +
  labs(x = "Relative dependency, adjusted\n(arm altered vs. arm neutral)",
       y = expression("\u2212" * log[10] ~ italic(P))) +
  guides(x = guide_axis(cap = "both"), y = guide_axis(cap = "both")) +
  pub_theme
save_figure(p_D, "panel_D_volcano", size_D)
fwrite(d_data[, .(type, event, gene, ensembl_id, rel_dep, p, hit)], file.path(figure_dir, "panel_D_source_data.csv.gz"))

# ---- Panel E: focus-gene dependency in LOY vs Y-WT organoids ----
# Raw scores, one screen per organoid, pooled over cancer types. The p value over the boxes follows
# panelE_p_method; the second line gives the adjusted relative dependency from the regression.
model_row <- res[type == "Arm loss" & event == "Y" & ensembl_id == x_id]
# Two-sided Wilcoxon rank-sum test on the plotted scores, with the same settings as Yi's Figure 1 script
p_wilcox <- wilcox.test(y_data[altered == 1, dependency], y_data[altered == 0, dependency],
                        exact = FALSE, correct = FALSE)$p.value
e_label <- format_p(p_wilcox)
e_bracket_y <- max(y_data$dependency) + 0.07 * diff(range(y_data$dependency))
e_counts    <- y_data[, .N, by = group]

e_breaks <- pretty(c(y_data$dependency, e_bracket_y))
p_E <- ggplot(y_data, aes(group, dependency)) +
  geom_boxplot(aes(fill = group, colour = group), width = 0.6, outlier.shape = NA, linewidth = 0.3,
               alpha = 0.45, fatten = 1.6) +
  geom_point(aes(colour = group, fill = after_scale(alpha(colour, 0.35))), shape = 21, size = 0.8, stroke = 0.3,
             position = position_jitter(width = 0.15, height = 0, seed = 1)) +
  annotate("segment", x = 1, xend = 2, y = e_bracket_y, yend = e_bracket_y, linewidth = 0.3) +
  annotate("text", x = 1.5, y = e_bracket_y, label = e_label, vjust = -0.3, lineheight = 0.95,
           size = txt(fig_base_pt - 0.5), family = fig_font) +
  #geom_text(data = e_counts, aes(y = -Inf, label = paste0("n = ", N)), vjust = -0.6,
            #size = txt(fig_base_pt - 1.5), family = fig_font, colour = "grey35") +
  scale_fill_manual(values = group_fill, guide = "none") +
  scale_colour_manual(values = group_dark, guide = "none") +
  scale_y_continuous(breaks = e_breaks, expand = expansion(0)) +
  scale_x_discrete(labels = function(g) sprintf("%s\n(n = %d)", g, e_counts$N[match(g, e_counts$group)])) +
  coord_cartesian(ylim = range(e_breaks), clip = "off") +
  labs(x = NULL, y = dependency_title) +
  guides(x = guide_axis(cap = "both"), y = guide_axis(cap = "both")) + pub_theme + x_label_theme + theme(plot.margin = margin(5, 3, 2, 2, "mm"))
save_figure(p_E, paste0("panel_E_", tolower(focus_x_gene), "_dependency_LOY_vs_YWT"), size_E)
fwrite(y_data[, .(sample_ID, individual_ID, primary_tumour_type, library, group, dependency)],
       file.path(figure_dir, "panel_E_source_data.csv"))

# ---- Panel F: focus-gene dependency vs its Y paralog's expression ----
# Styled after Yi's scatter plots: points coloured along the first principal component of the two
# variables (MetBrewer "Hokusai3"), a fitted regression line with its 95% band, Pearson r and p.
# Expression is log2(TPM + 1) (Yi's DepMap axis says "TPM" but those values are also log2(TPM + 1)).
f_data <- copy(y_data)[sample_ID %in% rownames(marker_tpm)]
f_data[, expression := log2(marker_tpm[sample_ID, focus_y_gene] + 1)]
pca <- prcomp(f_data[, .(dependency, expression)])
f_data[, pc := pca$x[, 1]]
if (cor(f_data$pc, f_data$expression) < 0) f_data[, pc := -pc]   # low expression = start of the palette

if (requireNamespace("MetBrewer", quietly = TRUE)) {
  f_palette <- MetBrewer::met.brewer("Hokusai3", n = 7, type = "continuous")
} else {
  f_palette <- grDevices::colorRampPalette(c("#d8d97a", "#95c36e", "#74c8c3", "#5a97c1", "#295384", "#0a2e57"))(7)
}
pearson <- cor.test(f_data$expression, f_data$dependency)
f_label <- sprintf("r = %.4f\n%s", unname(pearson$estimate), format_p(pearson$p.value))
cat(sprintf("Panel F: Pearson r = %.3f, %s (n = %d)\n", pearson$estimate, format_p(pearson$p.value), nrow(f_data)))

f_breaks <- pretty(f_data$dependency)
f_xbreaks <- pretty(f_data$expression)
p_F <- ggplot(f_data, aes(expression, dependency)) +
  geom_point(aes(colour = pc, fill = after_scale(alpha(colour, 0.35))), shape = 21, size = 1.1, stroke = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#3366FF", fill = "grey60", alpha = 0.25, linewidth = 0.3) +
  annotate("text", x = max(f_xbreaks), y = min(f_breaks) + c(0.08, 0) * diff(range(f_breaks)),
           label = c(sprintf('italic(r)~"="~"%.2f"', unname(pearson$estimate)), fmt_P(pearson$p.value)),
           parse = TRUE, hjust = 1, vjust = 0, size = txt(fig_base_pt), family = fig_font) +
  scale_colour_gradient2(low = "#FCAE1E", mid = "#E2DDFC", high = "#543D7B", midpoint = 0, guide = "none") +
  scale_x_continuous(breaks = f_xbreaks) +
  scale_y_continuous(breaks = f_breaks, expand = expansion(0)) +
  coord_cartesian(xlim = range(f_xbreaks), ylim = range(f_breaks), clip = "off") +
  labs(x = bquote(italic(.(focus_y_gene)) ~ "expression (log"[2] * "[TPM + 1])"), y = dependency_title) +
  guides(x = guide_axis(cap = "both"), y = guide_axis(cap = "both")) +
  pub_theme
save_figure(p_F, paste0("panel_F_", tolower(focus_x_gene), "_vs_", tolower(focus_y_gene)), size_F)
fwrite(f_data[, .(sample_ID, individual_ID, primary_tumour_type, group, expression_log2TPM1 = expression, dependency)],
       file.path(figure_dir, "panel_F_source_data.csv"))

## ==== Panel G: benchmark DDX3X / LOY against clinical biomarker-target pairs (Yi's Figure 6C) ======

# Each point compares the dependency on a target gene between biomarker-positive and biomarker-negative
# screened organoids:
#   * KRAS, PIK3CA, BRAF: activating driver mutation (missense or in-frame; Supplementary Table 3, sheet
#     "genomic_driver_events") vs all other screened organoids, on that gene's own dependency
#   * PRMT5 / MTAP del: homozygous MTAP deletion (lowest copy number across MTAP below 0.5) vs all other
#     organoids with a clean copy-number call at MTAP, on PRMT5 dependency
#   * DDX3X: LOY vs Y-WT organoids (the LOY comparison cohort), on DDX3X dependency
# Statistics: one-sided Wilcoxon rank-sum test (biomarker-negative scores greater
# than biomarker-positive), at least 3 organoids per group for a p value, and
#   delta = median dependency (biomarker-positive) - median dependency (biomarker-negative)
# (negative = biomarker-positive organoids more dependent). 
# Two scopes are computed: all screened organoids (plotted) and a within-lineage check in each pair's
# dominant cancer type (tables only).

# ---- Settings for this panel ----
# alleles: NA = any activating mutation in the gene; otherwise a regular expression matched against the
# protein change (e.g. "^V600"), in which case organoids with OTHER activating mutations in the gene are
# left out of both groups.
bench_comparisons <- list(
  list(label = "KRAS",   target = "KRAS",   gene = "KRAS",   alleles = NA, check_type = "Colorectal"),
  list(label = "PIK3CA", target = "PIK3CA", gene = "PIK3CA", alleles = NA, check_type = "Colorectal"),
  list(label = "BRAF",   target = "BRAF",   gene = "BRAF",   alleles = NA, check_type = "Colorectal"))
include_mtap_prmt5   <- TRUE           # homozygous MTAP deletion -> PRMT5
mtap_check_type      <- "Oesophageal"  # within-lineage check for MTAP (most of its cases are oesophageal)
mtap_del_cn          <- 0.5            # homozygous deletion: lowest absolute copy number across MTAP below this
bench_mutation_types <- c("missense", "inframe")
bench_min_group      <- 3              # minimum organoids per group for a p value (Yi: 3)
bench_exclude_amp_wt <- FALSE          # TRUE = leave organoids with an amplification (but no mutation) out of WT
bench_p_cut     <- 0.01                # Yi's colour rule: p < 0.01, delta < -0.1, mean biomarker+ score < -0.3
bench_delta_cut <- -0.1
bench_mean_cut  <- -0.3
col_bench_loy <- "#c7615a"   # Yi's DDX3X / LOY red
col_bench_sig <- "#7e95bf"   # Yi's blue for pairs passing the rule
col_bench_ns  <- "grey60"
size_G        <- c(46, 50)   # mm (width, height)

# ---- Driver mutations (organoid rows only) ----
drivers <- read_sheet_table(molecular_xlsx, "genomic_driver_events", "sample_ID")
drivers <- drivers[model == "organoid"]
screened_ids <- rownames(fitness)
bench_info   <- meta[sample_ID %in% screened_ids, .(sample_ID, primary_tumour_type)]
cat(sprintf("Driver events for %d organoids (%d of %d screened organoids)\n",
            uniqueN(drivers$sample_ID), sum(screened_ids %in% drivers$sample_ID), length(screened_ids)))

target_id <- function(target) {   # Ensembl ID of a target gene in the dependency matrix, or NA
  id <- names(symbol)[symbol == target][1]
  if (is.na(id) || !(id %in% colnames(fitness))) return(NA_character_)
  id
}

# Per-organoid table for one comparison. assessable_ids = organoids that can be placed in either group;
# pos_ids = the biomarker-positive organoids among them.
make_comparison <- function(label, target, comparison, pos_ids, assessable_ids, check_type, alleles = NULL) {
  gene_id <- target_id(target)
  if (is.na(gene_id)) {
    cat("  ", target, "is not among the screened genes; skipped\n")
    return(NULL)
  }
  d <- copy(bench_info)[sample_ID %in% assessable_ids]
  d[, group := fifelse(sample_ID %in% pos_ids, "biomarker", "non_biomarker")]
  d[, dependency := fitness[sample_ID, gene_id]]
  if (is.null(alleles)) {
    d[, alleles := NA_character_]
  } else {
    d <- merge(d, alleles, by = "sample_ID", all.x = TRUE)
  }
  d[, `:=`(label = label, target = target, comparison = comparison, check_type = check_type)]
  cat(sprintf("  %-16s %3d biomarker+ vs %3d biomarker-\n", label, sum(d$group == "biomarker"),
              sum(d$group == "non_biomarker")))
  d
}

bench_rows <- list()

# ---- Mutation -> self-dependency ----
for (cmp in bench_comparisons) {
  activating_events <- drivers[gene == cmp$gene & Genomic_alterations %in% bench_mutation_types]
  if (is.na(cmp$alleles)) {
    matching_events <- activating_events
  } else {
    matching_events <- activating_events[grepl(cmp$alleles, variant)]
  }
  pos_ids       <- unique(matching_events$sample_ID)
  other_mut_ids <- setdiff(unique(activating_events$sample_ID), pos_ids)   # other alleles: left out
  assessable    <- setdiff(screened_ids, other_mut_ids)
  amplified_ids <- unique(drivers[gene == cmp$gene & Genomic_alterations == "AMP", sample_ID])
  n_amp_wt <- length(setdiff(intersect(assessable, amplified_ids), pos_ids))
  if (bench_exclude_amp_wt) assessable <- setdiff(assessable, setdiff(amplified_ids, pos_ids))
  if (n_amp_wt > 0) {
    cat(sprintf("  %-16s %d WT organoids carry a %s amplification (%s)\n", cmp$label, n_amp_wt, cmp$gene,
                if (bench_exclude_amp_wt) "excluded" else "kept in WT"))
  }
  alleles <- matching_events[, .(alleles = paste(sort(unique(variant)), collapse = ";")), by = sample_ID]
  bench_rows[[cmp$label]] <- make_comparison(cmp$label, cmp$target, paste0(cmp$label, " mutation -> ", cmp$target),
                                             pos_ids, assessable, cmp$check_type, alleles)
}

# ---- Homozygous MTAP deletion -> PRMT5 ----
# Lowest absolute copy number across the supported segments overlapping MTAP; organoids with a mixed
# baseline profile, or no supported segment at MTAP, have no call and are left out.
if (include_mtap_prmt5) {
  mtap_pos <- gene_pos[ensembl_id == "ENSG00000099810"]
  if (nrow(mtap_pos) == 1) {
    mtap_seg <- seg[chromosome == mtap_pos$chromosome & depthWindowCount > 0 &
                      start <= mtap_pos$end & end >= mtap_pos$start & !(sample_ID %in% mixed_ids)]
    mtap_cn  <- mtap_seg[, .(mtap_min_cn = min(copyNumber)), by = sample_ID]
    bench_rows[["PRMT5 / MTAP del"]] <- make_comparison(
      "PRMT5 / MTAP del", "PRMT5", "homozygous MTAP deletion -> PRMT5",
      mtap_cn[mtap_min_cn < mtap_del_cn, sample_ID], intersect(mtap_cn$sample_ID, screened_ids), mtap_check_type)
  } else {
    cat("   MTAP position not found in the Ensembl gene table; MTAP deletion comparison skipped\n")
  }
}

# ---- LOY -> DDX3X ----
bench_rows[["DDX3X"]] <- make_comparison("DDX3X", "DDX3X", "LOY -> DDX3X",
                                         y_cohort[altered == 1, sample_ID], y_cohort$sample_ID, "Colorectal")
bench_data <- rbindlist(bench_rows, use.names = TRUE)

cat("Alleles in each mutation biomarker-positive group:\n")
print(bench_data[group == "biomarker" & !is.na(alleles), .N, by = .(label, alleles)][order(label, -N)])

# ---- Test
bench_test <- function(dep, group) {
  pos <- dep[group == "biomarker"]
  neg <- dep[group == "non_biomarker"]
  result <- list(delta = NA_real_, median_pos = NA_real_, median_neg = NA_real_, mean_pos = NA_real_,
                 p = NA_real_, n_pos = length(pos), n_neg = length(neg))
  if (length(pos) == 0 || length(neg) == 0) return(result)
  result$median_pos <- median(pos)
  result$median_neg <- median(neg)
  result$delta      <- median(pos) - median(neg)
  result$mean_pos   <- mean(pos)
  if (length(pos) >= bench_min_group && length(neg) >= bench_min_group) {
    result$p <- suppressWarnings(wilcox.test(neg, pos, alternative = "greater")$p.value)
  }
  result
}

bench_pooled <- bench_data[, bench_test(dependency, group), by = .(comparison, label, target)]
bench_pooled[, scope := "All screened organoids"]
bench_within <- bench_data[primary_tumour_type == check_type, bench_test(dependency, group),
                           by = .(comparison, label, target, check_type)]
bench_within[, scope := paste("Within", tolower(check_type), "organoids")]
bench_within[, check_type := NULL]
bench_results <- rbind(bench_pooled, bench_within)
bench_results[, colour := fcase(target == "DDX3X", col_bench_loy,
                                !is.na(p) & p < bench_p_cut & delta < bench_delta_cut & mean_pos < bench_mean_cut, col_bench_sig,
                                default = col_bench_ns)]
setorder(bench_results, scope, p)
print(bench_results[, .(scope, comparison, n_pos, n_neg, delta, mean_pos, p)], digits = 3)

# ---- Figure: all screened organoids, Yi's Figure 6C style ----
bench_plot <- bench_results[scope == "All screened organoids" & !is.na(p)]
bench_x_breaks <- pretty(c(0, bench_plot$delta))          # ticks that bracket the data (x includes 0)
bench_y_breaks <- pretty(c(0, -log10(bench_plot$p)))

p_G <- ggplot(bench_plot, aes(delta, -log10(p))) +
  geom_point(aes(size = -log10(p), colour = colour), shape = 16, alpha = 0.7) +
  ggrepel::geom_text_repel(aes(label = label, colour = colour), size = txt(fig_base_pt + 1), family = fig_font,
                           box.padding = 0.5, min.segment.length = 0.3, seed = 1) +
  scale_colour_identity() +
  scale_size_continuous(range = c(1, 3), guide = "none") +
  scale_x_continuous(breaks = bench_x_breaks) +
  scale_y_continuous(breaks = bench_y_breaks, expand = expansion(0)) +
  coord_cartesian(xlim = range(bench_x_breaks), ylim = range(bench_y_breaks), clip = "off") +
  labs(x = "\u0394 Median Dependency", y = expression(-log[10]*" (p value)")) +
  guides(x = guide_axis(cap = "both"), y = guide_axis(cap = "both")) + pub_theme
save_figure(p_G, "panel_G_benchmark_clinical_biomarkers", size_G)

# D-G together
p_row3 <- ggpubr::ggarrange(p_D, p_E, p_F, p_G, nrow = 1, align = "h",
                            widths = c(size_D[1], size_E[1], size_F[1], size_G[1]))
save_figure(p_row3, "panels_DEFG_combined", c(176, size_D[2]))


# ---- Tables ----
fwrite(bench_results[, .(scope, comparison, label, target, n_pos, n_neg, median_pos, median_neg, delta, mean_pos, p)],
       file.path(table_dir, "benchmark_clinical_biomarkers.csv"))
fwrite(bench_data[, .(comparison, label, target, sample_ID, primary_tumour_type, group, alleles, dependency)],
       file.path(figure_dir, "panel_G_source_data.csv"))


## ==== 11. Supplementary tables ==================================================================
cat("\n==== 11. Supplementary tables ====\n")

# ---- Organoid calls: one row per organoid, one column per arm (status) plus Y status ----
calls_wide <- dcast(arm_status, sample_ID ~ event, value.var = "status")
organoid_info <- meta[, .(sample_ID, individual_ID, primary_tumour_type, clinical_sex, wgs_gender,
                          sex_annotation_conflict, purity, ploidy, msStatus, crispr_screened,
                          Y_status, Y_evidence_tier = evidence_tier)]
organoid_info <- merge(organoid_info, screens[, .(sample_ID, crispr_library = library)], by = "sample_ID", all.x = TRUE)
organoid_calls <- merge(organoid_info, calls_wide, by = "sample_ID", all.x = TRUE)
setcolorder(organoid_calls, c(names(organoid_info), all_events))

# ---- Arm call details: one row per organoid x arm, with the measurements behind each call ----
arm_details <- arm_status[, .(sample_ID, event = as.character(event), status, measured_state, in_population,
                              coverage = round(coverage, 4), relative_median = round(relative_median, 4),
                              fraction_low = round(fraction_low, 4), fraction_neutral = round(fraction_neutral, 4),
                              fraction_gain = round(fraction_gain, 4))]

# ---- Y status evidence: DNA and RNA measurements behind every LOY call ----
y_evidence_cols <- c("sample_ID", "individual_ID", "primary_tumour_type", "clinical_sex", "wgs_gender",
                     "sex_annotation_conflict", "crispr_screened", "Y_status", "evidence_tier",
                     "dna_measurement_state", "rna_measurement_state",
                     "dna_Y_coverage", "dna_Y_median", "dna_Y_fraction_lt05", "dna_Y_fraction_ge05",
                     "dna_loss_stable", "dna_QC_eligible",
                     "rna_Y_median_TPM", "rna_low_marker_count", "rna_loss_stable", "rna_QC_eligible",
                     paste0("rna_", names(y_markers), "_TPM"))
y_evidence <- meta[, ..y_evidence_cols]

# ---- Comparison cohorts: which organoids were compared, for every tested comparison ----
cohort_rows <- list()
for (this_key in names(cohorts)) {
  info <- comparisons[match(this_key, comparisons$key)]
  cohort_rows[[this_key]] <- cohorts[[this_key]][, .(type = info$type, event = info$event, sample_ID, individual_ID,
                                                     primary_tumour_type, library,
                                                     group = fifelse(altered == 1, "altered", "neutral"))]
}
cohort_table <- rbindlist(cohort_rows)

# ---- Vulnerability hits, with the median dependency score in each group ----
hits_table <- res[hit == TRUE, .(type, event, gene, ensembl_id, gene_chromosome = chromosome, gene_arm,
                                 gene_on_altered_arm = on_altered_arm, rel_dep, rel_dep_CI_low, rel_dep_CI_high,
                                 p, n_altered, n_neutral, n_donors)]
hits_table[, median_dependency_altered := NA_real_]
hits_table[, median_dependency_neutral := NA_real_]
for (i in seq_len(nrow(hits_table))) {
  frame  <- cohorts[[paste(hits_table$type[i], hits_table$event[i])]]
  scores <- fitness[frame$sample_ID, hits_table$ensembl_id[i]]
  set(hits_table, i, "median_dependency_altered", median(scores[frame$altered == 1]))
  set(hits_table, i, "median_dependency_neutral", median(scores[frame$altered == 0]))
}
setorder(hits_table, p)
cat(sprintf("Vulnerability hits: %d (%d arm-loss, %d arm-gain)\n", nrow(hits_table),
            sum(hits_table$type == "Arm loss"), sum(hits_table$type == "Arm gain")))

# ---- README sheet ----
readme <- data.table(
  item = c("Organoid_calls", "Arm_call_details", "Y_status_evidence", "Tested_comparisons", "Comparison_cohorts",
           "Vulnerability_hits", "Arm statuses", "Y column / Y_status", "rel_dep", "Hit rule", "Model", "Source"),
  description = c(
    "One row per organoid: annotations, Y status, and the copy-number status of every chromosome arm.",
    "One row per organoid x arm: the status and the measurements behind it (coverage of the masked arm; base-pair-weighted median and fractions of the arm below 0.75, within 0.85-1.15 and above 1.25, all relative to ploidy).",
    "DNA and RNA evidence behind every loss-of-Y call.",
    "Every arm loss / gain (and Y loss) considered, the organoid counts at each stage, and why untested comparisons were not tested.",
    "The organoids in each tested comparison (altered vs neutral).",
    "Genes more essential in organoids with the arm altered: p < 0.01 and rel_dep >= 0.2.",
    "loss / neutral / gain / partial_loss / intermediate_or_mixed (measured); QC_limited (DNA QC failed: source QC, purity < 0.8); unassessable (mixed baseline profile, acrocentric short arm, or < 90% coverage); out_of_population (X arms in male organoids, Y in female organoids).",
    "Y loss = high-confidence LOY (DNA loss and low expression of >= 6 of 7 Y genes, both assays passing QC, stable across 27 DNA and 9 RNA threshold settings); Y neutral = Y-retained (DNA retained, RNA not low) with both QCs passing.",
    "Relative dependency = neutral - altered, adjusted (regression coefficient x -1). Positive = more essential when the arm is altered. CI = 95% confidence interval.",
    sprintf("p < %s and rel_dep >= %s; no multiple-testing correction.", p_cut, effect_cut),
    "dependency ~ altered + cancer type + sex/Y group (not for Y) + library + MSI + ploidy + screen quality (AUC_PR); OLS with donor-clustered standard errors; t tests on donors - 1 degrees of freedom.",
    "Herranz-Ors et al., Nature 2026 (Supplementary Tables 2, 3, 4, 6 and the Figshare copy-number segments)."))

# ---- Write the workbook, and CSV copies of everything ----
supp_tables <- list(README = readme, Organoid_calls = organoid_calls, Arm_call_details = arm_details,
                    Y_status_evidence = y_evidence, Tested_comparisons = eligibility,
                    Comparison_cohorts = cohort_table, Vulnerability_hits = hits_table)
workbook_path <- file.path(supp_table_dir, "supplementary_table_aneuploidy_calls_and_vulnerabilities.xlsx")
writexl::write_xlsx(supp_tables, workbook_path)
cat("Saved", workbook_path, "\n")

fwrite(organoid_calls, file.path(table_dir, "organoid_calls.csv"))
fwrite(arm_details,    file.path(table_dir, "arm_call_details.csv.gz"))
fwrite(y_evidence,     file.path(table_dir, "y_status_evidence.csv"))
fwrite(eligibility,    file.path(table_dir, "tested_comparisons.csv"))
fwrite(cohort_table,   file.path(table_dir, "comparison_cohorts.csv"))
fwrite(hits_table,     file.path(table_dir, "vulnerability_hits.csv"))
fwrite(hits_by_comparison, file.path(table_dir, "hits_by_comparison.csv"))
fwrite(res[, .(type, event, gene, ensembl_id, chromosome, gene_arm, on_altered_arm, beta, SE, CI95_low, CI95_high,
               rel_dep, rel_dep_CI_low, rel_dep_CI_high, p, q_BH_reference_only, hit, n_altered, n_neutral, n_donors)],
       file.path(table_dir, "all_tests.csv.gz"))
fwrite(data.table(marker = names(y_markers), female_99th_percentile_TPM = female_q99,
                  threshold_TPM = pmax(y_rna_floor_tpm, female_q99)), file.path(table_dir, "y_marker_thresholds.csv"))

## ==== 13. Two-sheet supplementary table (copy of two sheets from the main workbook) ====
# Writes a NEW file; the seven-sheet workbook is left untouched. Uses the tables already in memory, or,
# in a fresh R session, reads them back from the existing workbook, so nothing is recomputed.
# (In a fresh session, first set: supp_table_dir <- "path/to/your/output/supplementary_tables")
main_workbook <- file.path(supp_table_dir, "supplementary_table_aneuploidy_calls_and_vulnerabilities.xlsx")
if (exists("organoid_calls") && exists("hits_table")) {
  two_sheets <- list(Aneuploidy_calls = organoid_calls, Dependency_analysis = hits_table)
} else {
  two_sheets <- list(Aneuploidy_calls    = readxl::read_xlsx(main_workbook, sheet = "Organoid_calls"),
                     Dependency_analysis = readxl::read_xlsx(main_workbook, sheet = "Vulnerability_hits"))
}
two_sheet_path <- file.path(supp_table_dir, "supplementary_table_organoid_aneuploidy_and_dependencies.xlsx")
writexl::write_xlsx(two_sheets, two_sheet_path)
cat("Saved", two_sheet_path, "\n")


## ==== 12. Provenance ============================================================================
reference_files <- c(file.path(reference_dir, "hg38_gap.txt.gz"), file.path(reference_dir, "hg38-blacklist.v2.bed.gz"),
                     file.path(reference_dir, "hg38_cytoBandIdeo.txt.gz"), gene_pos_file)
checked <- c(inputs, reference_files)
fwrite(data.table(file = checked, md5 = unname(tools::md5sum(checked)), bytes = file.size(checked)),
       file.path(output_dir, "input_checksums.csv"))
settings <- data.table(
  setting = c("y_dna_min_coverage", "y_dna_loss_median", "y_dna_loss_fraction", "y_dna_keep_median",
              "y_dna_keep_fraction", "y_rna_floor_tpm", "y_rna_min_low", "min_purity",
              "arm_min_coverage", "loss_ratio", "loss_fraction", "neutral_low", "neutral_high", "neutral_fraction",
              "gain_ratio", "gain_fraction", "partial_fraction",
              "min_per_lineage", "min_altered", "min_neutral", "min_donors", "min_residual_df",
              "p_cut", "effect_cut", "min_tpm_neutral", "include_Y_panelC", "focus_x_gene", "focus_y_gene"),
  value = c(y_dna_min_coverage, y_dna_loss_median, y_dna_loss_fraction, y_dna_keep_median,
            y_dna_keep_fraction, y_rna_floor_tpm, y_rna_min_low, min_purity,
            arm_min_coverage, loss_ratio, loss_fraction, neutral_range[1], neutral_range[2], neutral_fraction,
            gain_ratio, gain_fraction, partial_fraction,
            min_per_lineage, min_altered, min_neutral, min_donors, min_residual_df,
            p_cut, effect_cut, min_tpm_neutral, include_Y_panelC, focus_x_gene, focus_y_gene))
fwrite(settings, file.path(output_dir, "analysis_settings.csv"))
writeLines(c(paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "", capture.output(sessionInfo())),
           file.path(output_dir, "session_info.txt"))
cat("\nDone. Everything is in", output_dir, "\n")

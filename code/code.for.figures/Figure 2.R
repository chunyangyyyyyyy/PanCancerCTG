library(tidyverse)
library(data.table)
library(magrittr)
library(foreach)
library(patchwork)

##- Fig. 2a --------------------------------------------------------------------
CTG.colors <- c("X-linked CTGs" = "#D53E4F", "Autosomal CTGs"  = "#4575B4")
load("./Proj_CTAs/output/manuscript.revision/X.CTG.methylation.input.RData")
load("./Proj_CTAs/output/manuscript.revision/methy.expr.pcor.RData")
plot.df <- pcor.df |> left_join(common.genes, by = "ensembl") |> 
  dplyr::filter(CTG_group %in% c("X-linked CTGs (Identified)", "Autosomal CTGs (Identified)")) |>
  mutate(
    CTG_group = sub(" \\(.+\\)$", "", CTG_group),
    gender    = str_to_title(gender)
  )

plot.df1   <- plot.df |> group_by(cancer, gender, CTG_group) |> summarise(med.rho = median(rho, na.rm = TRUE), .groups = "drop")
p1         <- ggplot(plot.df1, aes(x = CTG_group, y = med.rho)) +
  geom_hline(yintercept = 0, color = "#B0892E", linetype = "dotted", size = 0.45) +
  geom_line(aes(group = cancer), color = "grey75", size = 0.35, alpha = 0.65) +
  geom_point(aes(fill = CTG_group), shape = 21, size = 2.6, color = "grey25", stroke = 0.25, alpha = 0.85) +

  ggpubr::stat_compare_means(
    comparisons = list(c("X-linked CTGs", "Autosomal CTGs")), paired = TRUE,
    label = "p.format", tip.length = 0.01, hide.ns = FALSE, size = 3.5
  ) +

  scale_fill_manual(values = CTG.colors, guide = "none") +
  scale_x_discrete(
    limits = c("Autosomal CTGs", "X-linked CTGs"),
    labels = c("Autosomal CTGs" = "Autosomal\nCTGs", "X-linked CTGs" = "X-linked\nCTGs"),
    expand = expansion(add = 0.45)
  ) +
  scale_y_continuous(
    breaks = scales::breaks_width(0.1),
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  facet_wrap(~ gender, ncol = 2, strip.position = "top") +
  labs(x = NULL, y = expression(paste("Median partial Spearman's ", rho))) +

  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.title.y = element_text(size = 10.5, margin = margin(r = 7)),
    axis.text.x = element_text(size = 10, color = "black", lineheight = 0.95, margin = margin(t = 5)),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.line = element_line(color = "black", size = 0.35),
    axis.ticks = element_line(color = "black", size = 0.30),

    panel.grid.major.y = element_line(color = "grey94", size = 0.25),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(1.2, "lines"),

    strip.background = element_rect(fill = "grey96", color = NA),
    strip.text = element_text(size = 10.5, face = "plain", color = "black", margin = margin(t = 6, b = 6)),

    legend.position = "none",
    plot.margin = margin(t = 8, r = 10, b = 6, l = 8)
  )
pdf("./Proj_CTAs/output/figure.panels/Fig. 2a.pdf", width = 7, height = 4)
print(p1)
dev.off()

##- Fig. 2b --------------------------------------------------------------------
tcga.meth <- fread("./Datas/UCSC_Xena/PANCAN/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz")
setnames(tcga.meth, old = 1, new = "probe")
tcga.meth <- tcga.meth[, !duplicated(colnames(tcga.meth)), with = FALSE]
global.meth <- tcga.meth[, -1] |> colMeans(na.rm = TRUE)
global.meth <- data.frame(sample = names(global.meth), global.meth = global.meth)

tcga.sample <- fread("./Datas/UCSC_Xena/PANCAN/TCGA_phenotype_denseDataOnlyDownload.tsv.gz")
setnames(tcga.sample, old = "_primary_disease", new = "primary_disease")
tcga.sample <- tcga.sample |> 
  dplyr::filter(!is.na(sample_type_id)) |>
  mutate(TCGA_code = case_match(
    primary_disease,
    "acute myeloid leukemia"                ~ "LAML",
    "adrenocortical cancer"                 ~ "ACC",
    "bladder urothelial carcinoma"          ~ "BLCA",
    "brain lower grade glioma"              ~ "LGG",
    "breast invasive carcinoma"             ~ "BRCA",
    "cervical & endocervical cancer"        ~ "CESC",
    "cholangiocarcinoma"                    ~ "CHOL",
    "colon adenocarcinoma"                  ~ "COAD",
    "diffuse large B-cell lymphoma"         ~ "DLBC",
    "esophageal carcinoma"                  ~ "ESCA",
    "glioblastoma multiforme"               ~ "GBM",
    "head & neck squamous cell carcinoma"   ~ "HNSC",
    "kidney chromophobe"                    ~ "KICH",
    "kidney clear cell carcinoma"           ~ "KIRC",
    "kidney papillary cell carcinoma"       ~ "KIRP",
    "liver hepatocellular carcinoma"        ~ "LIHC",
    "lung adenocarcinoma"                   ~ "LUAD",
    "lung squamous cell carcinoma"          ~ "LUSC",
    "mesothelioma"                          ~ "MESO",
    "ovarian serous cystadenocarcinoma"     ~ "OV",
    "pancreatic adenocarcinoma"             ~ "PAAD",
    "pheochromocytoma & paraganglioma"      ~ "PCPG",
    "prostate adenocarcinoma"               ~ "PRAD",
    "rectum adenocarcinoma"                 ~ "READ",
    "sarcoma"                               ~ "SARC",
    "skin cutaneous melanoma"               ~ "SKCM",
    "stomach adenocarcinoma"                ~ "STAD",
    "testicular germ cell tumor"            ~ "TGCT",
    "thymoma"                               ~ "THYM",
    "thyroid carcinoma"                     ~ "THCA",
    "uterine carcinosarcoma"                ~ "UCS",
    "uterine corpus endometrioid carcinoma" ~ "UCEC",
    "uveal melanoma"                        ~ "UVM"
  ))

plot.df          <- inner_join(global.meth, tcga.sample, by = "sample") |> dplyr::filter(TCGA_code %in% names(tumor.ls))
norm.global.meth <- plot.df |> dplyr::filter(sample_type_id == 11) |> group_by(TCGA_code) |> summarise(norm.global.meth = median(global.meth))
plot.df          <- plot.df |> left_join(norm.global.meth, by = "TCGA_code") |> 
  mutate(delta.global.meth = global.meth - norm.global.meth,
         sample.type       = ifelse(sample_type_id < 10, "Tumor", "Normal"))

cancer.order      <- plot.df |> dplyr::filter(sample.type == "Tumor") |> group_by(TCGA_code) |>summarise(x = median(global.meth)) |> arrange(x) |> pull(TCGA_code)
plot.df$TCGA_code <- factor(plot.df$TCGA_code, levels = cancer.order)
x.labels          <- setNames(levels(plot.df$TCGA_code), levels(plot.df$TCGA_code))
x.labels["LIHC"]  <- "<span style='color:#C44E52;'>LIHC</span>"

p <- ggplot(plot.df, aes(x = TCGA_code, y = global.meth, fill = sample.type)) +
  
  geom_boxplot(
    lwd = 0.3, 
    width = 0.7,
    position = position_dodge2(preserve = "single", padding = 0, width = 0.8),
    
    outlier.size = 0.1,      
    outlier.alpha = 0.2,     
    outlier.color = "grey30",
    outlier.shape = 19) +
  
  scale_fill_manual(values = c("Normal" = "#E8E8E8", "Tumor" = "#C44E52")) +
  scale_x_discrete(labels = x.labels, expand = expansion(add = 0.55)) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  
  labs(x = NULL, y = "Global methylation level", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.title.y = element_text(size = 11, margin = margin(r = 8)),
    axis.text.x = ggtext::element_markdown(
      angle = 45, hjust = 1, vjust = 1, size = 9.5, color = "black", margin = margin(t = 5)
    ),
    axis.text.y = element_text(size = 9.5, color = "black"),
    axis.line = element_line(color = "black", size = 0.35),
    axis.ticks = element_line(color = "black", size = 0.30),
    axis.ticks.length = unit(2, "mm"),

    panel.grid.major.y = element_line(color = "grey94", size = 0.25),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),

    legend.position = c(0.9, 0.9),
    legend.justification = "right",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    legend.key.width = unit(0.45, "cm"),
    legend.key.height = unit(0.35, "cm"),
    legend.background = element_blank(),
    legend.margin = margin(t = 0, r = 0, b = 3, l = 0),

    plot.margin = margin(t = 6, r = 10, b = 6, l = 8)
  )

pdf("./Proj_CTAs/output/figure.panels/TCGA.global.methylation.pdf", width = 6, height = 4)
print(p)
dev.off()

##- Fig. 2c and 2d -------------------------------------------------------------
# Fig. 2c
load("./Proj_CTAs/output/manuscript.revision/pan.sample.CTGs.RData")
CTGs.df.updated <- readRDS("./Proj_CTAs/input/CTGs/CTGs.df.updated.rds")  ## candidates
CTGs.df.updated <- CTGs.df.updated %>%
  mutate(
    CTG_group = case_when(
      grepl("^X", location) & symbol %in% pan.sample.CTGs                      ~ "X-linked CTGs (Identified)",
      grepl("^X", location) & is_CTdb == "yes" & !symbol %in% pan.sample.CTGs  ~ "X-linked CTGs (Unidentified)",
      
      grepl("^\\d", location) & symbol %in% pan.sample.CTGs                     ~ "Autosomal CTGs (Identified)",
      grepl("^\\d", location) & is_CTdb == "yes" & !symbol %in% pan.sample.CTGs ~ "Autosomal CTGs (Unidentified)"
    )
  ) |>
  dplyr::filter(!is.na(CTG_group))
library(GEOquery)
gse        <- getGEO(filename = "./Proj_CTAs/input/bulk RNA-seq datasets/GSE202560_DNMTi_HCC/GSE202560_series_matrix.txt.gz", getGPL = FALSE)
meta       <- pData(gse) |> dplyr::select(geo_accession, cell_line = source_name_ch1, treatment = `treatment:ch1`)
grch38.ann <- fread("./Proj_CTAs/input/bulk RNA-seq datasets/GSE202560_DNMTi_HCC/Human.GRCh38.p13.annot.tsv.gz")
tpm        <- fread("./Proj_CTAs/input/bulk RNA-seq datasets/GSE202560_DNMTi_HCC/GSE202560_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")

gene.map     <- CTGs.df.updated |> inner_join(grch38.ann, by = c("symbol" = "Symbol")) |> dplyr::select(symbol, CTG_group, GeneID)
compare.df   <- gene.map |> inner_join(tpm, by = "GeneID") |> pivot_longer(starts_with("GSM"), names_to = "geo_accession", values_to = "tpm") |> inner_join(meta, by = "geo_accession")
compare.df   <- compare.df |> group_by(CTG_group, symbol, cell_line, treatment) |> summarise(tpm = mean(tpm)) |>
  group_modify(~ {
    d  <- .x |> dplyr::filter(treatment == "DPBS") |> pull(tpm)
    .x |> dplyr::select(treatment, tpm) |> mutate(logfc = log2((tpm + 1) / (d + 1)))
  }) |>
  ungroup()

group.order <- c("Autosomal CTGs", "X-linked CTGs")
plot.df     <- compare.df |> 
  dplyr::filter(treatment == "5-Aza-CdR", CTG_group %in% c("X-linked CTGs (Identified)", "Autosomal CTGs (Identified)")) |>
  mutate(
    CTG_group = sub(" \\(.+\\)$", "", CTG_group) |> factor(levels = group.order),
    cell_line = sprintf("%s (Male)", cell_line)
  )
p <- ggplot(plot.df, aes(x = CTG_group, y = logfc)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.30) +
  geom_violin(
    aes(fill = CTG_group), trim = TRUE, scale = "width", width = 0.80,
    alpha = 0.85, color = "grey30", linewidth = 0.30
  ) +
  geom_boxplot(
    width = 0.14, fill = "white", color = "grey30", linewidth = 0.30,
    outlier.shape = NA, show.legend = FALSE
  ) +
  ggpubr::stat_compare_means(
    comparisons = list(c("X-linked CTGs", "Autosomal CTGs")), label = "p.signif",
    label.y = 4.15, step.increase = 0, tip.length = 0.01, bracket.size = 0.30, size = 3
  ) +
  scale_fill_manual(values = CTG.colors, breaks = group.order, name = NULL) +
  scale_x_discrete(expand = expansion(add = 0.55)) +
  scale_y_continuous(breaks = c(-2, 0, 2, 4)) +
  coord_cartesian(ylim = c(-2.5, 5), expand = FALSE) +
  facet_wrap(~ cell_line, nrow = 1) +
  labs(x = NULL, y = expression(log[2] ~ "fold change")) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1, colour = NA))) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 10.5, margin = margin(r = 7)),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.30),

    panel.grid = element_blank(),
    panel.spacing.x = grid::unit(4, "mm"),
    strip.background = element_blank(),
    strip.text = element_text(size = 9.5, face = "plain", margin = margin(t = 4, b = 4)),

    legend.position = "inside",
    legend.position.inside = c(0.98, 0.03),
    legend.justification.inside = c(1, 0),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.30),
    legend.margin = margin(t = 3, r = 4, b = 3, l = 4),

    plot.margin = margin(t = 5, r = 6, b = 5, l = 5)
  )

p <- p +
  guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1, colour = NA))) +
  theme(
    legend.direction = "horizontal",
    legend.text = element_text(size = 7.5),
    legend.key.width = grid::unit(2.5, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.key.spacing.x = grid::unit(1.5, "mm"),
    legend.margin = margin(t = 1, r = 2, b = 1, l = 2)
  )

p -> p1

# Fig. 2d
library(DESeq2)
load("./Proj_CTAs/output/manuscript.revision/5-Aza-CdR.dosage.treatment.RData")

## batch 1
coldata    <- sample.info |> dplyr::filter(batch == "Batch1", cancer.type == "panc0213") |> dplyr::select(unique.sample, group) |> column_to_rownames("unique.sample")
countdata  <- counts[, rownames(coldata)]
countdata  <- countdata[rowSums(countdata ) > 0, ]
dds        <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ group)
dds        <- DESeq(dds)
res1       <- results(dds, contrast = c("group", "panc0213_Batch1_T0.1", "panc0213_Batch1_C")) |> as.data.frame() |> rownames_to_column("ensembl")

## batch 2
coldata    <- sample.info |> dplyr::filter(batch == "Batch2", cancer.type == "panc0213") |> dplyr::select(unique.sample, group) |> column_to_rownames("unique.sample")
countdata  <- counts[, rownames(coldata)]
countdata  <- countdata[rowSums(countdata ) > 0, ]
dds        <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ group)
dds        <- DESeq(dds)
res2       <- results(dds, contrast = c("group", "panc0213_Batch2_T1", "panc0213_Batch2_C")) |> as.data.frame() |> rownames_to_column("ensembl")
res3       <- results(dds, contrast = c("group", "panc0213_Batch2_T5", "panc0213_Batch2_C")) |> as.data.frame() |> rownames_to_column("ensembl")

n.sig.CTGs <- bind_rows(
  res1 |> inner_join(CTGs.df.updated, by = "ensembl") |> 
    group_by(CTG_group) |> summarise(
      n    = n(),
      up   = sum(log2FoldChange > 0.5 & padj < 0.05, na.rm = TRUE),
      down = sum(log2FoldChange < -0.5 & padj < 0.05, na.rm = TRUE)
    ) |>
    mutate(dose = "0.1 uM"),

  res2 |> inner_join(CTGs.df.updated, by = "ensembl") |> 
    group_by(CTG_group) |> summarise(
      n    = n(),
      up   = sum(log2FoldChange > 0.5 & padj < 0.05, na.rm = TRUE),
      down = sum(log2FoldChange < -0.5 & padj < 0.05, na.rm = TRUE)
    ) |>
    mutate(dose = "1 uM"),

  res3 |> inner_join(CTGs.df.updated, by = "ensembl") |> 
    group_by(CTG_group) |> summarise(
      n    = n(),
      up   = sum(log2FoldChange > 0.5 & padj < 0.05, na.rm = TRUE),
      down = sum(log2FoldChange < -0.5 & padj < 0.05, na.rm = TRUE)
    ) |>
    mutate(dose = "5 uM")
)

CTG.colors <- c("Autosomal CTGs" = "#4575B4", "X-linked CTGs"  = "#D53E4F")
dose.order <- c("0.1 uM", "1 uM", "5 uM")

plot.df    <- n.sig.CTGs |> dplyr::filter(CTG_group %in% c("X-linked CTGs (Identified)", "Autosomal CTGs (Identified)")) |>
  mutate(CTG_group = sub(" \\(.+\\)$", "", CTG_group)) |>
  group_by(dose) |>
  group_modify(~ {
    .x  <- .x |> mutate(other = n - up)
    a   <- .x |> dplyr::filter(CTG_group == "Autosomal CTGs") |> pull(up)
    b   <- .x |> dplyr::filter(CTG_group == "Autosomal CTGs") |> pull(other)
    c   <- .x |> dplyr::filter(CTG_group == "X-linked CTGs") |> pull(up)
    d   <- .x |> dplyr::filter(CTG_group == "X-linked CTGs") |> pull(other)
    
    mat <- matrix(c(c, d, a, b), nrow = 2, byrow = TRUE)
    fit <- fisher.test(mat, alternative = "two.sided")

    data.frame(
      prop.up.autosomal = a / (a + b),
      prop.up.xlinked   = c / (c + d),
      autosomal.label   = sprintf("%d/%d", a, a + b),
      xlinked.label     = sprintf("%d/%d", c, c + d),
      odds.ratio        = if (any(colSums(mat) == 0)) NA else fit$estimate, 
      p.val             = fit$p.value
    )

  }) |> 
  ungroup() |>
  mutate(
    dose      = factor(dose, levels = dose.order),
    dose.x    = as.integer(dose)
  )

bar.df     <- bind_rows(
  plot.df |> transmute(dose, CTG_group = "Autosomal CTGs", prop.up = prop.up.autosomal, count.label = autosomal.label, dose.x),
  plot.df |> transmute(dose, CTG_group = "X-linked CTGs", prop.up = prop.up.xlinked, count.label = xlinked.label, dose.x)
) |>
  mutate(
    dose = factor(dose, levels = dose.order),
    dose.x = as.integer(dose),
    CTG_group = factor(CTG_group, levels = c("Autosomal CTGs", "X-linked CTGs"))
  )
annot.df   <- plot.df |>
  mutate(
    dose = factor(dose, levels = dose.order),
    dose.x = as.integer(dose),
    bracket.y = pmax(prop.up.autosomal, prop.up.xlinked) + 0.085,

    p.label = case_when(
      !is.finite(p.val) ~ "italic(P) ~ 'not available'",
      p.val == 0       ~ "italic(P) ~ 'below numerical precision'",
      p.val < 0.001    ~ paste0(
        "italic(P) == ",
        sub("e", " %*% 10^", formatC(p.val, format = "e", digits = 2), fixed = TRUE)
      ),
      TRUE ~ paste0("italic(P) == ", formatC(p.val, format = "fg", digits = 3))
    ),

    or.label = case_when(
      is.na(odds.ratio)       ~ "'OR not estimable'",
      is.infinite(odds.ratio) ~ "'OR' == infinity",
      TRUE                   ~ paste0("'OR' == ", sprintf("%.1f", odds.ratio))
    )
  )

p <- ggplot(bar.df, aes(x = dose.x, y = prop.up, fill = CTG_group)) +
  geom_col(
    position = position_dodge(width = 0.72), width = 0.64,
    color = "grey30", size = 0.25, alpha = 0.85
  ) +
  geom_text(
    aes(label = count.label), position = position_dodge(width = 0.72),
    vjust = -0.4, size = 3.2, color = "black", show.legend = FALSE
  ) +

  geom_segment(
    data = annot.df,
    aes(x = dose.x - 0.18, xend = dose.x + 0.18, y = bracket.y, yend = bracket.y),
    inherit.aes = FALSE, color = "grey25", size = 0.3
  ) +
  geom_segment(
    data = annot.df,
    aes(x = dose.x - 0.18, xend = dose.x - 0.18, y = bracket.y - 0.015, yend = bracket.y),
    inherit.aes = FALSE, color = "grey25", size = 0.3
  ) +
  geom_segment(
    data = annot.df,
    aes(x = dose.x + 0.18, xend = dose.x + 0.18, y = bracket.y - 0.015, yend = bracket.y),
    inherit.aes = FALSE, color = "grey25", size = 0.3
  ) +

  geom_text(
    data = annot.df, aes(x = dose.x, y = bracket.y + 0.028, label = or.label),
    inherit.aes = FALSE, parse = TRUE, vjust = 0, size = 3.1, color = "black"
  ) +
  geom_text(
    data = annot.df, aes(x = dose.x, y = bracket.y + 0.092, label = p.label),
    inherit.aes = FALSE, parse = TRUE, vjust = 0, size = 3.2, color = "black"
  ) +

  scale_fill_manual(values = CTG.colors, name = NULL, breaks = c("Autosomal CTGs", "X-linked CTGs")) +
  scale_x_continuous(breaks = 1:3, labels = c("0.1", "1", "5"), expand = expansion(add = 0.5)) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25), labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, max(1.08, annot.df$bracket.y + 0.16, na.rm = TRUE))) +
  labs(
    title = "PANC 02.13 (Female)",
    x = expression(paste("Decitabine concentration (", mu, "M)")),
    y = "Upregulated CTGs (%)"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(alpha = 0.85))) +

  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "sans", color = "black"),
    plot.title = element_text(size = 12, face = "plain", hjust = 0.5, margin = margin(b = 8)),
    axis.title.x = element_text(size = 10.5, margin = margin(t = 7)),
    axis.title.y = element_text(size = 10.5, margin = margin(r = 7)),
    axis.text = element_text(size = 9.5, color = "black"),
    axis.line = element_line(color = "black", size = 0.35),
    axis.ticks = element_line(color = "black", size = 0.30),

    panel.grid = element_blank(),

    legend.position = "none",
    legend.text = element_text(size = 9.5),
    legend.key.width = unit(0.45, "cm"),
    legend.key.height = unit(0.35, "cm"),
    legend.background = element_blank(),
    legend.margin = margin(t = 4, b = 0),

    plot.margin = margin(t = 8, r = 10, b = 6, l = 8)
  )

p -> p2
# Fig. 2c and Fig. 2d
p1 <- p1 +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box.spacing = grid::unit(1, "mm")
  )

bar.df <- bar.df |> dplyr::mutate(cell.line = "PANC 02.13 (Female)")
p2     <- (p2 %+% bar.df) +
  labs(title = NULL) +
  facet_wrap(~ cell.line, nrow = 1)

p <- patchwork::wrap_plots(p1, p2, nrow = 1, widths = c(1, 1), guides = "keep") &
  theme(
    strip.background = element_blank(),
    strip.text.x.top = element_text(
      family = "sans", size = 9.5, face = "plain", color = "black",
      hjust = 0.5, margin = margin(t = 4, b = 4)
    ),
    plot.margin = margin(t = 5, r = 6, b = 5, l = 5)
  )

pdf("./Proj_CTAs/output/figure.panels/5-Aza-CdR.combined.pdf", width = 9, height = 4)
print(p)
dev.off()

##- Fig. 2e --------------------------------------------------------------------
load("./Proj_CTAs/output/3. define CTGs/pan.sample.CTGs.RData")
CTGs.df <- readRDS("./Proj_CTAs/input/CTGs/CTGs.df.rds")
CTGs.df <- CTGs.df |> 
  mutate(is_CTG        = ifelse(symbol %in% pan.sample.CTGs, "yes", "no"),
         chr_group     = case_when(
           grepl("^X", location) ~ "X", 
           grepl("^Y", location) ~ "Y", 
           .default = "Autosomes"
         ))

stat.df    <- CTGs.df |> dplyr::filter(chr_group %in% c("Autosomes", "X"))
fisher.mat <- table(stat.df$chr_group, stat.df$is_CTG)
p.val      <- fisher.test(fisher.mat)$p.value

p.sci      <- formatC(p.val, format = "e", digits = 2)
p.base     <- as.numeric(sub("e.*", "", p.sci))
p.exp      <- as.integer(sub(".*e", "", p.sci))
p.label    <- bquote(italic(P) == .(p.base) ~ "×" ~ 10^.(p.exp))

df.initial <- CTGs.df |> mutate(Stage = "Initial Candidates")
df.final   <- CTGs.df |> dplyr::filter(is_CTG == "yes") |> mutate(Stage = "Final CTGs")
plot.df    <- rbind(df.initial, df.final) |>
  group_by(Stage, chr_group) |>
  summarise(count = n()) |>
  group_by(Stage) |>
  mutate(
    total = sum(count),
    freq = count / total,
    stage_label = paste0(Stage, "\n(n=",scales::comma(total), ")")
  )

plot.df$Stage     <- factor(plot.df$Stage, levels = c("Initial Candidates", "Final CTGs"))
plot.df$chr_group <- factor(plot.df$chr_group, levels = c("X", "Autosomes", "Y"))

x.labels  <- plot.df |> distinct(Stage, stage_label) |> deframe()
my.colors <- c("Autosomes" = "#E8E8E8", "X" = "#C44E52", "Y" = "#7F7F7F")
p         <- ggplot(plot.df, aes(x = Stage, y = freq, fill = chr_group)) +
  
  geom_col(position = "fill", width = 0.6, color = "black", linewidth = 0.3) +
  
  scale_fill_manual(values = my.colors, name = NULL, breaks = c("Autosomes", "X", "Y")) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.025)),
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_x_discrete(labels = x.labels) +
  
  ## add p value
  annotate("segment", x = 1, xend = 2, y = 1.05, yend = 1.05, linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = 1.02, yend = 1.05, linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = 1.02, yend = 1.05, linewidth = 0.5) +
  annotate("text", x = 1.5, y = 1.10, label = p.label, size = 4.5, parse = FALSE) +
  
  geom_text(
    aes(label = ifelse(chr_group == "X", sprintf("%.1f%%", freq*100), "")),
    position = position_stack(vjust = 0.5), 
    color = "white", size = 3.5
  ) +
  
  labs( x = NULL, y = "Chromosomal Distribution") +
  
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    
    legend.position = "top",
    legend.justification = "center",
    legend.margin = margin(b = -5),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 10.5),
    
    aspect.ratio = 1.2
  )

pdf("./Fig. 2D.pdf", width = 5, height = 5)
p
dev.off()

##- Fig. 2f --------------------------------------------------------------------
load("./Proj_CTAs/output/manuscript.revision/X.CTG.methylation.input.RData")
load("./Proj_CTAs/output/manuscript.revision/methy.expr.specificity.RData")
CTG.colors <- c("X-linked CTGs" = "#D53E4F", "Autosomal CTGs"  = "#4575B4")

plot.df <- specificity.df |> left_join(common.genes, by = "ensembl") |> 
  dplyr::filter(CTG_group %in% c("X-linked CTGs (Identified)", "Autosomal CTGs (Identified)")) |>
  mutate(
    CTG_group = sub(" \\(.+\\)$", "", CTG_group),
    gender    = str_to_title(gender)
  ) |>

  group_by(CTG_group, ensembl, gender) |> 
  summarise(
    expr.gini = median(expr.gini, na.rm = TRUE),
    expr.skew = median(expr.skew, na.rm = TRUE),
    meth.gini = median(meth.gini, na.rm = TRUE),
    meth.skew = median(meth.skew, na.rm = TRUE), 
    median.meth = median(median.meth, na.rm = TRUE)
  ) |>
  ungroup()

p.main1 <- ggplot(plot.df |> dplyr::filter(gender == "Female"), aes(y = meth.skew, x = median.meth)) +
  
  geom_smooth(method = "lm", color = "grey35", se = TRUE, fill = "grey88", size = 0.55, alpha = 0.45) +
  
  geom_point(aes(fill = CTG_group), shape = 21, size = 2.3, color = "white", stroke = 0.25, alpha = 0.85) +
  
  ggpubr::stat_cor(
    aes(group = 1, label = paste(after_stat(r.label), gsub("italic\\(p\\)", "italic(P)", after_stat(p.label)), sep = "*`,`~")),
    method = "spearman", cor.coef.name = "rho", label.x.npc = "left", label.y.npc = "top",
    r.accuracy = 0.01, size = 3.5, hjust = 0, vjust = 0.5
  ) +
  
  scale_fill_manual(values = CTG.colors, name = NULL) +
  
  labs(x = "Median methylation level (MML)", y = "Hypomethylation specificity index (HSI)") +
  
  guides(
    fill = guide_legend(override.aes = list(size = 4, alpha = 1, stroke = 0.3))
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    aspect.ratio = 1,
    
    axis.line = element_line(size = 0.6),
    axis.text = element_text(color = "black", size = 11),
    
    legend.position = "none", 
    legend.background = element_blank(),
    legend.text = element_text(size = 10)
  )
p.top1   <- ggplot(plot.df |> dplyr::filter(gender == "Female"), aes(x = median.meth, fill = CTG_group)) +
  geom_density(alpha = 0.85, color = "white", size = 0.5) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = CTG.colors) +
  
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(t = 0, r = 2, b = 0, l = 0))
p.right1 <- p.top1 %+% aes(x = meth.skew) + coord_flip() + theme(plot.margin = margin(t = 2, r = 4, b = 0, l = 1))

p.main2  <- p.main1 %+% dplyr::filter(plot.df, gender == "Male")
p.top2   <- p.top1 %+% dplyr::filter(plot.df, gender == "Male")
p.right2 <- p.top2 %+% aes(x = meth.skew) + coord_flip()

p.main1 <- p.main1 +
  annotate(
    "label", x = Inf, y = -Inf, label = "Female", hjust = 1.15, vjust = -0.65,
    family = "sans", fontface = "plain", size = 4.5, color = "black", fill = "#F0F2F4", label.size = 0.35,
    label.padding = grid::unit(0.35, "lines"), label.r = grid::unit(0, "lines")
  )
p.main2 <- p.main2 +
  annotate(
    "label", x = Inf, y = -Inf, label = "Male", hjust = 1.15, vjust = -0.65,
    family = "sans", fontface = "plain", size = 4.5, color = "black", fill = "#F0F2F4", label.size = 0.35,
    label.padding = grid::unit(0.35, "lines"), label.r = grid::unit(0, "lines")
  )

library(aplot)
pdf("./Proj_CTAs/output/figure.panels/Fig. 2e.pdf", width = 7, height = 7)
p.main1 |> 
  insert_top(p.top1, height = 0.18) |>
  insert_right(p.right1, width = 0.22)

p.main2 |> 
  insert_top(p.top2, height = 0.18) |>
  insert_right(p.right2, width = 0.22)

dev.off()

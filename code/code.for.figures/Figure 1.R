library(tidyverse)
library(data.table)
library(readxl)
library(GSVA)
library(BiocParallel)
library(foreach)
library(magrittr)

##- Fig. 1a --------------------------------------------------------------------

##- Fig. 1b --------------------------------------------------------------------
CTGs.df.updated <- readRDS("./Proj_CTAs/input/CTGs/CTGs.df.updated.rds")  ## candidates
df              <- CTGs.df.updated |> mutate(value = ensembl, CTdatabase = is_CTdb == "yes", HPA = is_HPA == "yes")

p <- ggplot(df) +
  ggvenn::geom_venn(
    aes(A = CTdatabase, B = HPA), show_percentage = FALSE, auto_scale = TRUE,
    fill_color = c("#5E81B5", "#E19C25"), fill_alpha = 0.50,
    stroke_color = "grey80", stroke_size = 0.1,
    set_name_color = "black", set_name_size = 3.8,
    text_color = "black", text_size = 3.5
  ) +
  coord_fixed(clip = "off") +
  theme_void(base_size = 11, base_family = "sans") +
  theme(
    text = element_text(family = "sans", color = "black"),
    legend.position = "none",
    plot.margin = margin(t = 8, r = 10, b = 8, l = 10)
  )

pdf("./Proj_CTAs/output/1. data wrangling/Figures/venn plot.pdf", width = 4, height = 4)
print(p)
dev.off()

##- Fig. 1c --------------------------------------------------------------------

##- Fig. 1d, 1e and 1f ---------------------------------------------------------
source("./practice/customized.fuctions/add_sample_size.R")
hgnc.protein <- fread("./Datas/HGNC/gene_with_protein_product.txt")
load("./Proj_CTAs/output/manuscript.revision/CTDb.genes.RData")
CTdb.ensembl <- CTDb.hgnc.map |> filter(!is.na(ensembl_gene_id), ensembl_gene_id != "") |> pull(ensembl_gene_id) |> unique() ## 243

CTGs.Wang    <- read_xlsx("./Proj_CTAs/input/CTGs/Others/CTGs_Wang_NC_2016.xlsx", 2, skip = 1)
CTGs.Carter  <- read.csv("./Proj_CTAs/input/CTGs/Others/CTGs_Carter_JITC_2023.csv", header = TRUE)
CTGs.Carter  <- CTGs.Carter |> filter(!thymus_expression, tumor_expression)

load("./Proj_CTAs/output/manuscript.revision/pan.sample.CTGs.RData")
CTGs.df.updated    <- readRDS("./Proj_CTAs/input/CTGs/CTGs.df.updated.rds")  ## candidates
candidates.ensembl <- CTGs.df.updated$ensembl
CTGs.us.ensembl    <- CTGs.df.updated |> filter(symbol %in% pan.sample.CTGs) |> pull(ensembl)

CTGset       <- list(
  CTDatabase            = CTdb.ensembl,
  CTGs_Wang_NC_2016     = CTGs.Wang$ENname |> unique(),
  CTGs_Carter_JITC_2023 = CTGs.Carter$ensembl_gene_id |> unique(),
  candidates            = candidates.ensembl,
  CTGs_us               = CTGs.us.ensembl,
  CTDatabase.included   = intersect(CTdb.ensembl, CTGs.us.ensembl),
  CTDatabase.excluded   = setdiff(CTdb.ensembl, CTGs.us.ensembl)
)

CTGset.unique <- CTGset[c("CTDatabase", "CTGs_Wang_NC_2016", "CTGs_Carter_JITC_2023", "CTGs_us")]
unique.genes  <- unlist(CTGset.unique) |> table() %>% .[. == 1] |> names()
CTGset.unique <- map(CTGset.unique, \(x) intersect(x, unique.genes))

## TCGA cancer colors
cancer.info <- TCGAbiolinks::getGDCprojects() |> filter(grepl("^TCGA", id)) |> select(id, cancer = name, abbr = tumor)
cancers     <- cancer.info$abbr
colors      <- setNames(colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "Set3"))(length(cancers)), cancers)
# Fig. 1d
source("./Proj_CTAs/code/visualization.R")
load("./Proj_CTAs/input/bulk RNA-seq datasets/tcga.target.gtex.RData")  ## c("sample.ls", "tumor.sample.ls", "gene.info", "log2.tpm")
disc.sample.df <- foreach(i = names(sample.ls), .combine = 'rbind') %do% {sample.ls[[i]] |> mutate(cancer = i)}

load("./Datas/UCSC_Xena/normal.tumor.RData")
protein.genes <- gene.info |> filter(ensembl %in% hgnc.protein$ensembl_gene_id, id %in% rownames(log2.tpm)) |> pull(id)
ssgsea        <- gsva(expr          = log2.tpm[protein.genes, ],
                      gset.idx.list = lapply(CTGset, \(x) gene.info |> filter(ensembl %in% x) |> pull(id)),
                      method        = "ssgsea",
                      abs.ranking   = FALSE,
                      ssgsea.norm   = TRUE,
                      BPPARAM       = MulticoreParam(40))

## all TCGA cancer types
auc.df <- foreach(c = names(sample.ls), .combine = 'rbind') %do% {
  sample.info  <- sample.ls[[c]]
  y            <- factor(sample.info$group, levels = c("normal", "tumor"))
  ssgsea.score <- t(ssgsea)[sample.info$sample, , drop = FALSE]
  auc.df       <- foreach(i = colnames(ssgsea.score), .combine = 'rbind') %do% {
    auc   <- pROC::roc(y, ssgsea.score[, i], direction = "<")$auc |> as.numeric()
    p.val <- wilcox.test(ssgsea.score[, i] ~ y)$p.value

    data.frame(abbr = c, CTGset = i, auc, p.val)
  }
}
levels <- c("CTdatabase", "Wang et al. (2016)", "Carter et al. (2023)", "This Study", "Candidates", "CTDatabase (included)", "CTDatabase (excluded)")
auc.df <- auc.df |>
  mutate(CTGset = case_match(
    CTGset,
    "CTDatabase"            ~ "CTdatabase",
    "CTGs_Wang_NC_2016"     ~ "Wang et al. (2016)",
    "CTGs_Carter_JITC_2023" ~ "Carter et al. (2023)",
    "candidates"            ~ "Candidates",
    "CTGs_us"               ~ "This Study",
    "CTDatabase.included"   ~ "CTDatabase (included)",
    "CTDatabase.excluded"   ~ "CTDatabase (excluded)"
  )) |>
  mutate(CTGset = factor(CTGset, levels = levels))
ssgsea.df <- data.frame(sample = colnames(ssgsea), ssgsea = ssgsea["CTGs_us", ], row.names = NULL)
df        <- inner_join(disc.sample.df, ssgsea.df, by = "sample")
df$group  <- factor(df$group, levels = c("normal", "tumour"))
auc.disc  <- df |>
  group_by(cancer) |>
  group_modify(~ {
    auc <- pROC::roc(.x$group, .x$ssgsea, direction = "<")$auc |> as.numeric()
    data.frame(auc)
  }) |>
  left_join(cancer.abbr) |>
  mutate(abbr = factor(abbr, levels = cancer.abbr$abbr)) |>
  ungroup() |>
  select(abbr, auc) |>
  mutate(cohort = "Discovery Cohort")
drop.cancers <- c("GBM", "BRCA", "COAD", "READ", "HNSC", "KIRC", "LIHC", "LUAD", "OV", "PAAD", "PRAD", "SKCM")
auc.other    <- auc.df |> filter(CTGset == "This Study") |> select(abbr, auc) |> filter(abbr %in% drop.cancers == FALSE) |> mutate(cohort = "Other Cancer Types")
plot.df <- rbind(auc.disc |> filter(abbr != "NB"), auc.other)
p       <- ggplot(plot.df, aes(x = cohort, y = auc)) +
  gghalves::geom_half_boxplot(aes(fill = cohort), side = "l", position = position_nudge(x = -0.1), outlier.alpha = 0,
                              fill = "gray80", alpha = 0.7, lwd = 0.2) +
  ggbeeswarm::geom_beeswarm(aes(fill = abbr), side = 1,
                            cex = 3.5, shape = 21, size = 2.5, stroke = 0.3, alpha = 0.9) +
  
  scale_fill_manual(values = c(cancer.colours, colors[match(as.character(plot.df$abbr[plot.df$cohort == "Other Cancer Types"]), names(colors))])) +
  
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "black", linewidth = 0.15) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1)) +
  labs(x = NULL, y = "AUROC", title = NULL, subtitle = "TCGA and GTEx samples") +
  
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(hjust = 1, angle = 30),
        plot.subtitle = element_text(size = 10, face = "plain", hjust = 0, margin = margin(b = 8)),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.25))

fig.1d <- add_sample_size(p, equal_n_pos = "x_title")

pdf("./Proj_CTAs/output/figure.panels/Fig. 1d.pdf", width = 3.5, height = 5)
fig.1d
dev.off()

# Fig. 1e
tcga.purity     <- readRDS("./Datas/TCGA/pan-cancer/tcga.purity.rds")

gdc.log2.fpkm   <- fread("./Datas/UCSC_Xena/GDC-PANCAN/GDC-PANCAN.htseq_fpkm.tsv") |> column_to_rownames("Ensembl_ID") |> as.matrix()  ## log2(x+1)
gdc.sample.info <- fread("./Datas/UCSC_Xena/GDC-PANCAN/GDC-PANCAN.basic_phenotype.tsv")
gene.info.v22   <- fread("./Datas/UCSC_Xena/GDC-PANCAN/gencode.v22.annotation.gene.probeMap") |> mutate(ensembl_gene_id = str_sub(id, end = 15))

purity      <- tcga.purity |>
  select(sample = `Sample ID`, purity = CPE) |>
  mutate(purity = as.numeric(purity)) |>
  filter(!is.nan(purity))
sample.info <- gdc.sample.info |>
  filter(program == "TCGA", sample %in% intersect(purity$sample, colnames(gdc.log2.fpkm))) |>
  mutate(abbr = str_extract(project_id, "(?<=TCGA-).*"))

## all TCGA cancer types
protein.genes <- gene.info.v22 |> filter(ensembl_gene_id %in% hgnc.protein$ensembl_gene_id, id %in% rownames(gdc.log2.fpkm)) |> pull(id)
ssgsea        <- gsva(expr          = gdc.log2.fpkm[protein.genes, sample.info$sample],
                      gset.idx.list = lapply(CTGset, \(i) gene.info.v22 |> filter(ensembl_gene_id %in% i) |> pull(id)),
                      method        = "ssgsea",
                      abs.ranking   = FALSE,
                      ssgsea.norm   = TRUE,
                      BPPARAM       = MulticoreParam(40))

df     <- t(ssgsea) |> as.data.frame() |> rownames_to_column("sample") |> inner_join(sample.info) |> inner_join(purity)
rho.df <- foreach(i = rownames(ssgsea), .combine = 'rbind') %do% {
  fit.ls <- split(df, df$abbr) |> lapply(\(x) {cor.test(x[[i]], x$purity, method = "spearman")})
  rho    <- sapply(fit.ls, `[[`, "estimate")
  p.val  <- sapply(fit.ls, `[[`, "p.value")
  data.frame(rho, p.val, row.names = NULL) |> mutate(abbr = names(fit.ls), CTGset = i, .before = everything())
}
levels <- c("CTdatabase", "Wang et al. (2016)", "Carter et al. (2023)", "This Study", "Candidates", "CTDatabase (included)", "CTDatabase (excluded)")
rho.df <- rho.df |>
  mutate(CTGset = case_match(
    CTGset,
    "CTDatabase"            ~ "CTdatabase",
    "CTGs_Wang_NC_2016"     ~ "Wang et al. (2016)",
    "CTGs_Carter_JITC_2023" ~ "Carter et al. (2023)",
    "candidates"            ~ "Candidates",
    "CTGs_us"               ~ "This Study",
    "CTDatabase.included"   ~ "CTDatabase (included)",
    "CTDatabase.excluded"   ~ "CTDatabase (excluded)")) |>
  mutate(CTGset = factor(CTGset, levels = levels))
plot.df <- rho.df |> filter(!CTGset %in% c("Candidates", "CTDatabase (included)", "CTDatabase (excluded)"))
ls      <- lapply(c("Carter et al. (2023)", "Wang et al. (2016)", "CTdatabase"), \(x) c(x, "This Study"))
p       <- ggplot(plot.df, aes(x = reorder(CTGset, rho, median), y = rho)) +
  
  geom_hline(yintercept = seq(0, 0.8, by = 0.2), color = "grey95", linewidth = 0.4) +
  
  geom_line(aes(group = abbr), color = "grey80", linewidth = 0.4, alpha = 0.6) +
  
  stat_summary(
    fun.data = function(x) data.frame(y = median(x), ymin = quantile(x, 0.25), ymax = quantile(x, 0.75)),
    geom = "errorbar", width = 0.25, linewidth = 0.8, color = "#D53E4F"
  ) +
  stat_summary(
    fun.data = function(x) data.frame(y = median(x), ymin = median(x), ymax = median(x)),
    geom = "errorbar", width = 0.4, linewidth = 1.5, color = "#D53E4F" 
  ) +
  
  geom_point(aes(fill = abbr), shape = 21, size = 2.8, stroke = 0.2, alpha = 0.9, color = "black") +
  
  ggpubr::stat_compare_means(comparisons = ls, paired = TRUE, tip.length = 0.01, step.increase = 0.08, 
                             size = 3.5, label = "p.signif", vjust = 0.5) +
  
  scale_fill_manual(values = colors, name = "Cancer Type") +
  coord_cartesian(ylim = c(-0.15, 1)) +
  labs(x = NULL, y = expression(paste("Spearman's ", rho)), subtitle = "TCGA tumor samples") + 
  
  theme_classic(base_size = 12) +
  theme(
    aspect.ratio = 1.2,
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 11),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    legend.position = "none",
    plot.margin = margin(t = 15, r = 15, b = 10, l = 10),
    plot.subtitle = element_text(size = 11, face = "plain", hjust = 0, margin = margin(b = 8))
  )

fig.1e <- add_sample_size(p, merge_equal_n = TRUE, equal_n_pos = "x_title")

pdf("./Proj_CTAs/output/figure.panels/Fig. 1e.pdf", width = 5, height = 5)
fig.1e
dev.off()

# Fig. 1f
load("./Datas/UCSC_Xena/normal.tumor.RData")  ## c("gene.info", "sample.ls", "n.samples", "log2.tpm")
auc.cutoff        <- 0.70
p.adj.cutoff      <- 0.05
cancers.included  <- setdiff(names(sample.ls), "TGCT")
load("./Proj_CTAs/output/manuscript.revision/bulk.tumor.normal.auc.RData")
n.all.genes    <- sapply(CTGset, \(x) gene.info |> filter(ensembl %in% x, id %in% rownames(log2.tpm)) |> nrow())
pct.sig.all.df <- map_dfr(cancers.included, \(c) {
    auc.df.sub    <- bulk.tumor.normal.auc.df |> filter(cancer == c, auc > auc.cutoff, p.adj < p.adj.cutoff)
    n.sig.genes   <- sapply(CTGset, \(x) auc.df.sub |> filter(ensembl %in% x) |> nrow())
    pct.sig.genes <- n.sig.genes / n.all.genes

    data.frame(cancer = c, CTGset = names(CTGset), pct = pct.sig.genes, row.names = NULL)
})
n.all.genes       <- sapply(CTGset.unique, \(x) gene.info |> filter(ensembl %in% x, id %in% rownames(log2.tpm)) |> nrow())
pct.sig.unique.df <- map_dfr(cancers.included, \(c) {
    auc.df.sub    <- bulk.tumor.normal.auc.df |> filter(cancer == c, auc > auc.cutoff, p.adj < p.adj.cutoff)
    n.sig.genes   <- sapply(CTGset.unique, \(x) auc.df.sub |> filter(ensembl %in% x) |> nrow())
    pct.sig.genes <- n.sig.genes / n.all.genes

    data.frame(cancer = c, CTGset = names(CTGset.unique), pct = pct.sig.genes, row.names = NULL)
})
fill.colors  <- c("#FFFFFF", "#FEE8C8", "#F4A582", "#D47264", "#9E3D3A")
cancer.order <- cancers.included
n.cancers    <- length(cancer.order)

ctg.order    <- c("CTDatabase", "CTGs_Wang_NC_2016", "CTGs_Carter_JITC_2023", "CTGs_us")
subset.order <- c("CTDatabase.included", "CTDatabase.excluded")
ctg.labels   <- c(
  CTDatabase = "CTDatabase",
  CTGs_Wang_NC_2016 = "Wang et al. (2016)",
  CTGs_Carter_JITC_2023 = "Carter et al. (2023)",
  CTGs_us = "This study",
  CTDatabase.included = "CTDatabase (retained)",
  CTDatabase.excluded = "CTDatabase (excluded)"
)

row1.df    <- tibble(CTGset = ctg.order, y = rev(seq_along(ctg.order)), row.label = unname(ctg.labels[ctg.order]), panel = "all")
row2.order <- c(ctg.order, subset.order)
row2.df    <- tibble(CTGset = row2.order, y = rev(seq_along(row2.order))) |>
  mutate(y = y + if_else(CTGset %in% ctg.order, 0.8, 0), row.label = unname(ctg.labels[CTGset]), panel = "unique")

plot1.df <- pct.sig.all.df |> filter(CTGset %in% ctg.order) |> left_join(row1.df, by = "CTGset")
## keep the CTDatabase subset values identical to the full-set summary.
plot2.df <- bind_rows(
  pct.sig.unique.df |> filter(CTGset %in% ctg.order),
  pct.sig.all.df |> filter(CTGset %in% subset.order)
) |> left_join(row2.df, by = "CTGset")

plot.df  <- bind_rows(plot1.df, plot2.df) |>
  mutate(
    x            = match(cancer, cancer.order),
    fill.pct     = if_else(is.finite(pct), pct, NA_real_),
    value.label  = if_else(is.finite(pct), sprintf("%.1f", 100 * pct), "NA"),
    number.color = if_else(is.finite(pct) & pct >= 0.85, "white", "black")
  )

plot1.df <- plot.df |> filter(panel == "all")
plot2.df <- plot.df |> filter(panel == "unique")
row.df   <- bind_rows(row1.df, row2.df)

## statistical annotations
comparison.order <- c("CTDatabase", "CTGs_Wang_NC_2016", "CTGs_Carter_JITC_2023")
comparison.df    <- expand_grid(panel = c("all", "unique"), group2 = comparison.order) |>
  mutate(group1 = "CTGs_us", track = length(comparison.order) + 1L - match(group2, comparison.order)) |>
  bind_rows(
    tibble(panel = "unique", group1 = "CTDatabase.included", group2 = "CTDatabase.excluded", track = 1L)
  )

## align paired observations by cancer type before testing.
stat.df <- map_dfr(seq_len(nrow(comparison.df)), \(i) {
  comparison.row <- comparison.df[i, ]
  group.keys     <- c(comparison.row$group1, comparison.row$group2)

  pair.df <- plot.df |>
    filter(panel == comparison.row$panel, CTGset %in% group.keys) |>
    dplyr::select(cancer, CTGset, pct) |>
    pivot_wider(id_cols = cancer, names_from = CTGset, values_from = pct) |>
    filter(if_all(all_of(group.keys), is.finite))

  p.raw <- if (nrow(pair.df) > 0L) {
    wilcox.test(pair.df[[group.keys[1]]], pair.df[[group.keys[2]]], paired = TRUE, alternative = "two.sided")$p.value
  } else {
    NA_real_
  }

  bind_cols(comparison.row, tibble(n.pairs = nrow(pair.df), p.raw = p.raw))
}) |>
  left_join(row.df |> dplyr::select(panel, group1 = CTGset, y1 = y), by = c("panel", "group1")) |>
  left_join(row.df |> dplyr::select(panel, group2 = CTGset, y2 = y), by = c("panel", "group2")) |>
  mutate(
    y.low = pmin(y1, y2), y.high = pmax(y1, y2),
    x.bracket = n.cancers + 0.7 + 0.65 * (track - 1),
    x.label = x.bracket + 0.27,
    y.label = (y.low + y.high) / 2,
    p.label = case_when(
      !is.finite(p.raw) ~ "NA",
      p.raw <= 0.0001 ~ "****",
      p.raw <= 0.001 ~ "***",
      p.raw <= 0.01 ~ "**",
      p.raw <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

bracket.df <- bind_rows(
  stat.df |> transmute(panel, x = x.bracket, xend = x.bracket, y = y.low, yend = y.high),
  stat.df |> transmute(panel, x = x.bracket - 0.15, xend = x.bracket, y = y.low, yend = y.low),
  stat.df |> transmute(panel, x = x.bracket - 0.15, xend = x.bracket, y = y.high, yend = y.high)
)

stat2.df <- stat.df |> filter(panel == "unique")

p <- ggplot(mapping = aes(x = x, y = y)) +
  geom_tile(aes(fill = fill.pct), width = 0.96, height = 0.90, color = "white", linewidth = 0.20) +
  geom_text(aes(label = value.label, color = number.color), size = 2.7) +
  scale_color_identity() +
  scale_fill_gradientn(
    colors = fill.colors,
    values = c(0, 0.20, 0.40, 0.70, 1),
    limits = c(0, 1), breaks = seq(0, 1, 0.25),
    labels = scales::label_number(scale = 100, accuracy = 1), na.value = "grey90",
    guide = guide_colorbar(
      direction = "horizontal", title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(50, "mm"), barheight = grid::unit(3, "mm")
    )
  ) +
  scale_x_continuous(breaks = seq_along(cancer.order), labels = cancer.order) +
  labs(x = NULL, y = NULL, fill = "Tumor-enriched genes (%)") +
  theme_minimal(base_size = 11) +
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.text.x = element_text(size = 8.5, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9.5, color = "black"),
    axis.ticks = element_blank(),

    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "plain"),
    legend.text = element_text(size = 8.5),
    legend.margin = margin(t = 5),

    plot.margin = margin(t = 7, r = 8, b = 6, l = 7)
  )
p2 <- p %+% plot2.df +
  geom_segment(
    data = filter(bracket.df, panel == "unique"), aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE, color = "black", linewidth = 0.30
  ) +
  geom_text(
    data = stat2.df, aes(x = x.label, y = y.label, label = p.label),
    inherit.aes = FALSE, angle = 90, hjust = 0.5, vjust = 0.5, size = 3.2, color = "black"
  ) +
  scale_y_continuous(breaks = row2.df$y, labels = row2.df$row.label) +
  coord_cartesian(
    xlim = c(0.5, n.cancers + 2.6), ylim = c(0.5, max(row2.df$y) + 0.5),
    expand = FALSE, clip = "off"
  )
row.order    <- cancer.order[cancer.order %in% as.character(p2$data$cancer)]
column.order <- c(ctg.order, subset.order)
n.cancers    <- length(row.order)

column.df <- tibble(CTGset = column.order, x = seq_along(column.order)) |>
  mutate(x = x + if_else(CTGset %in% subset.order, 0.12, 0), label = unname(ctg.labels[CTGset]))

plot.df <- p2$data |>
  dplyr::select(cancer, CTGset, fill.pct, value.label, number.color) |>
  mutate(cancer = as.character(cancer), CTGset = as.character(CTGset)) |>
  left_join(column.df, by = "CTGset") |>
  mutate(y = n.cancers + 1L - match(cancer, row.order))

label.df <- stat2.df |>
  dplyr::select(group1, group2, p.label) |>
  left_join(column.df |> transmute(group1 = CTGset, x1 = x), by = "group1") |>
  left_join(column.df |> transmute(group2 = CTGset, x2 = x), by = "group2") |>
  mutate(
    x.low = pmin(x1, x2), x.high = pmax(x1, x2),
    track = dense_rank(x.high - x.low),
    track = if_else(group1 %in% subset.order & group2 %in% subset.order, max(track), track),
    y.bracket = n.cancers + 0.70 + 0.50 * (track - 1),
    x.label = (x.low + x.high) / 2, y.label = y.bracket + 0.20
  )

segment.df <- bind_rows(
  label.df |> transmute(x = x.low, xend = x.high, y = y.bracket, yend = y.bracket),
  label.df |> transmute(x = x.low, xend = x.low, y = y.bracket, yend = y.bracket - 0.10),
  label.df |> transmute(x = x.high, xend = x.high, y = y.bracket, yend = y.bracket - 0.10)
)

p2 <- ggplot(plot.df, aes(x = x, y = y)) +
  geom_tile(aes(fill = fill.pct), width = 0.96, height = 0.92, color = "white", linewidth = 0.20) +
  geom_text(aes(label = value.label, color = number.color), size = 2.7) +
  geom_segment(
    data = segment.df, aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE, color = "grey25", linewidth = 0.28
  ) +
  geom_text(
    data = label.df, aes(x = x.label, y = y.label, label = p.label),
    inherit.aes = FALSE, angle = 0, size = 3.2, color = "black"
  ) +
  scale_color_identity() +
  scale_fill_gradientn(
    colors = fill.colors, values = c(0, 0.20, 0.40, 0.70, 1),
    limits = c(0, 1), breaks = seq(0, 1, 0.25),
    labels = scales::label_number(scale = 100, accuracy = 1), na.value = "grey90",
    guide = guide_colorbar(
      direction = "vertical", title.position = "left",
      barwidth = grid::unit(3, "mm"), barheight = grid::unit(38, "mm")
    )
  ) +
  scale_x_continuous(breaks = column.df$x, labels = column.df$label, position = "bottom") +
  scale_y_continuous(breaks = seq_len(n.cancers), labels = rev(row.order), position = "left") +
  coord_cartesian(
    xlim = c(0.5, max(column.df$x) + 0.5), ylim = c(0.5, max(label.df$y.label) + 0.45),
    expand = FALSE, clip = "off"
  ) +
  labs(x = NULL, y = NULL, fill = "Tumor-enriched genes (%)") +
  theme_minimal(base_size = 11) +
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9, color = "black", angle = 0, hjust = 1),
    axis.ticks = element_blank(),
    
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "plain", angle = 90, hjust = 0.5, vjust = 0.5, margin = margin(r = 5)),
    legend.text = element_text(size = 8.5, angle = 0),
    legend.box.spacing = grid::unit(3, "mm"),
    plot.margin = margin(t = 7, r = 8, b = 7, l = 7)
  )

figure.width  <- max(9, 0.34 * n.cancers + 3.0) / 2
figure.height <- max(7, 0.28 * n.cancers + 2.4)

pdf("./Proj_CTAs/output/figure.panels/bulk.CTG.discrimination.unique.pdf", width = figure.width, height = figure.height)
print(p2)
dev.off()

##- Fig. 1g --------------------------------------------------------------------
load("./Proj_CTAs/output/manuscript.revision/CTG.annotation.RData")
criterion.order  <- c("Meets", "Does not meet", "Not evaluable")
criterion.colors <- c("Meets" = "#D47264", "Does not meet" = "grey75", "Not evaluable" = "white")

plot.df1   <- CTGs.annotated |> mutate(x = log2(max.non.testis.TPM + 1), y = log2(testis.TPM + 1), tier1 = factor(tier1, levels = criterion.order))
plot.df2   <- CTGs.annotated |> mutate(tier2 = factor(tier2, levels = criterion.order))

overlap.df <- CTGs.annotated |>
  transmute(
    tier1 = factor(tier1, levels = rev(criterion.order)),
    tier2 = factor(tier2, levels = criterion.order)
  ) |>
  count(tier1, tier2, name = "n", .drop = FALSE) |>
  mutate(
    group = case_when(
      tier1 == "Not evaluable" | tier2 == "Not evaluable" ~ "Not evaluable",
      tier1 == "Meets" & tier2 == "Meets"                 ~ "Both criteria",
      tier1 == "Meets"                                    ~ "Minimal expression only",
      tier2 == "Meets"                                    ~ "Tumor enrichment only",
      TRUE                                                ~ "Neither criterion"
    ),
    label = sprintf("%d\n(%.1f%%)", n, 100 * n / nrow(CTGs.annotated))
  )
overlap.colors <- c(
  "Both criteria"           = "#9970AB",
  "Minimal expression only" = "#6B9AC4",
  "Tumor enrichment only"   = "#D47264",
  "Neither criterion"       = "#BFBFBF",
  "Not evaluable"           = "#E7E6E6"
)

n.row <- nlevels(overlap.df$tier1)
n.col <- nlevels(overlap.df$tier2)

p1 <- ggplot(plot.df1, aes(x = x, y = y)) +
  geom_abline(intercept = 0, slope = 1, color = "grey75", linetype = "dashed", linewidth = 0.35) +
  geom_vline(xintercept = log2(3 + 1), color = "grey40", linetype = "dashed", linewidth = 0.35) +
  geom_point(
    aes(fill = tier1),
    shape = 21, size = 1.65, color = "grey35", stroke = 0.15, alpha = 0.75,
    na.rm = TRUE, show.legend = TRUE
  ) +
  annotate(
    "text", x = log2(3 + 1), y = Inf, label = "3 TPM",
    hjust = 1.15, vjust = 1.5, size = 2.8, color = "grey30"
  ) +
  scale_fill_manual(
    values = criterion.colors, limits = criterion.order, breaks = criterion.order,
    labels = criterion.order, name = NULL, drop = FALSE
  ) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 4), expand = expansion(mult = 0.05)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 4), expand = expansion(mult = 0.05)) +
  expand_limits(x = 0, y = 0) +
  labs(
    x = expression(paste("Max. non-testis ", log[2], "(TPM + 1)")),
    y = expression(paste("Testis ", log[2], "(TPM + 1)"))
  ) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2.6, alpha = 1))) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10.5, margin = margin(t = 7)),
    axis.title.y = element_text(size = 10.5, margin = margin(r = 7)),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.30),

    panel.grid = element_blank(),
    aspect.ratio = 1
  )

p2 <- ggplot(plot.df2, aes(x = n.cancers, y = max.auc)) +
  geom_vline(xintercept = 0.5, color = "grey70", linetype = "dashed", linewidth = 0.35) +
  geom_point(
    aes(fill = tier2),
    shape = 21, size = 1.65, color = "grey35", stroke = 0.15, alpha = 0.75,
    na.rm = TRUE, show.legend = TRUE
  ) +
  scale_fill_manual(
    values = criterion.colors, limits = criterion.order, breaks = criterion.order,
    labels = criterion.order, name = NULL, drop = FALSE
  ) +
  scale_x_continuous(breaks = scales::breaks_width(4), expand = expansion(add = 0.65)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = expansion(mult = c(0.05, 0.07))) +
  labs(x = "Cancer types with tumor enrichment", y = "Maximum AUROC") +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2.6, alpha = 1))) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10.5, margin = margin(t = 7)),
    axis.title.y = element_text(size = 10.5, margin = margin(r = 7)),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.30),

    panel.grid = element_blank(),
    aspect.ratio = 1
  )

p3 <- ggplot(overlap.df, aes(x = tier2, y = tier1)) +
  geom_tile(
    aes(fill = group),
    width = 0.96, height = 0.96, color = "white", linewidth = 0.5, alpha = 0.35
  ) +
  geom_text(aes(label = label), size = 3.6, lineheight = 1.15, color = "black") +

  ## draw column symbols below the matrix
  annotation_custom(
    grob = grid::pointsGrob(
      x = grid::unit((seq_len(n.col) - 0.5) / n.col, "npc"),
      y = grid::unit(rep(0, n.col), "npc") - grid::unit(rep(3.2, n.col), "mm"),
      pch = 21, size = grid::unit(3, "mm"),
      gp = grid::gpar(
        fill = unname(criterion.colors[levels(overlap.df$tier2)]),
        col = "grey35", lwd = 0.6
      )
    ),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +

  ## match row symbols to factor levels from bottom to top
  annotation_custom(
    grob = grid::pointsGrob(
      x = grid::unit(rep(-3.2, n.row), "mm"),
      y = grid::unit((seq_len(n.row) - 0.5) / n.row, "npc"),
      pch = 21, size = grid::unit(3, "mm"),
      gp = grid::gpar(
        fill = unname(criterion.colors[levels(overlap.df$tier1)]),
        col = "grey35", lwd = 0.6
      )
    ),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +

  scale_fill_manual(values = overlap.colors, guide = "none") +
  scale_x_discrete(labels = NULL, position = "bottom", drop = FALSE) +
  scale_y_discrete(labels = NULL, drop = FALSE) +
  coord_cartesian(
    xlim = c(0.5, n.col + 0.5), ylim = c(0.5, n.row + 0.5),
    expand = FALSE, clip = "off"
  ) +
  labs(x = "Tumor enrichment", y = "Minimal non-testicular expression") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x.bottom = element_text(size = 10.5, margin = margin(t = 21)),
    axis.title.y.left = element_text(size = 10.5, margin = margin(r = 17)),

    panel.grid = element_blank(),
    aspect.ratio = 1
  )

p <- patchwork::wrap_plots(p1, p2, p3, nrow = 1, widths = c(1, 1, 1), guides = "collect") &
  theme(
    text = element_text(family = "sans", color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = 9),
    legend.key = element_blank(),
    legend.key.width = grid::unit(5, "mm"),
    legend.key.height = grid::unit(4, "mm"),
    legend.margin = margin(t = 5, r = 0, b = 0, l = 0),

    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 6, r = 4, b = 6, l = 4)
  )

pdf("./Proj_CTAs/output/figure.panels/CTG.annotation.pdf", width = 10.8, height = 4.4)
print(p)
dev.off()

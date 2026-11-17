library(tidyverse)
library(magrittr)
library(foreach)
library(doParallel)
library(patchwork)

##- Fig. 4a --------------------------------------------------------------------

##- Fig. 4b --------------------------------------------------------------------
load("./Proj_CTAs/output/manuscript.revision/CTG.metrics.ls.RData")

registerDoParallel(cores = 40)
residuals.df <- foreach(s = names(CTG.metrics.ls), .combine = 'rbind') %dopar% {
  metrics.df <- CTG.metrics.ls[[s]]
  if (sum(metrics.df$malignant == "no") < 10) return(NULL)

  r1 <- lm(n.CTGs ~ n.genes, data = metrics.df)$residuals
  r2 <- lm(n.ctrl.genes ~ n.genes, data = metrics.df)$residuals

  x1 <- median(r1[metrics.df$malignant == "yes"])
  x2 <- median(r1[metrics.df$malignant == "no"])
  x3 <- median(r2[metrics.df$malignant == "yes"])
  x4 <- median(r2[metrics.df$malignant == "no"])
  
  data.frame(
    dataset_sample = s, 
    malignant      = c("malignant", "non-malignant", "malignant", "non-malignant"),
    gene           = c("n.CTGs", "n.CTGs", "n.ctrl.genes", "n.ctrl.genes"),
    residuals      = c(x1, x2, x3, x4)
  )
}
stopImplicitCluster()

malignant.colors <- c("malignant" = "#E64B35", "non-malignant" = "#4DBBD5")
facet.labels     <- c(
  "n.CTGs"       = "italic(N)[CTG]",
  "n.ctrl.genes" = "italic(N)[control]"
)

load("./Proj_CTAs/output/manuscript.revision/tumor.sample.info.updated.RData")
plot.df <- residuals.df |>
  inner_join(tumor.sample.info, by = "dataset_sample") |>
  mutate(
    malignant = factor(malignant, levels = c("malignant", "non-malignant")),
    gene = factor(gene, levels = c("n.CTGs", "n.ctrl.genes"))
  ) |>

  mutate(dataset_patient = paste0(data_folder, "_", patient)) |>
  group_by(dataset_group, dataset_patient, gene, malignant) |>
  summarise(residuals = median(residuals, na.rm = TRUE), .groups = "drop")

medians <- plot.df |>
  dplyr::filter(malignant == "malignant") |>
  group_by(dataset_group, gene) |>
  summarise(median = median(residuals), .groups = "drop")

label.df <- plot.df |>
  group_by(gene) |>
  summarise(
    y.min = min(residuals),
    y.max = max(residuals),
    .groups = "drop"
  ) |>
  mutate(bracket.y = y.max + 0.04 * pmax(y.max - y.min, 1))

p <- ggplot(plot.df, aes(x = malignant, y = residuals)) +
  geom_hline(
    data = medians, aes(yintercept = median),
    inherit.aes = FALSE,
    linetype = "dashed", color = "grey50", linewidth = 0.35
  ) +
  geom_point(
    aes(color = malignant),
    position = position_jitter(width = 0.12, height = 0, seed = 2026),
    shape = 21, fill = "white", size = 1.4, stroke = 0.25, alpha = 0.85
  ) +
  geom_boxplot(
    aes(fill = malignant),
    width = 0.5, outlier.shape = NA,
    color = "grey25", linewidth = 0.35, alpha = 0.2
  ) +
  ggpubr::stat_compare_means(
    data = dplyr::filter(plot.df, gene == "n.CTGs"),
    comparisons = list(c("malignant", "non-malignant")),
    method = "wilcox.test", paired = TRUE,
    label = "p.signif", hide.ns = FALSE,
    label.y = label.df$bracket.y[label.df$gene == "n.CTGs"],
    tip.length = 0.01, bracket.size = 0.3,
    size = 3.4, vjust = 0
  ) +
  ggpubr::stat_compare_means(
    data = dplyr::filter(plot.df, gene == "n.ctrl.genes"),
    comparisons = list(c("malignant", "non-malignant")),
    method = "wilcox.test", paired = TRUE,
    label = "p.signif", hide.ns = FALSE,
    label.y = label.df$bracket.y[label.df$gene == "n.ctrl.genes"],
    tip.length = 0.01, bracket.size = 0.3,
    size = 3.4, vjust = 0
  ) +
  scale_color_manual(values = malignant.colors, guide = "none") +
  scale_fill_manual(values = malignant.colors, guide = "none") +
  scale_x_discrete(
    labels = c("malignant" = "Malignant", "non-malignant" = "Non-malignant"),
    expand = expansion(add = 0.5)
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 5),
    expand = expansion(mult = c(0.04, 0.12))
  ) +
  facet_grid(
    gene ~ dataset_group, scales = "free_y",
    labeller = labeller(gene = as_labeller(facet.labels, default = label_parsed))
  ) +
  labs(x = NULL, y = "Median residual") +
  theme_minimal(base_size = 11) +
  theme(
    text = element_text(family = "sans", color = "black"),
    aspect.ratio = 1.2,

    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey95", linewidth = 0.25),

    axis.title.y = element_text(size = 11, margin = margin(r = 7)),
    axis.text.x = element_text(
      size = 9.5, color = "black",
      angle = 30, hjust = 1, margin = margin(t = 5)
    ),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.3),

    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.35),
    strip.text = element_text(size = 10.5, face = "plain", margin = margin(5, 5, 5, 5)),
    panel.spacing.x = grid::unit(0.75, "lines"),
    panel.spacing.y = grid::unit(0.9, "lines"),

    legend.position = "none",
    plot.margin = margin(t = 7, r = 7, b = 7, l = 7)
  )

pdf("./Proj_CTAs/output/figure.panels/adjusted n_CTG and n_ctrl.pdf", width = 6.5, height = 5.5)
print(p)
dev.off()

##- Fig. 4c --------------------------------------------------------------------
load("./Proj_CTAs/output/manuscript.revision/CTG.metrics.ls.RData")
load("./Proj_CTAs/output/manuscript.revision/tumor.sample.info.updated.RData")
tumor.sample.info  <- tumor.sample.info |> mutate(dataset_patient = paste(data_folder, patient, sep = "_"))
all.samples        <- intersect(tumor.sample.info$dataset_sample, names(CTG.metrics.ls))
sample.patient.map <- tumor.sample.info |> distinct(dataset_sample, dataset_patient, cancer, dataset_group)

## calculate the median of CTG burden for each cell group cross individual samples
med.n.CTGs.df <- foreach(s = all.samples, .combine = 'rbind') %do% {
  df <- CTG.metrics.ls[[s]]
  bind_rows(
    df |> group_by(normal.ref.cell) |> summarise(med.n.CTGs = median(n.CTGs, na.rm = TRUE), n.cell = n()) |> 
      rename(cell.type = normal.ref.cell) |> 
      mutate(cell.group = "group 2")
  ) |>
  mutate(dataset_sample = s)
}
source("./Proj_CTAs/code/visualization.R")
plot.df2 <- med.n.CTGs.df |> filter(cell.group == "group 2", cell.type %in% c("tumour", "normal"), n.cell >= 10, dataset_sample %in% all.samples) |>
  pivot_wider(id_cols = c("dataset_sample"), names_from = cell.type, values_from = med.n.CTGs) |>
  drop_na() |>

  left_join(sample.patient.map, by = "dataset_sample") |>
  group_by(dataset_group, cancer, dataset_patient) |>
  summarise(tumour = median(tumour), normal = median(normal), .groups = "drop") |>
  left_join(cancer.abbr, by = "cancer") |>
  mutate(abbr = ifelse(abbr %in% names(cancer.colours), abbr, "Other"))   # cancer types of Validation Cohort 2

my.colors     <- c(cancer.colours[intersect(names(cancer.colours), plot.df2$abbr)], Other = "#5F9EA0")
plot.df2$abbr <- factor(plot.df2$abbr, levels = names(my.colors))
axis.limit    <- max(plot.df2$tumour, plot.df2$normal) * 1.05

## per-cohort statistics; patients are the units of the paired test
label.df <- plot.df2 |>
  group_by(dataset_group) |>
  summarise(
    n = n(), pct = round(mean(tumour > normal) * 100, 1),
    p = if (n() >= 2) wilcox.test(tumour, normal, paired = TRUE)$p.value else NA_real_, .groups = "drop"
  ) |>
  mutate(
    p.sci = formatC(p, format = "e", digits = 2),
    line1 = sprintf("%.1f*'%%'~(italic(n)==%d)", pct, n),
    line2 = ifelse(
      p < 0.001,
      sprintf("italic(P)==%s%%*%%10^%d", sub("e.*", "", p.sci), as.integer(sub(".*e", "", p.sci))),
      sprintf("italic(P)=='%.2f'", p)
    ),
    x = axis.limit * 0.05, y1 = axis.limit * 0.95, y2 = axis.limit * 0.87
  )

p <- ggplot(plot.df2, aes(x = normal, y = tumour)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#B2182B", linewidth = 0.6) +
  geom_point(aes(color = abbr, fill = abbr), shape = 21, size = 1.8, stroke = 0.2, alpha = 0.5) +
  geom_text(data = label.df, aes(x = x, y = y1, label = line1), parse = TRUE, hjust = 0, vjust = 1, size = 3.3, inherit.aes = FALSE) +
  geom_text(data = label.df, aes(x = x, y = y2, label = line2), parse = TRUE, hjust = 0, vjust = 1, size = 3.3, inherit.aes = FALSE) +
  facet_wrap(~ dataset_group, nrow = 1) +
  coord_fixed(ratio = 1, xlim = c(0, axis.limit), ylim = c(0, axis.limit), expand = FALSE) +
  scale_color_manual(values = my.colors, name = NULL) +
  scale_fill_manual(values = my.colors, name = NULL) +
  labs(x = bquote("Normal reference cells (median" ~ italic(N)[CTG] * ")"), y = bquote("Malignant cells (median" ~ italic(N)[CTG] * ")")) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),

    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.35),
    strip.text = element_text(size = 11),
    panel.spacing = grid::unit(1, "lines"),

    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1, stroke = 0.4)))

pdf("./Proj_CTAs/output/figure.panels/nCTGs.comparison.with.normal.ref.cells.pdf", width = 12, height = 4.5)
p
dev.off()

##- Fig. 4d --------------------------------------------------------------------

##- Fig. 4e --------------------------------------------------------------------
registerDoParallel(cores = 50)

## metrics used to predict malignant cells
load("./Proj_CTAs/output/manuscript.revision/CTG.metrics.ls.RData")
cnv.score.ls       <- readRDS("./Proj_CTAs/output/6. n_CTGs/scRNA-seq/CTG.SingleR.cnv.score.ls.rds")
scCancer2.score.ls <- readRDS("./Proj_CTAs/output/6. n_CTGs/scRNA-seq/scCancer2.score.ls.rds")
ikarus.score.ls    <- readRDS("./Proj_CTAs/output/6. n_CTGs/scRNA-seq/ikarus.score.ls.rds")

## AUROCs and AUPRCs for each metric in each tumor sample
auc.df <- foreach(s = names(CTG.metrics.ls), .combine = 'rbind') %dopar% {
  metrics.df      <- CTG.metrics.ls[[s]]
  cnv.score       <- if (s %in% names(cnv.score.ls)) {
    cnv.score.ls[[s]] |> select(barcode, cnv.score = cnv_score)
  } else {
    data.frame(barcode = metrics.df$barcode, cnv.score = NA_real_)
  }
  scCancer2.score <- if (s %in% names(scCancer2.score.ls)) {
    scCancer2.score.ls[[s]] |> select(barcode, scCancer2.score = scCancer2_malignant_score)
  } else {
    data.frame(barcode = metrics.df$barcode, scCancer2.score = NA_real_)
  }
  ikarus.score    <- if (s %in% names(ikarus.score.ls)) {
    ikarus.score.ls[[s]] |> select(barcode, ikarus.score = ikarus_malignant_score)
  } else {
    data.frame(barcode = metrics.df$barcode, ikarus.score = NA_real_)
  }

  metrics.df <- list(metrics.df, cnv.score, scCancer2.score, ikarus.score) |>
    lapply(\(x) {x |> mutate(barcode = as.character(barcode))}) |>
    reduce(inner_join, by = "barcode")

  vars <- c("n.CTGs", "MaligCTG.malignant.score", "cnv.score", "scCancer2.score", "ikarus.score")
  df1  <- metrics.df
  y    <- ifelse(df1$malignant == "yes", 1, 0)
  if (min(sum(y == 0), sum(y == 1)) >= 10) {

    auroc <- vapply(df1[, vars], \(x) {
      if(all(is.na(x))) return(NA_real_)
      pROC::roc(y, x, direction = "<")$auc
    }, 0)

      auprc <- vapply(df1[, vars], \(x) {
        if(all(is.na(x))) return(NA_real_)
        PRROC::pr.curve(scores.class0 = x, weights.class0 = y)$auc.integral
      }, 0)

  } else {
    auroc <- NA_real_
    auprc <- NA_real_
  }
  res1 <- data.frame(dataset_sample = s, var = vars, auroc, auprc, ctrl.cells = "non-malignant cells", row.names = NULL)

  df2  <- metrics.df |> filter(normal.ref.cell %in% c("tumour", "normal"))
  y    <- ifelse(df2$normal.ref.cell == "tumour", 1, 0)
  if (min(sum(y == 0), sum(y == 1)) >= 10) {

    auroc <- vapply(df2[, vars], \(x) {
      if(all(is.na(x))) return(NA_real_)
      pROC::roc(y, x, direction = "<")$auc
    }, 0)

    auprc <- vapply(df2[, vars], \(x) {
      if(all(is.na(x))) return(NA_real_)
      PRROC::pr.curve(scores.class0 = x, weights.class0 = y)$auc.integral
    }, 0)

  } else {
    auroc <- NA_real_
    auprc <- NA_real_
  }
  res2 <- data.frame(dataset_sample = s, var = vars, auroc, auprc, ctrl.cells = "normal reference cells", row.names = NULL)

  rbind(res1, res2) |>
    mutate(
      n.mal        = sum(df1$malignant == "yes"),
      n.non.mal    = sum(df1$malignant == "no"),
      n.normal.ref = sum(df2$normal.ref.cell == "normal")
    )
}

stopImplicitCluster()
load("./Proj_CTAs/output/manuscript.revision/tumor.sample.info.updated.RData")
method.colors   <- c(
  "n.CTG"     = "#7F7F7F",
  "MaligCTG"  = "#DC0000B2",
  
  "inferCNV"  = "#3C5488B2",
  "CopyKAT"   = "#00A087B2",
  "SCEVAN"    = "#4DBBD5B2",
  "scCancer2" = "#8491B4B2",
  "ikarus"    = "#91D1C2B2"
)
x.axis.labels <- c(
  "n.CTG"     = expression(italic(N)[CTG]),
  "MaligCTG"  = "MaligCTG",
  "inferCNV"  = "inferCNV",
  "CopyKAT"   = "CopyKAT",
  "SCEVAN"    = "SCEVAN",
  "scCancer2" = "scCancer2",
  "ikarus"    = "ikarus"
)

plot.df1 <- auc.df |>
  inner_join(tumor.sample.info, by = "dataset_sample") |>
  filter(ctrl.cells == "non-malignant cells", n.non.mal >= 10) |> 

  mutate(method = case_match(
    var,
    "n.CTGs"                   ~ "n.CTG",
    "MaligCTG.malignant.score" ~ "MaligCTG",
    "cnv.score"                ~ "inferCNV",
    "ikarus.score"             ~ "ikarus",
    "scCancer2.score"          ~ "scCancer2",
    .default = as.character(var)
  )) |>
  mutate(method = factor(method, levels = c("n.CTG", "MaligCTG", "inferCNV", "CopyKAT", "SCEVAN", "scCancer2", "ikarus"))) |>
  
  pivot_longer(cols = c("auroc", "auprc"), names_to = "metric", values_to = "auc_value") |>
  mutate(metric = factor(metric, levels = c("auroc", "auprc"), labels = c("AUROC", "AUPRC"))) |>

  mutate(dataset_patient = paste0(data_folder, "_", patient)) |>
  group_by(dataset_group, metric, method, dataset_patient) |>
  summarise(auc_value = median(auc_value, na.rm = TRUE), .groups = "drop")
plot.df2 <- auc.df |>
  inner_join(tumor.sample.info, by = "dataset_sample") |>
  filter(ctrl.cells == "normal reference cells", n.normal.ref >= 10) |>

  mutate(method = case_match(
    var,
    "n.CTGs"                   ~ "n.CTG",
    "MaligCTG.malignant.score" ~ "MaligCTG",
    "cnv.score"                ~ "inferCNV",
    "ikarus.score"             ~ "ikarus",
    "scCancer2.score"          ~ "scCancer2",
    .default = as.character(var)
  )) |>
  mutate(method = factor(method, levels = c("n.CTG", "MaligCTG", "inferCNV", "CopyKAT", "SCEVAN", "scCancer2", "ikarus"))) |>
  
  pivot_longer(cols = c("auroc", "auprc"), names_to = "metric", values_to = "auc_value") |>
  mutate(metric = factor(metric, levels = c("auroc", "auprc"), labels = c("AUROC", "AUPRC"))) |>

  mutate(dataset_patient = paste0(data_folder, "_", patient)) |>
  group_by(dataset_group, metric, method, dataset_patient) |>
  summarise(auc_value = median(auc_value, na.rm = TRUE), .groups = "drop")
compa.methods <- c("n.CTG", "inferCNV", "scCancer2", "ikarus")
cont.df       <- bind_rows("non-malignant cells" = plot.df1, "normal reference cells" = plot.df2, .id = "ctrl.cells")

stat.df <- cont.df |>
  filter(method == "MaligCTG") |>
  select(ctrl.cells, dataset_group, metric, dataset_patient, value1 = auc_value) |>
  inner_join(
    cont.df |> filter(method %in% compa.methods) |> select(ctrl.cells, dataset_group, metric, dataset_patient, group2 = method, value2 = auc_value),
    by = c("ctrl.cells", "dataset_group", "metric", "dataset_patient")
  ) |>
  filter(!is.na(value1), !is.na(value2)) |>
  group_by(ctrl.cells, dataset_group, metric, group2) |>
  summarise(n.pairs = n(), p = suppressWarnings(wilcox.test(value1, value2, paired = TRUE)$p.value), .groups = "drop") |>
  mutate(group1 = "MaligCTG", group2 = factor(group2, levels = compa.methods)) |>
  arrange(ctrl.cells, dataset_group, metric, group2) |>
  group_by(ctrl.cells, dataset_group, metric) |>
  mutate(p.adj = p.adjust(p, method = "BH"), y.position = 1.02 + 0.06 * (row_number() - 1)) |>
  ungroup() |>
  mutate(
    group2   = as.character(group2),
    p.signif = cut(p, breaks = c(-Inf, 1e-4, 1e-3, 0.01, 0.05, Inf), labels = c("****", "***", "**", "*", "ns"))
  )
median.lines2 <- plot.df2 |> filter(method == "MaligCTG") |> group_by(dataset_group, metric) |>
  summarise(median_val = median(auc_value, na.rm = TRUE), .groups = "drop")
p1 <- ggplot(plot.df1, aes(x = method, y = auc_value)) +
  geom_boxplot(aes(color = method, fill = method), outlier.shape = NA, alpha = 0.2, width = 0.6, linewidth = 0.3) +
  geom_point(aes(color = method), size = 0.1, alpha = 0.4, position = position_jitter(width = 0.2, height = 0, seed = 1)) +
  facet_grid(metric ~ dataset_group, scales = "free_y", switch = "y") +
  scale_color_manual(values = method.colors) +
  scale_fill_manual(values = method.colors) +
  scale_x_discrete(labels = x.axis.labels) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = expansion(mult = c(0.05, 0.04))) +
  theme_bw(base_size = 14) +
  theme(
    aspect.ratio = 1,
    
    panel.grid = element_blank(),
    
    strip.background = element_rect(fill = "grey97", color = NA),
    strip.text = element_text(face = "plain", size = 12),
    strip.placement = "outside",
    
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_blank(),
    
    legend.position = "none"
  )
p2 <- p1 %+% plot.df2
p2 <- p2 +
  geom_hline(data = median.lines2, aes(yintercept = median_val), linetype = "dashed", color = "#DC0000B2", linewidth = 0.4, alpha = 0.8) +
  ggpubr::stat_pvalue_manual(
    stat.df |> filter(ctrl.cells == "normal reference cells", !is.na(p.signif)),
    label = "p.signif", xmin = "group1", xmax = "group2", y.position = "y.position", tip.length = 0.01, size = 3
  )

pdf("./Proj_CTAs/output/figure.panels/continuous.performance.comparison.pdf", width = 10.5, height = 10)
p2
dev.off()

##- Fig. 4f --------------------------------------------------------------------
registerDoParallel(cores = 50)
get.malignant.cells.based.on.posterior.p <- function(malignant_score, normal_score, p_threshold = 0.9) {
  library(mclust)
  
  set.seed(1234)
  mix.mal      <- Mclust(malignant_score, G = 2, verbose = FALSE)
  mal.comp.idx <- which.max(mix.mal$parameters$mean)
  mal.p        <- mix.mal$z[, mal.comp.idx]
  cells.A      <- mal.p > p_threshold
  
  set.seed(1234)
  mix.nor      <- Mclust(normal_score, G = 2, verbose = FALSE)
  mal.comp.idx <- which.min(mix.nor$parameters$mean)
  mal.p        <- mix.nor$z[, mal.comp.idx]
  cells.B      <- mal.p > p_threshold
  
  return(cells.A & cells.B)
}

## metrics used to predict malignant cells
load("./Proj_CTAs/output/manuscript.revision/CTG.metrics.ls.RData")
load("./Proj_CTAs/output/6. n_CTGs/scRNA-seq/malignant.cell.prediction.results.RData")  ## contained "SCEVAN.pred.ls" and "copykat.pred.ls"

## Precision, Recall, Specificity, and F1.score for each metric in each tumor sample
binary.performance.df <- foreach(s = names(CTG.metrics.ls), .combine = 'rbind') %dopar% {
  metrics.df    <- CTG.metrics.ls[[s]] |> 
    mutate(MaligCTG.label = get.malignant.cells.based.on.posterior.p(MaligCTG.malignant.score, MaligCTG.normal.score, p_threshold = 0.9))
  SCEVAN.label  <- if (s %in% names(SCEVAN.pred.ls)) {
    SCEVAN.pred.ls[[s]] |> select(barcode, SCEVAN.label = pred)
  } else {
    data.frame(barcode = metrics.df$barcode, SCEVAN.label = NA_real_)
  }
  copykat.label <- if (s %in% names(copykat.pred.ls)) {
    copykat.pred.ls[[s]] |> select(barcode, copykat.label = pred)
  } else {
    data.frame(barcode = metrics.df$barcode, copykat.label = NA_real_)
  }

  metrics.df <- list(metrics.df, SCEVAN.label, copykat.label) |>
    lapply(\(x) {x |> mutate(barcode = as.character(barcode))}) |>
    reduce(inner_join, by = "barcode")

  vars  <- c("MaligCTG.label", "SCEVAN.label", "copykat.label")
  df1   <- metrics.df
  truth <- ifelse(df1$malignant == "yes", 1, 0)
  if (min(sum(truth == 0), sum(truth == 1)) >= 10) {

    tmp <- map_dfr(vars, \(x) {
      pred <- df1[x]
      if(all(is.na(pred))) return(data.frame(Precision = NA_real_, Recall = NA_real_, Specificity = NA_real_, F1.score = NA_real_))
    
      TP <- sum(pred == 1 & truth == 1)
      FP <- sum(pred == 1 & truth == 0)
      TN <- sum(pred == 0 & truth == 0)
      FN <- sum(pred == 0 & truth == 1)
    
      Precision   <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
      Recall      <- TP / (TP + FN)
      Specificity <- TN / (TN + FP)
      F1.score    <- 2 * (Precision * Recall) / (Precision + Recall)

      data.frame(Precision, Recall, Specificity, F1.score)
    })
    
    Precision   <- tmp[, 1]
    Recall      <- tmp[, 2]
    Specificity <- tmp[, 3]
    F1.score    <- tmp[, 4]

  } else {
    Precision   <- NA_real_
    Recall      <- NA_real_
    Specificity <- NA_real_
    F1.score    <- NA_real_
  }
  res1 <- data.frame(dataset_sample = s, var = vars, Precision, Recall, Specificity, F1.score, ctrl.cells = "non-malignant cells", row.names = NULL)

  df2   <- metrics.df |> filter(normal.ref.cell %in% c("tumour", "normal"))
  truth <- ifelse(df2$normal.ref.cell == "tumour", 1, 0)
  if (min(sum(truth == 0), sum(truth == 1)) >= 10) {

    tmp <- map_dfr(vars, \(x) {
      pred <- df2[x]
      if(all(is.na(pred))) return(data.frame(Precision = NA_real_, Recall = NA_real_, Specificity = NA_real_, F1.score = NA_real_))
    
      TP <- sum(pred == 1 & truth == 1)
      FP <- sum(pred == 1 & truth == 0)
      TN <- sum(pred == 0 & truth == 0)
      FN <- sum(pred == 0 & truth == 1)
    
      Precision   <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
      Recall      <- TP / (TP + FN)
      Specificity <- TN / (TN + FP)
      F1.score    <- 2 * (Precision * Recall) / (Precision + Recall)

      data.frame(Precision, Recall, Specificity, F1.score)
    })
    
    Precision   <- tmp[, 1]
    Recall      <- tmp[, 2]
    Specificity <- tmp[, 3]
    F1.score    <- tmp[, 4]

  } else {
    Precision   <- NA_real_
    Recall      <- NA_real_
    Specificity <- NA_real_
    F1.score    <- NA_real_
  }
  res2 <- data.frame(dataset_sample = s, var = vars, Precision, Recall, Specificity, F1.score, ctrl.cells = "normal reference cells", row.names = NULL)

  rbind(res1, res2) |>
    mutate(
      n.mal        = sum(df1$malignant == "yes"),
      n.non.mal    = sum(df1$malignant == "no"),
      n.normal.ref = sum(df2$normal.ref.cell == "normal")
    )
}

stopImplicitCluster()
load("./Proj_CTAs/output/manuscript.revision/tumor.sample.info.updated.RData")
method.colors <- c(
  "MaligCTG"  = "#DC0000B2",
  "CopyKAT"   = "#00A087B2",
  "SCEVAN"    = "#4DBBD5B2"
)

plot.df1 <- binary.performance.df |>
  inner_join(tumor.sample.info, by = "dataset_sample") |>
  filter(ctrl.cells == "non-malignant cells", n.non.mal >= 10) |>

  mutate(method = case_match(
    var,
    "MaligCTG.label" ~ "MaligCTG",
    "SCEVAN.label"   ~ "SCEVAN",
    "copykat.label"  ~ "CopyKAT"
  )) |>
  mutate(method = factor(method, levels = c("MaligCTG", "CopyKAT", "SCEVAN"))) |>
  
  pivot_longer(cols = c("Precision", "Recall"), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric, levels = c("Precision", "Recall"), labels = c("Precision", "Recall"))) |>

  mutate(dataset_patient = paste0(data_folder, "_", patient)) |>
  group_by(dataset_group, metric, method, dataset_patient) |>
  summarise(value = median(value, na.rm = TRUE), .groups = "drop")
plot.df2 <- binary.performance.df |>
  inner_join(tumor.sample.info, by = "dataset_sample") |>
  filter(ctrl.cells == "normal reference cells", n.normal.ref >= 10) |>

  mutate(method = case_match(
    var,
    "MaligCTG.label" ~ "MaligCTG",
    "SCEVAN.label"   ~ "SCEVAN",
    "copykat.label"  ~ "CopyKAT"
  )) |>
  mutate(method = factor(method, levels = c("MaligCTG", "CopyKAT", "SCEVAN"))) |>
  
  pivot_longer(cols = c("Precision", "Recall"), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric, levels = c("Precision", "Recall"), labels = c("Precision", "Recall"))) |>

  mutate(dataset_patient = paste0(data_folder, "_", patient)) |>
  group_by(dataset_group, metric, method, dataset_patient) |>
  summarise(value = median(value, na.rm = TRUE), .groups = "drop")
compa.methods <- c("CopyKAT", "SCEVAN")
bin.df        <- bind_rows("non-malignant cells" = plot.df1, "normal reference cells" = plot.df2, .id = "ctrl.cells")

stat.df <- bin.df |>
  filter(method == "MaligCTG") |>
  select(ctrl.cells, dataset_group, metric, dataset_patient, value1 = value) |>
  inner_join(
    bin.df |> filter(method %in% compa.methods) |> select(ctrl.cells, dataset_group, metric, dataset_patient, group2 = method, value2 = value),
    by = c("ctrl.cells", "dataset_group", "metric", "dataset_patient")
  ) |>
  filter(!is.na(value1), !is.na(value2)) |>
  group_by(ctrl.cells, dataset_group, metric, group2) |>
  summarise(n.pairs = n(), p = suppressWarnings(wilcox.test(value1, value2, paired = TRUE)$p.value), .groups = "drop") |>
  mutate(group1 = "MaligCTG", group2 = factor(group2, levels = compa.methods)) |>
  arrange(ctrl.cells, dataset_group, metric, group2) |>
  group_by(ctrl.cells, dataset_group, metric) |>
  mutate(p.adj = p.adjust(p, method = "BH"), y.position = 1.02 + 0.06 * (row_number() - 1)) |>
  ungroup() |>
  mutate(
    group2   = as.character(group2),
    p.signif = cut(p, breaks = c(-Inf, 1e-4, 1e-3, 0.01, 0.05, Inf), labels = c("****", "***", "**", "*", "ns"))
  )
median.lines2 <- plot.df2 |> filter(method == "MaligCTG") |> group_by(dataset_group, metric) |>
  summarise(median_val = median(value, na.rm = TRUE), .groups = "drop")
p1 <- ggplot(plot.df1, aes(x = method, y = value)) +
  geom_boxplot(aes(color = method, fill = method), outlier.shape = NA, alpha = 0.2, width = 0.6, linewidth = 0.3) +
  geom_point(aes(color = method), size = 0.1, alpha = 0.4, position = position_jitter(width = 0.2, height = 0, seed = 1)) +
  facet_grid(metric ~ dataset_group, scales = "free_y", switch = "y") +
  scale_color_manual(values = method.colors) +
  scale_fill_manual(values = method.colors) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = expansion(mult = c(0.05, 0.04))) +
  theme_bw(base_size = 14) +
  theme(
    aspect.ratio = 1,
    
    panel.grid = element_blank(),
    
    strip.background = element_rect(fill = "grey97", color = NA),
    strip.text = element_text(face = "plain", size = 12),
    strip.placement = "outside",
    
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_blank(),
    
    legend.position = "none"
  )
p2 <- p1 %+% plot.df2
p2 <- p2 +
  geom_hline(data = median.lines2, aes(yintercept = median_val), linetype = "dashed", color = "#DC0000B2", linewidth = 0.4, alpha = 0.8) +
  ggpubr::stat_pvalue_manual(
    stat.df |> filter(ctrl.cells == "normal reference cells", !is.na(p.signif)),
    label = "p.signif", xmin = "group1", xmax = "group2", y.position = "y.position", tip.length = 0.01, size = 3
  )

pdf("./Proj_CTAs/output/figure.panels/binary.performance.comparison.pdf", width = 10.5, height = 10)
p2
dev.off()

##- Fig. 4g --------------------------------------------------------------------
load("./Proj_CTAs/output/6. n_CTGs/scRNA-seq/running.time.input.RData")
time.df <- readRDS("./Proj_CTAs/output/manuscript.revision/time.df.rds")
plot.df <- time.df |> group_by(dataset_sample, method) |> summarise(time = median(elapsed.sec), .groups = "drop") |>
  left_join(args.df |> select(dataset_sample, n_cell, n_cell_group), by = "dataset_sample") |>
  mutate(n_cell_group = factor(paste0("~", n_cell_group), levels = c("~1k", "~3k", "~5k", "~8k", "~10k")))
colors     <- c("MaligCTG" = "#DC0000B2", "infercnv" = "#3C5488B2")
malig.df   <- plot.df |> filter(method == "MaligCTG")
infer.df   <- plot.df |> filter(method == "infercnv") |> mutate(time = as.numeric(time) / 60)

malig.summ <- malig.df |> group_by(n_cell_group) |> summarise(mean = mean(time), .groups = "drop")
infer.summ <- infer.df |> group_by(n_cell_group) |> summarise(mean = mean(time), .groups = "drop")

p1 <- ggplot(malig.df, aes(x = n_cell_group, y = time)) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 1234),
    shape = 21, size = 1.7, fill = colors[["MaligCTG"]], color = "white", stroke = 0.20, alpha = 0.45
  ) +
  stat_summary(
    fun.data = mean_se, geom = "errorbar",
    width = 0.15, color = colors[["MaligCTG"]], linewidth = 0.40
  ) +
  geom_line(
    data = malig.summ, aes(y = mean, group = 1),
    color = colors[["MaligCTG"]], linewidth = 0.65
  ) +
  geom_point(
    data = malig.summ, aes(y = mean),
    shape = 21, size = 2.7, fill = colors[["MaligCTG"]], color = "white", stroke = 0.55
  ) +
  scale_x_discrete(drop = FALSE, expand = expansion(add = 0.40)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 5), expand = expansion(mult = c(0, 0.07))) +
  expand_limits(y = 0) +
  labs(x = "Number of cells", y = "Runtime (s)", title = "MaligCTG") +
  theme_classic(base_size = 11)

p2 <- ggplot(infer.df, aes(x = n_cell_group, y = time)) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 1234),
    shape = 21, size = 1.7, fill = colors[["infercnv"]], color = "white", stroke = 0.20, alpha = 0.45
  ) +
  stat_summary(
    fun.data = mean_se, geom = "errorbar",
    width = 0.15, color = colors[["infercnv"]], linewidth = 0.40
  ) +
  geom_line(
    data = infer.summ, aes(y = mean, group = 1),
    color = colors[["infercnv"]], linewidth = 0.65
  ) +
  geom_point(
    data = infer.summ, aes(y = mean),
    shape = 21, size = 2.7, fill = colors[["infercnv"]], color = "white", stroke = 0.55
  ) +
  scale_x_discrete(drop = FALSE, expand = expansion(add = 0.40)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 5), expand = expansion(mult = c(0, 0.07))) +
  expand_limits(y = 0) +
  labs(x = "Number of cells", y = "Runtime (min)", title = "inferCNV") +
  theme_classic(base_size = 11)

p <- (p1 | p2) &
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.title.x = element_text(size = 10.5, margin = margin(t = 7)),
    axis.title.y = element_text(size = 10.5, margin = margin(r = 7)),
    axis.text = element_text(size = 9, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.30),

    panel.grid.major.y = element_line(color = "grey94", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),

    plot.title = element_text(size = 11, face = "plain", hjust = 0.5, margin = margin(b = 7)),
    legend.position = "none",
    plot.margin = margin(t = 7, r = 8, b = 6, l = 6)
  )

pdf("./Proj_CTAs/output/figure.panels/runtime.pdf", width = 8, height = 4)
print(p)
dev.off()

##- Fig. 4h --------------------------------------------------------------------
library(Seurat)
prefix  <- "./Proj_CTAs/input/MaligCTG-validation/scONE_GSE185269"
counts  <- data.table::fread(file.path(prefix, "GSE185269_scONEseq_GBM_RNA_gene_counts_matrix.csv.gz")) |> column_to_rownames("V1") |> as.matrix()

cell.df <- read.csv(file.path(prefix, "GSE185269_scONEseq_GBM_sample_annotation.csv.gz")) |>
  transmute(
    barcode       = cellid,
    rna.cell.type = Neftel_type,
    dna.clone     = Clone_type,
    cell.type     = case_when(
      Clone_type %in% c("GBM clone 1", "GBM clone 2", "GBM clone 3") ~ "malignant",
      Clone_type %in% "Normal cell" & Neftel_type %in% "Astrocyte"   ~ "astrocytes",
      TRUE ~ NA_character_
    )
  ) |>
  remove_rownames() |>
  column_to_rownames("barcode")

kept <- intersect(rownames(cell.df), colnames(counts))
seu  <- CreateSeuratObject(counts = counts[, kept], meta.data = cell.df[kept, ])
seu  <- seu |> NormalizeData()

load("./Proj_CTAs/output/manuscript.revision/pan.sample.CTGs.RData")
CTGs   <- intersect(rownames(seu), pan.sample.CTGs)
n.CTGs <- Matrix::colSums(seu[["RNA"]]@data[CTGs, , drop = FALSE] > 0)

top       <- n.CTGs >= quantile(n.CTGs, 0.95)
tail      <- n.CTGs <= quantile(n.CTGs, 0.05)
top.expr  <- seu[["RNA"]]@counts[, top, drop = FALSE] |> Matrix::rowSums() |> (\(x) log2(1e6 * x / sum(x) + 1))()
tail.expr <- seu[["RNA"]]@counts[, tail, drop = FALSE] |> Matrix::rowSums() |> (\(x) log2(1e6 * x / sum(x) + 1))()
ref       <- cbind(top.expr, tail.expr)

pred      <- SingleR::SingleR(test = seu[["RNA"]]@data, ref = ref, labels = c("high", "low"))
seu$x     <- pred$scores[, "high"]

idx     <- !is.na(seu$cell.type)
truth   <- factor(seu$cell.type[idx], levels = c("astrocytes", "malignant"))
roc.fit <- pROC::roc(truth, seu$x[idx], direction = "<", quiet = TRUE)
auc.ci  <- pROC::ci.auc(roc.fit, method = "delong")

auc.df <- data.frame(
  n.malignant = sum(truth == "malignant"), n.normal = sum(truth == "astrocytes"),
  auc = as.numeric(roc.fit$auc), ci.low = as.numeric(auc.ci[1]), ci.high = as.numeric(auc.ci[3])
)
cell.colors <- c("astrocytes" = "#6B9AC4", "malignant" = "#D47264")
cell.labels <- c("astrocytes" = "Astrocytes", "malignant" = "Malignant cells")

seu.plot       <- seu[, idx] |> FindVariableFeatures() |> ScaleData() |> RunPCA() |> RunUMAP(dims = 1:30)
umap.df        <- Embeddings(seu.plot, reduction = "umap")[, 1:2, drop = FALSE] |> as.data.frame()
names(umap.df) <- c("umap.1", "umap.2")
umap.df        <- umap.df |> rownames_to_column("barcode") |>
  left_join(
    seu.plot[[]] |> rownames_to_column("barcode") |> dplyr::select(barcode, malignant.score = x, cell.type),
    by = "barcode"
  ) |>
  dplyr::mutate(cell.type = factor(cell.type, levels = names(cell.colors)))

## use identical square windows without changing the UMAP coordinates.
umap.span <- max(diff(range(umap.df$umap.1)), diff(range(umap.df$umap.2))) * 1.08
umap.xlim <- mean(range(umap.df$umap.1)) + c(-0.5, 0.5) * umap.span
umap.ylim <- mean(range(umap.df$umap.2)) + c(-0.5, 0.5) * umap.span

label <- sprintf(
  "AUROC = %.3f\n95%% CI [%.3f, %.3f]",
  auc.df$auc, auc.df$ci.low, auc.df$ci.high
)

p1 <- ggplot(umap.df, aes(x = umap.1, y = umap.2)) +
  geom_point(aes(color = malignant.score), size = 1.05) +
  scale_color_viridis_c(
    option = "D", name = "Malignant score",
    breaks = scales::breaks_pretty(n = 3), labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
  coord_equal(xlim = umap.xlim, ylim = umap.ylim, expand = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(
    color = guide_colorbar(
      direction = "vertical",
      theme = theme(
        legend.title.position = "left",
        legend.title = element_text(size = 9, angle = 90, hjust = 0.5, vjust = 0.5, margin = margin(r = 4)),
        legend.key.width = grid::unit(3, "mm"),
        legend.key.height = grid::unit(22, "mm")
      )
    )
  )

p2 <- ggplot(umap.df, aes(x = umap.1, y = umap.2)) +
  geom_point(aes(color = cell.type), size = 1.05) +
  scale_color_manual(
    values = cell.colors, breaks = names(cell.colors), labels = cell.labels,
    name = "", drop = FALSE
  ) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
  coord_equal(xlim = umap.xlim, ylim = umap.ylim, expand = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2.5, alpha = 1)))

## ROC coordinates from the same pROC fit that gives the AUROC label
roc.df <- data.frame(specificity = roc.fit$specificities, sensitivity = roc.fit$sensitivities)

p3 <- ggplot(roc.df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_abline(intercept = 0, slope = 1, color = "grey70", linetype = "dashed", linewidth = 0.35) +
  geom_path(color = "#D47264", linewidth = 0.80) +
  annotate(
    "text", x = 0.96, y = 0.10, label = label,
    hjust = 1, vjust = 0, size = 3.2, lineheight = 1.25, color = "black"
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(x = "1 - specificity", y = "Sensitivity")

p <- patchwork::wrap_plots(p1, p2, p3, nrow = 1, widths = c(1, 1, 1), guides = "keep") &
  theme_classic(base_size = 11) &
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.text = element_text(size = 8.5, color = "black"),
    axis.title.x = element_text(size = 10.5, margin = margin(t = 7)),
    axis.title.y = element_text(size = 10.5, margin = margin(r = 7)),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.30),
    
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    
    legend.position = "inside",
    legend.position.inside = c(0.97, 0.97),
    legend.justification.inside = c(1, 1),
    legend.title = element_text(size = 9, face = "plain"),
    legend.text = element_text(size = 8.5),
    legend.key = element_blank(),
    legend.key.size = grid::unit(3.5, "mm"),
    legend.background = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

pdf("./Proj_CTAs/output/figure.panels/Fig. 4g.pdf", width = 9.6, height = 3.5)
print(p)
dev.off()

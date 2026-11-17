library(foreach)
library(doParallel)
library(tidyverse)
library(magrittr)
library(Seurat)
library(patchwork)

##- Fig. 6a and 6d -------------------------------------------------------------
# Fig. 6a
source("./Proj_CTAs/code/visualization.R")
load("./Proj_CTAs/output/manuscript.revision/tumor.sample.info.updated.RData")
samples.ls <- tumor.sample.info |> 
  mutate(dataset_patient = paste0(data_folder, "_", patient)) |>
  filter(dataset_group == "Discovery Cohort") |>
  left_join(cancer.abbr, by = "cancer") |>
  mutate(abbr = factor(abbr, levels = cancer.abbr$abbr)) %$% 
  split(., abbr)

load("./Proj_CTAs/output/manuscript.revision/CTG.DE.auc.sc.RData")
med.auc.df.sc <- foreach(c = names(samples.ls), .combine = "rbind") %do% {
  sample.patient.map <- samples.ls[[c]] |> select(dataset_patient, dataset_sample) |> distinct()
  sample.patient.ls  <- split(sample.patient.map$dataset_sample, sample.patient.map$dataset_patient)
  samples            <- sample.patient.map$dataset_sample

  TME.samples        <- intersect(colnames(auc.mat.TME), samples)
  normal.ref.samples <- intersect(colnames(auc.mat.normal.ref), samples)

  auc.TME     <- sample.patient.ls |> lapply(intersect, TME.samples) %>% Filter(\(x) length(x) > 0, .) |> 
    sapply(\ (s) { apply(auc.mat.TME[, s, drop = FALSE], 1, median, na.rm = TRUE) })
  med.auc.TME <- apply(auc.TME, 1, \(x) {
    if ( sum(!is.na(x)) >= 3 ) median(x, na.rm = TRUE) else NA_real_
  })
  
  auc.normal.ref <- sample.patient.ls |> lapply(intersect, normal.ref.samples) %>% Filter(\(x) length(x) > 0, .) |> 
    sapply(\ (s) { apply(auc.mat.normal.ref[, s, drop = FALSE], 1, median, na.rm = TRUE) })
  if (length(auc.normal.ref) > 0) {
    med.auc.normal.ref <- apply(auc.normal.ref, 1, \(x) {
      if ( sum(!is.na(x)) >= 3 ) median(x, na.rm = TRUE) else NA_real_
    })
    lack.normal.ref    <- all(is.na(med.auc.normal.ref))
  } else {
    med.auc.normal.ref <- NA_real_
    lack.normal.ref    <- TRUE
  }

  data.frame(cancer = c, med.auc.TME, med.auc.normal.ref, lack.normal.ref) |> rownames_to_column("symbol")
}

load("./Proj_CTAs/output/manuscript.revision/CTG.annotation.RData")
load("./Proj_CTAs/output/4. characteristics of CTGs/cell.cycle.CTGs.RData")
CTG1    <- CTGs.annotated |> filter(tier1 == "Meets", tier2 == "Meets") |> pull(symbol)
CTG2    <- cc.CTGs

CTG3.df <- med.auc.df.sc |> 
  filter(med.auc.TME > 0.55 & (med.auc.normal.ref > 0.55 | lack.normal.ref)) |>
  filter(symbol %in% CTG1, !symbol %in% CTG2)

gene.type.df <- tibble::tribble(
  ~symbol,     ~approved.symbol, ~gene.type,
  "HMX1",      "HMX1",          "protein_coding",
  "CT83",      "CT83",          "protein_coding",
  "ROPN1",     "ROPN1",         "protein_coding",
  "DCAF4L2",   "DCAF4L2",       "protein_coding",
  "GAGE1",     "GAGE1",         "protein_coding",
  "MAGEA3",    "MAGEA3",        "protein_coding",
  "MAGEA6",    "MAGEA6",        "protein_coding",
  "DSCR8",     "DSCR8",         "lncRNA",
  "MAGEA4",    "MAGEA4",        "protein_coding",
  "MDGA2",     "MDGA2",         "protein_coding",
  "PRH1",      "PRH1",          "protein_coding",
  "C2CD6",     "CATSPERT",      "protein_coding",
  "LINC01194", "LINC01194",     "lncRNA"
)

sc.df       <- med.auc.df.sc |> filter(symbol %in% CTG1, !lack.normal.ref, med.auc.TME > 0.55, med.auc.normal.ref > 0.55)
main.pre.df <- sc.df |> filter(!symbol %in% CTG2)
supp.pre.df <- CTG3.df |> filter(lack.normal.ref)

audit.df    <- bind_rows(main.pre.df, supp.pre.df) |> distinct(symbol) |> left_join(gene.type.df, by = "symbol")
stopifnot(!anyNA(audit.df$gene.type))

main.df     <- main.pre.df |> left_join(gene.type.df, by = "symbol") |> filter(gene.type == "protein_coding")
supp.df     <- supp.pre.df |> left_join(gene.type.df, by = "symbol") |> filter(gene.type == "protein_coding")
excluded.df <- audit.df |> filter(gene.type != "protein_coding")

load("./Proj_CTAs/output/manuscript.revision/pan.sample.CTGs.RData")
n.catalog    <- n_distinct(pan.sample.CTGs)
n.annotation <- n_distinct(CTG1)
n.sc         <- n_distinct(sc.df$symbol)
n.cc         <- n.sc - n_distinct(main.pre.df$symbol)
n.pre        <- n_distinct(main.pre.df$symbol)
n.main       <- n_distinct(main.df$symbol)
n.noncoding  <- n.pre - n.main

gene.df <- main.df |> mutate(cancer = factor(cancer, levels = names(samples.ls))) |>
  group_by(cancer) |>
  summarise(genes = paste(unique(approved.symbol), collapse = ", "), .groups = "drop") |>
  mutate(y = 1.42 - (row_number() - (n() + 1) / 2) * 0.43)

box.df <- data.frame(
  xmin = 0, xmax = 6.2,
  ymin = c(8.30, 6.05, 4.20, 0.35), ymax = c(9.30, 7.45, 5.20, 3.35),
  title = c(
    "CTG catalog", "CTGs meeting both\nannotation criteria",
    "Malignancy-specific CTGs", "Candidate TCR-T targets"
  ),
  y.title = c(8.97, 7.00, 4.87, 3.03),
  n.label = paste0("n = ", c(n.catalog, n.annotation, n.sc, n.main)),
  y.n     = c(8.63, 6.45, 4.53, 2.69),
  fill    = c(rep("grey98", 3), "white"),
  border  = c(rep("grey60", 3), "grey25"),
  lwd     = c(rep(0.40, 3), 0.55)
)

arrow.df <- data.frame(x = 3.1, xend = 3.1, y = c(8.30, 6.05, 4.20), yend = c(7.45, 5.20, 3.35))
step.df  <- data.frame(
  x = 6.55, y = c(7.875, 5.625, 3.775),
  label = c(
    "Minimal non-testicular expression and\ntumor enrichment (bulk RNA-seq)",
    "Median AUROC > 0.55 vs. normal reference\nand TME cells (scRNA-seq)",
    paste0(
      "CTGs enriched in cycling malignant cells excluded (n = ", n.cc, ")\n",
      "Noncoding genes excluded (n = ", n.noncoding, ")"
    )
  )
)
branch.df <- data.frame(x = c(3.1, -2.35), xend = c(-2.35, -2.35), y = 5.625, yend = c(5.625, 5.25))

p <- ggplot() +
  geom_rect(
    data = box.df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = box.df$fill, color = box.df$border, linewidth = box.df$lwd
  ) +
  annotate("rect", xmin = -4.30, xmax = -0.40, ymin = 4.10, ymax = 5.25,
           fill = "white", color = "grey65", linewidth = 0.35, linetype = "dashed") +
  geom_segment(
    data = arrow.df, aes(x = x, xend = xend, y = y, yend = yend),
    color = "grey45", linewidth = 0.40, arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed")
  ) +
  geom_segment(
    data = branch.df[1, ], aes(x = x, xend = xend, y = y, yend = yend),
    color = "grey70", linewidth = 0.35, linetype = "dashed"
  ) +
  geom_segment(
    data = branch.df[2, ], aes(x = x, xend = xend, y = y, yend = yend),
    color = "grey70", linewidth = 0.35, linetype = "dashed",
    arrow = grid::arrow(length = grid::unit(1.4, "mm"), type = "closed")
  ) +
  geom_text(
    data = box.df, aes(x = 3.1, y = y.title, label = title),
    size = 3.0, color = "grey10", family = "sans", lineheight = 1.05
  ) +
  geom_text(data = box.df, aes(x = 3.1, y = y.n, label = n.label), size = 2.8, color = "grey35", family = "sans") +
  geom_text(
    data = step.df, aes(x = x, y = y, label = label),
    size = 2.6, color = "grey35", family = "sans", hjust = 0, lineheight = 1.15
  ) +
  annotate("text", x = -2.35, y = 4.675, size = 2.6, color = "grey35", family = "sans", lineheight = 1.15,
           label = "No normal reference cells:\nHGSOC, NB and SKCM assessed\nagainst TME cells only (Table S8)") +
  annotate("segment", x = 0.35, xend = 5.85, y = 2.38, yend = 2.38, color = "grey88", linewidth = 0.30) +
  geom_point(
    data = gene.df, aes(x = 0.40, y = y, fill = cancer),
    shape = 21, size = 2.1, color = "grey35", stroke = 0.25, show.legend = FALSE
  ) +
  geom_text(data = gene.df, aes(x = 0.65, y = y, label = cancer), size = 2.7, color = "grey30", family = "sans", hjust = 0) +
  geom_text(
    data = gene.df, aes(x = 1.65, y = y, label = genes),
    size = 2.7, color = "grey10", family = "sans", fontface = "italic", hjust = 0
  ) +
  scale_fill_manual(values = cancer.colours, guide = "none") +
  coord_cartesian(xlim = c(-4.70, 11.20), ylim = c(0.15, 9.45), expand = FALSE) +
  theme_void(base_size = 10) +
  theme(text = element_text(family = "sans", color = "black"), plot.margin = margin(t = 4, r = 5, b = 4, l = 5))

pdf("./Proj_CTAs/output/figure.panels/Fig. 6a.pdf", width = 8.7, height = 4.8)
print(p)
dev.off()

# Fig. 6d
HCC.auc.df        <- med.auc.df.sc |> filter(cancer == "HCC")

load("./Proj_CTAs/output/manuscript.revision/prop.mat.mal.RData")
sample.patient.ls <- samples.ls[["HCC"]] |> select(dataset_patient, dataset_sample) |> distinct() %$% split(dataset_sample, dataset_patient)
HCC.prop.mat      <- prop.mat.mal[HCC.auc.df$symbol, unlist(sample.patient.ls)]
HCC.prop.df       <- sapply(sample.patient.ls, \(s) {
  HCC.prop.mat[, s, drop = FALSE] |> apply(1, median, na.rm = TRUE)
}) |> apply(1, median, na.rm = TRUE) %>%
  data.frame(symbol = names(.), pct_median = .)

plot.df  <- inner_join(HCC.prop.df, HCC.auc.df, by = "symbol") |> mutate(is_target = med.auc.TME >= 0.55 & med.auc.normal.ref >= 0.55)
label.df <- plot.df |> filter(symbol == "DCAF4L2")

p <- ggplot(plot.df, aes(x = med.auc.normal.ref, y = med.auc.TME)) +
  annotate("rect", xmin = 0.55, xmax = Inf, ymin = 0.55, ymax = Inf, fill = "#D47264", alpha = 0.08) +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "grey80", linewidth = 0.25) +
  geom_vline(xintercept = 0.50, linetype = "dashed", color = "grey80", linewidth = 0.25) +
  geom_hline(yintercept = 0.55, linetype = "dashed", color = "grey50", linewidth = 0.35) +
  geom_vline(xintercept = 0.55, linetype = "dashed", color = "grey50", linewidth = 0.35) +
  geom_point(
    data = filter(plot.df, symbol != "DCAF4L2"), aes(size = pct_median, fill = is_target),
    shape = 21, color = "grey35", stroke = 0.25, alpha = 0.80, na.rm = TRUE
  ) +
  geom_point(
    data = label.df, aes(size = pct_median, fill = is_target),
    shape = 21, color = "grey15", stroke = 0.65, alpha = 1, show.legend = FALSE, na.rm = TRUE
  ) +
  ggrepel::geom_text_repel(
    data = label.df, aes(label = symbol), family = "sans", fontface = "italic", size = 3.5, seed = 1234,
    nudge_x = 0.009, nudge_y = 0.015, box.padding = 0.45, point.padding = 0.40,
    min.segment.length = 0, segment.color = "grey45", segment.size = 0.30, show.legend = FALSE
  ) +
  scale_fill_manual(values = c("TRUE" = "#D47264", "FALSE" = "grey75"), guide = "none") +
  scale_size_continuous(
    name = "Median detection (%)", range = c(1, 5),
    breaks = scales::breaks_pretty(n = 4), labels = scales::label_number(scale = 100, accuracy = 1)
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.05), labels = scales::label_number(accuracy = 0.01)) +
  scale_y_continuous(breaks = seq(0, 1, 0.05), labels = scales::label_number(accuracy = 0.01)) +
  labs(x = "AUROC (malignant vs normal reference)", y = "AUROC (malignant vs TME cells)") +
  guides(size = guide_legend(override.aes = list(fill = "grey75", color = "grey35", alpha = 1, stroke = 0.25))) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "sans", color = "black"),
    aspect.ratio = 1,
    axis.text = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10.5, margin = margin(t = 7)),
    axis.title.y = element_text(size = 10.5, margin = margin(r = 7)),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.30),

    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.40),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title.position = "left",
    legend.title = element_text(size = 9, face = "plain", angle = 90, hjust = 0.5, vjust = 0.5, margin = margin(r = 6)),
    legend.text = element_text(size = 8.5),
    legend.key.height = grid::unit(5, "mm"),
    plot.margin = margin(t = 7, r = 9, b = 7, l = 7)
  )

pdf("./Proj_CTAs/output/figure.panels/DCAF4L2.AUROCs.pdf", width = 6.5, height = 6.5)
print(p)
dev.off()

##- Fig. 6b --------------------------------------------------------------------
load("./Proj_CTAs/output/manuscript.revision/tumor.sample.info.updated.RData")
TNBC.sample.info <- tumor.sample.info |> 
  mutate(dataset_patient = paste0(data_folder, "_", patient)) |>
  filter(dataset_group == "Discovery Cohort", cancer == "Breast Cancer (TNBC)")
TNBC.datasets    <- unique(TNBC.sample.info$dataset)

registerDoParallel(cores = length(TNBC.datasets))
combined.df <- foreach(d = TNBC.datasets, .combine = 'rbind') %dopar% {
  prefix <- "./Proj_CTAs/input/scRNA-seq datasets/normalized Seurat object/"
  file   <- paste0(d, ".rds")
  seu    <- readRDS(paste0(prefix, file))
  seu    <- subset(seu, malignant == "yes")
  
  meta <- data.frame(dataset = d, sample = seu$sample) |> rownames_to_column("barcode") |> mutate(dataset_sample = paste(dataset, sample, sep = "_"))
  expr <- seu[["RNA"]]@data
  expr <- data.frame(ROPN1  = expr["ROPN1", ] > 0, 
                     ROPN1B = expr["ROPN1B", ] > 0, 
                     CT83   = expr["CT83", ] > 0)
  
  cbind(meta, expr)
}
stopImplicitCluster()

## malignant cells of all samples from one patient are pooled, so every patient contributes one bar
patient.map <- TNBC.sample.info |> distinct(dataset, sample, dataset_patient)
percent.df  <- combined.df |>
  left_join(patient.map, by = c("dataset", "sample")) |>
  group_by(dataset_patient) |>
  group_modify(~ {
    x <- .x$ROPN1
    y <- .x$CT83

    data.frame(
      percent = c(mean(x & y), mean(x & !y), mean(!x & y), mean(!x & !y)),
      group   = c("ROPN1+CT83+", "ROPN1+CT83-", "ROPN1-CT83+", "ROPN1-CT83-")
    )
  }) |>
  ungroup()

## patients ordered by the proportion of malignant cells positive for either gene (1 - double negative)
patient.orders <- percent.df |> filter(group == "ROPN1-CT83-") |> arrange(percent) |> pull(dataset_patient)
n.patients     <- length(patient.orders)

plot.df <- percent.df |>
  mutate(
    dataset_patient = factor(dataset_patient, levels = patient.orders),
    group = factor(group, levels = c("ROPN1-CT83-", "ROPN1-CT83+", "ROPN1+CT83+", "ROPN1+CT83-"))
  )

my.colors <- c(
  "ROPN1-CT83-" = "#E5E5E5",
  "ROPN1-CT83+" = "#D53E4F",
  "ROPN1+CT83+" = "#9970AB",
  "ROPN1+CT83-" = "#4575B4"
)
my.labels <- c(
  "ROPN1-CT83+" = expression(italic("CT83")^"+"~italic("ROPN1")^"-"),
  "ROPN1+CT83+" = expression(italic("CT83")^"+"~italic("ROPN1")^"+"),
  "ROPN1+CT83-" = expression(italic("CT83")^"-"~italic("ROPN1")^"+"),
  "ROPN1-CT83-" = expression(italic("CT83")^"-"~italic("ROPN1")^"-")
)
group.breaks <- c("ROPN1-CT83+", "ROPN1+CT83+", "ROPN1+CT83-", "ROPN1-CT83-")

p <- ggplot(plot.df, aes(x = dataset_patient, y = percent, fill = group)) +
  geom_col(width = 0.85, color = "grey45", linewidth = 0.15, alpha = 0.90) +
  scale_fill_manual(
    values = my.colors, labels = my.labels, breaks = group.breaks,
    name = NULL, drop = FALSE
  ) +
  scale_x_discrete(limits = patient.orders, drop = FALSE, expand = expansion(add = 0.85)) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25), labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = bquote("TNBC patients (" * italic(n) == .(n.patients) * ")"),
    y = "Proportion of malignant cells"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(alpha = 1))) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 9.5, color = "black"),
    axis.title.x = element_text(size = 10.5, margin = margin(t = 7)),
    axis.title.y = element_text(size = 10.5, margin = margin(r = 7)),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.30),

    panel.grid = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 9.5),
    legend.key.width = grid::unit(3.5, "mm"),
    legend.key.height = grid::unit(3.5, "mm"),
    legend.margin = margin(t = 0, r = 0, b = 3, l = 0),

    plot.margin = margin(t = 6, r = 8, b = 6, l = 6)
  )

pdf("./Proj_CTAs/output/figure.panels/Fig. 6B.pdf", width = 10, height = 5)
print(p)
dev.off()

##- Fig. 6c --------------------------------------------------------------------
load("./Proj_CTAs/output/manuscript.revision/ST.data.RData")

i          <- "HTAN_HT206B1-U1_ST_Bn1"
seu        <- readRDS(ST.data$processed_path[ST.data$unique_sample_id == i])
spatial.theme <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
)
seu.sub               <- subset(seu, idents = c("1", "3", "16"))
DefaultAssay(seu.sub) <- "Spatial"
seu.sub               <- NormalizeData(seu.sub)  ## LogNormalize is more suitable for visualizing CTG expression than SCTransform

cluster.colors <- c("1" = "#E69F00", "3" = "#56B4E9", "16" = "#009E73") 
CT83.colors    <- scale_fill_gradientn(
  colors = c("#F0F0F0", "#FEE08B", "#FDAE61", "#F46D43", "#D73027", "#A50026"),
  values = c(0, 0.1, 0.3, 0.6, 0.8, 1),
  name = ""
)
ROPN1.colors   <- scale_fill_gradientn(
  colors = c("#F0F0F0", "#ABD9E9", "#74ADD1", "#4575B4", "#313695", "#023858"), 
  values = c(0, 0.1, 0.3, 0.6, 0.8, 1),
  name = ""
)

p1 <- SpatialDimPlot(seu.sub, cols = cluster.colors) + 
  spatial.theme +
  labs(title = "Malignant Clusters", fill = "Cluster") +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 14))
p2 <- SpatialFeaturePlot(seu.sub, features = "ROPN1", stroke = 0) + 
  spatial.theme + ROPN1.colors +
  labs(title = expression(italic("ROPN1")~"Expression")) +
  theme(legend.position = "right")
p3 <- SpatialFeaturePlot(seu.sub, features = "CT83", stroke = 0) + 
  spatial.theme + CT83.colors +
  labs(title = expression(italic("CT83")~"Expression")) +
  theme(legend.position = "right")
p4 <- VlnPlot(seu.sub, features = c("ROPN1", "CT83"), pt.size = 0, ncol = 1, cols = cluster.colors) & 
  theme_classic(base_size = 12) &
  theme(
    axis.text.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 13),
    legend.position = "none"
  )

pdf("./Proj_CTAs/output/figure.panels/Fig. 6C.pdf", width = 15, height = 4.5)
p1 + p2 + p3 + p4 + plot_layout(widths = c(1, 1, 1, 0.8))
dev.off()

##- Fig. 6e --------------------------------------------------------------------
out            <- "./TCR.T.targets/"
mhcflurry.pred <- data.table::fread(paste0(out, "mhcflurry_predictions_unrolled.csv"))
hit.df         <- mhcflurry.pred |> filter(affinity < 500, affinity_percentile < 2)
plot.df        <- hit.df |> 
  filter(sequence_name == "DCAF4L2") |> 
  count(best_allele, name = "peptide_count") |> 
  ungroup() |> 
  mutate(
    best_allele = factor(best_allele), 
    HLA_group   = ifelse(grepl("^HLA-A", best_allele), "HLA-A", "HLA-B") 
  )
p <- ggplot(plot.df, aes(x = best_allele, y = peptide_count, fill = HLA_group)) +
  
  geom_bar(stat = "identity", width = 0.75, color = "grey30", linewidth = 0.3, alpha = 0.9) + 
  scale_fill_manual(values = c("HLA-A" = "#D47264", "HLA-B" = "#6B9AC4")) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  
  labs(x = "", y = "Number of High-Affinity Peptides", fill = "") +
  
  theme_classic(base_size = 12) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(color = "black"),
    
    legend.position = c(0.1, 0.9),
    legend.background = element_blank(),
    legend.key.size = unit(0.5, "cm")
  )

pdf("./Fig. 6F.pdf", width = 8, height = 4.5)
p
dev.off()

##- Fig. 6f --------------------------------------------------------------------
tabula.expr.df <- read.csv("./practice/cellxgene-census/output/Tabula_Sapiens_CT83_DCAF4L2.csv.gz")

gene.order     <- c("CT83", "DCAF4L2")
group.order    <- c("Testis-only cell types", "Other cell types")
group.colors   <- c("Testis-only cell types" = "#D47264", "Other cell types" = "#6B9AC4")
n.types        <- n_distinct(tabula.expr.df$cell_type)

plot.df <- tabula.expr.df |>
  transmute(
    cell_type, cell_type_ontology_term_id, tissue_general,
    CT83.umi    = ENSG00000204019, CT83.cp10k    = ENSG00000204019_CP10K,
    DCAF4L2.umi = ENSG00000176566, DCAF4L2.cp10k = ENSG00000176566_CP10K
  ) |>
  pivot_longer(
    cols = c(CT83.umi, CT83.cp10k, DCAF4L2.umi, DCAF4L2.cp10k),
    names_to = c("gene", ".value"), names_pattern = "^(CT83|DCAF4L2)\\.(umi|cp10k)$"
  ) |>
  group_by(gene, cell_type, cell_type_ontology_term_id) |>
  summarise(
    n.cells = n(), n.positive = sum(umi > 0),
    max.umi = max(umi), mean.cp10k = mean(cp10k),
    n.testis = sum(tissue_general == "testis"),
    .groups = "drop"
  ) |>
  mutate(
    pct.positive = 100 * n.positive / n.cells,
    tissue.group = if.else(n.testis == n.cells, "Testis-only cell types", "Other cell types"),
    gene         = factor(gene, levels = gene.order),
    tissue.group = factor(tissue.group, levels = group.order)
  )

label.df     <- plot.df |> filter(max.umi > 1) |> mutate(label = if.else(max.umi > 1, cell_type, ""))
group.labels <- c(
  "Testis-only cell types" = "Testis-restricted",
  "Other cell types" = "Other cell types"
)

ann.df      <- plot.df |> group_by(gene) |>
  summarise(n.detected = sum(n.positive > 0), .groups = "drop") |>
  dplyr::mutate(label = sprintf("Detected in %d/%d cell types", n.detected, n.types))

size.limits <- c(0, max(plot.df$mean.cp10k))
x.limits    <- c(0, max(plot.df$max.umi) * 1.12)

p <- ggplot(plot.df, aes(x = max.umi, y = pct.positive)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey65", linewidth = 0.30) +
  geom_point(
    aes(size = mean.cp10k, fill = tissue.group),
    shape = 21, color = "grey35", stroke = 0.25, alpha = 0.85
  ) +
  ggrepel::geom_text_repel(
    data = label.df, aes(label = label),
    size = 2.8, color = "grey20", seed = 1234, lineheight = 0.95,
    box.padding = 0.45, point.padding = 0.30,
    min.segment.length = 0, segment.color = "grey60", segment.size = 0.25,
    max.overlaps = Inf, max.time = 1, show.legend = FALSE
  ) +
  geom_text(
    data = ann.df, aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = 1.03, vjust = 1.20, size = 2.9, color = "grey30"
  ) +
  scale_fill_manual(
    values = group.colors, limits = group.order, labels = group.labels,
    name = NULL, drop = FALSE
  ) +
  scale_size_continuous(
    name = "Mean expression (CP10K)", limits = size.limits, range = c(1.8, 6.5),
    breaks = c(0, 0.5, 1, 1.5), labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_x_continuous(
    trans = scales::log1p_trans(), breaks = c(0, 1, 2, 4, 8, 16, 32, 64),
    limits = x.limits, expand = expansion(mult = c(0.04, 0.04)), oob = scales::oob_keep
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.001, base = 10),
    breaks = c(0, 0.001, 0.01, 0.1, 1, 10, 100),
    labels = c("0", "0.001", "0.01", "0.1", "1", "10", "100"),
    limits = c(0, 100), expand = expansion(mult = c(0.03, 0.08)), oob = scales::oob_keep
  ) +
  facet_wrap(~ gene, nrow = 1) +
  labs(x = "Maximum UMI count per cell", y = "Expressing cells (%)") +
  guides(
    fill = guide_legend(order = 1, nrow = 1, override.aes = list(size = 3, alpha = 1)),
    size = guide_legend(
      order = 2, nrow = 1,
      override.aes = list(fill = "grey70", color = "grey35", alpha = 1)
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10.5, margin = margin(t = 7)),
    axis.title.y = element_text(size = 10.5, margin = margin(r = 7)),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.30),
    
    panel.grid.major.y = element_line(color = "grey94", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = grid::unit(8, "mm"),
    
    strip.background = element_blank(),
    strip.text = element_text(face = "italic", size = 12, margin = margin(t = 3, b = 6)),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 9, face = "plain"),
    legend.text = element_text(size = 8.5),
    legend.key.width = grid::unit(5, "mm"),
    legend.key.height = grid::unit(5, "mm"),
    legend.spacing.x = grid::unit(3, "mm"),
    legend.box.spacing = grid::unit(2, "mm"),
    
    plot.margin = margin(t = 5, r = 8, b = 5, l = 6)
  )

pdf("./Proj_CTAs/output/figure.panels/Tabula.Sapiens.CT83.DCAF4L2.pdf", width = 10.4, height = 5)
print(p)
dev.off()

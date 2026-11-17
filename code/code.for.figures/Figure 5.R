library(tidyverse)
library(foreach)
library(doParallel)
library(Seurat)
library(patchwork)
library(magrittr)
library(ggridges)

##- Fig. 5a and 5b -------------------------------------------------------------
load("./Proj_CTAs/output/manuscript.revision/ST.data.RData")

rho.df <- foreach(i = names(ST.meta.ls), .combine = 'rbind') %do% {
  meta <- ST.meta.ls[[i]]
  vars <- c("e.CTGs", "n.CTGs", "n.genes", "percent.mt")
  
  map_dfr(vars, \(var) {
    x   <- meta[[var]]
    fit <- cor.test(x, meta$ESTIMATE_purity, method = "spearman")
    data.frame(unique_sample_id = i, var = var, rho = fit$estimate, p.val = fit$p.value, row.names = NULL)
  })

} |> inner_join(ST.data, by = "unique_sample_id")
# Fig. 5a
plot.df <- rho.df |> filter(var == "n.CTGs") |>
  group_by(cancer, unique_patient_id) |>
  summarise(rho = median(rho, na.rm = TRUE), .groups = "drop")

p       <- ggplot(plot.df, aes(x = reorder(cancer, rho, median), y = rho)) + 
  geom_hline(yintercept = c(0, 0.5), linetype = "dashed", color = "#606060", linewidth = 0.4) +
  
  geom_violin(fill = "#E9EEF2", color = NA, alpha = 1, scale = "width", width = 0.7, trim = TRUE) +
  
  geom_boxplot(width = 0.12, fill = "white", color = "grey20", 
               linewidth = 0.4, outlier.shape = NA, alpha = 0.9) +
  
  geom_jitter(aes(fill = rho), size = 1.8, alpha = 0.85, shape = 21, color = "grey30", stroke = 0.25,
              width = 0.12, height = 0) +
  
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",                  
                       midpoint = 0, limits = c(-1, 1),
                       name = expression(rho),
                       ) +
  
  scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1), limits = c(-0.5, 1)) +
  labs(x = NULL, y = expression(paste("Spearman's ", rho))) +
  
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(hjust = 1, angle = 30), 
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.4),
        panel.grid.major.x = element_blank(),
        legend.position = c(0.8, 0.08),
        legend.direction = "horizontal")

pdf("./Proj_CTAs/output/figure.panels/Fig. 5a.pdf", width = 7, height = 5)
p
dev.off()

# Fig. 5b
registerDoParallel(cores = 50)
samples <- c(
    "HTAN_HT242P1H2-S1Fc1U1Z1B1", "HTAN_S18-25657-Fs1U1Bp1", "SCAR_ST_000099", 
    "SCAR_ST_000091", "HTAN_TB-4951U1Bp1", "HTAN_S16-38794-A3U1Bp1"
)
seu.ls <- foreach(i = samples) %dopar% {
  path  <- ST.data$processed_path[ST.data$unique_sample_id == i]
  seu   <- readRDS(path)
  meta  <- ST.meta.ls[[i]]
  identical(colnames(seu), rownames(meta)) |> stopifnot()
  seu$x <- log10(meta$n.CTGs + 1)
  seu
} |> setNames(samples)

plots <- lapply(samples, \(i) {
  seu   <- seu.ls[[i]]
  
  Min   <- min(seu$ESTIMATE_purity)
  Max   <- max(seu$ESTIMATE_purity)
  p1    <- SpatialFeaturePlot(seu, features = "ESTIMATE_purity", pt.size.factor = 1.6, stroke = 0) + 
    viridis::scale_fill_viridis(
      option = "viridis", direction = 1, breaks = c(Min, Max), labels = c("Min", "Max"), name = "Purity",
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5, frame.colour = "black", ticks = FALSE)
    ) +
    labs(fill = "Tumor Purity", subtitle = paste("Sample ID:", i), title = ST.data$cancer[ST.data$unique_sample_id == i]) +
    theme_void() + 
    theme(aspect.ratio = 1, legend.position = "right", legend.title = element_text(size = 12, vjust = 1),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, color = "#7f8c8d", size = 5))
  
  Min   <- min(seu$x)
  Max   <- max(seu$x)
  p2    <- SpatialFeaturePlot(seu, features = "x", pt.size.factor = 1.6, stroke = 0) + 
    viridis::scale_fill_viridis(
      option = "magma", direction = 1, breaks = c(Min, Max), labels = c("Min", "Max"), name = expression(italic(N)[CTG]),
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5, frame.colour = "black", ticks = FALSE)
    ) +
    theme_void() + 
    theme(aspect.ratio = 1, legend.position = "right", legend.title = element_text(size = 12, vjust = 1))
  
  plot <- p1 / p2 
})

## remove redundant legends
for(i in c(1, 2, 4, 5)) {
  plots[[i]] <- plots[[i]] & theme(legend.position = "none")
}

pdf("./Proj_CTAs/output/figure.panels/Fig. 5b.pdf", width = 6.5, height = 8)
plots |> wrap_plots(ncol = 3)
dev.off()

##- Fig. 5c --------------------------------------------------------------------

##- Fig. 5d --------------------------------------------------------------------
load("./Proj_CTAs/output/manuscript.revision/ST.data.RData")
load("./Proj_CTAs/output/manuscript.revision/ST.tumor.layer.RData")

## correlation between tumor layers and purity adjusted n_CTGs for each sample
rho.df <- foreach(s = names(tumor.layer.ls), .combine = 'rbind') %do% {
  ST.meta  <- ST.meta.ls[[s]]
  layer.ls <- tumor.layer.ls[[s]]

  map_dfr(names(layer.ls), \(j) {
    x     <- layer.ls[[j]]

    if (!inherits(x, "data.frame")) return( data.frame(unique_sample_id = s, args = j, evaluable = FALSE, reason = as.character(x)[1], row.names = NULL) )

    idx   <- x$propagated.labels == "Malignant"
    layer <- x$propagated.layers[idx]
    bc    <- rownames(x)[idx]

    base  <- data.frame(unique_sample_id = s, args = j,
                        n.mal.spots = sum(idx),
                        n.layers    = dplyr::n_distinct(layer),
                        max.layer   = if (sum(idx)) max(layer, na.rm = TRUE) else NA_real_,
                        row.names   = NULL)

    if (base$n.mal.spots < 50 || base$n.layers < 3) return(cbind(base, evaluable = FALSE, reason = "too few spots or layers", slope = NA_real_, rho = NA_real_))

    burden <- ST.meta[bc, "n.CTGs"]
    purity <- ST.meta[bc, "ESTIMATE_purity"]
    adj    <- lm(burden ~ purity) |> residuals()

    cbind(base, evaluable = TRUE, reason = NA_character_,
          slope = unname(coef(lm(adj ~ layer))[2]),
          rho   = cor(adj, layer, method = "spearman"))
  })

} |> left_join(ST.data, by = "unique_sample_id")
patient.rho <- rho.df |> filter(evaluable) |>
  group_by(cancer, unique_patient_id, args) |>
  summarise(rho = median(rho), n.slides = n(), .groups = "drop")

args.chose   <- "k=6, seed=c(2,1)"
min.patients <- 5

conf.plot.df <- patient.rho |>
  filter(args == args.chose, is.finite(rho)) |>
  group_by(cancer) |>
  group_modify(~ {
    x         <- sort(.x$rho)
    n.patient <- length(x)

    averages  <- outer(x, x, "+") / 2
    hl.rho    <- median(averages[upper.tri(averages, diag = TRUE)])

    ## with n < 6, a finite two-sided exact 95% signed-rank CI is unavailable.
    fit      <- if (n.patient >= 6) {wilcox.test(x, mu = 0, conf.int = TRUE, conf.level = 0.95)} else NULL

    ci.level <- if (is.null(fit)) NA_real_ else attr(fit$conf.int, "conf.level")
    ci.ok    <- !is.null(fit) && isTRUE(ci.level >= 0.95) && all(is.finite(fit$conf.int))

    data.frame(
      n.patient = n.patient, med.rho = median(x), hl.rho = hl.rho,
      conf.low  = if (ci.ok) unname(fit$conf.int[1]) else NA_real_,
      conf.high = if (ci.ok) unname(fit$conf.int[2]) else NA_real_,
      ci.level  = ci.level
    )
  }) |>
  ungroup() |>
  filter(n.patient >= min.patients) |>
  arrange(med.rho)

cancer.levels <- as.character(conf.plot.df$cancer)
conf.plot.df  <- conf.plot.df |>
  mutate(
    cancer       = factor(cancer, levels = cancer.levels),
    effect.label = ifelse(
      is.na(conf.low), sprintf("%.2f [CI not shown]", hl.rho),
      sprintf("%.2f [%.2f, %.2f]", hl.rho, conf.low, conf.high)
    ),
    lab = ifelse(
      is.finite(conf.low) & is.finite(conf.high),
      sprintf("paste('HL(', rho, ') = %.2f; 95%% CI [%.2f, %.2f]')", hl.rho, conf.low, conf.high),
      sprintf("paste('HL(', rho, ') = %.2f')", hl.rho)
    )
  )

NBIN              <- 4
MIN.SPOTS.PER.BIN <- 10
min.bin.patients  <- min.patients
bin.labs          <- sprintf("%.0f-%.0f%%", (0:(NBIN - 1)) / NBIN * 100, (1:NBIN) / NBIN * 100)
bin.colors        <- setNames(
  colorRampPalette(c("#6B9AC4", "white", "#D47264"))(NBIN),
  as.character(seq_len(NBIN))
)

sample.df <- ST.data |> select(unique_sample_id, cancer, unique_patient_id) |> distinct()
slide.ids <- rho.df |> filter(args == args.chose, evaluable, is.finite(rho), cancer %in% cancer.levels) |>
  pull(unique_sample_id)

bin.df <- map_dfr(slide.ids, function(s) {
  x       <- tumor.layer.ls[[s]][[args.chose]]
  idx     <- x$propagated.labels %in% "Malignant"
  bc      <- rownames(x)[idx]
  layer   <- x$propagated.layers[idx]
  ST.meta <- ST.meta.ls[[s]]
  burden  <- ST.meta[bc, "n.CTGs"]
  purity  <- ST.meta[bc, "ESTIMATE_purity"]
  ok      <- is.finite(burden) & is.finite(purity) & is.finite(layer)
  burden  <- burden[ok]
  purity  <- purity[ok]
  layer   <- layer[ok]
  stopifnot(all(layer >= 0), max(layer) > 0)

  adj <- residuals(lm(burden ~ purity))
  bin <- pmin(floor(layer / max(layer) * NBIN) + 1L, NBIN)
  tibble(unique_sample_id = s, bin = bin, adj = adj)
}) |>
  left_join(sample.df, by = "unique_sample_id")

slide.bin  <- bin.df |>
  group_by(cancer, unique_patient_id, unique_sample_id, bin) |>
  summarise(adj = median(adj), n.spots = n(), .groups = "drop") |>
  filter(n.spots >= MIN.SPOTS.PER.BIN)

patient.bin <- slide.bin |>
  group_by(cancer, unique_patient_id, bin) |>
  summarise(adj = median(adj), n.slides = n(), .groups = "drop")

plot.df     <- patient.bin |> mutate(cancer = factor(cancer, levels = cancer.levels))
summary.df  <- plot.df |> group_by(cancer, bin) |>
  summarise(
    median.adj = median(adj),
    n.patient  = n_distinct(unique_patient_id),
    .groups    = "drop"
  ) |>
  complete(cancer, bin = seq_len(NBIN), fill = list(n.patient = 0L)) |>
  mutate(median.adj = ifelse(n.patient >= min.bin.patients, median.adj, NA_real_))

median.line.df <- summary.df |> group_by(cancer) |> arrange(bin, .by_group = TRUE) |>
  mutate(next.bin = lead(bin), next.adj = lead(median.adj)) |> ungroup() |>
  filter(next.bin == bin + 1L, is.finite(median.adj), is.finite(next.adj))
box.df         <- plot.df |> semi_join(filter(summary.df, n.patient >= min.bin.patients), by = c("cancer", "bin"))

p <- ggplot(plot.df, aes(x = bin, y = adj)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.35) +
  geom_boxplot(
    data = box.df, aes(group = bin, fill = factor(bin, levels = seq_len(NBIN))),
    color = "grey35", outlier.shape = NA, width = 0.65, linewidth = 0.35
  ) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 1234),
    color = "grey35", size = 0.85, alpha = 0.50
  ) +
  geom_segment(
    data = median.line.df, aes(x = bin, xend = next.bin, y = median.adj, yend = next.adj),
    color = "#B5423A", linewidth = 0.65
  ) +
  geom_point(data = summary.df, aes(y = median.adj), color = "#B5423A", size = 1.6, na.rm = TRUE) +

  geom_text(
    data = conf.plot.df, aes(x = -Inf, y = Inf, label = lab),
    inherit.aes = FALSE, parse = TRUE,
    hjust = -0.025, vjust = 1.8, size = 2.7, color = "grey20", na.rm = TRUE
  ) +
  geom_text(
    data = summary.df, aes(x = bin, y = Inf, label = paste0("n = ", n.patient)),
    inherit.aes = FALSE, vjust = 5, size = 2.5, color = "grey30"
  ) +

  scale_fill_manual(values = bin.colors, guide = "none", drop = FALSE) +
  scale_x_continuous(breaks = seq_len(NBIN), labels = bin.labs, expand = expansion(add = 0.45)) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.42))) +
  facet_wrap(
    ~ cancer, scales = "free_y", ncol = 6, drop = FALSE
  ) +
  labs(
    x = "Normalized layer index (%)",
    y = expression(paste("Median ", italic(N)[CTG], " (purity-adjusted)"))
  ) +
  theme_bw(base_size = 11) +
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.text = element_text(size = 8.5, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(size = 11, margin = margin(t = 8)),
    axis.title.y = element_text(size = 11, margin = margin(r = 8)),
    axis.ticks = element_line(color = "black", linewidth = 0.30),

    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.35),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),

    strip.background = element_rect(fill = "grey94", color = NA),
    strip.text = element_text(face = "plain", size = 10, margin = margin(t = 5, b = 5)),
    legend.position = "none",
    plot.margin = margin(t = 7, r = 8, b = 7, l = 8)
  )

pdf("./Proj_CTAs/output/figure.panels/nCTG.layer.line.plots.pdf", width = 14, height = 5)
print(p)
dev.off()

##- Fig. 5e --------------------------------------------------------------------
load("./Proj_CTAs/output/manuscript.revision/ST.data.RData")
Keratin.samples <- c("GSE208253_OSCC_sample_3", "GSE208253_OSCC_sample_4", "GSE208253_OSCC_sample_5", "GSE208253_OSCC_sample_6")
my.cols         <- c("Tumor" = "#DC0000B2", "Keratin" = "#00A087B2")
addMedian       <- function(...) {
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, 
               fill = "white", color = "black", stroke = 0.5, 
               position = position_dodge(width = 0.9), ...)
}

plots <- foreach(i = Keratin.samples) %do% {
  seu       <- readRDS(ST.data$processed_path[ST.data$unique_sample_id == i])

  ST.meta   <- ST.meta.ls[[i]] |> select(n.CTGs, ESTIMATE_purity)
  stopifnot(identical(colnames(seu), rownames(ST.meta)))
  seu$n.CTGs <- ST.meta$n.CTGs

  seu           <- subset(seu, pathologist_anno.x %in% c("SCC", "Keratin"))
  seu@meta.data <- seu@meta.data |> mutate(pathologist_anno.x = ifelse(pathologist_anno.x == "SCC", "Tumor", "Keratin"))
  
  p1    <- SpatialDimPlot(seu, group.by = "pathologist_anno.x", pt.size.factor = 1.6, stroke = 0, cols = my.cols, alpha = 0.9) + 
    labs(fill = "Histology", subtitle = paste("Sample ID:", i)) +
    theme_void() + 
    theme(aspect.ratio = 1, 
          legend.position = "top", legend.title = element_text(size = 12, vjust = 1), 
          plot.subtitle = element_text(hjust = 0, color = "#7f8c8d", size = 8)) +
    guides(fill = guide_legend(override.aes = list(size = 4)))
  
  Min   <- min(seu$n.CTGs)
  Max   <- max(seu$n.CTGs)
  p2    <- SpatialFeaturePlot(seu, features = "n.CTGs", pt.size.factor = 1.6, stroke = 0) + 
    viridis::scale_fill_viridis(
      option = "magma", direction = 1, limits = c(Min, Max), breaks = c(Min, Max), labels = c("Min", "Max"), name = expression(italic(N)[CTG]),
      guide = guide_colorbar(title.position = "left", title.hjust = 0.5, frame.colour = "black", ticks = FALSE)
    ) +
    theme_void() + 
    theme(aspect.ratio = 1, legend.position = "top", legend.title = element_text(size = 12, vjust = 1))
  
  y.max <- max(seu$n.CTGs, na.rm = TRUE)
  p3    <- VlnPlot(seu, group.by = "pathologist_anno.x", features = "n.CTGs", pt.size = 0, cols = my.cols) + 
    addMedian() + 
    ggpubr::stat_compare_means(
      comparisons = list(c("Tumor", "Keratin")), 
      label = "p.format", label.y = y.max * 1.05, tip.length = 0.01,
      size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = NULL, y = expression(italic(N)[CTG]), title = NULL) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12, color = "black", angle = 0, hjust = 0.5),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.y = element_text(size = 12, face = "italic", margin = margin(r = 5)))
  
  p <- p1 + p2 + p3 + plot_layout(widths = c(4, 4, 2))
}

p1 <- plots[[4]][[1]]
p2 <- plots[[4]][[2]]
p3 <- plots[[4]][[3]]
p3 <- p3 + 
  geom_boxplot(width = 0.25, fill = "white", color = "black", outliers = FALSE) +
  theme(plot.margin = margin(0, 0, 0, 10))
p3 <- wrap_plots(plot_spacer(), p3, plot_spacer(), 
                 ncol = 1, heights = c(0.4, 2.5, 0.3))

pdf("./Proj_CTAs/output/figure.panels/HNSCC.n-CTGs.violin.one.sample.pdf", width = 9, height = 5)
wrap_plots(p1, p2, p3, ncol = 3, widths = c(1, 1, 0.6))
dev.off()

##- Fig. 5f --------------------------------------------------------------------
load("./Proj_CTAs/output/manuscript.revision/HNSCC.GSE268014.RData")
merged.metadata.Extend <- readRDS("./Proj_CTAs/input/ST/GSE268014/metadata/merged_metadata_Extend.rds")
compare.df             <- score.df |> inner_join(merged.metadata.Extend, by = c("Sample", "Key")) |>
  filter(binCNAstatus == "Malignant", Zone %in% paste0("Zone_", 1:3)) |>
  mutate(Zone = sub("_", " ", as.character(Zone))) |>
  group_by(Sample) |>
  group_modify(~ {
    adj <- lm(.x$n.CTGs ~ .x$ESTIMATE.purity) |> residuals()
    data.frame(value = tapply(adj, .x$Zone, median, na.rm = TRUE)) |> rownames_to_column("Zone")
  })

plot.df     <- compare.df
zone.order  <- paste0("Zone ", 1:3)
zone.colors <- c("Zone 1" = "#D47264", "Zone 2" = "#D9B382", "Zone 3" = "#6B9AC4")
y.range     <- range(plot.df$value)
y.span      <- max(diff(y.range), 1)

## keep the widest comparison above the adjacent-zone comparisons.
label.y     <- y.range[2] + y.span * c(0.22, 0.06, 0.14)
y.limits    <- c(y.range[1] - 0.04 * y.span, y.range[2] + 0.28 * y.span)

ls <- combn(paste0("Zone ", c(3, 1, 2)), 2, simplify = FALSE)
p  <- ggplot(plot.df, aes(x = Zone, y = value)) +
  geom_line(aes(group = Sample), color = "grey75", linewidth = 0.35, alpha = 0.65) +
  
  geom_boxplot(
    aes(fill = Zone), width = 0.58, outlier.shape = NA,
    color = "grey30", linewidth = 0.35, alpha = 0.20
  ) +
  geom_point(aes(fill = Zone), shape = 21, size = 2.3, color = "white", stroke = 0.20, alpha = 0.85) +
  
  ggpubr::stat_compare_means(
    comparisons = ls, method = "wilcox.test", paired = TRUE,
    label = "p.format", label.y = label.y, step.increase = 0,
    tip.length = 0.01, bracket.size = 0.30, size = 3.3
  ) +
  
  scale_fill_manual(values = zone.colors) +
  scale_x_discrete(labels = zone.order) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 5)) +
  coord_cartesian(
    xlim = c(0.4, 3.6),
    ylim = c(y.range[1] - 0.08 * y.span, y.range[2] + 0.4 * y.span), expand = FALSE
  ) +
  
  labs(x = NULL, y = expression(paste("Median ", italic(N)[CTG], " (purity-adjusted)"))) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 9.5, color = "black"),
    axis.title.y = element_text(size = 10.5, margin = margin(r = 8)),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.30),
    
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 8, r = 8, b = 7, l = 7)
  )

pdf("./Proj_CTAs/output/figure.panels/HNSCC.GSE268014.nCTG.pdf", width = 3.8, height = 4.2)
print(p)
dev.off()

##- Fig. 5g --------------------------------------------------------------------
load("./Proj_CTAs/output/manuscript.revision/tumor.sample.info.updated.RData")
load("./Proj_CTAs/output/manuscript.revision/hypoxia.score.RData")
hypoxia.auc <- hypoxia.score |> group_by(dataset_sample) |>
  group_modify(~ {
    x <- .x$HALLMARK_HYPOXIA
    y <- factor(.x$CTG_group, levels = c("CTG-low", "CTG-high"))

    if (min(table(y)) >= 10) {
      auc <- pROC::roc(y, x, direction = "<")$auc |> as.numeric()
      data.frame(auc)
    } else {
      data.frame(auc = NA)
    }
  }) |> ungroup() |>
  left_join(tumor.sample.info, by = "dataset_sample")

source("./Proj_CTAs/code/visualization.R")
plot.df     <- hypoxia.auc |> filter(dataset_group %in% c("Discovery Cohort"), !is.na(auc)) |> 
  mutate(unique_patient = paste0(dataset, "_", patient)) |>
  group_by(cancer, unique_patient) |>
  group_modify(~ { data.frame(auc = median(.x$auc, na.rm = TRUE)) }) |>
  ungroup() |>

  left_join(cancer.abbr, by = "cancer") |>
  mutate(abbr = fct_reorder(abbr, auc, .fun = median, .desc = TRUE))
sig.df      <- plot.df |>
  group_by(abbr) |>
  summarise(p.val = wilcox.test(auc, mu = 0.5)$p.value, n.patients = n()) |>
  mutate(
    sig.star = case_when(
      p.val <= 0.0001 ~ "****",
      p.val <= 0.001  ~ "***",
      p.val <= 0.01   ~ "**",
      p.val <= 0.05   ~ "*",
      TRUE            ~ "ns"
    ))

ridge.scale <- 1.1
p           <- ggplot(plot.df, aes(x = auc, y = abbr, fill = abbr)) +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey60", size = 0.5) +
  
  stat_density_ridges(
    geom = "density_ridges_gradient",
    scale = ridge.scale,
    alpha = 0.85,
    quantile_lines = TRUE,  ## median lines
    quantiles = 2, 
    quantile_fun = median,
    linetype = "dashed",
    color = "grey30",
    size = 0.3) +
  
  ## draw a white solid border to cover the dashed edge of the upper layer
  geom_density_ridges(
    scale = ridge.scale,
    fill = NA,
    color = "white",
    size = 0.7,
    quantile_lines = FALSE) +
  
  ## add significance
  geom_text(
    data = sig.df, 
    aes(x = 0.95, y = abbr, label = sig.star), 
    inherit.aes = FALSE, 
    hjust = 0,
    vjust = -0.2,
    size = 4.5, 
    fontface = "plain", 
    color = "black"
  ) +
  
  scale_fill_manual(values = cancer.colours) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.1), expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "AUROC", y = NULL, title = "Hypoxia (Hallmark)") +
  
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.y = element_text(size = 11), 
        axis.text.x = element_text(size = 11),
        legend.position = "none",
        plot.margin = margin(t = 15, r = 15, b = 15, l = 10))

pdf("./Proj_CTAs/output/figure.panels/hypoxia.ridge.plot.pdf", height = 6, width = 4)
p
dev.off()

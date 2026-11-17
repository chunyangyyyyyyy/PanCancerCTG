library(tidyverse)
library(Seurat)
library(foreach)
library(doParallel)
library(clusterProfiler)
library(magrittr)

##- Fig. 3a and 3b -------------------------------------------------------------
source("./Proj_CTAs/code/visualization.R")
load("./Proj_CTAs/output/manuscript.revision/pan.sample.CTGs.RData")
load("./Proj_CTAs/output/manuscript.revision/gene.CTG.cor.ls.RData")

rho.df <- foreach(i = cancer.abbr$abbr, .combine = 'rbind') %do% {
  gene.CTG.cor.ls[[i]] %>%
    mutate(gene_group = case_match(
      gene,
      pan.sample.CTGs.ls[[i]]          ~ "CTGs",
      setdiff(.$gene, pan.sample.CTGs) ~ "Non-CTGs"
    )) |>
    dplyr::filter(!is.na(gene_group), CTG %in% pan.sample.CTGs.ls[[i]], gene != CTG, !is.na(rho)) |>  ## !CTG %in% cc.CTGs, 
    group_by(CTG) |>
    group_modify(~ {
      CTG.median     <- .x |> dplyr::filter(gene_group == "CTGs")     |> pull(rho) |> median()
      non.CTG.median <- .x |> dplyr::filter(gene_group == "Non-CTGs") |> pull(rho) |> median()
      data.frame(abbr = i, rho = c(CTG.median, non.CTG.median), gene_group = c("CTGs", "Non-CTGs"))
    }) |>
    ungroup()
}
# Fig. 3a
group.colors <- c("CTGs" = "#D07B5D", "Non-CTGs" = "#687698")
plot.df      <- rho.df |>
  mutate(abbr = factor(abbr, levels = cancer.abbr$abbr), gene_group = factor(gene_group, levels = c("CTGs", "Non-CTGs"))) |>
  arrange(abbr, CTG, gene_group)

rho.range <- range(plot.df$rho)
rho.span  <- max(diff(rho.range), 0.1)

p <- ggplot(plot.df, aes(x = gene_group, y = rho)) +
  geom_hline(
    yintercept = 0, linetype = "dashed",
    color = "grey80", linewidth = 0.25
  ) +
  geom_violin(
    aes(fill = gene_group, color = gene_group),
    trim = TRUE, scale = "width",
    width = 0.85, alpha = 0.85, linewidth = 0.35,
    key_glyph = "rect"
  ) +
  geom_boxplot(
    width = 0.15, fill = "white", color = "grey20",
    outlier.shape = NA, linewidth = 0.3,
    show.legend = FALSE
  ) +
  ggpubr::stat_compare_means(
    comparisons = list(c("CTGs", "Non-CTGs")),
    method = "wilcox.test", paired = TRUE,
    label = "p.signif", hide.ns = FALSE,
    label.y = rho.range[2] + 0.015 * rho.span,
    tip.length = 0.012, bracket.size = 0.25,
    size = 2.8, vjust = 0
  ) +
  facet_wrap(~ abbr, ncol = 7, scales = "fixed") +
  scale_color_manual(values = group.colors, guide = "none") +
  scale_fill_manual(
    values = group.colors,
    breaks = c("CTGs", "Non-CTGs"),
    name = NULL
  ) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_continuous(
    breaks = seq(-1, 1, 0.5),
    expand = expansion(mult = c(0.03, 0.10))
  ) +
  labs(x = NULL, y = expression(paste("Median Spearman's ", rho))) +
  guides(
    fill = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(
        alpha = 0.85,
        colour = unname(group.colors[c("CTGs", "Non-CTGs")])
      )
    )
  ) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(family = "sans", color = "black"),

    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(
      face = "plain", size = 9, color = "black",
      margin = margin(t = 4, b = 4)
    ),

    axis.title.y = element_text(size = 10, margin = margin(r = 6)),
    axis.text.y = element_text(size = 8.5, color = "black"),
    axis.line.y = element_line(linewidth = 0.35, color = "black"),
    axis.ticks.y = element_line(linewidth = 0.3, color = "black"),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),

    panel.grid = element_blank(),
    panel.spacing.x = grid::unit(0.45, "lines"),
    panel.spacing.y = grid::unit(0.75, "lines"),

    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = grid::unit(0.32, "cm"),
    legend.key.spacing.x = grid::unit(8, "pt"),
    legend.text = element_text(size = 9, color = "black"),
    legend.margin = margin(t = 2, r = 0, b = 0, l = 0),

    plot.margin = margin(t = 6, r = 6, b = 6, l = 6)
  )

pdf("./Proj_CTAs/output/figure.panels/Fig. 3A.pdf", width = 6, height = 5)
print(p)
dev.off()

# Fig. 3b
x.ls        <- list(x.all, x.non.cc, x.cc)
names(x.ls) <- c("Total CTGs", "Non-Cell Cycle CTGs", "Cell Cycle-CTGs")

seq.palette     <- c("#EAEAEA", "#FEE08B", "#F46D43", "#D73027", "#82001A")
target.pathways <- c("Chromatin modifying enzymes", "Chromatin organization")
plots           <- foreach(i = names(x.ls)) %do% {
  x                      <- x.ls[[i]]
  p                      <- enrichplot::dotplot(x, font.size = 10, label_format = Inf, showCategory = 3)

  df                     <- p$data |> mutate(
    abbr           = str_extract(Cluster, ".+(?=\\n)"),
    FoldEnrichment = GeneRatio / (as.numeric(sub("/.*$", "", BgRatio)) / as.numeric(sub("^.*/", "", BgRatio)))
  )
  levels(df$Description) <- str_wrap(levels(df$Description), width = 55)
  cancer.omitted         <- setdiff(cancer.abbr$abbr, df$abbr)
  df.add                 <- df[rep(1, length(cancer.omitted)), ] |> mutate(abbr = cancer.omitted, GeneRatio = NA, FoldEnrichment = NA, p.adjust = NA)
  df                     <- rbind(df, df.add) |> 
    mutate(Cluster = factor(abbr, levels = cancer.abbr$abbr), log10.p = -log10(p.adjust))
  
  y.levels    <- levels(df$Description)
  text.colors <- ifelse(y.levels %in% target.pathways, "#82001A", "black")
  text.faces  <- ifelse(y.levels %in% target.pathways, "bold", "plain")
  
  p1 <- ggplot(df, aes(x = Cluster, y = Description)) + 
    
    geom_point(aes(size = FoldEnrichment, color = log10.p), alpha = 0.85) +
    
    labs(title = i, x = NULL, y = NULL) +
    
    scale_size_continuous(
      range = c(2, 5.5), 
      name = "Fold\nenrichment",
      guide = guide_legend(order = 1, override.aes = list(color = "grey60"))
    ) + 
    
    scale_color_gradientn(
      colors = seq.palette,
      name = bquote(-log[10]("p.adj")),
      na.value = "transparent",
      limits = c(0, 12),
      oob = scales::squish,
      guide = guide_colorbar(
        order = 2,
        ticks.linewidth = 0.5, 
        frame.colour = "black", 
        frame.linewidth = 0.5
      )
    ) +
    
    theme_classic(base_size = 12) +
    theme(
      plot.title          = element_text(hjust = 0.5, face = "plain", size = 12, color = "black"),
      
      axis.line           = element_line(linewidth = 0.5, color = "black"),
      axis.ticks          = element_line(linewidth = 0.5, color = "black"),
      
      axis.text.x         = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 11),
      
      axis.text.y         = element_text(color = text.colors, face = text.faces, size = 10),
      
      panel.grid.major.y  = element_line(linewidth = 0.25, color = "grey85", linetype = "dashed"),
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_blank(),
      
      legend.position     = "right",
      legend.background   = element_blank(),
      legend.title        = element_text(face = "plain", size = 8, color = "black"),
      legend.text         = element_text(size = 7, color = "black"),
      legend.key.size     = unit(0.35, "cm")
    )
}

pdf("./Proj_CTAs/output/figure.panels/Fig. 3B.pdf", width = 10, height = 5)
plots[[1]]
dev.off()

##- Fig. 3c --------------------------------------------------------------------
load("./Proj_CTAs/output/manuscript.revision/regulator.ranking.RData")
gene.swi <- c("ARID1A", "ARID1B", "ARID2", "PBRM1", "SMARCC1")
gene.mod <- NULL  ## c("EP300", "KMT2D", "KDM6A")
plot.df  <- regulator.ranking |> arrange(desc(recurrence)) |>
  mutate(
    rank.val = row_number(),
    grp.hl   = case_when(
      gene %in% gene.swi ~ "SWI/SNF",
      gene %in% gene.mod ~ "Modifiers",
      TRUE ~ "Other"
    )) |>
  mutate(grp.hl = factor(grp.hl, levels = c("Other", "Modifiers", "SWI/SNF")))

col.hl       <- c("Other" = "#E8E8E8", "Modifiers" = "#4575B4", "SWI/SNF" = "#D53E4F")
highlight.df <- plot.df |> filter(grp.hl != "Other")

p <- ggplot(plot.df, aes(x = rank.val, y = recurrence)) +
  geom_hline(yintercept = 2.5, linetype = "dashed", color = "grey45", linewidth = 0.35) +

  geom_point(
    data = filter(plot.df, grp.hl == "Other"),
    color = "grey75", shape = 16, size = 2.2, alpha = 0.65
  ) +
  geom_point(
    data = highlight.df, aes(fill = grp.hl),
    shape = 21, size = 3.2, color = "grey20", stroke = 0.35, alpha = 0.90
  ) +

  annotate(
    "text", x = max(plot.df$rank.val, na.rm = TRUE), y = 2.25,
    label = "'Recurrence threshold:' >= 3~'cancer types'",
    parse = TRUE, hjust = 1, vjust = 0.5, size = 3, color = "grey35"
  ) +
  ggrepel::geom_text_repel(
    data = highlight.df, aes(label = gene, color = grp.hl),
    fontface = "italic", size = 3.4, seed = 1234, max.overlaps = Inf,
    box.padding = 0.5, point.padding = 0.35, force = 4, nudge_x = 2, nudge_y = 0.55,
    segment.color = "grey50", segment.size = 0.3, min.segment.length = 0,
    show.legend = FALSE
  ) +

  scale_fill_manual(values = col.hl) +
  scale_color_manual(values = col.hl) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 4), expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(
    breaks = seq(0, max(plot.df$recurrence, na.rm = TRUE) + 2, 2),
    expand = expansion(mult = c(0.06, 0.10))
  ) +

  labs(x = "Rank of chromatin regulators", y = "Recurrence (number of cancer types)") +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11, face = "plain", margin = margin(t = 7)),
    axis.title.y = element_text(size = 11, face = "plain", margin = margin(r = 7)),
    axis.line = element_line(color = "black", linewidth = 0.40),
    axis.ticks = element_line(color = "black", linewidth = 0.35),

    panel.grid.major.y = element_line(color = "grey94", linewidth = 0.25, linetype = "dashed"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),

    legend.position = "none",
    plot.margin = margin(t = 7, r = 9, b = 6, l = 7)
  )

pdf("./Proj_CTAs/output/figure.panels/Fig. 3C.pdf", width = 5, height = 4.5)
print(p)
dev.off()

##- Fig. 3d --------------------------------------------------------------------
library(GEOquery)
library(DESeq2)
library(enrichplot)
library(ggtext)
gse        <- getGEO(filename = "./Proj_CTAs/input/bulk RNA-seq datasets/GSE69567_siARID1B_HepG2/GSE69567_series_matrix.txt.gz", getGPL = FALSE)
meta       <- pData(gse)
grch38.ann <- data.table::fread("./Proj_CTAs/input/bulk RNA-seq datasets/GSE69567_siARID1B_HepG2/Human.GRCh38.p13.annot.tsv.gz") |>
  mutate(GeneID = as.character(GeneID))
counts     <- data.table::fread("./Proj_CTAs/input/bulk RNA-seq datasets/GSE69567_siARID1B_HepG2/GSE69567_raw_counts_GRCh38.p13_NCBI.tsv.gz") |>
  column_to_rownames("GeneID") |> as.matrix()

coldata       <- meta |> dplyr::select(sample = geo_accession, group = source_name_ch1) |>
  mutate(group = gsub("HepG2_|_siRNA", "", group)) |>
  filter(group %in% c("Arid1a", "Arid1b", "Arid2", "non-targeting"))
coldata$group <- factor(coldata$group)

protein.genes <- grch38.ann |> filter(GeneType == "protein-coding") |> pull(GeneID) |> intersect(rownames(counts))
counts        <- counts[protein.genes, rownames(coldata)]
counts        <- counts[rowSums(counts) > 0, ]

dds           <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ group)
dds           <- DESeq(dds)

## CTG terms
load("./Proj_CTAs/output/manuscript.revision/pan.sample.CTGs.RData")
term2gene       <- data.frame(term = "Total CTGs", gene = pan.sample.CTGs)

## GSEA
gsea.ls        <- foreach(i = c("Arid1a", "Arid1b", "Arid2")) %do% {
  res      <- results(dds, contrast = c("group", i, "non-targeting")) |> as.data.frame() |> 
    rownames_to_column("GeneID") |>
    left_join(grch38.ann, by = "GeneID")
  
  geneList <- setNames(res$log2FoldChange, res$Symbol) |> sort(decreasing = TRUE)
  set.seed(1234)
  gsea     <- clusterProfiler::GSEA(geneList = geneList, TERM2GENE = term2gene, pvalueCutoff = 1, minGSSize = 1, maxGSSize = 500)
}
names(gsea.ls) <- c("Arid1a", "Arid1b", "Arid2")
curve.ls <- foreach(i = c("Arid1a", "Arid1b", "Arid2")) %do% {
  stat.df <- gsea.ls[[i]]@result |> filter(ID == "Total CTGs")
  p.plot  <- enrichplot::gseaplot2(gsea.ls[[i]], geneSetID = "Total CTGs")

  p.label <- ifelse(
    stat.df$pvalue < 0.001,
    "<i>P</i> &lt; 0.001",
    sprintf("<i>P</i> = %.3f", stat.df$pvalue)
  )

  p.plot[[1]]$data |>
    mutate(
      condition.name = toupper(i),
      legend.label   = sprintf("<i>%s</i> KD (NES = %.2f, %s)", toupper(i), stat.df$NES, p.label)
    )
}

group.order      <- c("ARID1A", "ARID1B", "ARID2")
merged.gsea.data <- bind_rows(curve.ls) |>
  mutate(condition.name = factor(condition.name, levels = group.order))

label.df      <- merged.gsea.data |> distinct(condition.name, legend.label) |> arrange(condition.name)
legend.labels <- setNames(label.df$legend.label, as.character(label.df$condition.name))
gsea.colors   <- c("ARID1A" = "#D53E4F", "ARID1B" = "#4575B4", "ARID2" = "#009E73")

p <- ggplot(merged.gsea.data, aes(x = x, y = runningScore, color = condition.name)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  geom_line(linewidth = 0.85, alpha = 0.85) +

  scale_color_manual(values = gsea.colors, breaks = group.order, labels = legend.labels, name = NULL) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 5), expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 5), expand = expansion(mult = c(0.09, 0.07))) +

  labs(
    x = expression(paste("Gene rank (decreasing ", log[2], "FC)")),
    y = "Running enrichment score (ES)"
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(alpha = 1, linewidth = 0.9))) +

  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11, margin = margin(t = 7)),
    axis.title.y = element_text(size = 11, margin = margin(r = 7)),
    axis.line = element_line(color = "black", linewidth = 0.40),
    axis.ticks = element_line(color = "black", linewidth = 0.35),

    panel.grid = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.04, 0.07),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.width = grid::unit(0.65, "cm"),
    legend.key.height = grid::unit(0.44, "cm"),
    legend.text = ggtext::element_markdown(size = 10, color = "black", lineheight = 1.15),
    legend.margin = margin(0, 0, 0, 0),

    plot.margin = margin(t = 8, r = 10, b = 7, l = 8)
  )

pdf("./Proj_CTAs/output/figure.panels/Fig. 3D.pdf", width = 6.5, height = 5)
print(p)
dev.off()

##- Fig. 3e --------------------------------------------------------------------
hepg2.essential.log2FC <- data.table::fread("./Proj_CTAs/input/Perturb-seq_Nadig_2025/HepG2Essential_Log2FoldChange.csv.gz")
log2FC.mat             <- hepg2.essential.log2FC |> column_to_rownames("V1") |> as.matrix()
load("./Proj_CTAs/output/manuscript.revision/pan.sample.CTGs.RData")
CTGs.df.updated <- readRDS("./Proj_CTAs/input/CTGs/CTGs.df.updated.rds")  ## candidates
auto.CTGs       <- CTGs.df.updated |> 
  dplyr::filter(grepl("^\\d", location), symbol %in% pan.sample.CTGs) |> 
  pull(ensembl) |> intersect(rownames(log2FC.mat))  ## 85
load("./Proj_CTAs/input/Perturb-seq_Nadig_2025/hepg2.seu.ctrl.RData")
auto.genes <- seu@assays[["RNA"]]@meta.features |> filter(chr %in% paste0("chr", 1:22)) |> rownames()

## the code was inspired from `AddModuleScore {Seurat}`
CreateRandomGenes <- function (x, feature, nbin, seed = 1) {
  
  set.seed(seed)
  
  data.avg <- Matrix::rowMeans(x)
  data.avg <- data.avg[order(data.avg)]
  data.cut <- ggplot2::cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                                  n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  
  feature   <- intersect(feature, names(data.cut))
  ref.genes <- setdiff(names(data.cut), feature)
  
  nfeature.each.bin  <- as.list(table(data.cut[feature]))
  ref.genes.each.bin <- split(names(data.cut[ref.genes]), data.cut[ref.genes])
  
  random.genes <- lapply(names(nfeature.each.bin), \(i) {
    
    sample(x       = ref.genes.each.bin[[i]], 
           size    = nfeature.each.bin[[i]], 
           replace = FALSE)
  }
  )
  random.genes <- as.character(unlist(random.genes))
  
  return(random.genes)
}
set.seed(1234)
seeds           <- sample.int(1e6, size = 1000, replace = FALSE)
random.genesets <- lapply(seeds, \(i) CreateRandomGenes(seu[["RNA"]]@data[auto.genes, ], feature = auto.CTGs, nbin = 24, seed = i))

load("./Proj_CTAs/output/manuscript.revision/regulator.ranking.RData")
regulator.order <- intersect(colnames(log2FC.mat), top.regulators)
main.regulators <- c("SMARCC1", "KANSL1", "TRRAP", "SETD2")
n.permutations  <- length(random.genesets)

res.ls <- foreach(r = regulator.order) %do% {
  observed.values <- log2FC.mat[auto.CTGs, r]
  median.obs      <- median(observed.values, na.rm = TRUE)
  null.median     <- vapply(
    random.genesets, \(genes) median(log2FC.mat[genes, r], na.rm = TRUE), numeric(1)
  )

  if (!is.finite(median.obs) || any(!is.finite(null.median))) {
    stop("Non-finite median log2FC values for regulator: ", r)
  }

  p.empirical <- (1 + sum(null.median <= median.obs)) / (length(null.median) + 1)

  list(
    null.median = null.median,
    stat = data.frame(
      regulator = r, median.observed = median.obs, median.null = median(null.median, na.rm = TRUE),
      p.empirical = p.empirical, n.genes = length(auto.CTGs),
      n.genes.used = sum(is.finite(observed.values)), n.permutations = n.permutations
    )
  )
} |> setNames(regulator.order)

null.ls <- lapply(res.ls, "[[", "null.median")
stat.df <- bind_rows(lapply(res.ls, "[[", "stat")) |>
  mutate(
    delta  = median.observed - median.null,
    fdr.bh = p.adjust(p.empirical, method = "BH"),
    panel  = ifelse(regulator %in% main.regulators, "Main", "Supplementary")
  ) |>
  arrange(fdr.bh, delta, match(regulator, regulator.order))
permutation.colors <- c("Null fill" = "#DDE5ED", "Null line" = "#6F8192", "Observed" = "#D09A38")
legend.labels      <- c("Expression-matched non-CTG sets", "Autosomal CTGs")
legend.colors      <- setNames(unname(permutation.colors[c("Null line", "Observed")]), legend.labels)

main.x.range  <- range(c(
  unlist(null.ls[main.regulators], use.names = FALSE),
  stat.df$median.observed[match(main.regulators, stat.df$regulator)], 0
))
main.x.span   <- max(diff(main.x.range), 0.001)
main.x.limits <- main.x.range + c(-0.12, 0.08) * main.x.span
main.x.breaks <- scales::breaks_pretty(n = 5)(main.x.limits)

plot.ls <- foreach(r = regulator.order) %do% {
  stat.row <- stat.df |> filter(regulator == r)
  plot.df  <- data.frame(x = null.ls[[r]])

  null.density <- stats::density(
    plot.df$x, adjust = 0.85, from = min(plot.df$x), to = max(plot.df$x)
  )
  density.max <- max(null.density$y)
  delta.label <- sprintf("Delta == %.3f", stat.row$delta)
  fdr.label   <- sprintf("'FDR' == %.3f", stat.row$fdr.bh)

  if (r %in% main.regulators) {
    x.limits <- main.x.limits
    x.breaks <- main.x.breaks

    annotation.labels <- c(
      sprintf("'CTG median' == %.3f", stat.row$median.observed), delta.label, fdr.label
    )
    annotation.y <- density.max * c(1.50, 1.31, 1.12)
    y.upper      <- density.max * 1.66
  } else {
    x.range  <- range(c(plot.df$x, stat.row$median.observed, 0))
    x.span   <- max(diff(x.range), 0.001)
    x.limits <- x.range + c(-0.12, 0.08) * x.span
    x.breaks <- scales::breaks_pretty(n = 4)(x.limits)

    annotation.labels <- c(delta.label, fdr.label)
    annotation.y      <- density.max * c(1.32, 1.13)
    y.upper           <- density.max * 1.49
  }

  annotation.df <- data.frame(
    x = x.limits[2] - 0.025 * diff(x.limits),
    y = annotation.y, label = annotation.labels
  )

  ggplot(plot.df, aes(x = x)) +
    geom_density(
      aes(color = "Expression-matched non-CTG sets"),
      fill = permutation.colors["Null fill"], linewidth = 0.50, alpha = 0.85,
      adjust = 0.85, trim = TRUE, key_glyph = "rect"
    ) +

    ## draw the zero reference above the density fill.
    annotate(
      "segment", x = 0, xend = 0, y = 0, yend = density.max * 1.04,
      color = "black", linetype = "dotted", linewidth = 0.35
    ) +
    geom_segment(
      data = stat.row, aes(x = median.observed, xend = median.observed, color = "Autosomal CTGs"),
      y = 0, yend = density.max * 1.04, inherit.aes = FALSE,
      linewidth = 0.75, linetype = "dashed"
    ) +
    geom_text(
      data = annotation.df, aes(x = x, y = y, label = label), inherit.aes = FALSE,
      parse = TRUE, hjust = 1, vjust = 0.5, size = 3.1, color = "black"
    ) +

    scale_color_manual(values = legend.colors, breaks = legend.labels, name = NULL) +
    scale_x_continuous(breaks = x.breaks, expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 4), expand = expansion(mult = 0)) +
    coord_cartesian(xlim = x.limits, ylim = c(0, y.upper), expand = FALSE) +

    labs(
      x = expression(paste("Median ", log[2], "FC (CRISPRi vs control)")),
      y = "Density", title = bquote(italic(.(r)) ~ "CRISPRi")
    ) +
    guides(color = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +

    theme_classic(base_size = 11) +
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(size = 11.5, face = "plain", hjust = 0.5, margin = margin(b = 5)),
      axis.title.x = element_text(size = 10, margin = margin(t = 6)),
      axis.title.y = element_text(size = 10, margin = margin(r = 6)),
      axis.text = element_text(size = 9, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.30),
      axis.line = element_line(color = "black", linewidth = 0.35),

      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      legend.key.width = grid::unit(0.55, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.margin = margin(t = 3, r = 0, b = 0, l = 0),
      plot.margin = margin(t = 6, r = 7, b = 6, l = 6)
    )
}
names(plot.ls) <- regulator.order

fig.main <- patchwork::wrap_plots(
  plot.ls[main.regulators], nrow = 1, ncol = 4, byrow = TRUE, guides = "collect"
) & theme(legend.position = "bottom")
pdf("./Proj_CTAs/output/figure.panels/Perturb-seq.HepG2.main.pdf", width = 12, height = 3.5)
print(fig.main)
dev.off()

##- Fig. 3f --------------------------------------------------------------------
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
load("./Proj_CTAs/output/manuscript.revision/regulator.ranking.RData")   ## top.regulators
load("./Proj_CTAs/output/4. characteristics of CTGs/peak.anno.ls.RData") ## peak.anno.ls, cistromeDB.data.kept
txdb     <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene.gr  <- GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE)
gene.chr <- data.frame(entrez = names(gene.gr), chr = as.character(GenomicRanges::seqnames(gene.gr)))
gene.chr$symbol <- suppressMessages(
  AnnotationDbi::mapIds(org.Hs.eg.db, gene.chr$entrez, "SYMBOL", "ENTREZID")
)
gene.chr <- gene.chr |> dplyr::filter(!is.na(symbol)) |> dplyr::distinct(symbol, .keep_all = TRUE)

auto.universe  <- gene.chr$symbol[gene.chr$chr %in% paste0("chr", 1:22)]
X.universe     <- gene.chr$symbol[gene.chr$chr == "chrX"]

load("./Proj_CTAs/output/manuscript.revision/pan.sample.CTGs.RData")
CTGs.df.updated <- readRDS("./Proj_CTAs/input/CTGs/CTGs.df.updated.rds")
auto.CTGs       <- CTGs.df.updated |>
  dplyr::filter(grepl("^\\d", location), symbol %in% pan.sample.CTGs) |>
  pull(symbol) |> intersect(auto.universe)
X.CTGs          <- CTGs.df.updated |>
  dplyr::filter(grepl("^X", location),  symbol %in% pan.sample.CTGs) |>
  pull(symbol) |> intersect(X.universe)

auto.bg <- setdiff(auto.universe, pan.sample.CTGs)
X.bg    <- setdiff(X.universe,    pan.sample.CTGs)

n.auto  <- length(auto.CTGs)  ## 387
n.X     <- length(X.CTGs)  ## 28

## (3) per-dataset occupancy

## log2 odds ratio of "promoter bound", autosomal vs X-linked, with a
## Haldane-Anscombe continuity correction. The identical function is applied to
## the observed CTG sets and to every null draw.
log2.OR <- function(a1, n.a, b1, n.b) {
  log2(((a1 + 0.5) / (n.a - a1 + 0.5)) / ((b1 + 0.5) / (n.b - b1 + 0.5)))
}

## hypogeometric null draws of "promoter bound" counts, size-matched to the observed CTG sets
draw.null <- function(B, n.pick, n.bound, n.pool) {
  rhyper(B, m = n.bound, n = n.pool - n.bound, k = n.pick)
}
n.perm      <- 1000
seed.base   <- 1234
std.chr     <- paste0("chr", c(1:22, "X", "Y"))
registerDoParallel(cores = 40)

occupancy.df <- foreach(k = seq_along(peak.anno.ls), .combine = "rbind") %dopar% {

  id  <- names(peak.anno.ls)[k]
  obj <- peak.anno.ls[[k]]
  if (!inherits(obj, "csAnno")) return(NULL)  ## tryCatch() stored an error message

  pk <- as.data.frame(obj)

  ## (a) canonical chromosomes only (the peak files contain scaffolds such as
  ##     chrKI270466.1)
  pk <- pk[as.character(pk$seqnames) %in% std.chr, , drop = FALSE]
  if (nrow(pk) == 0) return(NULL)

  ## (c) genes whose promoter carries at least one peak
  bound <- unique(pk$SYMBOL[grepl("Promoter", pk$annotation, ignore.case = TRUE)])
  bound <- bound[!is.na(bound)]

  A.bound <- sum(auto.CTGs %in% bound)
  B.bound <- sum(X.CTGs    %in% bound)
  obs     <- log2.OR(A.bound, n.auto, B.bound, n.X)

  ## size-matched null
  set.seed(seed.base + k)
  null <- log2.OR(draw.null(n.perm, n.auto, sum(auto.bg %in% bound), length(auto.bg)),
                  n.auto,
                  draw.null(n.perm, n.X,    sum(X.bg    %in% bound), length(X.bg)),
                  n.X)

  p.hi <- (1 + sum(null >= obs)) / (n.perm + 1)
  p.lo <- (1 + sum(null <= obs)) / (n.perm + 1)

  data.frame(
    DCid          = id,
    A.bound       = A.bound,
    B.bound       = B.bound,
    obs.log2.OR   = obs,
    null.median   = median(null),
    delta.log2.OR = obs - median(null),
    p.emp         = min(1, 2 * min(p.hi, p.lo)),
    row.names     = NULL
  )
}

stopImplicitCluster()
solid.tumour.lines <- c(
  "MCF-7", "T47D", "SUM159PT", "MDA-MB-231", "ZR-75-30",
  "breast triple negative tumor cell line BCSC1",
  "A549", "Lu130", "HepG2", "HeLa", "HeLa-S3",
  "LNCaP", "LNCaP-abl", "PC-3", "VCaP", "DU145",
  "LoVo", "LS180", "RKO", "HCT-116",
  "786-O", "G-401",
  "TTC1240", "TTC-549", "BT-16", "U2OS", "SK-N-MC", "HS-ES-2M", "Rh4",
  "T98G", "U87MG", "SF8628", "SF7761",
  "SK-N-SH", "SK-N-SH_RA", "SH-SY5Y", "KELLY", "LAN-1",
  "501-Mel", "WM262", "WM451-lu", "WM853.2",
  "OECM1", "BICR31", "KYSE-150", "ECC-1",
  "MKL-1", "TC-797", "ZNF532-NUT 24335"
)
lost.DCids        <- c(
  c("81895", "81896", "81899", "88642", "88644"),  ## VA-ES-BJ (epithelioid sarcoma)
  c("74196", "74197", "74198", "74199"),  ## SU-DIPG-IV (diffuse intrinsic pontine glioma),
  c("82732", "86643")  ## colorectal adenocarcinoma
)

solid.tumour.idx     <- cistromeDB.data.kept$Cell_line %in% solid.tumour.lines | cistromeDB.data.kept$DCid %in% lost.DCids
cistromeDB.data.used <- cistromeDB.data.kept |> 
  dplyr::filter(Factor %in% top.regulators, solid.tumour.idx) |>
  dplyr::select(DCid, Factor, Cell_line)

## dataset level
ds.df     <- occupancy.df |>
  inner_join(cistromeDB.data.used, by = "DCid") |>
  mutate(
    auto.CTG  = A.bound / n.auto,
    X.CTG     = B.bound / n.X,
    dominance = case_when(auto.CTG > X.CTG ~ "auto",
                          auto.CTG < X.CTG ~ "X",
                          TRUE             ~ "tie")
  )
## factor level
factor.df <- ds.df |>
  group_by(Factor) |>
  group_modify(~ {
    x  <- .x$delta.log2.OR[is.finite(.x$delta.log2.OR)]
    ci <- if (length(x) >= 3) {
      suppressWarnings(wilcox.test(x, mu = 0, conf.int = TRUE, exact = FALSE))$conf.int
    } else c(NA_real_, NA_real_)

    data.frame(delta = median(x), n.datasets = length(x), ci.lo = ci[1], ci.hi = ci[2])
  }) |>
  ungroup() |> arrange(desc(delta))

plot.df1 <- factor.df |> mutate(Factor = factor(Factor, levels = rev(factor.df$Factor)))
point.df <- ds.df |>
  mutate(
    Factor    = factor(Factor, levels = levels(plot.df1$Factor)),
    dominance = factor(dominance, levels = c("auto", "tie", "X"))
  ) |>
  left_join(dplyr::select(factor.df, Factor, n.datasets), by = "Factor")

x.range  <- range(c(point.df$delta.log2.OR, plot.df1$delta, 0), na.rm = TRUE)
x.limits <- x.range + c(-0.08, 0.10) * diff(x.range)
count.x  <- x.limits[2] + 0.045 * diff(x.limits)

## Hodges-Lehmann estimate
ht           <- suppressWarnings(wilcox.test(factor.df$delta, mu = 0, conf.int = TRUE, exact = FALSE))
p.label      <- if (ht$p.value < 0.001) {
  "italic(P) < 0.001"
} else {
  paste0("italic(P) == ", format(signif(ht$p.value, 2), scientific = FALSE, trim = TRUE))
}
overall.label <- sprintf(
  "'Overall HL = %.2f (95%% CI [%.2f, %.2f]);'~%s",
  unname(ht$estimate), ht$conf.int[1], ht$conf.int[2], p.label
)

p <- ggplot(plot.df1, aes(x = delta, y = Factor)) +
  ## overall null-adjusted HL estimate and its 95% confidence interval
  annotate(
    "rect", xmin = ht$conf.int[1], xmax = ht$conf.int[2],
    ymin = -Inf, ymax = Inf, fill = "#D09A38", alpha = 0.12
  ) +
  geom_vline(xintercept = unname(ht$estimate), color = "#B8862E", linewidth = 0.55) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +

  geom_boxplot(
    data = filter(point.df, n.datasets >= 10), aes(x = delta.log2.OR, y = Factor),
    orientation = "y", width = 0.55, fill = "grey94", color = "grey70",
    outlier.shape = NA, linewidth = 0.30, inherit.aes = FALSE
  ) +
  ## dataset colors encode the direction of the raw occupancy difference
  geom_point(
    data = point.df, aes(x = delta.log2.OR, y = Factor, color = dominance),
    position = position_jitter(width = 0, height = 0.14, seed = 1234),
    shape = 16, size = 1.35, alpha = 0.85, inherit.aes = FALSE
  ) +

  ## open circles show factor-level medians of the null-adjusted effects
  geom_point(shape = 21, size = 2.1, fill = "white", color = "black", stroke = 0.35) +
  geom_text(aes(label = n.datasets), x = count.x, hjust = 0.5, size = 2.8, color = "grey35") +
  annotate("text", x = count.x, y = Inf, label = "n", vjust = 0, size = 3, color = "black") +

  scale_color_manual(
    values = c("auto" = "#4575B4", "tie" = "grey70", "X" = "#D53E4F"),
    breaks = c("auto", "X"), labels = c("Autosomal CTGs", "X-linked CTGs"),
    name = "Higher raw occupancy (before null adjustment)"
  ) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 6)) +
  coord_cartesian(xlim = x.limits, clip = "off") +

  labs(x = expression(Delta ~ log[2] ~ "OR"), y = NULL, subtitle = parse(text = overall.label)) +
  guides(
    color = guide_legend(
      nrow = 1, title.position = "top",
      override.aes = list(size = 2.4, alpha = 0.85)
    )
  ) +

  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "sans", color = "black"),
    plot.subtitle = element_text(size = 9.5, hjust = 0, color = "grey25", margin = margin(b = 9)),
    axis.text.y = element_text(size = 9, face = "italic", color = "black", margin = margin(r = 4)),
    axis.text.x = element_text(size = 9.5, color = "black"),
    axis.title.x = element_text(size = 10.5, margin = margin(t = 7)),
    axis.line.x = element_line(color = "black", linewidth = 0.35),
    axis.line.y = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 0.30),
    axis.ticks.y = element_blank(),

    panel.grid.major.y = element_line(color = "grey96", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),

    legend.position = "bottom",
    legend.justification = "left",
    legend.title = element_text(size = 9, face = "plain"),
    legend.text = element_text(size = 9),
    legend.key.width = grid::unit(0.45, "cm"),
    legend.key.height = grid::unit(0.35, "cm"),
    legend.margin = margin(t = 4, r = 0, b = 0, l = 0),

    plot.margin = margin(t = 8, r = 28, b = 6, l = 6)
  )

pdf("./Proj_CTAs/output/figure.panels/Fig. 3E.pdf", width = 6.4, height = max(4.2, 0.19 * nrow(plot.df1) + 1.8))
print(p)
dev.off()

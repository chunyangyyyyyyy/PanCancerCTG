# MaligCTG

MaligCTG scores cells in a Seurat object against malignant-like and normal-like references defined by cancer-testis gene (CTG) burden. It optionally assigns binary malignant labels from both scores.

## Contents

- `MaligCTG/`: the R package, exposing `run.MaligCTG()` and the bundled CTG catalog `CTG.df`.

## Installation

```r
install.packages(c("remotes", "BiocManager"))
BiocManager::install(c("SingleR", "BiocParallel"))
remotes::install_github("chunyangyyyyyyy/PanCancerCTG", subdir = "MaligCTG")
```

Other declared dependencies are Matrix, SeuratObject, mclust, methods, and stats. Seurat is needed separately if you use it to prepare your input object.

## Run

Prepare `seu` (a Seurat object with raw counts) yourself. The 443-gene CTG catalog from the paper is bundled as `CTG.df` and is the default value of `CTGs`; pass another data frame to `CTGs` to use a different gene set.

```r
library(MaligCTG)
result <- run.MaligCTG(seu = seu, CTGs = CTG.df, gene.id = "symbol", assay = "RNA")
head(result)

# Continuous scores only
scores <- run.MaligCTG(seu = seu, CTGs = CTG.df, binary = FALSE)

# Otsu splitting instead of the default GMM
otsu <- run.MaligCTG(seu = seu, CTGs = CTG.df, binary.method = "otsu")
```

## Output

Each row is an analyzed cell. With the default `subset.cells = NULL`, row names and `barcode` follow the raw-count matrix cell barcodes in their original order. Character subsets preserve that order; numeric subsets follow the supplied indices.

| Column | Meaning |
|---|---|
| `barcode` | Cell barcode. |
| `n.CTG` | Number of matching CTG matrix rows with raw counts greater than zero. |
| `reference` | `"high"` or `"low"` for malignant-like or normal-like reference cells; otherwise `NA`. |
| `malignant.score` | SingleR score against the high-CTG reference. |
| `normal.score` | SingleR score against the low-CTG reference. |
| `malignant` | Logical intersection of both malignant-side decisions; included only when `binary = TRUE`. |

SingleR scores are similarity scores, not malignant probabilities. GMM uses `p.threshold = 0.9` for component posterior probabilities. Ties at the reference quantiles are all included. Degenerate scores or failed fits can yield missing labels or errors; two nonempty groups are not guaranteed.

## Input notes

- Use raw counts from the selected assay (`assay = "RNA"` by default).
- Match gene identifiers to the selected CTG-table column (`gene.id = "symbol"` by default, or `"ensembl"`).
- For a Seurat v5 assay with multiple counts layers, join those layers with `SeuratObject::JoinLayers()` first.

## GMM diagnostics

For the default GMM classification, the result also contains a two-row diagnostic table:

```r
gmm.diagnostics <- attr(result, "gmm.diagnostics")
print(gmm.diagnostics)
```

It reports component proportions, means and standard deviations, Ashman's D, the Bimodality Index (BI), the BIC difference between two- and one-component models, and the empirical-versus-fitted CDF distance. These are descriptive diagnostics; they do not apply QC cutoffs or change cell labels. One-component diagnostic-fit failures are recorded in the table. Errors from the original two-component classification fit still propagate.

BI is computed as `sqrt(pro.low * pro.high) * ashman.D`. For the equal-variance E model, it uses the common component standard deviation. For the unequal-variance V model, it is an extension using the square root of the mean component variance. The package does not apply a BI cutoff.

Diagnostics describe all analyzed cells in that call. For per-sample diagnostics, call the function separately for each sample. Subsetting a returned result does not refit its diagnostics. Use `saveRDS(result, "MaligCTG.result.rds")` to retain the attribute, or export the diagnostic table separately; CSV output of the cell table does not retain it. Otsu and continuous-score-only results do not have this attribute.

## Practical interpretation

Use the continuous malignant score together with expression-based clustering and cell-identity evidence to annotate candidate malignant populations. Inspect score distributions within clusters and their sizes relative to the analyzed cell set. Small clusters with consistently high scores may represent rare malignant populations and warrant closer review of GMM or Otsu labels. Binary labels are auxiliary outputs and should be interpreted alongside these observations.

## License

MIT. See [`LICENSE`](LICENSE).

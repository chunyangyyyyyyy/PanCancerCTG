# MaligCTG

MaligCTG scores cells in a Seurat object against malignant-like and normal-like references defined by cancer-testis gene (CTG) burden. It optionally assigns binary malignant labels from both scores.

## Contents

- `MaligCTG/`: the R package, exposing `run.MaligCTG()` and the bundled CTG catalog `CTG.df`.
- `code/code.for.analysis/FindTumorLayers.r`: a separate spatial analysis function.
- `code/code.for.figures/`: figure scripts retain the original project data paths; required input data are not bundled.

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
- Match gene identifiers to the selected CTG-table column (`gene.id = "symbol"` by default, or `"ensembl"`). Ensembl version suffixes are stripped from matrix row names only; supply unversioned IDs in the CTG table. Matching is case-sensitive, and duplicate matching rows are counted separately.
- For a Seurat v5 assay with multiple counts layers, join those layers with `SeuratObject::JoinLayers()` first.

## License

MIT. See [`LICENSE`](LICENSE).

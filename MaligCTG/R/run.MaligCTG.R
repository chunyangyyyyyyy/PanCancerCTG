#' Annotate malignant cells by CTG burden
#'
#' Score cells against CTG-derived malignant-like and normal-like references, and optionally assign binary malignant labels.
#'
#' @param seu A Seurat object with raw counts in the assay selected by \code{assay}. For Seurat v5, join multiple counts layers before calling this function.
#' @param CTGs A data.frame containing the column selected by \code{gene.id}, normally both \code{symbol} and \code{ensembl}. Default: \code{\link{CTG.df}}, the 443-gene CTG catalog bundled with this package. Ensembl IDs in this table must be unversioned.
#' @param gene.id Gene identifier type: \code{"symbol"} (default) or \code{"ensembl"}. For Ensembl matching, numeric version suffixes are removed from matrix row names only.
#' @param assay Assay containing raw counts. Default: \code{"RNA"}.
#' @param subset.cells Character barcodes, numeric indices, or a logical vector selecting cells to analyze. Default: \code{NULL} (all cells). Character selections retain matrix order; numeric selections follow index order.
#' @param normal.cells Character barcodes, numeric indices, or a logical vector selecting normal reference cells, interpreted against all cells of \code{seu}. These cells must be contained in \code{subset.cells}. Default: \code{NULL} (infer the low-CTG reference).
#' @param prop Quantile-tail proportion in \eqn{(0, 0.5)}. Default: \code{0.05}. All cells tied at a threshold are included, so a reference need not contain exactly 5 percent of cells.
#' @param binary Whether to append binary malignant labels. Default: \code{TRUE}.
#' @param binary.method Binary splitting method: \code{"gmm"} (default) or \code{"otsu"}.
#' @param p.threshold GMM component posterior probability threshold, compared strictly using \code{>}. Default: \code{0.9}. It is validated for every call but used for splitting only by GMM.
#' @param seed Seed passed to \code{set.seed()} before each GMM fit. Default: \code{1234}. With \code{NULL}, no seed is set. Otsu does not use this parameter.
#' @param BPPARAM Parallelization parameter passed to SingleR. Default: \code{BiocParallel::SerialParam()}.
#'
#' @details
#' CTG burden counts matching matrix rows whose raw count is greater than zero. Matching is case-sensitive. In Ensembl mode, a trailing dot followed by digits is stripped from matrix row names, but not from the CTG table. Matching rows are not aggregated after suffix removal.
#'
#' With \code{normal.cells = NULL}, reference cells have burdens at or below the \code{prop} quantile and at or above the \code{1 - prop} quantile. All threshold ties are included. When normal cells are supplied, they form the low reference, and the high-reference quantile is computed from the remaining analyzed cells. Overlapping high and low references cause an error; fewer than 10 cells in either reference cause a warning.
#'
#' Each pseudo-bulk profile is \code{log2(1e6 * x / sum(x) + 1)}, where \code{x} is the gene-wise sum of raw reference counts. SingleR receives raw test counts and these two profiles, with \code{fine.tune = FALSE}, \code{prune = FALSE}, and the supplied \code{BPPARAM}. The returned scores are extracted directly from \code{pred$scores}.
#'
#' GMM fits two components separately to each score. The high-mean malignant-score component and low-mean normal-score component must each exceed \code{p.threshold}. If Mclust returns \code{NULL}, the helper warns and returns missing side labels. Errors thrown by Mclust are not caught, including possible failures for degenerate scores.
#'
#' Otsu searches midpoints between consecutive distinct sorted scores, maximizing between-class variance. The first maximizing split is used. A constant score causes a warning and missing side labels. Side labels are combined with R logical \code{&}, so a missing side combined with \code{FALSE} gives \code{FALSE}. Neither method guarantees two nonempty output groups.
#'
#' @return A data.frame with one row per analyzed cell. Row names and
#' \code{barcode} follow the selected raw-count matrix columns (all cells in
#' matrix order by default). Columns:
#' \describe{
#' \item{\code{barcode}}{Character cell barcode.}
#' \item{\code{n.CTG}}{Integer number of matching CTG matrix rows with counts greater than zero.}
#' \item{\code{reference}}{\code{"high"} for malignant-like reference cells, \code{"low"} for normal-like reference cells, otherwise \code{NA}.}
#' \item{\code{malignant.score}}{SingleR score from \code{pred$scores[, "high"]}; not a malignant probability.}
#' \item{\code{normal.score}}{SingleR score from \code{pred$scores[, "low"]}; not a normal probability.}
#' \item{\code{malignant}}{Logical intersection of the malignant sides of both scores; present only when \code{binary = TRUE}. May contain \code{NA}.}
#' }
#'
#' @author Chunyang Fu
#' @examples
#' \dontrun{
#' # seu: a Seurat object with raw RNA counts, prepared by the user.
#' # CTG.df: the CTG catalog bundled with this package (also the default of CTGs).
#' result <- run.MaligCTG(seu = seu, CTGs = CTG.df,
#'                       gene.id = "symbol", assay = "RNA")
#' head(result)
#' scores <- run.MaligCTG(seu = seu, CTGs = CTG.df, binary = FALSE)
#' otsu <- run.MaligCTG(seu = seu, CTGs = CTG.df, binary.method = "otsu")
#' }
#' @seealso \code{\link[SingleR]{SingleR}}, \code{\link[mclust]{Mclust}}
#' @export
#' @importFrom Matrix colSums rowSums
#' @importFrom BiocParallel SerialParam
#' @importFrom SingleR SingleR
#' @importFrom stats quantile
run.MaligCTG <- function(seu, CTGs = MaligCTG::CTG.df, gene.id = c("symbol", "ensembl"), assay = "RNA",
                         subset.cells = NULL, normal.cells = NULL, prop = 0.05,
                         binary = TRUE, binary.method = c("gmm", "otsu"), p.threshold = 0.9,
                         seed = 1234, BPPARAM = BiocParallel::SerialParam())
{
    ## checking inputs
    if (!inherits(seu, "Seurat")) {
        stop("'seu' must be a Seurat object")
    }
    gene.id       <- match.arg(gene.id)
    binary.method <- match.arg(binary.method)
    if (!is.data.frame(CTGs) || !gene.id %in% colnames(CTGs)) {
        stop(sprintf("'CTGs' must be a data.frame with a '%s' column", gene.id))
    }
    if (!is.numeric(prop) || length(prop) != 1L || prop <= 0 || prop >= 0.5) {
        stop("'prop' must be a number in (0, 0.5)")
    }
    if (!isTRUE(binary) && !isFALSE(binary)) {
        stop("'binary' must be TRUE or FALSE")
    }
    if (!is.numeric(p.threshold) || length(p.threshold) != 1L || p.threshold <= 0 || p.threshold >= 1) {
        stop("'p.threshold' must be a number in (0, 1)")
    }

    ## raw counts of the analyzed cells
    counts <- .get_counts(seu, assay)
    use    <- .to_cells(subset.cells, colnames(counts), "subset.cells")
    counts <- counts[, use, drop = FALSE]

    ## CTGs present in the matrix
    ids   <- rownames(counts)
    if (gene.id == "ensembl") {
        ids <- sub("\\.[0-9]+$", "", ids)
    }
    keep  <- ids %in% CTGs[[gene.id]]
    n.hit <- length(unique(ids[keep]))
    n.all <- length(unique(CTGs[[gene.id]]))
    if (n.hit == 0L) {
        stop(sprintf("no CTGs found in the row names of 'seu'; is 'gene.id = \"%s\"' correct?", gene.id))
    }
    if (n.hit < 0.5 * n.all) {
        warning(sprintf("only %i of %i CTGs found in the row names of 'seu'; is 'gene.id = \"%s\"' correct?", n.hit, n.all, gene.id))
    }
    n.CTG <- as.integer(Matrix::colSums(counts[keep, , drop = FALSE] > 0))

    ## internal references
    if (is.null(normal.cells)) {
        low  <- n.CTG <= stats::quantile(n.CTG, prop)
        high <- n.CTG >= stats::quantile(n.CTG, 1 - prop)
    } else {
        normal <- .to_cells(normal.cells, colnames(seu), "normal.cells")
        if (!all(normal %in% use)) {
            stop("'normal.cells' must be contained in 'subset.cells'")
        }
        low  <- use %in% normal
        high <- !low & n.CTG >= stats::quantile(n.CTG[!low], 1 - prop)
    }
    if (any(high & low)) {
        stop("CTG counts do not vary enough to separate the high and low references")
    }
    if (sum(high) < 10L || sum(low) < 10L) {
        warning(sprintf("few reference cells (high: %i, low: %i)", sum(high), sum(low)))
    }

    ref <- cbind(high = .pseudo_bulk(counts[, high, drop = FALSE]),
                 low  = .pseudo_bulk(counts[, low, drop = FALSE]))
    pred <- SingleR::SingleR(test = counts, ref = ref, labels = c("high", "low"),
                             fine.tune = FALSE, prune = FALSE, BPPARAM = BPPARAM)

    output <- data.frame(
        barcode         = use,
        n.CTG           = n.CTG,
        reference       = ifelse(high, "high", ifelse(low, "low", NA_character_)),
        malignant.score = unname(pred$scores[, "high"]),
        normal.score    = unname(pred$scores[, "low"]),
        row.names       = use
    )

    ## binary label: malignant side of both scores
    if (binary) {
        side.FUN <- switch(binary.method, gmm = .gmm_side, otsu = .otsu_side)
        output$malignant <- side.FUN(output$malignant.score, upper = TRUE, p.threshold = p.threshold, seed = seed) &
                            side.FUN(output$normal.score, upper = FALSE, p.threshold = p.threshold, seed = seed)
    }

    output
}


## internal functions -----------------------------------------------------------------------
.get_counts <- function(seu, assay) {
    if (!assay %in% names(seu@assays)) {
        stop(sprintf("assay '%s' not found in 'seu'", assay))
    }
    assay.obj <- seu[[assay]]
    if (inherits(assay.obj, "Assay5")) {
        if (length(SeuratObject::Layers(assay.obj, search = "^counts")) > 1L) {
            stop("multiple 'counts' layers found; run 'JoinLayers()' first")
        }
        SeuratObject::LayerData(assay.obj, layer = "counts")
    } else {
        methods::slot(assay.obj, "counts")
    }
}

.to_cells <- function(i, cells, arg) {
    if (is.null(i)) {
        return(cells)
    }
    if (is.logical(i)) {
        if (length(i) != length(cells) || anyNA(i)) {
            stop(sprintf("logical '%s' must have one non-NA value per cell", arg))
        }
        return(cells[i])
    }
    if (is.numeric(i)) {
        if (anyNA(i) || any(i < 1 | i > length(cells))) {
            stop(sprintf("integer '%s' out of range", arg))
        }
        return(cells[i])
    }
    if (is.character(i)) {
        if (!all(i %in% cells)) {
            stop(sprintf("some cells in '%s' not found in 'seu'", arg))
        }
        return(cells[cells %in% i])
    }
    stop(sprintf("'%s' must be a character, integer or logical vector", arg))
}

.pseudo_bulk <- function(counts) {
    x <- Matrix::rowSums(counts)
    log2(1e6 * x / sum(x) + 1)
}

## TRUE for cells on the malignant side: high-mean component if 'upper', low-mean component otherwise
.gmm_side <- function(x, upper, p.threshold, seed) {
    mclustBIC <- mclust::mclustBIC    # Mclust() looks up mclustBIC() in the calling frame
    if (!is.null(seed)) {
        set.seed(seed)
    }
    fit <- mclust::Mclust(x, G = 2, verbose = FALSE)
    if (is.null(fit)) {
        warning("two-component mixture model could not be fitted; returning NA labels")
        return(rep(NA, length(x)))
    }
    comp <- if (upper) which.max(fit$parameters$mean) else which.min(fit$parameters$mean)
    unname(fit$z[, comp] > p.threshold)
}

## Otsu's threshold: maximize the between-class variance w0 * w1 * (mu0 - mu1)^2 over all splits of the sorted scores
.otsu_side <- function(x, upper, ...) {
    s     <- sort(x)
    n     <- length(s)
    k     <- which(diff(s) > 0)            # split after position k, only between distinct values
    if (length(k) == 0L) {
        warning("score is constant and cannot be split; returning NA labels")
        return(rep(NA, length(x)))
    }
    cs     <- cumsum(s)
    mu0    <- cs[k] / k
    mu1    <- (cs[n] - cs[k]) / (n - k)
    sb     <- (k / n) * (1 - k / n) * (mu0 - mu1)^2
    best   <- k[which.max(sb)]
    cutoff <- (s[best] + s[best + 1L]) / 2
    if (upper) x > cutoff else x < cutoff
}
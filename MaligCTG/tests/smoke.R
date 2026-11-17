# Integration checks with synthetic raw counts; no downloads or external data.
args <- commandArgs(trailingOnly = TRUE)
mode <- match.arg(if (length(args)) args[1L] else "namespaced",
                  c("namespaced", "attached"))
stopifnot(!"package:MaligCTG" %in% search())
catalog <- MaligCTG::CTG.df
stopifnot(is.data.frame(catalog), nrow(catalog) == 443L,
          all(c("symbol", "ensembl") %in% names(catalog)),
          !anyDuplicated(catalog$symbol), !anyDuplicated(catalog$ensembl))
set.seed(20261117)
ng <- 1000L
nc <- 240L
nctg <- nrow(catalog)
genes <- c(catalog$symbol, paste0("BACKGROUND_", seq_len(ng - nctg)))
cells <- paste0("cell_", seq_len(nc))
counts <- matrix(rpois(ng * nc, 2), nrow = ng,
                 dimnames = list(genes, cells))
p <- c(runif(nc / 2L, 0.02, 0.18), runif(nc / 2L, 0.55, 0.95))
counts[seq_len(nctg), ] <-
    matrix(rbinom(nctg * nc, 1L, rep(p, each = nctg)), nrow = nctg) *
    (1L + matrix(rpois(nctg * nc, 3), nrow = nctg))
counts[444:643, 121:240] <- matrix(rpois(200L * 120L, 9), nrow = 200L)
counts[644:843, 1:120] <- matrix(rpois(200L * 120L, 9), nrow = 200L)
seu <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    min.cells = 0, min.features = 0
)
# A workspace decoy must never replace the bundled default catalog.
CTG.df <- data.frame(symbol = "NOT_A_CTG", ensembl = "NOT_AN_ENSEMBL_ID")
if (mode == "attached") {
    library(MaligCTG)
    run <- run.MaligCTG
} else {
    run <- MaligCTG::run.MaligCTG
}
base.columns <- c("barcode", "n.CTG", "reference", "malignant.score", "normal.score")
settings <- list(
    gmm = list(binary = TRUE, binary.method = "gmm"),
    otsu = list(binary = TRUE, binary.method = "otsu"),
    scores = list(binary = FALSE)
)
for (setting in names(settings)) {
    call.args <- c(list(seu = seu, seed = 1234,
                        BPPARAM = BiocParallel::SerialParam()), settings[[setting]])
    a <- do.call(run, call.args)
    b <- do.call(run, call.args)
    explicit <- do.call(run, c(call.args, list(CTGs = catalog)))
    expected <- if (setting == "scores") base.columns else c(base.columns, "malignant")
    stopifnot(identical(names(a), expected),
              identical(rownames(a), colnames(seu)),
              identical(a$barcode, colnames(seu)),
              identical(a, b), identical(a, explicit),
              is.integer(a$n.CTG),
              identical(a$n.CTG, as.integer(colSums(counts[seq_len(nctg), ] > 0))),
              all(is.finite(a$malignant.score)), all(is.finite(a$normal.score)))
    if (setting != "scores") stopifnot(is.logical(a$malignant))
    cat("PASS:", mode, setting, "columns, order, repeatability and default catalog\n")
}
if (mode == "namespaced") stopifnot(!"package:MaligCTG" %in% search())
stopifnot(identical(getNamespaceExports("MaligCTG"), "run.MaligCTG"))
cat("PASS: bundled catalog, public exports and workspace isolation\n")

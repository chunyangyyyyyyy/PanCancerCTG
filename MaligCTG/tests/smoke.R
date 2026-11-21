# Integration checks with synthetic raw counts; no downloads or external data.
args <- commandArgs(trailingOnly = TRUE)
mode <- match.arg(if (length(args)) args[1L] else "namespaced",
                  c("namespaced", "attached"))
stopifnot(!"package:MaligCTG" %in% search())
baseline.file <- if (file.exists("helpers/gmm-baseline.R")) {
    "helpers/gmm-baseline.R"
} else {
    file.path("MaligCTG", "tests", "helpers", "gmm-baseline.R")
}
source(baseline.file, local = TRUE)
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
    d <- attr(a, "gmm.diagnostics")
    if (setting == "gmm") {
        expected.diagnostics <- c(
            "score", "n.total", "n.valid", "fit.available", "model",
            "pro.low", "pro.high", "mean.low", "mean.high", "sd.low", "sd.high",
            "ashman.D", "BIC.1", "BIC.2", "delta.BIC", "ks.D",
            "fit1.available", "fit1.error", "fit1.warning",
            "diagnostics.error", "diagnostics.warning"
        )
        stopifnot(is.data.frame(d), nrow(d) == 2L,
                  identical(names(d), expected.diagnostics),
                  identical(d$score, c("malignant.score", "normal.score")),
                  all(d$n.total == nc), all(d$n.valid == nc),
                  all(d$fit.available), all(d$fit1.available),
                  all(is.na(d$diagnostics.error)),
                  all(is.finite(d$ashman.D)), all(d$ks.D >= 0 & d$ks.D <= 1),
                  isTRUE(all.equal(d$delta.BIC, d$BIC.2 - d$BIC.1)))
    } else {
        stopifnot(is.null(d))
    }
    cat("PASS:", mode, setting, "columns, order, repeatability and default catalog\n")
}
if (mode == "namespaced") stopifnot(!"package:MaligCTG" %in% search())
stopifnot(identical(getNamespaceExports("MaligCTG"), "run.MaligCTG"))
cat("PASS: bundled catalog, public exports and workspace isolation\n")

# Reuse the original classifier with the same scoring/reference pipeline.
# Compare columns, not the new diagnostic attribute.
baseline.run <- run
baseline.env <- new.env(parent = environment(baseline.run))
baseline.env$.gmm_side <- .baseline_gmm_side
environment(baseline.run) <- baseline.env
for (classification.seed in list(1234, NULL)) {
    set.seed(9981)
    previous <- baseline.run(seu = seu, seed = classification.seed)
    previous.rng <- .Random.seed
    set.seed(9981)
    current <- run(seu = seu, seed = classification.seed)
    current.rng <- .Random.seed
    for (column in names(previous)) {
        stopifnot(identical(previous[[column]], current[[column]]))
    }
    stopifnot(identical(previous.rng, current.rng))
}
cat("PASS:", mode, "all original columns and post-call RNG match the original classifier\n")

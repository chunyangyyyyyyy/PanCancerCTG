baseline.file <- if (file.exists("helpers/gmm-baseline.R")) {
    "helpers/gmm-baseline.R"
} else {
    file.path("MaligCTG", "tests", "helpers", "gmm-baseline.R")
}
source(baseline.file, local = TRUE)

library(MaligCTG)
current.side <- getFromNamespace(".gmm_side", "MaligCTG")
diagnose <- getFromNamespace(".gmm_diagnostics", "MaligCTG")
fit.one <- getFromNamespace(".gmm_fit_one", "MaligCTG")
fit.model <- function(x, G, modelNames = NULL) {
    mclustBIC <- mclust::mclustBIC
    mclust::Mclust(x, G = G, modelNames = modelNames, verbose = FALSE)
}
rng.state <- function() {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else NULL
}
capture.side <- function(fun, x, seed, upper) {
    set.seed(717)
    warnings <- character()
    value <- tryCatch(
        withCallingHandlers(
            fun(x, upper = upper, p.threshold = 0.9, seed = seed),
            warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        ), error = function(e) e
    )
    list(value = value, warnings = warnings, rng = rng.state())
}
assert.same <- function(old, new) {
    stopifnot(identical(old$warnings, new$warnings), identical(old$rng, new$rng))
    if (inherits(old$value, "error")) {
        stopifnot(inherits(new$value, "error"),
                  identical(class(old$value), class(new$value)),
                  identical(conditionMessage(old$value), conditionMessage(new$value)))
    } else {
        stopifnot(!inherits(new$value, "error"),
                  identical(as.logical(old$value), as.logical(new$value)))
    }
}
set.seed(28)
x <- c(rnorm(300, -0.4, 0.08), rnorm(300, 0.4, 0.16))
cases <- list(
    separated = x,
    unbalanced = c(rnorm(970, -0.4, 0.1), rnorm(30, 0.5, 0.1)),
    constant = rep(0.25, 50),
    nonfinite = c(x, NA_real_, Inf),
    large = c(rnorm(1200, -0.4, 0.1), rnorm(1200, 0.4, 0.1))
)
for (case in names(cases)) {
    for (seed in list(1234, NULL)) {
        for (upper in c(TRUE, FALSE)) {
            old <- capture.side(.baseline_gmm_side, cases[[case]], seed, upper)
            new <- capture.side(current.side, cases[[case]], seed, upper)
            assert.same(old, new)
        }
    }
}
cat("PASS: original side labels, warnings, errors and RNG preserved across boundary cases\n")

# Both supported univariate variance structures and independent KS statistic.
for (model in c("E", "V")) {
    fit <- fit.model(x, G = 2, modelNames = model)
    d <- diagnose(x, fit, seed = 1234)
    stopifnot(d$fit.available, d$fit1.available, identical(d$model, model),
              is.na(d$diagnostics.error), d$mean.low <= d$mean.high,
              d$sd.low > 0, d$sd.high > 0,
              isTRUE(all.equal(d$pro.low + d$pro.high, 1)),
              isTRUE(all.equal(d$BIC.2, unname(fit$bic))),
              isTRUE(all.equal(d$BIC.1, unname(fit.one(x)$bic))))
    mixture.cdf <- function(q) {
        d$pro.low * pnorm(q, d$mean.low, d$sd.low) +
            d$pro.high * pnorm(q, d$mean.high, d$sd.high)
    }
    # Only compare D, not the uncalibrated P value for fitted parameters.
    ks <- suppressWarnings(stats::ks.test(x, mixture.cdf))
    stopifnot(isTRUE(all.equal(d$ks.D, unname(ks$statistic), tolerance = 1e-12)))
}
cat("PASS: E/V component diagnostics, BIC values and KS distance\n")

# Inject failures through a private-function copy; never change a namespace.
fit <- fit.model(x, G = 2)
clone.diagnostic <- function(replacement) {
    fun <- diagnose
    env <- new.env(parent = environment(fun))
    env$.gmm_fit_one <- replacement
    environment(fun) <- env
    fun
}
broken <- clone.diagnostic(function(x) {
    runif(3)
    stop("injected one-component failure")
})
set.seed(19)
before <- rng.state()
d <- broken(x, fit, seed = NULL)
stopifnot(!d$fit1.available, is.na(d$delta.BIC),
          identical(d$fit1.error, "injected one-component failure"),
          identical(before, rng.state()))
# Also preserve absence of .Random.seed, including an error after random draws.
rm(".Random.seed", envir = .GlobalEnv)
d <- broken(x, fit, seed = 1234)
stopifnot(is.null(rng.state()), !d$fit1.available)
one <- fit.one(x)
warning.fit <- local({
    saved <- one
    function(x) {
        runif(2)
        warning("injected diagnostic warning")
        saved
    }
})
warned <- clone.diagnostic(warning.fit)
set.seed(31)
before <- rng.state()
outer.warnings <- character()
d <- withCallingHandlers(warned(x, fit, seed = NULL), warning = function(w) {
    outer.warnings <<- c(outer.warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
})
stopifnot(d$fit1.available, identical(d$fit1.warning, "injected diagnostic warning"),
          !length(outer.warnings), identical(before, rng.state()))
for (unusable in list(NULL, list(bic = NA_real_))) {
    fake <- local({saved <- unusable; function(x) saved})
    d <- clone.diagnostic(fake)(x, fit, seed = 1234)
    stopifnot(!d$fit1.available, is.na(d$delta.BIC), !is.na(d$fit1.error))
}
cat("PASS: diagnostic failure/warning recording and exact RNG restoration\n")

# Bad diagnostic metadata must not gate otherwise usable classification.
corrupt.diagnostics <- function(x, fit, seed) {
    fit$parameters$variance$sigmasq <- NULL
    diagnose(x, fit, seed)
}
for (replacement in list(broken, corrupt.diagnostics)) {
    side <- current.side
    env <- new.env(parent = environment(side))
    env$.gmm_diagnostics <- replacement
    environment(side) <- env
    old <- capture.side(.baseline_gmm_side, x, seed = NULL, upper = TRUE)
    new <- capture.side(side, x, seed = NULL, upper = TRUE)
    assert.same(old, new)
    d <- attr(new$value, "gmm.diagnostics")
    stopifnot(!is.null(d),
              !is.na(d$fit1.error) || !is.na(d$diagnostics.error))
}
d <- diagnose(c(1, NA, Inf), NULL, seed = NULL)
stopifnot(!d$fit.available, d$n.total == 3L, d$n.valid == 1L,
          is.na(d$fit1.available), !is.na(d$diagnostics.error))
stopifnot(!any(c("BI", "qc.pass") %in% names(d)))
cat("PASS: diagnostic availability is separate from cell classification\n")

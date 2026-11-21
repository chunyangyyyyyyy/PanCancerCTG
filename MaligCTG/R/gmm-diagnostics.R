# Descriptive diagnostics only: classification is completed before this runs.
.gmm_diagnostics <- function(x, fit, seed) {
    result <- data.frame(
        n.total = length(x), n.valid = sum(is.finite(x)),
        fit.available = !is.null(fit), model = NA_character_,
        pro.low = NA_real_, pro.high = NA_real_,
        mean.low = NA_real_, mean.high = NA_real_,
        sd.low = NA_real_, sd.high = NA_real_, ashman.D = NA_real_,
        BI = NA_real_,
        BIC.1 = NA_real_, BIC.2 = NA_real_, delta.BIC = NA_real_,
        ks.D = NA_real_, fit1.available = NA,
        fit1.error = NA_character_, fit1.warning = NA_character_,
        diagnostics.error = NA_character_, diagnostics.warning = NA_character_,
        stringsAsFactors = FALSE
    )

    # Restore exactly the state after the original G = 2 classification fit.
    rng.exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (rng.exists) {
        rng.state <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
        if (rng.exists) {
            assign(".Random.seed", rng.state, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)

    diagnostic.warnings <- character()
    tryCatch(
        withCallingHandlers({
            if (is.null(fit)) {
                result$diagnostics.error <- "Two-component fit is unavailable."
                return(result)
            }
            result$model <- as.character(fit$modelName)
            if (length(fit$bic) == 1L && is.finite(fit$bic)) {
                result$BIC.2 <- as.numeric(fit$bic)
            }
            mu <- fit$parameters$mean
            pro <- fit$parameters$pro
            vars <- fit$parameters$variance$sigmasq
            valid <- length(x) > 0L && all(is.finite(x)) &&
                length(mu) == 2L && all(is.finite(mu)) &&
                length(pro) == 2L && all(is.finite(pro)) && all(pro > 0) &&
                length(vars) %in% c(1L, 2L) &&
                all(is.finite(vars)) && all(vars > 0)
            if (!valid) {
                stop("Two-component parameters or scores are unusable for diagnostics.")
            }
            if (length(vars) == 1L) vars <- rep(vars, 2L)
            ord <- order(mu)
            result$pro.low <- pro[ord[1L]]
            result$pro.high <- pro[ord[2L]]
            result$mean.low <- mu[ord[1L]]
            result$mean.high <- mu[ord[2L]]
            result$sd.low <- sqrt(vars[ord[1L]])
            result$sd.high <- sqrt(vars[ord[2L]])
            result$ashman.D <- abs(mu[1L] - mu[2L]) / sqrt(mean(vars))
            result$BI <- sqrt(pro[1L] * pro[2L]) * result$ashman.D
            sorted <- sort(x)
            n <- length(sorted)
            cdf <- pro[1L] * stats::pnorm(sorted, mu[1L], sqrt(vars[1L])) +
                pro[2L] * stats::pnorm(sorted, mu[2L], sqrt(vars[2L]))
            result$ks.D <- max(abs(seq_len(n) / n - cdf),
                               abs((seq_len(n) - 1L) / n - cdf))

            # A failed diagnostic fit must never discard existing labels.
            fit1.error <- NULL
            fit1.warnings <- character()
            if (!is.null(seed)) set.seed(seed)
            fit1 <- tryCatch(
                withCallingHandlers(
                    .gmm_fit_one(x),
                    warning = function(w) {
                        fit1.warnings <<- c(fit1.warnings, conditionMessage(w))
                        invokeRestart("muffleWarning")
                    }
                ),
                error = function(e) {
                    fit1.error <<- conditionMessage(e)
                    NULL
                }
            )
            if (length(fit1.warnings)) {
                result$fit1.warning <- paste(unique(fit1.warnings), collapse = "; ")
            }
            result$fit1.available <- !is.null(fit1) &&
                length(fit1$bic) == 1L && is.finite(fit1$bic)
            if (result$fit1.available) {
                result$BIC.1 <- as.numeric(fit1$bic)
                if (is.finite(result$BIC.2)) {
                    result$delta.BIC <- result$BIC.2 - result$BIC.1
                }
            } else {
                result$fit1.error <- if (is.null(fit1.error)) {
                    "One-component fit did not supply a finite BIC."
                } else fit1.error
            }
            result
        }, warning = function(w) {
            diagnostic.warnings <<- c(diagnostic.warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }),
        error = function(e) {
            result$diagnostics.error <- conditionMessage(e)
            result
        }
    ) -> result
    if (length(diagnostic.warnings)) {
        result$diagnostics.warning <- paste(unique(diagnostic.warnings), collapse = "; ")
    }
    result
}

# Kept separate so diagnostic failures can be tested without mocking mclust.
.gmm_fit_one <- function(x) {
    mclustBIC <- mclust::mclustBIC
    mclust::Mclust(x, G = 1, verbose = FALSE)
}

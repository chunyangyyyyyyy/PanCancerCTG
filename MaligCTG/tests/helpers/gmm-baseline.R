# Frozen classifier from f3bc10746fc3b558c77923c4fc8c8fa1f82452ed.
.baseline_gmm_side <- function(x, upper, p.threshold, seed) {
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

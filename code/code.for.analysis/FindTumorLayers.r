library(Seurat)
library(FNN)

FindTumorLayers <- function(object,
                            k                 = 6,
                            n.seed            = 2,  ## c(2, 1) will only use the normal refrence cluster as a seed (similar to "Cottrazm")
                            min.seed.size     = 1,
                            exclude.reference = FALSE,
                            reference.level   = "normal",
                            max.nn.dist       = "auto",  ## "auto" will drop over-distant neighbours according to the value of k
                            vote.margin       = 1,
                            max.iter          = 20) {

  ## ---- seed clusters ------------------------------------------------------
  cnv.score.by.cluster <- tapply(object$cnv_score, object$cnv_cluster, median, na.rm = TRUE)
  ord     <- sort(cnv.score.by.cluster, decreasing = TRUE, na.last = NA)
  cl.size <- table(object$cnv_cluster)
  ord     <- ord[cl.size[names(ord)] >= min.seed.size]
  if (exclude.reference) ord <- ord[names(ord) != reference.level]

  n.seed <- rep_len(n.seed, 2)            # c(malignant, normal)
  if (sum(n.seed) > length(ord))
    stop("Only ", length(ord), " usable CNV clusters: n.seed = c(",
         n.seed[1], ", ", n.seed[2], ") cannot be disjoint.")
  malignant.clusters <- names(ord)[seq_len(n.seed[1])]
  normal.clusters    <- rev(names(ord))[seq_len(n.seed[2])]

  current.labels <- object$cnv_cluster
  states         <- rep("Unknown", length(current.labels))
  names(states)  <- rownames(object@meta.data)
  states[current.labels %in% malignant.clusters] <- "Malignant"
  states[current.labels %in% normal.clusters]    <- "Normal"

  ## ---- initialize layer ---------------------------------------------------
  layers        <- rep(NA, length(states))
  names(layers) <- names(states)
  layers[states == "Malignant"] <- 0
  layers[states == "Normal"]    <- 0

  ## ---- position matrix ----------------------------------------------------
  images  <- Images(object)
  idx     <- if (length(images) == 1) 1 else {sapply(images, \(x) sum(rownames(object@images[[x]]@coordinates) %in% colnames(object))) |> which.max()}
  coords  <- GetTissueCoordinates(object, image = images[idx])
  if (ncol(coords) >= 2) coords <- coords[, 1:2]
  if (setequal(rownames(coords), rownames(object@meta.data))) {
    coords <- coords[rownames(object@meta.data), ]
  } else {
    stop("Error: Row names of coordinates do not match Seurat object. Please ensure they contain the exact same set of cells.")
  }

  ## Visium spot centres form a triangular lattice, so neighbours occur in
  ## discrete coordination shells. In units of the spot spacing d (estimated
  ## per slide as the median first-nearest-neighbour distance):
  ##
  ##   shell  radius      spots   cumulative k   gap to next shell
  ##     1    1.000 d       6           6        (1.000, 1.732)
  ##     2    1.732 d       6          12        (1.732, 2.000)   sqrt(3) d
  ##     3    2.000 d       6          18        (2.000, 2.646)
  ##     4    2.646 d      12          30        (2.646, 3.000)   sqrt(7) d
  ##     5    3.000 d       6          36        (3.000, 3.464)

  SHELL.CUM.K  <- c(6, 12, 18, 30, 36)
  SHELL.CUTOFF <- c(1.5, 1.8, 2.2, 2.8, 3.2)
  if (identical(max.nn.dist, "auto")) {
    s <- match(k, SHELL.CUM.K)
    if (is.na(s))
      stop("k = ", k, " does not complete a coordination shell of the Visium ",
           "lattice. Use k in {", paste(SHELL.CUM.K, collapse = ", "),
           "} or set max.nn.dist to a number.")
    max.nn.dist <- SHELL.CUTOFF[s]
  }
  stopifnot(is.numeric(max.nn.dist), length(max.nn.dist) == 1)

  knn     <- get.knn(as.matrix(coords), k = k)
  knn.idx <- knn$nn.index
  nn1.med <- median(knn$nn.dist[, 1])                 ## ~= the spot spacing d
  knn.idx[knn$nn.dist > max.nn.dist * nn1.med] <- NA  ## Drop over-distant neighbours.

  ## ---- propagate layer ----------------------------------------------------
  exit.reason <- "max.iter.reached"
  for (i in 1:max.iter) {
    unknown.idx <- which(states == "Unknown")

    if (length(unknown.idx) == 0) {
      message("All cells have been classified.")
      exit.reason <- "all.classified"
      break
    }

    new.states   <- states
    n.mal.added  <- 0
    n.norm.added <- 0

    for (j in unknown.idx) {
      neighbor.idx    <- knn.idx[j, ]
      neighbor.states <- states[neighbor.idx]

      n.mal  <- sum(neighbor.states == "Malignant", na.rm = TRUE)
      n.norm <- sum(neighbor.states == "Normal",    na.rm = TRUE)

      ## voting; vote.margin = 1 is the simple majority used originally
      if ((n.mal - n.norm) >= vote.margin) {
        new.states[j] <- "Malignant"
        layers[j]     <- i
        n.mal.added   <- n.mal.added + 1
      } else if ((n.norm - n.mal) >= vote.margin) {
        new.states[j] <- "Normal"
        layers[j]     <- -i
        n.norm.added  <- n.norm.added + 1
      }
    }

    ## update
    states <- new.states
    if (n.mal.added == 0 && n.norm.added == 0) {
      exit.reason <- "converged"
      break
    }
  }

  object$propagated.labels <- states
  object$propagated.layers <- layers


  object@misc$FindTumorLayers <- list(
    params = list(k = k, n.seed = n.seed, min.seed.size = min.seed.size,
                  exclude.reference = exclude.reference,
                  max.nn.dist = max.nn.dist, vote.margin = vote.margin,
                  max.iter = max.iter),
    seed.clusters      = list(malignant = malignant.clusters, normal = normal.clusters),
    cluster.median.cnv = cnv.score.by.cluster,
    nn1.median         = nn1.med,
    mean.degree        = mean(rowSums(!is.na(knn.idx))),
    n.spots            = length(states),
    n.malignant        = sum(states == "Malignant"),
    n.normal           = sum(states == "Normal"),
    n.unknown          = sum(states == "Unknown"),
    pct.unknown        = 100 * mean(states == "Unknown"),
    exit.reason        = exit.reason,
    converged          = exit.reason != "max.iter.reached",
    last.iter          = suppressWarnings(max(abs(layers), na.rm = TRUE)),
    n.layer.levels     = length(unique(na.omit(layers[states == "Malignant"]))))

  return(object)
}
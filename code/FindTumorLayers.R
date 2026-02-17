library(Seurat)
library(FNN)

FindTumorLayers <- function(object, max.iter = 20) {
  cnv.score.by.cluster <- tapply(object$cnv_score, object$cnv_cluster, median)
  malignant.clusters   <- names(sort(cnv.score.by.cluster, decreasing = TRUE)[1:2])
  normal.clusters      <- names(sort(cnv.score.by.cluster)[1:2])
  
  
  current.labels <- object$cnv_cluster
  states         <- rep("Unknown", length(current.labels))
  names(states)  <- rownames(object@meta.data)
  
  
  states[current.labels %in% malignant.clusters] <- "Malignant"
  states[current.labels %in% normal.clusters]    <- "Normal"
  
  ## Initialize Layer
  layers        <- rep(NA, length(states))
  names(layers) <- names(states)
  layers[states == "Malignant"] <- 0
  layers[states == "Normal"]    <- 0
  
  
  ## position matrix
  images  <- Images(object)
  idx     <- if (length(images) == 1) 1 else {sapply(images, \(x) sum(rownames(object@images[[x]]@coordinates) %in% colnames(object))) |> which.max()}
  coords  <- GetTissueCoordinates(object, image = images[idx])
  if(ncol(coords) >= 2) coords <- coords[, 1:2]
  
  
  if (setequal(rownames(coords), rownames(object@meta.data))) {
    coords <- coords[rownames(object@meta.data), ]
  } else {
    stop("Error: Row names of coordinates do not match Seurat object. Please ensure they contain the exact same set of cells.")
  }
  
  
  knn     <- get.knn(as.matrix(coords), k = 6)
  knn.idx <- knn$nn.index
  
  
  ## propogate layer
  for (i in 1:max.iter) {
    
    unknown.idx <- which(states == "Unknown")
    
    if (length(unknown.idx) == 0) {
      message("All cells have been classified.")
      break
    }
    
    new.states   <- states
    n.mal.added  <- 0
    n.norm.added <- 0
    
    
    for (idx in unknown.idx) {
      neighbor.idx    <- knn.idx[idx, ]
      neighbor.states <- states[neighbor.idx]
      
      n.mal  <- sum(neighbor.states == "Malignant")
      n.norm <- sum(neighbor.states == "Normal")
      
      ## voting
      if ((n.mal + n.norm) >= 1) { 
        if (n.mal > n.norm) {
          new.states[idx] <- "Malignant"
          layers[idx]     <- i
          n.mal.added     <- n.mal.added + 1
        } else if (n.norm > n.mal) {
          new.states[idx] <- "Normal"
          layers[idx]     <- -i
          n.norm.added    <- n.norm.added + 1
        }
      }
    }
    
    
    ## update
    states <- new.states
    if (n.mal.added == 0 && n.norm.added == 0) break
  }
  
  
  object$propagated.labels <- states
  object$propagated.layers <- layers
  return(object)
}
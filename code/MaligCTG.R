library(Seurat)


get.malignant.cells.based.on.posterior.p <- function(malignant_score, normal_score, p_threshold = 0.9) {
  library(mclust)
  
  
  set.seed(1234)
  mix.mal      <- Mclust(malignant_score, G = 2, verbose = FALSE)
  mal.comp.idx <- which.max(mix.mal$parameters$mean)
  mal.p        <- mix.mal$z[, mal.comp.idx]
  cells.A      <- mal.p > p_threshold
  
  set.seed(1234)
  mix.nor      <- Mclust(normal_score, G = 2, verbose = FALSE)
  mal.comp.idx <- which.min(mix.nor$parameters$mean)
  mal.p        <- mix.nor$z[, mal.comp.idx]
  cells.B      <- mal.p > p_threshold
  
  return(cells.A & cells.B)
}


CTGs   <- readRDS("../data/CTGs.rds")
seu  ## a seurat object


## (1) CTG burden --------------------------------------------------------------
genes  <- intersect(rownames(seu), pan.sample.CTGs) 
n.CTGs <- colSums(seu[["RNA"]]@data[genes, , drop = FALSE] > 0)


## (2) Similarity score --------------------------------------------------------
n.CTGs    <- colSums(x[["RNA"]]@data[genes, , drop = FALSE] > 0)
top       <- n.CTGs >= quantile(n.CTGs, 0.95)
tail      <- n.CTGs <= quantile(n.CTGs, 0.05)
top.expr  <- x[["RNA"]]@counts[, top, drop = FALSE] |> rowSums() |> (\(x) log2(1e6 * x / sum(x) + 1))()
tail.expr <- x[["RNA"]]@counts[, tail, drop = FALSE] |> rowSums() |> (\(x) log2(1e6 * x / sum(x) + 1))()
ref       <- cbind(top.expr, tail.expr)
pred      <- SingleR::SingleR(test = x[["RNA"]]@data, ref = ref, labels = c("high", "low"))

df        <- data.frame(barcode = colnames(x), 
                        n_CTGs  = n.CTGs, 
                        SingleR_malignant_score = pred$scores[, "high"], 
                        SingleR_normal_score    = pred$scores[, "low"],
                        row.names = NULL)


## (3) Classification ----------------------------------------------------------
pred      <- get.malignant.cells.based.on.posterior.p(df$SingleR_malignant_score, df$SingleR_normal_score, p_threshold = 0.9)



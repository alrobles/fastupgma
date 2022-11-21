#' fastupgma
#' @description  UPGMA and WPGMA clustering.
#' Just a wrapper function around hclust from fastcluster package
#' instead base hcluse.
#' @param D A distance matrix object to create the phylogenetic tree
#' @param method The agglomeration method to be used. This should be
#'(an unambiguous abbreviation of) one of "ward", "single", "complete",
#' "average", "mcquitty", "median" or "centroid". The default is "average".
#' @param ... Further arguments passed to or from other methods.
#' @importFrom stats as.dist reorder
#' @importFrom fastcluster hclust
#' @importFrom ape as.phylo
#'
#' @return A phylogenetic tree of class phylo.
#' @export
#'
#' @examples
#' fastupgma(CrociduraDistance)
#' fastupgma(rodentiaDistance, method = "mcquitty")
fastupgma <- function (D, method = "average", ...)
{
  DD <- stats::as.dist(D)
  hc <- fastcluster::hclust(DD, method = method, ...)
  result <- ape::as.phylo(hc)
  result <- stats::reorder(result, "postorder")
  result
}

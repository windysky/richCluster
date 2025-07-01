#' richCluster: Fast, Robust Clustering Algorithms for Gene Enrichment Data
#'
#' richCluster is a fast C++ agglomerative hierarchical clustering algorithm packaged 
#' into easily callable R functions, designed to help cluster biological 'terms' based on how similar 
#' of genes are expressed in their activation.
#'
#' Clustering:
#' \itemize{
#'   \item \code{\link{cluster}}: Clusters biological terms by gene content.
#' }
#'
#' @section Visualization:
#' 
#' We support two broad categories of cluster visualization:
#' 
#' Cluster-level: Aggregates pvalue across all terms in the cluster, to compare different clusters to each other
#' Term-level: Displays individual term pvalues within a cluster
#' 
#' Heatmaps:
#' \itemize{
#'   \item \code{\link{cluster_hmap}}: Heatmap of p-values across all clusters and datasets
#'   \item \code{\link{term_hmap}}: Heatmap of terms within the specified cluster(s)
#'   \item \code{\link{cluster_correlation_hmap}}: Correlation heatmap between terms within a cluster
#' }
#'
#' Bar Plots:
#' \itemize{
#'   \item \code{\link{cluster_bar}}: Bar plot of p-values across all clusters and datasets
#'   \item \code{\link{term_bar}}: Bar plot of term-level enrichment for given cluster(s).
#' }
#'
#' Dot Plots:
#' \itemize{
#'   \item \code{\link{cluster_dot}}: Dot plot summarizing cluster size and enrichment across datasets.
#'   \item \code{\link{term_dot}}: Dot plot of terms + gene size within the specific cluster(s).
#' }
#'
#' Network Graphs:
#' \itemize{
#'   \item \code{\link{cluster_network}}: Network graph of all term correlations within a cluster
#'   \item \code{\link{full_network}}: Displays the full network graph of all terms
#' }
#'
#' Export:
#' \itemize{
#'   \item \code{\link{export_df}}: Converts clustering results to a tidy data frame
#' }
#' @docType _PACKAGE
#' @name richCluster
#' @useDynLib richCluster
#' @include load_all.R
NULL

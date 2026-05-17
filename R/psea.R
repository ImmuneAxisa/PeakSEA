utils::globalVariables(c("seqnames", "start", "end", "term", "gene"))

#' Peak-Set Enrichment Analysis (PSEA)
#'
#' Rank-based enrichment analysis for genomic peak sets. Peaks in `gr` are
#' overlapped with annotated peak sets in `gr_set` and GSEA is applied to
#' identify enriched peak sets.
#'
#' @param gr      GRanges object with a peak ranking metric in its metadata
#'   columns.
#' @param rank_col Bare column name in `gr` metadata to use for ranking peaks.
#' @param gr_set  GRanges object annotating peak set membership. Must contain
#'   a metadata column named `"term"` that labels each range with its peak set
#'   identifier.
#' @param ...     Additional arguments passed to [plyranges::join_overlap_left()]
#'   (e.g. `maxgap`, `minoverlap`).
#'
#' @return A `gseaResult` object as returned by [clusterProfiler::GSEA()].
#'
#' @importFrom dplyr select filter arrange desc
#' @importFrom plyranges join_overlap_left
#' @importFrom tidyr unite
#' @importFrom tibble deframe
#' @importFrom clusterProfiler GSEA
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(plyranges)
#'
#' gr <- data.frame(
#'   seqnames = "chr1",
#'   start    = seq(1e4, by = 1e3, length.out = 100),
#'   width    = 500L,
#'   log2FC   = seq(-50, 49)
#' ) |>
#'   plyranges::as_granges()
#'
#' gr_set <- c(
#'   shift_right(gr[1:50],  100L) |> mutate(term = "peakSet1"),
#'   shift_right(gr[51:100], 100L) |> mutate(term = "peakSet2")
#' )
#'
#' psea(gr, log2FC, gr_set)
#' }
psea <- function(gr, rank_col, gr_set, ...) {
  rank_col_name <- deparse(substitute(rank_col))
  validate_psea_inputs(gr, rank_col_name, gr_set)

  peak_sets <- join_overlap_left(gr, gr_set, ...) |>
    as.data.frame() |>
    tidyr::unite("gene", seqnames, start, end) |> # needs "gene" column name for clusterProfiler
    dplyr::filter(!is.na(term)) |>
    dplyr::select(term, gene)

  ranked_vector <- gr |>
    as.data.frame() |>
    tidyr::unite("gene", seqnames, start, end) |>
    dplyr::select(gene, {{ rank_col }}) |>
    dplyr::arrange(dplyr::desc({{ rank_col }})) |>
    tibble::deframe()

  GSEA(ranked_vector, TERM2GENE = peak_sets)
}
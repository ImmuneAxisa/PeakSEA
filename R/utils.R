#' Validate inputs to psea()
#'
#' Stops with an informative message if any argument is malformed.
#'
#' @param gr           Input GRanges object.
#' @param rank_col_name Character; name of the ranking column in `gr`.
#' @param gr_set       GRanges object with a `term` metadata column.
#'
#' @noRd
validate_psea_inputs <- function(gr, rank_col_name, gr_set) {
  if (!inherits(gr, "GRanges")) {
    stop("`gr` must be a GRanges object, got: ", class(gr)[1L], call. = FALSE)
  }
  if (!inherits(gr_set, "GRanges")) {
    stop("`gr_set` must be a GRanges object, got: ", class(gr_set)[1L],
         call. = FALSE)
  }
  if (!"term" %in% names(as.data.frame(gr_set))) {
    stop("`gr_set` must contain a metadata column named 'term'.", call. = FALSE)
  }
  if (!rank_col_name %in% names(as.data.frame(gr))) {
    stop("Column '", rank_col_name, "' not found in `gr` metadata.", call. = FALSE)
  }
}

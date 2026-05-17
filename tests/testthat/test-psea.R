# Happy-path ------------------------------------------------------------------

test_that("psea returns a gseaResult object", {
  result <- suppressMessages(psea(gr_peaks, log2FC, gr_sets))
  expect_s4_class(result, "gseaResult")
})

test_that("psea result contains numeric NES scores", {
  result <- suppressMessages(psea(gr_peaks, log2FC, gr_sets))
  expect_true(is.numeric(result@result$NES))
})

test_that("psea geneList is sorted in descending order", {
  result <- suppressMessages(psea(gr_peaks, log2FC, gr_sets))
  expect_true(all(diff(result@geneList) <= 0))
})

# Input validation ------------------------------------------------------------

test_that("psea errors if gr is not a GRanges object", {
  expect_error(
    psea(data.frame(seqnames = "chr1", start = 1, end = 100, log2FC = 1.0),
         log2FC, gr_sets),
    regexp = "`gr` must be a GRanges object"
  )
})

test_that("psea errors if gr_set is not a GRanges object", {
  expect_error(
    psea(gr_peaks, log2FC,
         data.frame(seqnames = "chr1", start = 1, end = 100, term = "A")),
    regexp = "`gr_set` must be a GRanges object"
  )
})

test_that("psea errors if gr_set has no 'term' column", {
  gr_set_noterm <- plyranges::shift_right(gr_peaks[1:50], 100L)
  expect_error(
    psea(gr_peaks, log2FC, gr_set_noterm),
    regexp = "'term'"
  )
})

test_that("psea errors if rank_col is absent from gr metadata", {
  expect_error(
    psea(gr_peaks, nonexistent_col, gr_sets),
    regexp = "nonexistent_col"
  )
})

# Argument forwarding ---------------------------------------------------------

test_that("psea passes ... to join_overlap_left", {
  # maxgap = 600 extends the matching window beyond the 500 bp peak width,
  # drawing in additional matches compared to the default (maxgap = 0).
  result_gap <- suppressMessages(
    psea(gr_peaks, log2FC, gr_sets, maxgap = 600L)
  )
  expect_s4_class(result_gap, "gseaResult")
})

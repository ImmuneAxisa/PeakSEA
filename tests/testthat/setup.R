library(plyranges)

# Shared test fixtures ---------------------------------------------------------
# 100 peaks on chr1, 500 bp wide, spaced 1 kb apart, with a log2FC metric.
set.seed(123)
gr_peaks <- data.frame(
  seqnames = "chr1",
  start    = seq(1e4, by = 1e3, length.out = 100),
  width    = 500L,
  log2FC   = rnorm(100)
) |>
  plyranges::as_granges()

# Two non-overlapping peak sets: each is a 100 bp-shifted copy of a disjoint
# subset of gr_peaks. The 100 bp shift keeps peaks inside the original 500 bp
# window, guaranteeing genomic overlap with their source peaks.
gr_sets <- c(
  plyranges::shift_right(gr_peaks[1:50],  100L) |> plyranges::mutate(term = "setA"),
  plyranges::shift_right(gr_peaks[51:100], 100L) |> plyranges::mutate(term = "setB")
)

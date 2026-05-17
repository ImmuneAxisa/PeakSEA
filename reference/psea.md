# Peak-Set Enrichment Analysis (PSEA)

Rank-based enrichment analysis for genomic peak sets. Peaks in `gr` are
overlapped with annotated peak sets in `gr_set` and GSEA is applied to
identify enriched peak sets.

## Usage

``` r
psea(gr, rank_col, gr_set, ...)
```

## Arguments

- gr:

  GRanges object with a peak ranking metric in its metadata columns.

- rank_col:

  Bare column name in `gr` metadata to use for ranking peaks.

- gr_set:

  GRanges object annotating peak set membership. Must contain a metadata
  column named `"term"` that labels each range with its peak set
  identifier.

- ...:

  Additional arguments passed to
  [`plyranges::join_overlap_left()`](https://tidyomics.github.io/plyranges/reference/overlap-joins.html)
  (e.g. `maxgap`, `minoverlap`).

## Value

A `gseaResult` object as returned by
[`clusterProfiler::GSEA()`](https://rdrr.io/pkg/clusterProfiler/man/GSEA.html).

## Examples

``` r
# \donttest{
library(plyranges)
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: IRanges
#> Loading required package: S4Vectors
#> Loading required package: stats4
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: GenomicRanges
#> Loading required package: Seqinfo
#> Loading required package: dplyr
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:GenomicRanges’:
#> 
#>     intersect, setdiff, union
#> The following object is masked from ‘package:Seqinfo’:
#> 
#>     intersect
#> The following objects are masked from ‘package:IRanges’:
#> 
#>     collapse, desc, intersect, setdiff, slice, union
#> The following objects are masked from ‘package:S4Vectors’:
#> 
#>     first, intersect, rename, setdiff, setequal, union
#> The following objects are masked from ‘package:BiocGenerics’:
#> 
#>     combine, intersect, setdiff, setequal, union
#> The following object is masked from ‘package:generics’:
#> 
#>     explain
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
#> 
#> Attaching package: ‘plyranges’
#> The following objects are masked from ‘package:dplyr’:
#> 
#>     between, n, n_distinct

gr <- data.frame(
  seqnames = "chr1",
  start    = seq(1e4, by = 1e3, length.out = 100),
  width    = 500L,
  log2FC   = seq(-50, 49)
) |>
  plyranges::as_granges()

gr_set <- c(
  plyranges::shift_right(gr[1:50],  100L) |> plyranges::mutate(term = "peakSet1"),
  plyranges::shift_right(gr[51:100], 100L) |> plyranges::mutate(term = "peakSet2")
)
#> Error: 'mutate' is not an exported object from 'namespace:plyranges'

psea(gr, log2FC, gr_set)
#> Error: object 'gr_set' not found
# }
```

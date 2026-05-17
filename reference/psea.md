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
  shift_right(gr[1:50],  100L) |> mutate(term = "peakSet1"),
  shift_right(gr[51:100], 100L) |> mutate(term = "peakSet2")
)

psea(gr, log2FC, gr_set)
#> Warning: For some pathways, in reality P-values are less than 1e-10. You can set the eps argument to zero for better estimation.
#> Warning: qvalue::qvalue() failed, returning NA for qvalue. Error: missing values and NaN's not allowed if 'na.rm' is FALSE
#> #
#> # Gene Set Enrichment Analysis
#> #
#> #...@organism     unknown 
#> #...@setType      unknown 
#> #...@keytype      unknown 
#> #...@geneList     Named int [1:100] 49 48 47 46 45 44 43 42 41 40 ...
#>  - attr(*, "names")= chr [1:100] "chr1_109000_109499" "chr1_108000_108499" "chr1_107000_107499" "chr1_106000_106499" ...
#> #...nPerm     1000 
#> #...pvalues adjusted by 'BH' with cutoff < 0.05
#> #...2 enriched terms found
#> 'data.frame':    2 obs. of  12 variables:
#>  $ ID             : chr  "peakSet2" "peakSet1"
#>  $ Description    : chr  "peakSet2" "peakSet1"
#>  $ setSize        : int  50 50
#>  $ enrichmentScore: num  1 -1
#>  $ NES            : num  4.04 -4
#>  $ pvalue         : num  1e-10 1e-10
#>  $ p.adjust       : num  1e-10 1e-10
#>  $ qvalue         : num  NA NA
#>  $ rank           : int  49 51
#>  $ leading_edge   : chr  "tags=98%, list=49%, signal=100%" "tags=100%, list=51%, signal=98%"
#>  $ core_enrichment: chr  "chr1_109000_109499/chr1_108000_108499/chr1_107000_107499/chr1_106000_106499/chr1_105000_105499/chr1_104000_1044"| __truncated__ "chr1_59000_59499/chr1_58000_58499/chr1_57000_57499/chr1_56000_56499/chr1_55000_55499/chr1_54000_54499/chr1_5300"| __truncated__
#>  $ log2err        : num  NaN NaN
#> #...Citation
#> S Xu, E Hu, Y Cai, Z Xie, X Luo, L Zhan, W Tang, Q Wang, B Liu, R Wang, W Xie, T Wu, L Xie, G Yu. Using clusterProfiler to characterize multiomics data. Nature Protocols. 2024, 19(11):3292-3320 
#> 
# }
```

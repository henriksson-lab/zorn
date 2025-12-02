# Plot minimum number of dimensions needed to retain distance between samples

https://en.wikipedia.org/wiki/Johnson–Lindenstrauss_lemma

## Usage

``` r
PlotJohnsonLindenstraussMinDim(listEps, minCells = 10, maxCells = 1e+06)
```

## Arguments

- listEps:

  List of eps to plot for

- minCells:

  Plot range min cells to consider

- maxCells:

  Plot range max cells to consider

## Value

A ggplot object

## Details

(1 - eps) \|\|u - v\|\|^2 \< \|\|p(u) - p(v)\|\|^2 \< (1 + eps) \|\|u -
v\|\|^2

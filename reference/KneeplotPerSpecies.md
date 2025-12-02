# Produce a kneeplot

Produce a kneeplot

## Usage

``` r
KneeplotPerSpecies(adata, maxSpecies = NULL)
```

## Arguments

- adata:

  A Seurat object with the DefaultAssay having counts per species (or
  similar)

- maxSpecies:

  Maximum number of species to show. The most abundant species will be
  shown first

## Value

A ggplot object

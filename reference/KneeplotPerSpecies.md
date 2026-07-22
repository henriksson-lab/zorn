# Produce a kneeplot

Produce a kneeplot

## Usage

``` r
KneeplotPerSpecies(adata, maxSpecies = NULL, showTotal = FALSE)
```

## Arguments

- adata:

  A Seurat object with the DefaultAssay having counts per species (or
  similar)

- maxSpecies:

  Maximum number of species to show. The most abundant species will be
  shown first

- showTotal:

  Include a separate kneeplot line for the total count across all
  species

## Value

A ggplot object

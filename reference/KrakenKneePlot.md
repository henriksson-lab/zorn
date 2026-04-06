# Take a KRAKEN2 adata object and generate per-species kneeplots

Take a KRAKEN2 adata object and generate per-species kneeplots

## Usage

``` r
KrakenKneePlot(
  adata,
  groupby = c("phylum", "class", "order", "family", "genus", "species"),
  showNumSpecies = 15,
  sortByName = FALSE,
  includeTotal = TRUE
)
```

## Arguments

- adata:

  Seurat object

- groupby:

  Which taxonomic level to group cells by

- showNumSpecies:

  Max number of species to show

- sortByName:

  Sort by name

- includeTotal:

  Show also a total count kneeplot

## Value

A ggplot object

# Add metadata to a Seurat object, handling cell mismatches

Unlike Seurat's AddMetaData, this function handles the case where the
metadata data.frame has fewer or more rows than cells in the object.
Missing cells get NA values.

## Usage

``` r
BascetAddMetaData(
  adata,
  metadata,
  columns = NULL,
  subsetCommon = FALSE,
  excludeColumns = c("_index", "cell_index", "taxid_index", "cnt"),
  prefix = NULL
)
```

## Arguments

- adata:

  A Seurat object

- metadata:

  A data.frame with rownames matching cell names

- columns:

  Columns to add. Default: all except excludeColumns

- subsetCommon:

  If TRUE, subset the Seurat object to only cells present in both adata
  and metadata

- excludeColumns:

  Columns excluded when columns is NULL

- prefix:

  Optional prefix added to column names before appending

## Value

The Seurat object with added metadata columns

# From a fragments file, get a chromatin assay for Signac.

Note: signac is dirt slow at counting as of writing. It might be
scalable enough for certain inputs and tasks, but we still provide the
option of using it.

## Usage

``` r
FragmentsToSignac(fragpath)
```

## Arguments

- fragpath:

  List to fragment file

## Value

A ChromatinAssay

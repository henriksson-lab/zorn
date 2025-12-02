# Show the FASTQC HTML report for a cell, in the available web browser

Show the FASTQC HTML report for a cell, in the available web browser

## Usage

``` r
ShowFASTQCforCell(
  bascetRoot,
  inputName,
  cellID,
  readnum = c(1, 2),
  useBrowser = FALSE,
  verbose = FALSE
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- inputName:

  Name of input shard

- cellID:

  Name of the cell

- readnum:

  1 or 2, for R1 and R2

- useBrowser:

  Use operating system browser to open file

- verbose:

  Show debug output

## Value

Nothing

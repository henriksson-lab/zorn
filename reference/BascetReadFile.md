# Read one file from a Bascet

TODO This can be made faster by, e.g., once and for all reading the
location of all objects in the file

## Usage

``` r
BascetReadFile(
  bascetFile,
  cellID,
  filename,
  as = c("tempfile", "text"),
  bascetInstance = GetDefaultBascetInstance(),
  verbose = FALSE
)
```

## Arguments

- bascetFile:

  Bascet file instance

- cellID:

  Name of the cell

- filename:

  Name of the file

- as:

  Format requested; affects return value

- bascetInstance:

  A Bascet instance

- verbose:

  Print additional information, primarily to help troubleshooting

## Value

If as="tempfile": Name of a temporary file where the read content is
stored

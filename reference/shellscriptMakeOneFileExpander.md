# Helper function that takes content of a list and generates a BASH script that stores the content in a temporary file during execution

Helper function that takes content of a list and generates a BASH script
that stores the content in a temporary file during execution

## Usage

``` r
shellscriptMakeOneFileExpander(tmpname, list_lines)
```

## Arguments

- tmpname:

  Name of environment variable in which the name of the file will be
  stored

- list_lines:

  Content to write to the file

## Value

BASH script content for the expander

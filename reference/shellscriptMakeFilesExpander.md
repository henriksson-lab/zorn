# Generate a shell script command to produce a file of list of strings

The name of each file will be in VARIABLETASK_ID where TASK_ID starts
from 0

## Usage

``` r
shellscriptMakeFilesExpander(for_variable, list_content)
```

## Arguments

- for_variable:

  Environment variable with array of files having the content

- list_content:

  Content for each file

## Value

Part of a shell script

## Details

Note that there must be one file per task, or each task deletes its own
file. This will result in race conditions

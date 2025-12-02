# Call a MAP function for all cells

Call a MAP function for all cells

## Usage

``` r
BascetMapCell(
  bascetRoot,
  withfunction,
  inputName,
  outputName,
  args = list(),
  overwrite = FALSE,
  runner = GetDefaultBascetRunner(),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- withfunction:

  Function to apply

- inputName:

  Name of input shard

- outputName:

  Name of output shard

- args:

  List of arguments (key,value) to provide to the script

- overwrite:

  Force overwriting of existing files. The default is to do nothing
  files exist

- runner:

  The job manager, specifying how the command will be run (e.g. locally,
  or via SLURM)

- bascetInstance:

  A Bascet instance

## Value

A job to be executed, or being executed, depending on runner settings

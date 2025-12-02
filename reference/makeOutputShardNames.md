# Helper function: Generate suitable output filenames according to shard system i.e. root/name.##.ext

Helper function: Generate suitable output filenames according to shard
system i.e. root/name.##.ext

## Usage

``` r
makeOutputShardNames(bascetRoot, outputName, ext, num_shards)
```

## Arguments

- bascetRoot:

  The root folder where all Bascets are stored

- outputName:

  Name of output shard

- ext:

  File extension, not including the .

- num_shards:

  Number of output shards

## Value

List of files to write shards to

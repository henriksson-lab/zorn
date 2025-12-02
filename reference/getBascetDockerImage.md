# Get and install a Bascet docker image. It will be cached to avoid downloading it each the time the function is called

Get and install a Bascet docker image. It will be cached to avoid
downloading it each the time the function is called

## Usage

``` r
getBascetDockerImage(
  storeAt = getwd(),
  tempdir = NULL,
  forceInstall = FALSE,
  mapDirs = NULL,
  verbose = FALSE
)
```

## Arguments

- storeAt:

  Directory to store the container in. Default is current directory but
  it is likely better to provide a single systems level directory

- tempdir:

  Default is to create a directory for temporary files in the current
  directory. Place it on a fast disk if possible

- forceInstall:

  Set to true to overwrite any existing Docker image

- mapDirs:

  Directories to map through to the container

- verbose:

  Print additional information, primarily to help troubleshooting

## Value

A Bascet instance

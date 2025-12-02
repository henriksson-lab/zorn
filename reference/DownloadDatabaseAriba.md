# Download database for Ariba

Download database for Ariba

## Usage

``` r
DownloadDatabaseAriba(
  dbdir,
  ref = c("ncbi", "argannot", "card", "megares", "plasmidfinder", "resfinder",
    "srst2_argannot", "vfdb_core", "vfdb_full", "virulencefinder"),
  bascetInstance = GetDefaultBascetInstance()
)
```

## Arguments

- dbdir:

  Directory in which to download database

- ref:

  Which reference to download

- bascetInstance:

  A Bascet instance

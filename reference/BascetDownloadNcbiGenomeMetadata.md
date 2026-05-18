# Download and cache an NCBI genome assembly summary

Downloads the NCBI `assembly_summary.txt` metadata file for a
RefSeq/GenBank group. By default this caches RefSeq bacteria metadata
under
`file.path(tools::R_user_dir("Zorn", "cache"), "ncbi", "refseq", "bacteria", "assembly_summary.txt")`.
This is platform-specific and portable across Linux, macOS, and Windows.
If the cached file already exists, it is reused unless
`overwrite = TRUE`.

## Usage

``` r
BascetDownloadNcbiGenomeMetadata(
  db = "refseq",
  group = "bacteria",
  cacheDir = defaultNcbiGenomeCacheDir(),
  overwrite = FALSE,
  dest = NULL
)
```

## Arguments

- db:

  NCBI source database, usually `"refseq"` or `"genbank"`.

- group:

  NCBI genome group, default `"bacteria"`.

- cacheDir:

  Cache root directory. Defaults to the Zorn cache under the user cache
  directory.

- overwrite:

  Re-download when the cached file already exists.

- dest:

  Optional explicit output file path. When provided, this path is used
  instead of `file.path(cacheDir, db, group, "assembly_summary.txt")`.

## Value

Path to the cached `assembly_summary.txt`.

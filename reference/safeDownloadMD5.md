# Download a file, check MD5 to ensure success. This assumes a file.md5 is stored on the server

Download a file, check MD5 to ensure success. This assumes a file.md5 is
stored on the server

## Usage

``` r
safeDownloadMD5(url, file)
```

## Arguments

- url:

  URL to the file to download

- file:

  Name of the file to download content to

## Value

Nothing; panics if the download fails

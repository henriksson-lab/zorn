# Installation

## Requirements

We see two types of users and requirements: - **Users who do heavy
preprocessing** Our tools are designed to maximally benefit from a beefy
computer, ideally an HPC cluster. But a 16GB RAM laptop is technically
enough. - **Users who do postprocessing in R** Any computer should be
enough, but most single-cell users will likely want to have 32 or 64GB
installed to be able to run Seurat on large modern datasets. Bascet
integrates with Seurat.

We strive to support Linux, OSX and Windows.

> ⚠️ Bascet is compiled with the assumption that your CPU supports the
> BMI instruction set. Any computer past 2015 should have it. You will
> need to [compile Bascet
> yourself](https://henriksson-lab.github.io/zorn/articles/for_developers.md)
> if:
>
> - Your CPU does not support the BMI instruction set (uncommon for
>   machines past 2015)

## Step 1: Install Zorn

Zorn is under active development, so just get it straight from GitHub
and install from source:

``` r

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("henriksson-lab/zorn")
```

or newer method (above is becoming deprecated):

``` r

pak::pkg_install("henriksson-lab/zorn")
```

You can then load it:

``` r

library(Zorn)
```

## Step 2: Install Bascet

You also need Bascet, our command-line suite for single-cell analysis.

``` r

bascetInstance.default <- getBascetBinary()
```

or if you want to store the binary in a specific directory (it is
\<500mb):

``` r

bascetInstance.default <- getBascetBinary(storeAt="/somewhere/on/your/disk/")
```

If the file is already downloaded, it will not be downloaded again.

## Step 3: Test the installation

Run the following and see if you get “ok” before proceeding:

``` r

TestBascetInstance(bascetInstance.default)
```

## Using the Bascet instance

You need to pass the bascetInstance variable to each Zorn command that
requires it.

## For developers: Accessing a locally compiled Bascet

See [Bascet for
Developers](https://henriksson-lab.github.io/zorn/articles/for_developers.md)

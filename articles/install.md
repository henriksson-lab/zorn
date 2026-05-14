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
> if your CPU does not support the BMI instruction set (uncommon for
> machines past 2015)

## Step 1: Install upstream package

You will likely need to first install the following packages (might or
might not get pulled in by Zorn automatically):

- [Bioconductor](https://www.bioconductor.org/install/)
- [Seurat](https://satijalab.org/seurat/articles/install.html)
- [Signac](https://stuartlab.org/signac/articles/install)

## Step 2: Install Zorn

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

## Step 3a: Install Bascet: as a precompiled binary

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

## Step 3b: Install Bascet: from source, for developers or if your computer does not support out binaries

See section [for
developers](https://henriksson-lab.github.io/zorn/articles/for_developers.md).

## Step 4: Test the installation

Run the following and see if you get “ok” before proceeding:

``` r

TestBascetInstance(bascetInstance.default)
```

## Summary: Minimal Zorn/Bascet code for local use

This code is all you will need at the end. But if you want to run this
on [SLURM, see separate article
next](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r

library(Zorn)
bascetInstance.default <- getBascetBinary()
bascetRunner.default <- LocalRunner(mem="30g")
```

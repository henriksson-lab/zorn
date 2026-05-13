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

## Step 3a: Install Bascet, precompiled binary

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

## Step 3b: Install Bascet, developer or install from source

If you want to develop Bascet further, or compile from source, you can
follow these instructions to get Bascet, compile it, and use your own
build.

First install Rust from
[rust-lang.org](https://rust-lang.org/tools/install/):

``` bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Then clone the Bascet repository (adjust the directory to wherever you
keep your source checkouts):

``` bash
cd /home/me/git
git clone https://github.com/henriksson-lab/bascet.git
```

You have to make a release build before using Bascet from Zorn. Note
that there are few checks to ensure the code is properly compiled!

``` bash
make
```

You can now point to this build:

``` r

bascetInstance.default <- getBascetDevDir("/home/me/git/bascet/")
```

More [developer information
here](https://henriksson-lab.github.io/zorn/articles/for_developers.md).

## Step 4: Test the installation

Run the following and see if you get “ok” before proceeding:

``` r

TestBascetInstance(bascetInstance.default)
```

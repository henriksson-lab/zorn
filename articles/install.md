# Installation

## Step 1: Install Zorn

Zorn is under active development, so just get it straight from GitHub
and install from source:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("henriksson-lab/zorn")
```

You can then load it:

``` r
library(Zorn)
```

## Step 2: Install Bascet

You also need Bascet, our command-line suite for single-cell analysis.
Most users will find it easiest to run Bascet through Singularity or
Docker; instructions for each are described below. Singularity is for
Linux users, while OSX and Windows users have to use Docker.

### Install Bascet via Singularity (Linux only)

First [install
Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).

Note one possibly annoying default setting! We could not read Bascet
files on some file systems, because they were invisible to Bascet. To
avoid this, edit /etc/singularity/singularity.conf to below:

    mount hostfs = yes

Once installed, Zorn can then pull down the latest Singularity image for
you. We recommend keeping it in a separate directory from your
workspace, as multiple workspaces can share the same image

``` r
bascetInstance.default <- getBascetSingularityImage(storeAt="/somewhere/on/your/disk/")
```

### Install Bascet via Docker (Singularity is recommended for Linux)

First [install Docker](https://docs.docker.com/get-started/).

Once installed, Zorn can then pull down the latest Docker image for you.
We recommend keeping it in a separate directory from your workspace, as
multiple workspaces can share the same image

``` r
bascetInstance.default <- getBascetDockerImage(storeAt="/somewhere/on/your/disk/")
```

## Step 3: Test the installation

Run the following and see if you get “ok” before proceeding:

``` r
TestBascetInstance(bascetInstance.default)
```

## Using the Bascet instance

You need to pass the bascet_inst variable to each Zorn command that
requires it. The singularity image will be cached, so if you run this
again, it will instead use the previous image. You will need to delete
it manually if you wish to replace it.

## For developers: Accessing a locally compiled Bascet

See [Bascet for
Developers](https://henriksson-lab.github.io/zorn/articles/for_developers.md)

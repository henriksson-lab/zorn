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
> - You are running OSX, where our Docker/Podman/Singularity images
>   currently only run via Apple’s x86 emulation layer for Apple
>   Silicon. For native Apple Silicon support you will need to [compile
>   Bascet
>   yourself](https://henriksson-lab.github.io/zorn/articles/for_developers.md).

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
Most users will find it easiest to run Bascet through Singularity or
Docker; instructions for each are described below. Singularity is for
Linux users, while OSX and Windows users have to use Docker.

### Install Bascet via Singularity (Linux only; likely the best option for HPC environments)

First [install
Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).

> ⚠️ Note one possibly annoying default setting! We could not read
> Bascet files on some file systems, because they were invisible to
> Bascet. To avoid this, edit /etc/singularity/singularity.conf to
> below:

    mount hostfs = yes

Once installed, Zorn can then pull down the latest Singularity image for
you. We recommend keeping it in a separate directory from your
workspace, as multiple workspaces can share the same image

``` r
bascetInstance.default <- getBascetSingularityImage(storeAt="/somewhere/on/your/disk/")
```

### Install Bascet via Podman

First [install Podman](https://podman.io/getting-started/installation).

Once installed, Zorn can then pull down the latest Podman image for you.
It can be deleted once it has been loaded into Podman.

``` r
bascetInstance.default <- getBascetPodmanImage(storeAt="/somewhere/on/your/disk/")
```

Note that Podman only exposes certain directories, listed upon start.
You may have to add more directories depending on where your data is
located.

### Install Bascet via Docker (Singularity or Podman are recommended for Linux)

First [install Docker](https://docs.docker.com/get-started/).

Once installed, Zorn can then pull down the latest Docker image for you.
We recommend keeping it in a separate directory from your workspace, as
multiple workspaces can share the same image

``` r
bascetInstance.default <- getBascetDockerImage(storeAt="/somewhere/on/your/disk/")
```

Note that Docker only exposes certain directories, listed upon start.
You may have to add more directories depending on where your data is
located.

> ⚠️ Also note that you need to have Docker desktop running
> (OSX/Windows) whenever you use the container. If you don’t like this,
> please use Podman or Singularity instead.

## Step 3: Test the installation

Run the following and see if you get “ok” before proceeding:

``` r
TestBascetInstance(bascetInstance.default)
```

## Using the Bascet instance

You need to pass the bascetInstance variable to each Zorn command that
requires it. The singularity image will be cached, so if you run this
again, it will instead use the previous image. You will need to delete
it manually if you wish to replace it.

## For developers: Accessing a locally compiled Bascet

See [Bascet for
Developers](https://henriksson-lab.github.io/zorn/articles/for_developers.md)

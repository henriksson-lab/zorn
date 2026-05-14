# De novo assembly

First set up Zorn/Bascet according to the [install
instructions](https://henriksson-lab.github.io/zorn/articles/install.md).
This tutorial assumes that you have [debarcoded the
reads](https://henriksson-lab.github.io/zorn/articles/debarcoding.md).

## When to perform *de novo* assembly

*De novo* assembly is the process of reconstructing a genome (or
transcriptome) from sequencing reads *without using a reference genome
as a guide*. It is used primarily in microbiology, as research on
eukaryotes is commonly done on organisms with known genomes. Possible
exceptions are the assembly of cancer genomes, but assembly of large
genomes is difficult, and Zorn/Bascet is currently optimized for the
handling of smaller microbial genomes.

## De novo assembly using SKESA

Assembly of each cell’s genome can be performed using SKESA:

[(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r

### Assemble all genomes
BascetMapCellSKESA(
  bascetRoot,
  inputName = "filtered",
  outputName = "contigs"
)
```

This produce contigs for each cell, all grouped together in “contigs”.
Note that both Skesa and Spades only assembles cells for which there are
sufficient number of reads.

## Investigating assemblies

The output is zip-files with contigs for each cell, separately. You can
thus easily extract the contigs of interest and analyze with any tool of
choice (not just Zorn/Bascet).

You can obtain the contigs for one cell like this (this is the fastest
way):

``` r

b <- OpenBascet(bascetRoot, "contigs", bascetInstance)
contigs <- BascetReadFile(b, cellID = "<cell_id>", filename = "contigs.fa", as = "text")
CloseBascet(b)
```

We also provide convencience functions for getting a subset. The worked
example below will work for future file formats as well, not just ZIP:

``` r

bascetRoot <- "/path/to/your/bascet"
bascetInstance <- GetDefaultBascetInstance()

# Get cell names
cells <- BascetCellNames(bascetRoot, "contigs", bascetInstance = bascetInstance)

# Pick the cells you want — e.g. all of them, or a subset
listCells <- cells$cell[1:10]

# Make sure the output dir exists
outputDir <- "contigs_out"
dir.create(outputDir, showWarnings = FALSE)

# Dump contigs — one <cellid>.fa file per cell
BascetDumpContigs(
  bascetRoot     = bascetRoot,
  inputName      = "contigs",          # shard holding contigs.fa per cell
  listCells      = listCells,
  outputDir      = outputDir,
  bascetInstance = bascetInstance
)
```

## Further analysis of assembled contigs

To assess assemblies, have a look at our [Map
scripts](https://henriksson-lab.github.io/zorn/articles/map_scripts.md)
vignette, which shows how to run per-cell tools (e.g. QUAST for assembly
QC).

## Coassembly

If you have cells with highly similar genomes, you might be able to
generate “coassemblies” - consensus assemblies using reads from multiple
similar cells.

To do this, simply extract the debarcoded reads from the cells of
interest. There are many ways of picking them, where one method is to
use the clustering function in Seurat. But you can use any method you
wish to come up with a list of cell names. To learn more about
clustering, [see this Seurat tutorial
first](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html). Note
that it is an open research problem how to best pick cells for
coassembly!

After clustering, you can get the names of cells in cluster “0”:

``` r
### Get cell names
listCells <- rownames(adata[adata$cluster_id="0",])
```

Next extract a FASTQ with reads:

``` r

BascetDumpContigs(
  bascetRoot     = bascetRoot,
  inputName      = "contigs",          # shard holding contigs.fa per cell
  listCells      = listCells,
  outputDir      = outputDir,
  bascetInstance = bascetInstance
)
```

And finally, run [SKESA](https://github.com/ncbi/SKESA) or your
favourite software to assemble the contigs: (in BASH)

``` r

#Run SKESA


########################################################### TODO: run skesa
```

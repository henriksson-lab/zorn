# Map scripts

## Bascet has two types of MAP functions

Zorn/Bascet is designed to let you run all kinds of software that
operates on your cells. This can be either the raw reads, the contigs,
or any other data that you produce. Because these operations can be
computationally intense, all of this happens through the MAP framework.

Bascet currently hosts two type of MAP functions:

- MAP functions written in Rust for highest performance. If the original
  software is not Rust, we perform [static analysis-mediated LLM
  translation](https://github.com/henriksson-lab/rustification) to make
  it compatible
- MAP functions calling software through shell scripts. This is not as
  fast but avoids the creation of thousands (or millions) of small files

## Running first-class Rust MAP functions

These wrappers call native Bascet subcommands directly, so they are the
fastest way to process all cells. Zorn currently provides:

- `BascetMapCellSKESA` — *de novo* assembly with SKESA, producing
  contigs from filtered reads
- `BascetMapCellFASTQC` — read-level quality control with FastQC
- `BascetMapCellGECCO` — biosynthetic gene cluster prediction on
  assembled contigs

### GECCO — biosynthetic gene clusters (for assembled genomes)

[GECCO](https://gecco.embl.de/) scans assembled contigs for biosynthetic
gene clusters (BGCs).

``` r

BascetMapCellGECCO(
  bascetRoot,
  inputName  = "contigs",
  outputName = "gecco"
)
```

Aggregate the per-cell cluster tables into a single list of data.frames:

``` r

gecco_aggr <- BascetAggregateGECCO(
  bascetRoot,
  inputName = "gecco"
)
```

### SKESA — assembly (runs on reads as input)

``` r

BascetMapCellSKESA(
  bascetRoot,
  inputName  = "filtered",
  outputName = "contigs",
  kmer       = 21,
  minContig  = 50
)
```

### FastQC — QC of reads (runs on reads as input)

``` r

BascetMapCellFASTQC(
  bascetRoot,
  inputName  = "filtered",
  outputName = "fastqc"
)

fastqc_aggr <- BascetAggregateFASTQC(
  bascetRoot,
  inputName = "fastqc"
)
```

## Running a shell-script based MAP function

Shell-script wrappers are the easiest way to run external bioinformatics
tools without writing your own MAP function. Zorn includes wrappers for
the most common single-cell microbial workflows:

| Wrapper | Tool | Typical input | Typical output |
|----|----|----|----|
| `BascetMapCellQUAST` | QUAST assembly QC | `contigs` | `quast` |
| `BascetMapCellAbricate` | Abricate AMR / virulence | `contigs` | `abricate` |
| `BascetMapCellBakta` | Bakta genome annotation | `contigs` | `bakta` |
| `BascetMapCellAriba` | Ariba AMR from reads | `filtered` | `ariba` |
| `BascetMapCellAMRfinder` | NCBI AMRfinder | `contigs` | `AMRfinder` |

Each of these is a thin wrapper around `BascetMapCell` that knows the
right script name and arguments.

### QUAST — assembly quality

[(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r

BascetMapCellQUAST(
  bascetRoot,
  inputName  = "contigs",
  outputName = "quast"
)
```

If you prefer the raw form, this is equivalent to:

``` r

BascetMapCell(
  bascetRoot,
  withfunction = "_quast",
  inputName    = "contigs",
  outputName   = "quast"
)
```

### Abricate — AMR / virulence screening

``` r

BascetMapCellAbricate(
  bascetRoot,
  inputName  = "contigs",
  outputName = "abricate",
  db         = "ncbi"
)

abricate_mat <- BascetAggregateAbricate(
  bascetRoot,
  inputName = "abricate"
)
```

### Bakta — genome annotation

``` r

DownloadDatabaseBakta("/path/to/bakta_db", dbtype = "light")

BascetMapCellBakta(
  bascetRoot,
  inputName  = "contigs",
  outputName = "bakta",
  db         = "/path/to/bakta_db"
)
```

### Ariba — AMR identification from reads

``` r

BascetMapCellAriba(
  bascetRoot,
  inputName  = "filtered",
  outputName = "ariba",
  db         = "/path/to/ariba_db/out.prepareref"
)

ariba_mat <- BascetAggregateAriba(
  bascetRoot,
  inputName = "ariba"
)
```

### AMRfinder

``` r

DownloadDatabaseAMRfinder("/path/to/amrfinder_db")

BascetMapCellAMRfinder(
  bascetRoot,
  inputName  = "contigs",
  outputName = "AMRfinder",
  db         = "/path/to/amrfinder_db"
)

amr_df <- BascetAggregateAMRfinder(
  bascetRoot,
  inputName = "AMRfinder"
)
```

## Aggregating MAP results

Once you have run your map function, you most likely want to load the
results into R. We call this procedure “aggregate”. In case of QUAST,
this procedure loads all quality metrics into an R data.frame object:

``` r

quast_aggr <- MapListAsDataFrame(BascetAggregateMap(
  bascetRoot,
  inputName="quast",
  aggr.quast
))
```

## Arguments to MAP functions

Some scripts require additional arguments to be sent (such as a link to
a database file). This is done by setting the args argument. Below will
set two environment variables such that the contents can be picked up
the script:

``` r
BascetMapCell(
  ...
  args(DB="some/path",OTHERARG="hi")
  ...
)
```

## Custom MAP functions - introduction

It is easy to add new functions! Easiest way is to simply copy and
modify the code for an existing script. You can start from either \*
[QUAST](https://github.com/henriksson-lab/bascet/blob/main/src/mapcell_scripts/quast.sh),
which takes contigs as input \*
[SKESA](https://github.com/henriksson-lab/bascet/blob/main/src/mapcell_scripts/skesa.sh),
which takes FASTQ as input

Once you have written your script, you invoke it with a direct path:

``` r

BascetMapCell(
  bascetRoot,
  withfunction = "/path/to/your/script.sh",
  inputName = "...",
  outputName = "..."
)
```

In most cases you want to write your own aggregate function. This
function will take the output from your tool, parse it, and put in a
sensible R object. Have a look at [example and existing aggregate
functions](https://github.com/henriksson-lab/zorn/blob/main/R/aggr_functions.R)
for inspiration.

There is also a catch-all aggregate function that requires a bit of a
special way of calling. The example below takes “out.txt”, generated by
each tool, and stores the raw file content in a list. This is not pretty
but it may help you in debugging and development:

``` r

quast_aggr <- MapListAsDataFrame(BascetAggregateMap(
  bascetRoot,
  inputName="..",
  aggr.raw("out.txt")
))
```

## Custom MAP functions - details

If you look at any [of our MAP
functions](https://github.com/henriksson-lab/bascet/tree/main/crates/bascet-cli/src/mapcell_scripts),
you will find that it is a BASH script that conforms to a certain
pattern. It actually is just a script (in any language) that takes
certain command line arguments.

``` console
--bascet-api
```

The script returns the API version, also validating that it is a valid
script for MAP calls

``` console
--expect-files
```

The script returns a list of what files to extract from the Bascet, for
each cell. Here, “\*” means to get everything. Asking for less means
higher performance

``` console
--missing-file-mode
```

The Bascet what to do if the files are missing. “skip” means to just
proceed with the next cell

``` console
--compression-mode
```

How to compress the output files. “default” means to compress. However,
if your tool generates compressed files already, it is just a waste of
time trying to do it again, in which case the script can return
“uncompressed”.

``` console
--input-dir XXX
```

This is the directory where input files are located

``` console
--output-dir YYY
```

Where to store output to. This directory is already created

``` console
--num-threads ZZZ
```

How many threads to use for this particular process. Note that Bascet is
already calling multiple MAP scripts in parallel and there is thus
typically little benefit in making individual process multithreaded

``` console
--recommend-threads
```

Return how many threads (at least 1) the job should get. This is used if
the user runs mapcell but only specifies the total number of threads.
Bascet will then try to allocate workers accordingly. Return 1 if your
mapcell script does not support multithreading

``` console
--preflight-check
```

This is called once only, to check that the script has the needed
software dependencies. In such case, it returns “MAPCELL-CHECK”

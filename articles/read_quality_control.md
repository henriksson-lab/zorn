# Read-based quality control

First set up Zorn/Bascet according to the [install
instructions](https://henriksson-lab.github.io/zorn/articles/install.md).
This tutorial assumes that you have [debarcoded the
reads](https://henriksson-lab.github.io/zorn/articles/debarcoding.md).

## When to perform read-based quality control

We find that read-based quality control is not a very useful metric for
single-cell analysis. Adapters are already handled in the pipeline, and
the insert size distribution (distance from R1 to R2) is already checked
prior to sequencing. Instead we recommend to just generate a UMAP using,
e.g., the [KRAKEN2
workflow](https://henriksson-lab.github.io/zorn/articles/kraken.md),
will better tell you about common single-cell problems: - If you have
many doublets (mixing of cells) - If you have background RNA/DNA mixing
into cells - If you are able to identity cells

That said, we wrap FASTQC for read-based quality control, which might
help pinpoint certain problems that can arise in the library
preparation.

## FASTQC on library as whole

A common tool to investigate raw reads is
[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
You can always download the original software and run it on the FASTQ
files. But you also have a Rust-translation of FASTQC built into Bascet,
and you can invoke it under exttool (… means FASTQC CLI parameters; see
original tool documentation).

    ./bascet exttool fastqc ...

You will likely find a ton of repeated sequences for single-cell data as
adapters are kept until after debarcoding. You can transform TIRP files
to FASTQ to have adapters trimmed before FASTQC analysis, if relevant
(e.g., to check quality of adapter trimming, etc):

[(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r

### Convert debarcoded reads to FASTQ format
BascetMapTransform(
  bascetRoot,
  inputName="filtered",
  outputName="asfq",
  out_format="R1.fq.gz"
)
```

## FASTQC per cell

You can also run FASTQC on each individual cell (TIRP as input). This
takes a fair bit of time, but can help tell if, e.g., a cluster of cells
is caused by technical issues such as adapter content. You first run
FASTQC on each cell like this:

[(SLURM-compatible
step)](https://henriksson-lab.github.io/zorn/articles/slurm.md)

``` r

BascetMapCellFASTQC(
  bascetRoot,
  inputName = "filtered"  #or other source of reads
)
```

If you have an outlier cell in your dataset, you can investigate its
FASTQ HTML report in the follow manner (opening in the RShiny plot pane,
or separate browser):

``` r

ShowFASTQCforCell(
    bascetFile,
    cellID="xyz", #name of your cell 
    readnum="1", #for R1
)
```

You can also compare cells by aggregating the data. Note that FASTQC
creates rather complex statistics that need further extraction for
simple plotting

``` r

aggr_fastqc <- BascetAggregateFASTQC(
  bascetRoot
)
```

One relevant statistic is the adapter content across the read:

``` r
PlotFASTQCadapterContent <- function(
    aggr_fastqc,
    readnum="1" #for R1
)
```

You can also retrieve a table of pass/fail statistics:

``` r

fastqc_passfail <- GetFASTQCpassfailStats(
    aggr_fastqc,
    readnum="1" #for R1
)
```

Because there are so many things you can do with this statistics, we
provide a general interface to each table that FASTQC generates:

``` r

mystats <- GetFASTQCassembledDF(
    aggr_fastqc, 
    section="see below", 
    readnum="1"
)
```

Possible values of section are:

- “Basic Statistics”
- “Per base sequence quality” d
- “Per sequence quality scores”
- “Per base sequence content”
- “Per sequence GC content”
- “Per base N content”
- “Sequence Length Distribution”
- “Sequence Duplication Levels”
- “Overrepresented sequences”
- “Adapter Content”

Other statistics can also be extracted. Refer to the full R command
reference manual for a list

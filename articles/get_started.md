# Getting Started with Zorn

Zorn and Bascet can do many things:

- Preprocess reads for single-cell metagenomics, WGS and RNA-seq
- Perform de novo assembly and other tasks relevant for whole-genome
  sequencing
- Perform custom per-cell detailed analysis operations

The interface is designed to interoperate with Signac and Seurat for
high-level analysis of the data, after the preprocessing.

Depending on your use case, you will need different tutorials. However,
these are in common for all workflows:

- [Installing Zorn and
  Bascet](https://henriksson-lab.github.io/zorn/articles/install.md)
- [Using
  SLURM](https://henriksson-lab.github.io/zorn/articles/slurm.md) - If
  you want to run this in a high-performance environment
- [Debarcoding
  tutorial](https://henriksson-lab.github.io/zorn/articles/debarcoding.md) -
  Ingesting raw FASTQ, trimming, and preparing for later steps
- [Read-based quality
  control](https://henriksson-lab.github.io/zorn/articles/read_quality_control.md) -
  Quality control of reads

Available for early testing is also:

- [Isolate/SRA-workflow](https://henriksson-lab.github.io/zorn/articles/sra_fetching.md) -
  Load raw read isolate data for processing
- [Isolate/NCBI-workflow](https://henriksson-lab.github.io/zorn/articles/ncbi_genomes.md) -
  Load isolate assembled genome data for processing

The next step to get an overview is likely:

- [KRAKEN-based
  workflow](https://henriksson-lab.github.io/zorn/articles/kraken.md) -
  if you don’t know what your sample contains, and want to compare to
  known species
- [Alignment-based
  workflow](https://henriksson-lab.github.io/zorn/articles/alignment.md) -
  if you have some clue what your sample contains
- [Informative KMER-based
  workflow](https://henriksson-lab.github.io/zorn/articles/kmer.md) - if
  you don’t know what your sample contains, and want to find novel
  genome compositions - Option \#1
- [Count sketch KMER-based
  workflow](https://henriksson-lab.github.io/zorn/articles/countsketch.md) -
  if you don’t know what your sample contains, and want to find novel
  genome compositions - Option \#2

Then there are several options on how to continue with deeper analysis:

- [Doublet
  removal](https://henriksson-lab.github.io/zorn/articles/doublet_removal.md) -
  to flag and remove cells that are mixtures of two cells
- [*De novo*
  assembly](https://henriksson-lab.github.io/zorn/articles/assembly.md)
- [Genome
  annotation](https://henriksson-lab.github.io/zorn/articles/genome_annotation.md) -
  annotate your assemblies
- [SNP
  analysis](https://henriksson-lab.github.io/zorn/articles/snp_analysis.md) -
  to check genetic variation
- [Running software for each cell using
  MAP](https://henriksson-lab.github.io/zorn/articles/map_scripts.md) -
  the general method for performing calculations on each cell

If you are in a hurry, you can check our examples:

- [Atrandi short-read WGS with host
  removal](https://henriksson-lab.github.io/zorn/articles/ex_atrandi_host.md)

These may help you understand the design of Zorn and Bascet:

- [Overview of file
  formats](https://henriksson-lab.github.io/zorn/articles/fileformats.md)
- [The types of Bascet-ZIP
  files](https://henriksson-lab.github.io/zorn/articles/bascet_zip.md)
- [Working with
  Bascet-FASTQ](https://henriksson-lab.github.io/zorn/articles/bascet_fastq.md)
- [For
  developers](https://henriksson-lab.github.io/zorn/articles/for_developers.md) -
  if you wish to modify Bascet
- [Example
  data](https://henriksson-lab.github.io/zorn/articles/example_data.md) -
  you can run the workflow with some small datasets before you scale up

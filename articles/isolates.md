# Isolate analysis

If you want to analyze bulk samples (e.g. isolates) with Zorn/Bascet,
this is absolutely possible, and you then get the benefits of 1. a
workflow manager 2. a suite of analysis tools 3. fewer files to manage

Depending on your use case, there are two ways to do this:

## If you have a few input FASTQ or FASTA files

There are two special file formats, “list-fastq” and “list-fasta”, which
is a TSV that simply consists of \* A name for each sample \* Links to
raw FASTQ-files or FASTA-files

This file can be used in place of a Bascet TIRP or ZIP-file, e.g., as
the source of filtered reads. In effect, you skip debarcoding and
sharding.

TODO file format example

## If you have really many FASTA/FASTQ files

In this case we recommend you to produce TIRP files yourself, in a way
that makes sense based on your infrastructure. See the section on file
formats to learn how to produce it. You can easily append one input file
after another, and thus avoid the need to keep all files on disk at the
same time.

## Mixing single-cell and bulk samples

Bascet currently does not support multiple input file formats for a
given stage; i.e. if you have filtered reads, they must all be in TIRP
format, or another format.

The easiest way to mix samples is to 1. Perform debarcoding and sharding
of your single-cell data into TIRP 2. In a separate working directory
for bulk data, use Bascet transform to convert a list-FASTQ/FASTA to
TIRP 3. Copy this bulk TIRP as filtered.xxx.tirp.gz in your single-cell
directory, where xxx is the next number after your current single-cell
shards

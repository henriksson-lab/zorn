# Package index

## All functions

- [`AggregateMinhashes()`](https://henriksson-lab.github.io/zorn/reference/AggregateMinhashes.md)
  : Aggregate frequency of minhashes across cells
- [`AtrandiBarcodeStats()`](https://henriksson-lab.github.io/zorn/reference/AtrandiBarcodeStats.md)
  : Given a Bascet, produce a matrix showing for each combinatorial
  barcode, how many times it occurs across the cells. Presented as a
  96-well plate matrix
- [`BarnyardPlotMatrix()`](https://henriksson-lab.github.io/zorn/reference/BarnyardPlotMatrix.md)
  : Produce a matrix of Barnyard plots, i.e., counts for one species vs
  another, for all combinations of species.
- [`Bascet-class`](https://henriksson-lab.github.io/zorn/reference/Bascet-class.md)
  : A bascet, along with all the shards
- [`BascetAddMetaData()`](https://henriksson-lab.github.io/zorn/reference/BascetAddMetaData.md)
  : Add metadata to a Seurat object, handling cell mismatches
- [`BascetAggregateAMRfinder()`](https://henriksson-lab.github.io/zorn/reference/BascetAggregateAMRfinder.md)
  : Aggregate data from AMRfinder This is a thin wrapper around
  BascetAggregateMap
- [`BascetAggregateAbricate()`](https://henriksson-lab.github.io/zorn/reference/BascetAggregateAbricate.md)
  : Aggregate data from Abricate This is a thin wrapper around
  BascetAggregateMap
- [`BascetAggregateAriba()`](https://henriksson-lab.github.io/zorn/reference/BascetAggregateAriba.md)
  : Aggregate data from Ariba This is a thin wrapper around
  BascetAggregateMap
- [`BascetAggregateFASTQC()`](https://henriksson-lab.github.io/zorn/reference/BascetAggregateFASTQC.md)
  : Aggregate data from FASTQC This is a thin wrapper around
  BascetAggregateMap
- [`BascetAggregateGECCO()`](https://henriksson-lab.github.io/zorn/reference/BascetAggregateGECCO.md)
  : Aggregate data from GECCO This is a thin wrapper around
  BascetAggregateMap
- [`BascetAggregateMap()`](https://henriksson-lab.github.io/zorn/reference/BascetAggregateMap.md)
  : Aggregate data from previous Map call
- [`BascetAggregateQUAST()`](https://henriksson-lab.github.io/zorn/reference/BascetAggregateQUAST.md)
  : Aggregate data from QUAST This is a thin wrapper around
  BascetAggregateMap
- [`BascetAlignToReference()`](https://henriksson-lab.github.io/zorn/reference/BascetAlignToReference.md)
  : Align from FASTQ, generate sorted and indexed BAM file
- [`BascetAlignmentToBigwig()`](https://henriksson-lab.github.io/zorn/reference/BascetAlignmentToBigwig.md)
  : Generate a bigwig out of all reads in a sorted BAM. Note that the
  caller is responsible for sorting the BAM first
- [`BascetBam2Fragments()`](https://henriksson-lab.github.io/zorn/reference/BascetBam2Fragments.md)
  : Take aligned BAM file and produce Fragments.tsv.gz, compatible with
  Signac ATAC-seq style analysis
- [`BascetCacheComputation()`](https://henriksson-lab.github.io/zorn/reference/BascetCacheComputation.md)
  : A wrapper to cache a computation. Put your function in as an
  argument, as R will only compute its value if needed. If the cache
  file exist, it will not be run again
- [`BascetCellNames()`](https://henriksson-lab.github.io/zorn/reference/BascetCellNames.md)
  : Get list of cells in a Bascet
- [`BascetComputeMinhash()`](https://henriksson-lab.github.io/zorn/reference/BascetComputeMinhash.md)
  : Compute minhashes for each cell. This is a thin wrapper around
  BascetMapCell
- [`BascetCountChrom()`](https://henriksson-lab.github.io/zorn/reference/BascetCountChrom.md)
  : From aligned BAM file, compute counts per chromosome
- [`BascetCountFeature()`](https://henriksson-lab.github.io/zorn/reference/BascetCountFeature.md)
  : From aligned BAM file, compute counts per feature
- [`BascetCountMatrixToAssay()`](https://henriksson-lab.github.io/zorn/reference/BascetCountMatrixToAssay.md)
  : Convert a BascetCountMatrix to a Seurat Assay
- [`BascetDumpContigs()`](https://henriksson-lab.github.io/zorn/reference/BascetDumpContigs.md)
  : Store all contigs in an output directory, as cell_id.fa
- [`BascetFilterAlignment()`](https://henriksson-lab.github.io/zorn/reference/BascetFilterAlignment.md)
  : Filter an alignment (BAM-file).
- [`BascetGatherCountSketch()`](https://henriksson-lab.github.io/zorn/reference/BascetGatherCountSketch.md)
  : Gather all count sketches into a single count sketch matrix
- [`BascetGetRaw()`](https://henriksson-lab.github.io/zorn/reference/BascetGetRaw.md)
  : Extract barcodes and trim input raw FASTQ
- [`BascetIndexGenomeBWAMEM2()`](https://henriksson-lab.github.io/zorn/reference/BascetIndexGenomeBWAMEM2.md)
  : Index a genome using BWA-MEM2 such that it can be used for alignment
- [`BascetIndexGenomeBowtie2()`](https://henriksson-lab.github.io/zorn/reference/BascetIndexGenomeBowtie2.md)
  : Index a genome using Bowtie2 such that it can be used for alignment
- [`BascetIndexGenomeMinimap2()`](https://henriksson-lab.github.io/zorn/reference/BascetIndexGenomeMinimap2.md)
  : Index a genome using minimap2 such that it can be used for alignment
- [`BascetIndexGenomeSTAR()`](https://henriksson-lab.github.io/zorn/reference/BascetIndexGenomeSTAR.md)
  : Index a genome using STAR such that it can be used for alignment
- [`BascetInstance()`](https://henriksson-lab.github.io/zorn/reference/BascetInstance.md)
  : Create a new bascet instance. For advanced users only
- [`BascetListFilesForCell()`](https://henriksson-lab.github.io/zorn/reference/BascetListFilesForCell.md)
  : List files for a cell in a Bascet
- [`BascetLoadCountSketchMatrix()`](https://henriksson-lab.github.io/zorn/reference/BascetLoadCountSketchMatrix.md)
  : Load count sketch matrix as Seurat object
- [`BascetMakeMinhashHistogram()`](https://henriksson-lab.github.io/zorn/reference/BascetMakeMinhashHistogram.md)
  : Gather all minhashes into a single histogram file
- [`BascetMapCell()`](https://henriksson-lab.github.io/zorn/reference/BascetMapCell.md)
  : Call a MAP function for all cells
- [`BascetMapCellAMRfinder()`](https://henriksson-lab.github.io/zorn/reference/BascetMapCellAMRfinder.md)
  : Run AMRfinder on contigs of all cells. This is a thin wrapper around
  BascetMapCell
- [`BascetMapCellAbricate()`](https://henriksson-lab.github.io/zorn/reference/BascetMapCellAbricate.md)
  : Run Abricate on contigs of all cells. This is a thin wrapper around
  BascetMapCell
- [`BascetMapCellAriba()`](https://henriksson-lab.github.io/zorn/reference/BascetMapCellAriba.md)
  : Run Ariba on reads of all cells. This is a thin wrapper around
  BascetMapCell
- [`BascetMapCellBakta()`](https://henriksson-lab.github.io/zorn/reference/BascetMapCellBakta.md)
  : Run Bakta on contigs of all cells. This is a thin wrapper around
  BascetMapCell
- [`BascetMapCellFASTQC()`](https://henriksson-lab.github.io/zorn/reference/BascetMapCellFASTQC.md)
  : Run FASTQC on reads of all cells.
- [`BascetMapCellGECCO()`](https://henriksson-lab.github.io/zorn/reference/BascetMapCellGECCO.md)
  : Run GECCO on contigs of all cells.
- [`BascetMapCellQUAST()`](https://henriksson-lab.github.io/zorn/reference/BascetMapCellQUAST.md)
  : Run QUAST on reads of all cells. This is a thin wrapper around
  BascetMapCell
- [`BascetMapCellSKESA()`](https://henriksson-lab.github.io/zorn/reference/BascetMapCellSKESA.md)
  : Run SKESA on reads of all cells. This is a thin wrapper around
  BascetMapCell
- [`BascetMapCellSKESAintegrated()`](https://henriksson-lab.github.io/zorn/reference/BascetMapCellSKESAintegrated.md)
  : Run integrated SKESA on reads of all cells.
- [`BascetMapTransform()`](https://henriksson-lab.github.io/zorn/reference/BascetMapTransform.md)
  : Transform data
- [`BascetQueryFq()`](https://henriksson-lab.github.io/zorn/reference/BascetQueryFq.md)
  : Build count table from FASTQ reads and a list of selected kmers
- [`BascetReadFile()`](https://henriksson-lab.github.io/zorn/reference/BascetReadFile.md)
  : Read one file from a Bascet
- [`BascetReadMinhashHistogram()`](https://henriksson-lab.github.io/zorn/reference/BascetReadMinhashHistogram.md)
  : Read histogram of KMERs, the output of BascetMakeMinhashHistogram
- [`BascetRunCellSNP()`](https://henriksson-lab.github.io/zorn/reference/BascetRunCellSNP.md)
  : Align from FASTQ, generate sorted and indexed BAM file
- [`BascetRunFASTP()`](https://henriksson-lab.github.io/zorn/reference/BascetRunFASTP.md)
  : Run FASTP for each cell. Input must be in FASTQ file format
- [`BascetRunKraken()`](https://henriksson-lab.github.io/zorn/reference/BascetRunKraken.md)
  : Run KRAKEN2 for each cell. Then produce a count matrix of taxonomy
  IDs from the output
- [`BascetShardify()`](https://henriksson-lab.github.io/zorn/reference/BascetShardify.md)
  : Take debarcoded reads, merge them, and split them into suitable
  numbers of shards.
- [`BascetToFastq()`](https://henriksson-lab.github.io/zorn/reference/BascetToFastq.md)
  : Convert data to Bascet-FASTQ
- [`ChooseInformativeKMERs()`](https://henriksson-lab.github.io/zorn/reference/ChooseInformativeKMERs.md)
  : Pick random KMERs from KMC3 database. The choice is among KMERs
  within a frequency range
- [`ChromToSpeciesCount()`](https://henriksson-lab.github.io/zorn/reference/ChromToSpeciesCount.md)
  : Produce a count matrix on strain level
- [`CloseBascet()`](https://henriksson-lab.github.io/zorn/reference/CloseBascet.md)
  : Close a Bascet file.
- [`CountDataFrameToSparseMatrix()`](https://henriksson-lab.github.io/zorn/reference/CountDataFrameToSparseMatrix.md)
  : Count entries in long format data frame and return as a sparse
  matrix
- [`CountGrangeFeatures()`](https://henriksson-lab.github.io/zorn/reference/CountGrangeFeatures.md)
  : Obtain a feature matrix (as seurat object) given an seurat object
  having Fragments associated
- [`CountSketchUMAP()`](https://henriksson-lab.github.io/zorn/reference/CountSketchUMAP.md)
  : Run UMAP on a count sketch reduction
- [`CreateSeuratObject.BascetCountMatrix()`](https://henriksson-lab.github.io/zorn/reference/CreateSeuratObject.BascetCountMatrix.md)
  : Create a Seurat object from a BascetCountMatrix
- [`CreateSeuratObjectWithReduction()`](https://henriksson-lab.github.io/zorn/reference/CreateSeuratObjectWithReduction.md)
  : Create a seurat object from e.g. count sketch reduction
- [`DebarcodedKneePlot()`](https://henriksson-lab.github.io/zorn/reference/DebarcodedKneePlot.md)
  : Produce summary kneeplot given debarcoded statistics
- [`DetectRawFileMeta()`](https://henriksson-lab.github.io/zorn/reference/DetectRawFileMeta.md)
  : Detect metadata for raw input FASTQ files
- [`DownloadDatabaseAMRfinder()`](https://henriksson-lab.github.io/zorn/reference/DownloadDatabaseAMRfinder.md)
  : Download a database for AMRfinder
- [`DownloadDatabaseAriba()`](https://henriksson-lab.github.io/zorn/reference/DownloadDatabaseAriba.md)
  : Download database for Ariba
- [`DownloadDatabaseBakta()`](https://henriksson-lab.github.io/zorn/reference/DownloadDatabaseBakta.md)
  : Download a database for Bakta
- [`FragmentCountsPerChrom()`](https://henriksson-lab.github.io/zorn/reference/FragmentCountsPerChrom.md)
  : From a Signac chromatin assay with fragments, for each cell, count
  how many reads per chromosome
- [`FragmentCountsPerChromAssay()`](https://henriksson-lab.github.io/zorn/reference/FragmentCountsPerChromAssay.md)
  : From a Signac chromatin assay with fragments, for each cell, count
  how many reads per chromosome. This function directly returns an assay
  that can be added to a Seurat multimodal object
- [`FragmentsToSignac()`](https://henriksson-lab.github.io/zorn/reference/FragmentsToSignac.md)
  : From a fragments file, get a chromatin assay for Signac.
- [`GetBascetTempDir()`](https://henriksson-lab.github.io/zorn/reference/GetBascetTempDir.md)
  : Get a temp directory to use; need to be created
- [`GetDefaultBascetInstance()`](https://henriksson-lab.github.io/zorn/reference/GetDefaultBascetInstance.md)
  : Get default Bascet instance from global variable
  (bascetInstance.default)
- [`GetDefaultBascetRunner()`](https://henriksson-lab.github.io/zorn/reference/GetDefaultBascetRunner.md)
  : Get the current default runner
- [`GetFASTQCassembledDF()`](https://henriksson-lab.github.io/zorn/reference/GetFASTQCassembledDF.md)
  : Get a data frame for one type of FASTQ statistics across across all
  cells
- [`GetFASTQCbasicStats()`](https://henriksson-lab.github.io/zorn/reference/GetFASTQCbasicStats.md)
  : From aggregated FASTQC data, get basic statistics for overlay on
  UMAP etc
- [`GetFASTQCpassfailStats()`](https://henriksson-lab.github.io/zorn/reference/GetFASTQCpassfailStats.md)
  : From aggregated FASTQC data, get overall pass-fail statistics for
  overlay on UMAP etc
- [`KneeplotPerSpecies()`](https://henriksson-lab.github.io/zorn/reference/KneeplotPerSpecies.md)
  : Produce a kneeplot
- [`KrakenFindConsensusTaxonomy()`](https://henriksson-lab.github.io/zorn/reference/KrakenFindConsensusTaxonomy.md)
  : For a KRAKEN2 count matrix, return consensus taxID for each cell as
  metadata
- [`KrakenKneePlot()`](https://henriksson-lab.github.io/zorn/reference/KrakenKneePlot.md)
  : Take a KRAKEN2 adata object and generate per-species kneeplots
- [`KrakenSpeciesDistribution()`](https://henriksson-lab.github.io/zorn/reference/KrakenSpeciesDistribution.md)
  : Using a KRAKEN2 count matrix, produce a "kneeplot" of species
- [`ListDatabaseAbricate()`](https://henriksson-lab.github.io/zorn/reference/ListDatabaseAbricate.md)
  : List installed databases available for Abricate
- [`LocalRunner()`](https://henriksson-lab.github.io/zorn/reference/LocalRunner.md)
  : Create new local runner instance
- [`MapCellMultiListAsDataFrame()`](https://henriksson-lab.github.io/zorn/reference/MapCellMultiListAsDataFrame.md)
  : Convenience function; alternative is to somehow implement
  as.data.frame
- [`MapListAsDataFrame()`](https://henriksson-lab.github.io/zorn/reference/MapListAsDataFrame.md)
  : Convenience function; alternative is to somehow implement
  as.data.frame.
- [`MergeBascetCountMatrix()`](https://henriksson-lab.github.io/zorn/reference/MergeBascetCountMatrix.md)
  : Merge a list of count matrices as produced by Bascet
- [`NoRunner()`](https://henriksson-lab.github.io/zorn/reference/NoRunner.md)
  : Create new no-runner instance, used for debugging
- [`OpenBascet()`](https://henriksson-lab.github.io/zorn/reference/OpenBascet.md)
  : Open a Bascet, prepare it for reading individual files
- [`PlotFASTQCadapterContent()`](https://henriksson-lab.github.io/zorn/reference/PlotFASTQCadapterContent.md)
  : From aggregated FASTQC data, plot adapter content
- [`PlotHistogram()`](https://henriksson-lab.github.io/zorn/reference/PlotHistogram.md)
  : Plot a histogram, loaded by ReadHistogram
- [`PlotJohnsonLindenstraussMinDim()`](https://henriksson-lab.github.io/zorn/reference/PlotJohnsonLindenstraussMinDim.md)
  : Plot minimum number of dimensions needed to retain distance between
  samples
- [`PrepareSharding()`](https://henriksson-lab.github.io/zorn/reference/PrepareSharding.md)
  : Prepare to shard reads by collecting statistics about each barcode,
  and filtering out cells with few reads
- [`ReadBascetCountMatrix()`](https://henriksson-lab.github.io/zorn/reference/ReadBascetCountMatrix.md)
  : Read a count matrix as produced by Bascet (hdf5 format). This can be
  output from both BascetQueryFq and BascetCountChrom
- [`ReadCellSNPmatrix()`](https://henriksson-lab.github.io/zorn/reference/ReadCellSNPmatrix.md)
  : Read a count matrix as produced by CellSNP, but as shards
- [`ReadHistogram()`](https://henriksson-lab.github.io/zorn/reference/ReadHistogram.md)
  : Read the count histogram associated with a Bascet. Not all Bascets
  have one, but it is typically produced after debarcoding
- [`SetTaxonomyNamesFeatures()`](https://henriksson-lab.github.io/zorn/reference/SetTaxonomyNamesFeatures.md)
  : Take a KRAKEN2 count matrix where the column is the taxonomyID.
  Convert to a matrix where the columns instead are the names of each
  taxonomy. Unused taxonomyID columns will not be kept
- [`ShowFASTQCforCell()`](https://henriksson-lab.github.io/zorn/reference/ShowFASTQCforCell.md)
  : Show the FASTQC HTML report for a cell, in the available web browser
- [`SlurmRunner()`](https://henriksson-lab.github.io/zorn/reference/SlurmRunner.md)
  : Create a runner that submits jobs to SLURM
- [`TabixGetFragmentsSeqs()`](https://henriksson-lab.github.io/zorn/reference/TabixGetFragmentsSeqs.md)
  : Using Tabix, get list of sequences in a fragment file
- [`TestBascetInstance()`](https://henriksson-lab.github.io/zorn/reference/TestBascetInstance.md)
  : Check if a Bascet instance works
- [`aggr.abricate()`](https://henriksson-lab.github.io/zorn/reference/aggr.abricate.md)
  : Callback function for aggregating ABRicate data for each cell. To be
  called from BascetAggregateMap
- [`aggr.amrfinder()`](https://henriksson-lab.github.io/zorn/reference/aggr.amrfinder.md)
  : Callback function for aggregating ABRicate data for each cell. To be
  called from BascetAggregateMap
- [`aggr.ariba()`](https://henriksson-lab.github.io/zorn/reference/aggr.ariba.md)
  : Callback function for aggregating ARIBA data for each cell. To be
  called from BascetAggregateMap
- [`aggr.fastqc()`](https://henriksson-lab.github.io/zorn/reference/aggr.fastqc.md)
  : Callback function for aggregating FASTQC data for each cell. To be
  called from BascetAggregateMap
- [`aggr.gecco()`](https://henriksson-lab.github.io/zorn/reference/aggr.gecco.md)
  : Callback function for aggregating GECCO data for each cell. To be
  called from BascetAggregateMap
- [`aggr.minhash()`](https://henriksson-lab.github.io/zorn/reference/aggr.minhash.md)
  : Callback function for aggregating min-hashes for each cell. To be
  called from BascetAggregateMap
- [`aggr.quast()`](https://henriksson-lab.github.io/zorn/reference/aggr.quast.md)
  : Callback function for aggregating QUAST data. To be called from
  BascetAggregateMap
- [`aggr.rawtext()`](https://henriksson-lab.github.io/zorn/reference/aggr.rawtext.md)
  : Callback function for just getting raw file contents To be called
  from BascetAggregateMap
- [`colSums(`*`<BascetCountMatrix>`*`)`](https://henriksson-lab.github.io/zorn/reference/colSums-BascetCountMatrix-method.md)
  : Column sums of a BascetCountMatrix (counts per feature)
- [`colnames(`*`<BascetCountMatrix>`*`)`](https://henriksson-lab.github.io/zorn/reference/colnames-BascetCountMatrix-method.md)
  : Column names of a BascetCountMatrix (feature names)
- [`` `colnames<-`( ``*`<BascetCountMatrix>`*`)`](https://henriksson-lab.github.io/zorn/reference/colnames-set-BascetCountMatrix-method.md)
  : Set column names of a BascetCountMatrix (feature names)
- [`createSlurmJobFromExisting()`](https://henriksson-lab.github.io/zorn/reference/createSlurmJobFromExisting.md)
  : This creates a job object, linking to a running command. Mainly used
  for development but can be used in case a Zorn session died and you
  want to create a new monitor
- [`dim(`*`<BascetCountMatrix>`*`)`](https://henriksson-lab.github.io/zorn/reference/dim-BascetCountMatrix-method.md)
  : Dimensions of a BascetCountMatrix
- [`formatPlainNumber()`](https://henriksson-lab.github.io/zorn/reference/formatPlainNumber.md)
  : Parse a string with a size, such as 1g, 1m, 1k, or just 123 (bytes)
- [`getBascetBinary()`](https://henriksson-lab.github.io/zorn/reference/getBascetBinary.md)
  : Get a Bascet binary for the current platform It will be cached in
  the provided directory to avoid downloading it each the time the
  function is called
- [`getBascetDockerImage()`](https://henriksson-lab.github.io/zorn/reference/getBascetDockerImage.md)
  : Get and install a Bascet docker image. It will be cached to avoid
  downloading it each the time the function is called
- [`getBascetExecutable()`](https://henriksson-lab.github.io/zorn/reference/getBascetExecutable.md)
  : Get a Bascet executable It will be cached in the provided directory
  to avoid downloading it each the time the function is called
- [`getBascetPodmanImage()`](https://henriksson-lab.github.io/zorn/reference/getBascetPodmanImage.md)
  : Get and install a Bascet podman image. It will be cached to avoid
  downloading it each the time the function is called
- [`getBascetSingularityImage()`](https://henriksson-lab.github.io/zorn/reference/getBascetSingularityImage.md)
  : Get a Bascet image (singularity or docker). It will be cached in the
  provided directory to avoid downloading it each the time the function
  is called
- [`is.BascetCountMatrix()`](https://henriksson-lab.github.io/zorn/reference/is.BascetCountMatrix.md)
  : Check that parameter is a BascetCountMatrix
- [`isBamPairedAlignment()`](https://henriksson-lab.github.io/zorn/reference/isBamPairedAlignment.md)
  : Figure out if a BAM-file is a paired alignment or not
- [`lseq()`](https://henriksson-lab.github.io/zorn/reference/lseq.md) :
  logarithmic spaced sequence; taken from emdbook library
- [`removeBascetDockerImage()`](https://henriksson-lab.github.io/zorn/reference/removeBascetDockerImage.md)
  : Remove current Bascet docker image
- [`removeBascetPodmanImage()`](https://henriksson-lab.github.io/zorn/reference/removeBascetPodmanImage.md)
  : Remove current Bascet podman image
- [`rowSums(`*`<BascetCountMatrix>`*`)`](https://henriksson-lab.github.io/zorn/reference/rowSums-BascetCountMatrix-method.md)
  : Row sums of a BascetCountMatrix (counts per cell)
- [`rownames(`*`<BascetCountMatrix>`*`)`](https://henriksson-lab.github.io/zorn/reference/rownames-BascetCountMatrix-method.md)
  : Row names of a BascetCountMatrix (cell names)
- [`` `rownames<-`( ``*`<BascetCountMatrix>`*`)`](https://henriksson-lab.github.io/zorn/reference/rownames-set-BascetCountMatrix-method.md)
  : Set row names of a BascetCountMatrix (cell names)
- [`` `[`( ``*`<BascetCountMatrix>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://henriksson-lab.github.io/zorn/reference/sub-BascetCountMatrix-ANY-ANY-ANY-method.md)
  : Subset a BascetCountMatrix

The OCTOPUS pipeline performs the following steps:

# 1. Demultiplexing

Due to the nature of the iGenomX protocol, there are effectively two demultiplexing steps. The first, automatically performed by the sequencer, reads standard Illumina indices to split your experiments up into plates. The second, handled here, uses [Fulcrum Genomic's fgbio](https://github.com/fulcrumgenomics/fgbio) to demultiplex each plate into individual wells. Note that `fgbio` will parse [src/igenomx-meta.txt](src/igenomx-meta.txt) for iGenomX's pre-specified primer indices. If you are using custom primers, please modify `src/igenomx-meta.txt` and/or the `--metadata` flag in the `Makefile`.

# 2. Read Pre-processing

To ensure a high quality _de novo_ assembly, we perform a number of processing steps. This protocol is adopted from one included with the Joint Genome Institute's [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/), and is handled by [jgi-preproc.sh](src/jgi-preproc.sh). Broadly, it removes optical duplicates, trims Illumina adapters, filters contaminants, and error-corrects the remaining reads. For this application, we filter out PhiX, a list of known sequencing artifacts (included with BBTools), and the NEB 5a genome.

## Alternative Contaminants

Users with other applications can filter out a different set of contaminants (in the form of a fasta file) by replacing `src/background.fasta` with their own `src/background.fasta`. Alternatively you can modify the `Makefile` as follows:

```
old: pipeline/%/preproc: pipeline/%/demux src/background.fasta
new: pipeline/%/preproc: pipeline/%/demux path/to/your/fasta
```

If you would like to ignore PhiX reads or the list of known artifacts, update our preprocessing script [jgi-preproc.sh](src/jgi-preproc.sh) as follows:

```
old: ref=artifacts,phix,${CONTAM_REF} \
new: ref=${CONTAM_REF} \
```

# 3. _De Novo_-based Identification

With the reads processed, we then attempt to assemble each well using [SPAdes](http://cab.spbu.ru/software/spades/). Following the JGI protocol, we attempt to merge these reads, quality trim any overlaps, and feed those to SPAdes for assembly with [jgi-denovo.sh](src/jgi-denovo.sh). Depending on your application, you may need to modify the SPAdes settings (or try a different assembler).

We then align the _de novo_ assembly products to the user specified library (the `input.fasta` that you dropped into your sequencing run folder) to identify what's in each well using [minimap2](https://lh3.github.io/minimap2/). Since we will not know the orientation of the resulting assembly, we concatenate (or flatten) our reference library before the alignment with [flatten-fasta.py](src/flatten-fasta.py). This ensures we can align the entire assembly. For example, if our plasmid (`ref` below) starts at 0, and the _de novo_ assembly (`asm` below) happens to start at 6, the resulting alignment would look like

```
ref: 01234567890123456789
asm:       6789012345
```

as opposed to

```
ref: 0123456789
asm:       6789012345
```

It should be noted that the curent version of `SPAdes` (3.13.0) produces an assembly with the same starting and ending k-mer. This will not affect the alignment (see below) but users relying on these assemblies (found in `pipeline/your-run-id/spades-contigs.fasta`) should take this into account.

```
ref: 12345678901234567890
asm:      67890123456 <- 6 is repeated!
```

# 4. Variant Calling

While aligning the _de novo_ assembly will can reveal variants, we wanted a finer-grained control over the process. Thus, we use [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/) to align the processed reads in each well to their cognate plasmid identified in the previous stemp. We then use [freebayes](https://github.com/ekg/freebayes) to call variants at positions where at least one read makes up 50% of the reads (e.g. 4 reads total, 2 call A, 2 call T -> variant call). To reduce the possibility of a sequencing error introducing a false positive, we exclude basecalls with Q<20. If you would like to adjust these parameters, please edit [denovo-guided-assembly.sh](src/denovo-guided-assembly.sh).

# 5. Quality Control

OCTOPUS provides a number of different quality control metrics to ensure the plasmids you select are correct.

## Percent Contaminants

The percentage of reads in each well that are from "contaminants" as specified in the [preprocessing](1.-read-preprocessing) step (PhiX, Illumina artifacts, and DH5a).

## Coverage

We report the percent of bases with <10x and <3x coverage. We recommend inspecting plasmids with a high percentage of bases at <3x to ensure critical regions are adequately covered. We provide `.bam` files for each well to assist in this process. For example, in `/path/to/octopus/pipeline/your-run-id/` you could run

```
samtools index plate-id/well-id.map.bam
samtools tview plate-id/well-id.map.bam lib/<reference>.fasta
```

to view a pileup of all the reads. Note that you will have to specify what reference file the well aligns to.

## Barcode Filter

Adding a barcode to each of our plasmids enables a number of useful downstream applications. For cloning, the barcode enables us to detect colonies that have multiple plasmids. To specify a barcode, simply place a string of N's as long as the barcode in the reference fasta. OCTOPUS will automatically detect barcodes declared in this fashion use them to check for plasmid contamination. First, we generate a pileup of all the reads at the barcode and collapse them at a [Levenshtein distance](https://en.wikipedia.org/wiki/Levenshtein_distance) of 1 to eliminate potential sequencing errors

```
ATGC  ---> ATGC 4
ATGC   /   TTAA 1
ATGC  /
ATGA /
TTAA
```

We can use the relative frequencies of the barcodes to determine if the well is contaminated. Importantly, the iGenomX protocol will have a low level of template switching. This will result in a large amount of unique barcodes with few reads

```
AACCTTGGCCTTAA 50
ATGCATTACAGACA 5
TTACCATTCATGAT 2
AGGGACCGATTAGC 1
GGTATTAGGCCATA 1
CTATAGCATTGCAT 1
...
```

True contamination typically results in a secondary barcode with many reads

```
AACCTTGGCCTTAA 50
ATGCATTACAGACA 10
TTACCATTCATGAT 2
AGGGACCGATTAGC 1
GGTATTAGGCCATA 1
CTATAGCATTGCAT 1
...
```

With this in mind, we've developed a filtering heuristic that empirically eliminates wells that are actually contaminated without being too conservative. Specifically, if the second most common barcode is >10% of the total number of reads (10/65 ~15% in the above example), we call that well contaminated. If the second most common barcode is <4% of the total reads, it's most likely template swapping and we call the well clean. Lastly, if the second most common barcode is >4% and <10% of the reads (5/60 ~8% in the first example), we check the ratio of the second and third most common barcodes is < 2 (5/2 > 2). This ratio test is designed to capture our observation that true plasmid contaminations are a distinct population, while barcodes from template swapping are uniformly represented with low counts. In the case where there is only two barcodes, we will call the well contaminated if the second barcode makes up >4% of the reads.
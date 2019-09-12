# OCTOPUS

A scalable cloning system.

- [Getting Started](#getting-started)
- [Running OCTOPUS](#running-octopus)
- [Details](#details)
- [Contributing](#contributing)
- [License](#license)

# Getting Started

## Prerequisites

### Docker

The computational pipeline is contained in a [docker image](https://hub.docker.com/octantbio/octopus). We strongly recommend using it for your analyses. Follow the official docs for: [MacOS](https://docs.docker.com/docker-for-mac/install/), [Linux](https://docs.docker.com/install/linux/docker-ce/ubuntu/), or [Windows](https://docs.docker.com/docker-for-windows/install/).

Next, choose one of the following:

**A.** Pull from DockerHub

```
docker pull octant/octopus
```

**B.** Build it yourself

```
git clone https://github.com/octantbio/octopus.git
cd octopus/docker
docker build .
```

It is possible to run OCTOPUS without `docker`, but not recommended. If you insist, see the [dockerfile](docker/Dockerfile) for help.

### Hardware

For optimal performance, we recommend deploying on a machine with at least 64 GB of RAM and 16 cores. We use [GNU Parallel](https://www.gnu.org/software/parallel/) to distribute tasks where possible. If you run out of RAM, you can tune the number of jobs by setting the `--jobs` flag (e.g. 75% of cores `--jobs 75%` or 4 simultaneous jobs `--jobs 4`)

# Running OCTOPUS

## Set Up

Clone this repository

```
git clone https://github.com/octantbio/octopus.git
```

or unzip the [latest release](https://github.com/octantbio/octopus/releases) tarball.

## Linking Data

As OCTOPUS uses `make` to orchestrate everything, there are some conventions your data must adhere to. First, deposit the output folders from your sequencing run must be deposited in the [data](data) directory

```
cp -r /path/to/run-id /path/to/octopus/data
```

under a unique folder. We typically use the default folder name produced by the sequencer as an identifier. Second, many steps in the OCTOPUS pipeline will process the file name of the fastq's. To avoid issues, make sure all unique information is contained before the first underscore in your SampleSheet (most Illumina sequencers will automatically convert any \_'s in the `Sample_Name` column of the SampleSheet to -'s anyways). Importantly, the pipeline will trim out anything between the first underscore and the read specifier (e.g. `my-reads_foo_bar_baz_R1.fastq.gz -> my-reads_R1.fastq.gz`) to ensure everything behaves properly.

Alternatively, you can manually add fastq's to the [pipeline](pipeline) folder provided they are in their own folder and are not symlinks to outside of the `octopus` folder (if you are following our docker instructions).

### Reference Library

If you know the sequences of you plasmids _a priori_ you should include them as `./octopus/data/run-id/input.fasta`. The OCTOPUS pipeline will also automatically parse any barcodes in the form of N's for downstream analyses. If you do not know your input, run `make de-novo` instead to take the pipeline through the _de novo_ assembly step.

## The Pipeline

With the data in place, we can drop into our docker image with

```
docker run --rm -it -v /path/to/your/octopus:/root/octopus octant/octopus /bin/bash
```

This links your octopus folder (`/path/to/your/octopus`) to the docker image (`/root/octopus`). Note that `--rm` makes the image ephemeral so anything written outside of the octopus directory will be lost if you logout of the shell. Next, we can

```
cd /path/to/your/octopus
make all
```

to run the pipeline on every sequencing run in the `./data` directory and produce `octopus/pipeline/run-id/aggregated-stats.tsv`

### aggregated-stats.tsv

As the name suggests, results pertinent to an OCTOPUS run are aggregated into a `tsv` file for your analysis. The columns are:

- `Run`: Illumina run ID
- `Plate`: plate ID
- `Well`: well address
- `Plate_Well`: unique plate\_well identifier
- `DeNovo_Ref`: well identity based on aligning _de novo_ assembly to reference library
- `CIGAR`: CIGAR string from aligning the _de novo_ assembly to `DeNovo_Ref`
- `LT_10`: percentage of input reference with < 10x coverage (ideally close to 0)
- `LT_3`: percentage of input reference sequence with < 3x coverage (if not 0 inspect read pileup)
- `BC_Contam`: are there multiple plasmids in this well (TRUE/FALSE)? [More details](#barcode-filter)
- `n_vars`: number of variants detected by FreeBayes (note barcodes count as variants)
- `n_barcodes`: number of barcodes detected
- `expected_bcs`: expected number of barcodes based on the reference (in a perfect plasmid `n_vars = n_barcodes = expected_bcs`)
- `bc_1`: sequence of barcode 1 pulled from the variant caller (may be reverse compliment)
- `pos_1`: position of barcode 1 in _de novo_ assembly
- `bc_N`: sequence of barcode N pulled from the variant caller (may be reverse compliment; NA if missing)
- `pos_N`: position of barcode N in _de novo_ assembly (NA if missing)
- `Contaminants`: number of reads from "contaminants" ([more details](#percent-contaminants))
- `Leftover`: number of reads in well leftover after filtering out "contaminants"
- `Percent_aligned`: percentage of Leftover reads that align with the reference sequence.
- `Contig`: the _de novo_ assembly. Note the first and last N bases (often 55 or 125) are repeated
- `Ref_Seq`: sequence that _de novo_ assembly aligns to

# Details

This will perform the following steps:

## 1. Demultiplexing

Due to the nature of the iGenomX protocol, there are effectively two demultiplexing steps. The first is automatically performed by the sequencer and reads standard Illumina indices to split your experiments up (e.g. one index per plate of samples). The second, handled here, uses [Fulcrum Genomic's fgbio](https://github.com/fulcrumgenomics/fgbio) to demultiplex each plate into individual wells. `fgbio` will parse [src/igenomx-meta.txt](src/igenomx-meta.txt) for iGenomX's pre-specified primer indices. If you are using custom primers, please modify `src/igenomx-meta.txt` and/or the `--metadata` flag in the `Makefile`.

## 2. Read Pre-processing

To ensure a high quality _de novo_ assembly, we perform a number of processing steps. This protocol is adopted from one included with the Joint Genome Institute's [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/), and is handled by [jgi-preproc.sh](src/jgi-preproc.sh). Broadly, it removes optical duplicates, trims Illumina adapters, filters contaminants, and error-corrects the remaining reads. For this application, we filter out PhiX, a list of known sequencing artifacts (included with BBTools), and the NEB 5a genome. Users with other applications can filter out a different set of contaminants (in the form of a fasta file) by modifying the `Makefile` as follows

```
old: pipeline/%/preproc: pipeline/%/demux src/ecoli.fasta
new: pipeline/%/preproc: pipeline/%/demux path/to/your/fasta
```

If you would like to ignore PhiX reads or the list of known artifacts, update our preprocessing script [jgi-preproc.sh](src/jgi-preproc.sh) as follows:

```
old: ref=artifacts,phix,${CONTAM_REF} \
new: ref=${CONTAM_REF} \
```

## 3. _De Novo_-based Identification

With the reads processed, we then attempt to assemble each well using [SPAdes](http://cab.spbu.ru/software/spades/). Following the JGI protocol, we attempt to merge these reads, quality trim any overlaps, and feed those to SPAdes for assembly in [src/jgi-denovo.sh](src/jgi-denovo.sh). Depending on your application, you may need to modify the SPAdes settings (or try a different assembler).

We then align the _de novo_ assembly products to a user specified library to identify what's in each well using [minimap2](https://lh3.github.io/minimap2/). Since we will not know the orientation of the resulting assembly, we concatenate (or flatten) our reference library before the alignment. This ensures we can align the entire assembly.

```
ref: 1234567890#1234567890
asm:       6789012345
```

It should be noted that the curent version of `SPAdes` (3.13.0) produces an assembly with the same starting and ending k-mer. This will not affect the alignment (see below) but users relying on these assemblies (found in `pipeline/your-run-id/spades-contigs.fasta`) should take this into account.

```
ref: 1234567890#1234567890
asm:       67890123456 <- 6 is repeated!
```

## 4. Variant Calling

While aligning the _de novo_ assembly will reveal any variants, we want a finer-grained control over the process. Thus, we use [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/) to align the processed reads in each well to their cognate plasmid identified by the _de novo_ assembly. We then use [freebayes](https://github.com/ekg/freebayes) to call variants at positions where at least one read makes up 50% of the reads (e.g. 4 reads total, 2 call A, 2 call T -> variant call). To reduce the possibility of a sequencing error introducing a false positive, we exclude basecalls with Q<20. If you would like to adjust these parameters, please edit [src/denovo-guided-assembly.sh](src/denovo-guided-assembly.sh).

## 5. Quality Control

OCTOPUS provides a number of different quality control metrics to ensure the plasmids you select are correct.

### Percent Contaminants

The percentage of reads in each well that are from "contaminates" as specified in the [preprocessing](1.-read-preprocessing) step (PhiX, Illumina artifacts, and DH5a).

### Coverage

We report the percent of bases with <10x and <3x coverage. We reccomend inspecting plasmids with a high percentage of bases at <3x to ensure critical regions are adequately covered. We provide `.bam` files for each well to assist in this process. For example, in `/path/to/octopus/pipeline/your-run-id/` you could run

```
samtools index plate-id/well-id.map.bam
samtools tview plate-id/well-id.map.bam lib/<reference>.fasta
```

to view a pileup of all the reads. Note that you will have to specify what reference file the well aligns to.

### Barcode Filter

Adding a barcode to each of our plasmids enables a number of useful downstream applications. For cloning, the barcode enables us to detect colonies that have multiple plasmids. To specify a barcode, simply place a string of N's as long as the barcode in the reference fasta. OCTOPUS will automatically detect barcodes declared in this fashion use them to check for plasmid contamination. First, we generate a pileup of all the reads at the barcode and collapse them at a [Levenshtein distance](https://en.wikipedia.org/wiki/Levenshtein_distance) of 1 to eliminate potential sequencing errors

```
ATGC  ---> ATGC 4
ATGC   /   TTAA 1
ATGC  /
ATGA /
TTAA
```

We can use the relative frequencies of the barcodes to determine if the well is contaminated. Importantly, the IgenomX protocol will have a low level of template switching. This will result in a large amount of unique barcodes with few reads

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

With this in mind, we've developed a filtering heuristic that empirically eliminates wells that are actually contaminated without being too conservative. Specifically, if the second most common barcode is >10% of the total number of reads (10/65 ~15% above), we call that well contaminated. If the second most common barcode is <4% of the total reads, its most template swapping and we call the well clean. Lastly, if the second most common barcode is >4% and <10% (e.g. 5/60 ~8% in the first example), we check the ratio of the second and third most common barcodes is < 2 (e.g. 5/2 > 2). This ratio test is designed to capture our experience of true plasmid contaminations being represented as a distinct population, while barcodes from template swapping primarily being lowly represented. In the case where there is only two barcodes, we will call the well contaminated if the second barcode makes up >4% of the reads.

# Contributing

Please feel free to open an issue or pull request.

# License

This project is licensed under the Apache 2.0 License - see the [LICENSE](LICENSE) file for details. Additional licensing information:

- [BBTools](docker/bbtools-license)
- [mlr](docker/mlr-license)
- [starcode](docker/starcode-license)


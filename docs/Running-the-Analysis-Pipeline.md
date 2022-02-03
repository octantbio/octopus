# 1. Linking Data

As OCTOPUS uses `make` to orchestrate everything, there are some conventions your data must adhere to. First, deposit the output folder from your sequencing run in the `data` directory

```
cp -r /path/to/run-id /path/to/octopus/data
```

under a unique folder. We typically use the default folder name produced by the sequencer as an identifier. Second, many steps in the OCTOPUS pipeline will process the file name of the fastq's. To avoid issues, make sure all unique information is contained before the first underscore in your SampleSheet (most Illumina sequencers will automatically convert any \_'s in the `Sample_Name` column of the SampleSheet to -'s anyways). Importantly, the pipeline will trim out anything between the first underscore and the read specifier (e.g. `my-reads_foo_bar_baz_R1.fastq.gz -> my-reads_R1.fastq.gz`) to ensure everything behaves properly.

Alternatively, you can manually add fastq's under `./octopus/pipeline/*run-id*/fastqs` provided they are not symlinks to outside of the `octopus` folder (if you are following our docker instructions).

## Reference Library

Next, place a fasta file containing the sequences of the plasmids you are trying to sequence under `./octopus/data/*run-id*/input.fasta`. The OCTOPUS pipeline will also automatically parse any barcodes in the form of N's for downstream analyses.

Similar to the fastq's, you can manually place `input.fasta` at `./octopus/pipeline/*run-id*/input.fasta`.

### De Novo Assembly

If you do not know your input, run `make de-novo` instead to take the pipeline through the _de novo_ assembly step. If you forget, and run `make all` the pipeline will throw an error.

# 2. Running the Pipeline

After getting the data in place, make sure you `cd` into your octopus folder. From there we can drop into our docker image with

```
docker run --rm -it -v "$(pwd)":/root/octopus octant/octopus /bin/bash
```

This links your octopus folder (`/path/to/your/octopus`) to the docker image (`/root/octopus`). Note that Docker requires you to specify the *absolute* path to the folder (`$(pwd)` is a handy shortcut to do that for you). Also, be aware that that `--rm` makes the image ephemeral so anything written outside of the octopus directory will be lost if you logout of the shell. From the Docker image, we can

```
cd octopus
make all
```

to run the pipeline on every sequencing run in the `./data` directory and produce `octopus/pipeline/*run-id*/aggregated-stats.tsv`. You will get an error if you did not place the `input.fasta` file under `data/*run-id*/input.fasta`. If you don't have one try `make denovo`.

## aggregated-stats.tsv

As the name suggests, results pertinent to an OCTOPUS run are aggregated into a `tsv` file for your analysis. The columns are:

- `Run`: Illumina run ID
- `Plate`: plate ID
- `Well`: well address
- `Plate_Well`: unique plate\_well identifier
- `DeNovo_Ref`: well identity based on aligning _de novo_ assembly to reference library
- `CIGAR`: CIGAR string from aligning the _de novo_ assembly to `DeNovo_Ref`
- `LT_10`: percentage of input reference with < 10x coverage (ideally close to 0)
- `LT_3`: percentage of input reference sequence with < 3x coverage (if not 0 inspect read pileup)
- `BC_Contam`: are there multiple plasmids in this well (TRUE/FALSE)? ([more details](https://github.com/octantbio/octopus/wiki/Pipeline-Details#barcode-filter))
- `n_vars`: number of variants detected by FreeBayes (note barcodes count as variants)
- `n_barcodes`: number of barcodes detected
- `expected_bcs`: expected number of barcodes based on the reference (in a perfect plasmid `n_vars = n_barcodes = expected_bcs`)
- `bc_1`: sequence of barcode 1 pulled from the variant caller (may be reverse complement)
- `pos_1`: position of barcode 1 in _de novo_ assembly
- `bc_N`: sequence of barcode N pulled from the variant caller (may be reverse complement; NA if missing)
- `pos_N`: position of barcode N in _de novo_ assembly (NA if missing)
- `Contaminants`: number of reads from "contaminants" ([more details](https://github.com/octantbio/octopus/wiki/Pipeline-Details#alternative-contaminants))
- `Leftover`: number of reads in well leftover after filtering out "contaminants"
- `Percent_aligned`: percentage of Leftover reads that align with the reference sequence.
- `Contig`: the _de novo_ assembly. Note the first and last N bases (often 55 or 125) are repeated
- `Ref_Seq`: sequence that _de novo_ assembly aligns to

# 3. Picking perfects

One way you can analyze the results is by pasting the `aggregated-stats.tsv` into a spreadsheet

1. If applicable, filter out any "TRUE" values under `BC_Contam`
2. If applicable, flag or filter out any duplicate barcodes
3. Filter out any unexpected variants. The pipeline will automatically detect any strings of N's in the `input.fasta` and report the number of `expected_bcs` for that reference. Perfect clones should have `expected_bcs = n_barcodes = n_vars`
4. Ensure that there is adequate coverage by checking `LT_10` and `LT_3`. We recommend only picking wells with `LT_3 = 0`. You can be more conservative by using `LT_10` to specify your cutoff. (For a 10kb plasmid an `LT_3` of 0.001 means that 10 bp of the plasmid did not have a coverage of at least three).
5. If there happens to be a clone that does not have sufficient coverage (by `LT_10` or `LT_3` but is absolutely required, use `samtools tview` to manually inspect the read pileup in critical areas of your plasmid:
    1. In a new terminal, `cd` into your octopus directory
    2. Open up a new docker instance - `docker run --rm -it -v "$(pwd)":/root/octopus octant/octopus /bin/bash`
    3. Navigate into the folder that contains the analyzed data from that run - `cd octopus/pipeline/your_run_id`
    4. View the pileup - `samtools tview your_plate/your_well.map.bam lib/your_ref.fasta`
    - Note `your_ref` will be `DeNovo_Ref` in `aggregated-stats.tsv`
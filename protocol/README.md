# Equipment/Reagents:

  - [Riptide Kit](https://igenomx.com/product/riptide/) from iGenomX
  - ddH2O
  - 30% [Glycerol](https://www.fishersci.com/shop/products/glycerol-molecular-biology-fisher-bioreagents-2/BP2291)
  - 0.1 M [Sodium Hydroxide](https://www.fishersci.com/shop/products/sodium-hydroxide-pellets-certified-acs-fisher-chemical-7/S318100)
  - [2xYT](https://us.vwr.com/store/product/7437420/vwr-life-science-2xyt-medium-broth)
  - [Flat Bottom 96-well plates](https://www.sigmaaldrich.com/catalog/product/aldrich/br781602?lang=en&region=US)
  - [96-well PCR plates](https://www.thermofisher.com/order/catalog/product/AB0600)
  - [96-well thermal cycler](https://www.bio-rad.com/en-us/product/c1000-touch-thermal-cycler?ID=LGTW9415)
  - 384-well PCR plates (optional)
  - [384-well thermal cycler](https://www.thermofisher.com/order/catalog/product/4388444) (optional)
  - Foil Plate seal (Axygen PCR-AS-200)
  - Plastic Plate seal (Axygen PCR-SP)
  - Multichannel 10
  - Multichannel 200
  - [Liquidator 20](https://www.shoprainin.com/Products/Pipettes-and-Tips/Pipettes/High-throughput-Pipetting/Liquidator%E2%84%A2-96/Liquidator-96-2-20-%C2%B5L-LIQ-96-20/p/17014207) (optional)
  - [Liquidator 200](https://www.shoprainin.com/Products/Pipettes-and-Tips/Pipettes/High-throughput-Pipetting/Liquidator%E2%84%A2-96/Liquidator96%2C-5-200-%C2%B5L-LIQ-96-200/p/17010335) (optional)
  - 37˚C shaking Incubator
  - Nanodrop, Qubit, BioAnalyzer or Tapestation
  - [Illumina Miseq](https://www.illumina.com/systems/sequencing-platforms/miseq.html)
  - [Miseq sequencing kits](https://www.illumina.com/products/by-type/sequencing-kits/cluster-gen-sequencing-reagents/miseq-reagent-kit-v2.html)
  - [Dynamag](https://www.thermofisher.com/order/catalog/product/12321D) (or magnetic tube stand)
  - [Multi-tube vortexer](https://ru.vwr.com/store/product/596528/multi-tube-vortexers) (optional)

# Table of Contents

  - [Picking colonies](#picking-colonies) (1hr) (o/n growth)
  - [Lysate and glycerol stock generation](#lysate-and-glcerol) (1hr)
  - [Miniaturized Riptide](#"minaturized") (4-6hrs total)
      - Random primer extension and biotinylated termination: “A” reaction (1-2hrs)
      - DNA capture and library conversion: “B” reaction (1-2hrs)
      - Library Amplification (1-2hrs)
      - Size Selection (1hr)
  - Library quantification and sequencing (1hr) (16-28hr sequencing)
  - Running the computational pipeline (1-2hrs see github for processor times)

Interpreting the data and picking perfects

# Picking colonies \[~10-15min per 96 colonies\] – O/N growth

  - Pick single colonies into 150 µl of 2xYT with appropriate antibiotic (this is important for good plasmid yield in the lysate) in a flat bottom 96-well plate. When cloning in a pooled library format, we usually pick a minimum of 4 colonies per desired clone.
      - *Tip: We pick colonies using 200 µL pipette tips, put them back into the box, and then use a multichannel pipette to inoculate the 2xYT plate. This saves a lot of time and should allow you to pick a 96-well plate in about ~10 minutes with some practice (potential upgrade is a colony picker - we have yet to identify a cost effective and robust model*)
  - Tape the lid tightly to the plate to reduce evaporation and grow overnight at 37˚C in a shaking incubator (a 96-well plate specific incubator with a short throw and high rpm is great, but we just use a normal incubator/shaker meant for flasks at 250 rpm)

# Lysate and glycerol stock generation \[~1-2 hours\] (Day 2)

  - Pull 50 µL of resuspended culture from each well using a liquidator/multichannel into a 96-well PCR plate and pellet the bacteria by spinning at 4000 rcf for 10 min at RT in a swinging bucket rotor.
  - Store the remaining 100 µL as glycerol stock by adding 100 µL of 30% glycerol to the plates, shake in incubator for 15 min, foil seal, and store in -80˚C freezer.
  - Flick and blot the media from the pelleted bacteria (media in this step will inhibit proper lysis)
  - Add 40 µL of MQ water to each pellet, plastic seal, and resuspend pellet by pipetting or with a multi-tube vortexer (tip: pulse by going from medium to max speed, then leave at near max speed for 5-10 secs)
  - Lyse cells in 96-well thermocycler by heating at 95˚C for 3 min and then cooling to 4˚C
  - Clarify lysate by spinning at 4000+ rcf for 10 min at RT
  - Remove the top 20 µLof clarified lysate and store in a separate 96 or 384-well plate
  - *At this point the protocol can be paused by freezing the lysate at -80˚C*

# “Miniaturized” Riptide protocol (~8-12 hours)

This protocol is adapted from the iGenomX Riptide library prep protocol [here](https://igenomx.com/resources/)

## STEP 1: Random primer extension and biotinylated termination; “A” Reaction

  - For each 96-well plate, prepare 200 µL of Reaction A master-mix by mixing
      - 100 µL DNTP Mix 1
      - 50 µL 10X Enzyme 1 Buffer
      - 50 µL Enzyme 1
  - Fill a 96-well (or 384 if working with 4 96-well plates) plate with 2 µL each of the master-mix using a multichannel or liquidator
  - Pipette into each well 1 µL of Primer A, using liquidator or multichannel (make sure to properly mix primer A plate if thawing
  - Pipette into each well 2 µL of clarified lysate using liquidator or multichannel (make sure to properly mix lysate if thawing (If sequencing mini-preps, use 2 µL of 5 ng/µL per well)
  - Seal the plate with foil and run the following protocol on a thermocycler
      - 92˚C for 3 min
      - 16˚C for 5 min
      - Slow ramp (0.1˚C/sec to 68˚C)
      - 68˚C for 15 min
      - Hold at 4˚C
  - At this point the protocol can be paused by freezing the plate(s) at -20˚C
  - Warm the SPRI Beads I to RT.
  - This next step is stopping the reaction with EDTA and pooling the contents of each 96-well plate
      - We use a liquidator to add 1 µL of 75 mM EDTA to each well then use those same tips to aspirate and dispense the contents of each 96-well plate onto a clean liquidator tip-box lid and pool by tapping the lid at an angle, using the pooled liquid to wash and collect any drops that are still stuck, and recovering into a low-retention Eppendorf tube (should recover around 400-500 µL of sample).
      - *Alternatively: pool samples into a pcr strip tube already containing EDTA, using a multichannel*
      - After this step, each 96-well plate will have it’s own pooled tube
  - Add 1.8 volumes of well resuspended RT SPRI Beads I to each pooled tube (use a 1 mL pipette to measure the volume). Mix by pipetting, incubate 10 min at RT.
  - Place tube(s) on dynamag, allow the solution to clear (2-5 min), discard supernatant.
  - Keeping tube(s) on the dynamag, add 1300 µL of “freshly” prepared 80% ethanol to each tube, wait 30 secs, then aspirate ethanol.
      - Repeat this step, and ensure that all ethanol is removed.
  - Open the caps and allow the beads to air dry on the dynamag for 10 min, don’t overdry.
  - Add 50 µL of RT 10 mM Tris-HCl pH 8 to the beads,
  - Remove from dynamag, and resuspend the beads until homogenous. Incubate at RT for 10 min to elute.
  - Allow the beads to clear on dynamag for 2 min, then transfer clarified elution to new low retention PCR tube(s)
  - At this point the protocol can be paused by freezing tubes in the -20

## STEP 2: DNA Capture and Library Conversion; “B” Reaction

  - Heat-denature the elution at 95C for 3 min and hold at 4C in a thermocycler
  - While it’s heating, prep the Capture Beads
      - Warm the Capture Beads and HS-Buffer to RT, resuspend thoroughly, and transfer 20 µL of slurry for each sample into a new PCR tube.
      - Dynamag the capture beads and discard the supernatant.
      - Remove tubes from dynamag and add 100 µL of HS-Buffer
      - Dynamag and remove the wash.
      - Resuspend the beads in 20 µL of HS-Buffer.
  - Add all 50 µL of the heat-denatured elution to the washed Capture Beads, mix, and incubate at RT for 10 min.
  - Pipette mix beads again and incubate at RT for 10 min.
  - Dynamag the beads and discard the supernatant.
  - Resuspend beads with 50 µL of 0.1M NaOH, incubate 4 min at RT, dynamag and remove sup
  - Wash the beads 3 times (resuspend beads in 100 µL of RT Bead Wash Buffer, dynamag, remove supernatant). Make sure to remove any remaining liquid after final wash.
  - Prepare a mastermix for “Reaction B”: for every tube of beads, mix together on ice
      - 4 µLof 5x Enzyme II Buffer
      - 1.5 µL DNTP mix II
      - 2 µL Primer B
      - 12 µL Nuclease-Free Water
      - 0.5 µL Enzyme II
  - Quickly add 19.5 µL of “Reaction B” mastermix to each tube (try not to let mastermix sit around at RT)
  - Incubate the tubes in a thermocycler for 20 min at 24˚C and hold at 4˚C for at least 3 min
  - Pipette mix the beads, dynamag and discard supernatant
  - Wash the beads 3 times (resuspend beads in 100 µL of RT Bead Wash Buffer, dynamag, remove supernatant). Make sure to remove any remaining liquid after final wash.

## STEP 3: Amplification; “PCR”

  - Resuspend the beads in 21 µL of nuclease free water, then setup the PCR reaction by adding:
      - 2 µL Universal PCR primer
      - 2 µL Index PCR primer (1-12) (choose one barcoded primer per pool)
      - 25 µL 2X PCR Amplification Mix
      - 50 µL total
  - Input the following program into a thermocycler

```
1 cycle:   98˚C, 2 min
11 cycles: 98˚C, 20 sec
           60˚C, 30 sec
           72˚C, 30 sec
1 cycle:   72˚C, 5 min
           4˚C, hold
```

  - Record which samples received which Index PCR primer
  - Sample can be left in the thermocycler at 4˚C overnight.
  - Briefly spin the PCR tube in a picofuge, dynamag and transfer the supernatant to new low retention eppendorf 1.5 mL tubes. Discard the PCR tubes containing the Capture Beads.

## STEP 4: Size Selection; SPRI bead cleanup and gel-isolation

  - Add 70 µL of well resuspended RT SPRI Beads II to the samples. Mix well and incubate at RT for 10 min.
  - During this incubation, pour a 1.0% agarose gel with Sybr-Safe for gel-extraction
  - Dynamag the beads for at least 2 min, and discard supernatant.
  - Add 200 µL of 80% ethanol to each tube(s), wait 30 sec, then remove and discard the ethanol. (it is unnecessary to remove from magnet for the ethanol wash)
  - Repeat with another 200 µL of 80% ethanol, carefully remove all ethanol from tube without disturbing beads.
  - Open cap and allow to air dry for 10 min on dynamag (careful not to overdry).
  - Add 25 µL of RT 10 mM Tris-HCl (pH8) to beads. Remove from dynamag and fully resuspend the beads.
  - Incubate at RT for 10min to elute, place back on dynamag and transfer supernatant to new low-retention tubes.
  - Run sample(s) on a 1.0% agarose gel and gel extract. It should come out as a visible smear, isolate the 400-1200bp region, taking care to avoid the potential primer dimer band.
  - At this point the protocol can be paused by putting the tube(s) in the -20

# Library Quantification and Sequencing

  - Every quantification method is slightly different and each comes with its own pros and cons. These issues are exacerbated by the Riptide library product which is a heterogenous mixture of DNA fragment sizes (*see above*).
  - The way we reduce variability in quantification is to quantify a previously run OCTOPUS library and use that as a baseline to estimate the concentration of the current libraries. In our hands, a fluorescence based assay gives sufficiently accurate quantification alongside a previously run library. If no previously run libraries are available, quantification by qPCR has been the most accurate method.
  - Select the appropriate sequencing kit: as a rule of thumb, 10,000 paired-end 150 reads/well, or 1,000,000 reads per 96-well plate (nano-kit V2 for one 96-well plate, micro-kit V2 for a 384-well plate, and a standard V2 kit for anything more than 384-well, the most we do is 3 x 384-wells due to the limited number of index primers in the riptide kit)
  - The only information necessary to generate the sample sheet is to provide the plate index sequences
  - ***Example sample-sheet.csv:***

```
[Header],,,,,,
IEMFileVersion,5,,,,,
Investigator Name,Octonaut,,,,,
Experiment Name,20190826_OCTOPUS_plate001-012,,,,,
Date,8/26/19,,,,,
Workflow,GenerateFASTQ,,,,,
Application,FASTQ Only,,,,,
Instrument Type,MiSeq,,,,,
Assay,Nextera DNA,,,,,
Index Adapters,"Nextera Index Kit (24 Indexes, 96 Samples)",,,,,
Description,Test on known samples,,,,,
Chemistry,Amplicon,,,,,
,,,,,,
[Reads],,,,,,
151,,,,,,
151,,,,,,
[Settings],,,,,,
,,,,,,
[Data],,,,,,
Sample_ID,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2
OCTOPUS-plate001,plate001,,A001,ATCACG,,
OCTOPUS-plate002,plate002,,A002,CGATGT,,
OCTOPUS-plate003,plate003,,A003,TTAGGC,,
OCTOPUS-plate004,plate004,,A004,TGACCA,,
OCTOPUS-plate005,plate005,,A005,ACAGTG,,
OCTOPUS-plate006,plate006,,A006,GCCAAT,,
OCTOPUS-plate007,plate007,,A007,CAGATC,,
OCTOPUS-plate008,plate008,,A008,ACTTGA,,
OCTOPUS-plate009,plate009,,A009,GATCAG,,
OCTOPUS-plate010,plate010,,A010,TAGCTT,,
OCTOPUS-plate011,plate011,,A011,GCCTAC,,
OCTOPUS-plate012,plate012,,A012,CTTGTA,,
```

  - The sequencing run will take ~16-30 hours depending on the size of the kit (nano, micro, standard)

# Running the computational pipeline

To run the computational pipeline, you must first set up the OCTOPUS docker image and pull the source code. Detailed instructions are available on our [GitHub](https://github.com/octantbio/octopus). Once this is set up and your sequencing run is finished, transfer/copy the raw data folder into the `/octopus/data` folder. Next, you will need to generate a fasta file named `input.fasta` that contains all possible sequences that are expected to be retrieved from the sequencing run. This is not necessary if you only plan on performing *de novo* assembly. Note this is outside the scope of the standard pipeline and will not include the typical quality control metrics.

**Example input.fasta:**
```
>Plasmid1
CTTGAAGTAATGTATACGACAGAGTCCGTGCACCTACCAAACCTCTTTAGTCTAAGTTCAGACTAGTTGGAAGTTTGTCTAGATCTCAGATTTTGTCACTAGAGGACGCACGCTCTATTTTTATGATCCATTGATGTCCCTGACGCTGCAAAATTTGCAACCAGGCAGTCTTCGCGGTAGGTCCTAGTGCAATGGGGCTTTTTTTCCATAGTCCTCGAGAGGAGGAGACGTCAGTCCAGATATCTTTGATGTCGTGATTGGAAGGACCCTTGGCCCTCCACCCTTAGGCAGTGTATACTCTTCCATAAACGGGCTATTAGTTATGAGGTCCGAAGATTGAAAAAGGTGAGGGAACTCGGCCGAACGGGAAAGACGGACATCTAGGCAACCTGACCACGGTTGCGCGTCCGTATCAAGGTCCTCTTAATAGGCCCCCGTTACTGTTGGTCGTAGAGCCCAGAACGGGTTGGCCAGATGTGCGACAATA
>Plasmid2
CGCTTAGTCGCTCTTGGGCCGCGGTGCGCTACCTTGCAGGAATTGAGACCGTCCGTTAATTTCCCTTGCATATATATTGCGTTTCTTTGACCTTTTAACCGCTCTCTTAGAAGAGAGACAGATAGCTTCTTACCGGTGCGCCACCGTAGGCAGTACGATCGCACGCCCCATGTGAACGATTGGTAAACCCAGTGTCCTGTGAGCGACAAAAGCTTAAATGGGAAATACGCGCCCATAACTTGGTGCGAATACGGGTCGTAGCAATGTTCGTCTGACTATGATCTACATATTACAGGCGGTACGTCTGCTTTGGTCAGCCTCTAATGGCTCGTAAGATAGTGCAGCCGCTGGTGATCACTCGATGACCTCGGCTCCCCATTGCAACTACGGGGATTCTTGGAGAGCCAGCTGCGTTCGCTAATGTGAGGACAGTGTAGTATTAGCAAACGATAAGTCCCGAACTGGTTGTGACCTAACGAAAAGTGAACTTCATAATACGTGCTGTCCCACGC
>Plasmid3-with_1BC
ATACTCTCGTAGTTAACATCTAGCCCGGCCCTATCAGTACAGCAGTGCCTTGAATGACATACTCATCATTAAATTTTCTCTACAGCCAAACGACCAAGTGCATTTCCAGGGAGCGCGATGGAGATTCATTCTCTCGCCAGCACTGTAATAGGCACTAAAAGAGTGATGATAATCATGAGTGCCGCGCTAAGGTGGTGTCGGAACAAAGCGGTCTTACGGTCAGTCGTATTCCTTCTCGAGTTCCGTCCAGTTGAGCGTGTCACTCCCAGTGTACCTGCAAGCCGAGATGGCTGTGCTTGGAGTCAATCGCATGTAGGATGGTCTCCAGACACCGGGGCACCAGTTTTCACGNNNNNNNNNNNNNNNAATACCTGGAGCTGTACCGTTATTGCGCTGCATAGATGCAGTGCTGCTCTTATCACATTTGTTTCGACGACAGCCGCCTTCGCAGTTTCCTCAGACACTTAAGAATAAGCGCTTATTGTAGGCAGAGGCACGCCCTATTAGTGGCTGCGGCAAAATATCTTCGGATCCCCTTGTCCAACCAAATTGATCGAATTCTTTCATTTAAGACCCTAATATGTCATCATTAGTGATTAAATGCCACTCCGAAAATACCGCCTAGAAATGTCTAAGAT
```

After preparing the directory, start a docker image in the `octopus` folder and type `make all` in the command prompt to run the pipeline. This will take 1-2 hours depending on your computational resources. If you have multiple sequencing runs, `make` will run the pipeline on them sequentially. Users without an `input.fasta` must type `make denovo` to run the pipeline through the *de novo* assembly step.


# Interpreting the output

The computational pipeline will automatically generate a new folder with the name of the run-file in the /pipeline/ folder. This folder will contain the results of the analysis and outputs a ‘aggregated-stats.tsv’ file that contains contains the following columns:

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
- `bc_1`: sequence of barcode 1 pulled from the variant caller (may be reverse complement)
- `pos_1`: position of barcode 1 in _de novo_ assembly
- `bc_N`: sequence of barcode N pulled from the variant caller (may be reverse complement; NA if missing)
- `pos_N`: position of barcode N in _de novo_ assembly (NA if missing)
- `Contaminants`: number of reads from "contaminants" ([more details](#percent-contaminants))
- `Leftover`: number of reads in well leftover after filtering out "contaminants"
- `Percent_aligned`: percentage of Leftover reads that align with the reference sequence.
- `Contig`: the _de novo_ assembly. Note the first and last N bases (often 55 or 125) are repeated
- `Ref_Seq`: sequence that _de novo_ assembly aligns to

# Picking perfects

We analyze our data by pasting the aggregated-stats.tsv into a spreadsheet.

1. If applicable, filter out any "TRUE" values under `BC_Contam`
2. If applicable, flag or filter out any duplicate barcodes
3. Filter out any non expected variants. The pipeline will automatically detect any strings of N’s in the input.fasta of the `DeNovo_ref` and identify those as an `expected bcs`. The number of expected barcodes in the aligning sequence will be listed in `expected_bcs`. Assuming that the expected BC(s) was detected, perfect sequences will have the same value in `n_vars` as `expected_bcs` and `n_barcodes`.
4. You are now ready to select wells to prep. The `plate_well` column will identify which plates/wells from the glycerol plate will correspond to the desired sequences. If possible, select wells with the best coverage as represented by a larger `Leftover` value and with `LT_10` of 0 and `LT_3` of 0. `LT_10` is less critical but `LT_3` should never be larger than 0 (for a 10kb plasmid an `LT_3` of 0.001 means that 10 bp of the plasmid did not have a coverage of at least three).
5. If there happens to be a clone that does not have sufficient coverage but is absolutely required, use the pileup from `samtools` to manually inspect aligned reads to the `DeNovo_ref` if the missing coverage is in a critical part of the plasmid.
      - In a new terminal, `cd` into your octopus directory
      - Open up a new docker instance - `docker run --rm -it -v "$(pwd)":/root/octopus octant/octopus /bin/bash`
      - Navigate into the folder that contains the analyzed data from that run - `cd octopus/pipeline/your_run_id`
      - Index the well of interest - `samtools index your_plate/your_well.map.bam`
      - View the pileup - `samtools tview your_plate/your_well.map.bam lib/your_ref.fasta`
	  - Note `your_ref` can be determined by `samtools view your_plate/your_well.map.bam | head` and looking for the reference field

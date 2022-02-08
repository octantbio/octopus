**UPDATE:** There is now an update to the protocol we want to share that uses RCA for OCTOPUS sample plate generation as an alternative to the lysate generation steps below. Our motivation for this was to minimize genomic DNA contaminant reads in the final NGS pool which increases plasmid/genomic DNA ratio as well as to offload some of the labor required in Day 2 to Day 1. Additionally, this adds further flexibility to scaling OCTOPUS runs, as lysate generation was the least scalable step in the older protocol. The updated protocol can be [read here](Bench-Protocol-with-RCA.md).

# Table of Contents

- [Reagents](#reagents)
- [Consumables](#consumables)
- [Equipment](#equipment)
- [Day 1 - Picking colonies](#day-1---picking-colonies-10-15-min-per-96-colonies)
- [Day 2 - Lysate and glycerol stock generation](#day-2---lysate-and-glycerol-stock-generation-1-2-hours)
- [Day 2 - "Miniaturized" Riptide protocol](#day-2---miniaturized-riptide-protocol-4-6-hours)
    - [Random primer extension and biotinylated termination](#a-reaction-random-primer-extension-and-biotinylated-termination-1-hour)
    - [DNA capture and library conversion](#b-reaction-dna-capture-and-library-conversion-1-2-hours)
    - [Library Amplification](#pcr-amplification-45-mins)
    - [Size Selection](#spri-bead-cleanup-and-gel-isolation-1-hour)
- [Day 2 - Library quantification and sequencing](#day-2---library-quantification-and-sequencing-setup-45-mins)

# Reagents:

- ddH2O
- 30% [Glycerol](https://www.fishersci.com/shop/products/glycerol-molecular-biology-fisher-bioreagents-2/BP2291)
- 0.1 M [Sodium Hydroxide](https://www.fishersci.com/shop/products/sodium-hydroxide-pellets-certified-acs-fisher-chemical-7/S318100)
- [2xYT](https://us.vwr.com/store/product/7437420/vwr-life-science-2xyt-medium-broth)
- 75 mM EDTA
- Ethanol
- Agarose
- TAE

# Consumables:

- [Riptide Kit](https://igenomx.com/product/riptide/) from iGenomX
- [Flat Bottom 96-well plates](https://www.sigmaaldrich.com/catalog/product/aldrich/br781602?lang=en&region=US)
- [96-well PCR plates](https://www.thermofisher.com/order/catalog/product/AB0600)
- [384-well PCR plates](https://www.thermofisher.com/order/catalog/product/4483317?SID=srch-hj-4483317) (optional)
- Foil Plate seal (Axygen PCR-AS-200)
- Plastic Plate seal (Axygen PCR-SP)
- [Gel Extraction Kit](https://www.zymoresearch.com/collections/zymoclean-gel-dna-recovery-kits/products/zymoclean-gel-dna-recovery-kit)
- [MiSeq sequencing kits](https://www.illumina.com/products/by-type/sequencing-kits/cluster-gen-sequencing-reagents/miseq-reagent-kit-v2.html)

# Equipment:

- [96-well thermal cycler](https://www.bio-rad.com/en-us/product/c1000-touch-thermal-cycler?ID=LGTW9415)
- [384-well thermal cycler](https://www.thermofisher.com/order/catalog/product/4388444) (optional)
- Multichannel 10
- Multichannel 200
- [Liquidator 20](https://www.shoprainin.com/Products/Pipettes-and-Tips/Pipettes/High-throughput-Pipetting/Liquidator%E2%84%A2-96/Liquidator-96-2-20-%C2%B5L-LIQ-96-20/p/17014207) (optional)
- [Liquidator 200](https://www.shoprainin.com/Products/Pipettes-and-Tips/Pipettes/High-throughput-Pipetting/Liquidator%E2%84%A2-96/Liquidator96%2C-5-200-%C2%B5L-LIQ-96-200/p/17010335) (optional)
- 37˚C shaking Incubator
- Nanodrop, Qubit, BioAnalyzer or Tapestation
- [Illumina Miseq](https://www.illumina.com/systems/sequencing-platforms/miseq.html)
- [Dynamag](https://www.thermofisher.com/order/catalog/product/12321D) (or magnetic tube stand)
- [Multi-tube vortexer](https://ru.vwr.com/store/product/596528/multi-tube-vortexers) (optional)

# Day 1 - Picking colonies (~10-15 min per 96 colonies)

1. Pick single colonies into 150 µl of 2xYT with appropriate antibiotic in a flat-bottom 96-well plate
    - Different antibiotics have different plasmid yield in the lysate
    - When cloning in a pooled library format, we usually pick a minimum of 4 colonies per desired clone.
    - *Tip: We pick colonies using 200 µL pipette tips, put them back into the box, and then use a multichannel pipette to inoculate the 2xYT plate. This saves a lot of time and should allow you to pick a 96-well plate in about ~10 minutes with some practice (potential upgrade is a colony picker - we have yet to identify a cost effective and robust model)*
2. Tape the lid tightly to the plate to reduce evaporation and grow overnight at 37˚C in a shaking incubator
    - a 96-well plate specific incubator with a short throw and high rpm is great, but we just use a normal incubator/shaker meant for flasks at 250 rpm

# Day 2 - Lysate and glycerol stock generation (~1-2 hours)

1. Transfer 50 µL of resuspended culture from each well using a liquidator/multichannel into a 96-well PCR plate
2. Pellet the bacteria by spinning at 4000 rcf for 10 min at RT in a swinging bucket rotor
3. Store the remaining 100 µL as glycerol stock by adding 100 µL of 30% glycerol to the plates, shake in incubator for 15 min, foil seal, and store in -80˚C freezer.
4. Flick and blot the media from the pelleted bacteria (media in this step will inhibit proper lysis)
5. Add 40 µL of MillQ water to each pellet, plastic seal, and resuspend pellet by pipetting or with a multi-tube vortexer
    - *tip: pulse by going from medium to max speed, then leave at near max speed for 5-10 secs*
6. Lyse cells in 96-well thermocycler by heating at 95˚C for 3 min and then cooling to 4˚C
7. Clarify lysate by spinning at 4000+ rcf for 10 min at RT
8. Remove the top 20 µLof clarified lysate and store in a separate 96 or 384-well plate
9. *At this point the protocol can be paused by freezing the lysate at -80˚C*

# Day 2 - "Miniaturized" Riptide protocol (~4-6 hours)

This protocol is adapted from the iGenomX Riptide library prep protocol [here](https://igenomx.com/resources/)

## A Reaction: Random primer extension and biotinylated termination (~1 hour)

Consumables needed for "A Reaction":

```
dNTP Mix 1
10X Enzyme 1 Buffer
Primer A
Enzyme 1
SPRI Beads 1
75 mM EDTA
80% ethanol (freshly prepared)
10 mM Tris-HCl pH 8
```

1. For each 96-well plate, prepare 200 µL of Reaction A master-mix by mixing
    - 100 µL **dNTP Mix 1**
    - 50 µL **10X Enzyme 1 Buffer**
    - 50 µL **Enzyme 1**
2. Fill a 96-well plate with 2 µL each of the master-mix using a multichannel or liquidator
    - Consider 384-well plates if working with 4 96-well plates
3. Pipette into each well 1 µL of **Primer A**, using liquidator or multichannel
    - Make sure to properly mix primer A plate if thawing
4. Pipette into each well 2 µL of clarified lysate using liquidator or multichannel
    - Make sure to properly mix lysate if thawing
    - If sequencing mini-preps, use 2 µL of 5 ng/µL per well
5. Seal the plate with foil and run the following protocol on a thermocycler
    - 92˚C for 3 min
    - 16˚C for 5 min
    - Slow ramp (0.1˚C/sec to 68˚C)
    - 68˚C for 15 min
    - Hold at 4˚C
6. **At this point the protocol can be paused by freezing the plate(s) at -20˚C**
7. Warm the **SPRI Beads I** to RT.
8. Stop the reaction with EDTA and pool the contents of each 96-well plate
    - Option A:
    1. Use a liquidator to add 1 µL of **75 mM EDTA** to each well
    2. Use those same tips to aspirate and dispense the contents of each 96-well plate onto a clean liquidator tip-box lid and pool everything by tapping the lid at an angle. Make sure all drops are collected.
    3. Dispense pooled liquid into a low-retention Eppendorf tube (should recover around 400-500 µL of sample).
    - Option B:
    1. Pool samples into a PCR strip tube already containing **75 mM EDTA**, using a multichannel
9. Ensure each 96-well plate has it's own pooled tube
10. Add 1.8 volumes of thoroughly-mixed room temperature **SPRI Beads I** to each pooled tube (use a 1 mL pipette to measure the volume). Mix by pipetting, incubate 10 min at room temperature.
11. Place tube(s) on dynamag, allow the solution to clear (2-5 min), and discard supernatant.
12. Keeping tube(s) on the dynamag, add 1300 µL of freshly prepared **80% ethanol** to each tube, wait 30 secs, then aspirate ethanol.
      - Repeat this step, and ensure that all ethanol is removed.
13. Open the caps and allow the beads to air dry on the dynamag for 10 min, don't overdry.
14. Add 50 µL of RT **10 mM Tris-HCl pH 8** to the beads
15. Remove tube(s) from dynamag, and resuspend the beads until homogenous. Incubate at RT for 10 min to elute.
16. Allow the beads to clear on dynamag for 2 min, then transfer clarified elution to new low retention PCR tube(s)
17. **At this point the protocol can be paused by freezing tubes in the -20˚C**

## B Reaction: DNA Capture and Library Conversion (~1-2 hours)

Consumables needed for "B Reaction":

```
HS-Buffer
Capture Beads
Bead Wash Buffer
0.1M NaOH
Enzyme II
Enzyme II Buffer
dNTP Mix II
Primer B
Nuclease Free Water
```

1. Heat-denature the elution at 95˚C for 3 min and hold at 4˚C in a thermocycler
2. While it's heating, prep the **Capture Beads**
    1. Warm the **Capture Beads** and **HS-Buffer** to RT, resuspend thoroughly, and transfer 20 µL of slurry for each sample into a new PCR tube.
    2. Dynamag the **Capture Beads** and discard the supernatant.
    3. Remove tubes from dynamag and add 100 µL of **HS-Buffer**
    4. Dynamag and remove the wash.
    5. Resuspend the beads in 20 µL of **HS-Buffer**.
3. Add all 50 µL of the heat-denatured elution to the washed **Capture Beads**, mix, and incubate at RT for 10 min.
4. Pipette mix beads again and incubate at RT for 10 min.
5. Dynamag the beads and discard the supernatant.
6. Resuspend beads with 50 µL of **0.1M NaOH**, incubate 4 min at RT, dynamag and remove supernatant
7. Wash the beads 3 times (resuspend beads in 100 µL of RT **Bead Wash Buffer**, dynamag, remove supernatant). Make sure to remove any remaining liquid after final wash.
8. Prepare a mastermix for **Reaction B**: for every tube of beads, mix together on ice
    - 4 µL 5x **Enzyme II Buffer**
    - 1.5 µL **dNTP Mix II**
    - 2 µL **Primer B**
    - 12 µL **Nuclease-Free Water**
    - 0.5 µL **Enzyme II**
9. Quickly add 19.5 µL of **Reaction B** mastermix to each tube (try not to let mastermix sit around at RT)
10. Incubate the tubes in a thermocycler for 20 min at 24˚C and hold at 4˚C for at least 3 min
11. Pipette mix the beads, dynamag and discard supernatant
12. Wash the beads 3 times (resuspend beads in 100 µL of RT **Bead Wash Buffer**, dynamag, remove supernatant). Make sure to remove any remaining liquid after final wash.

## PCR Amplification (45 mins)

Consumables needed for PCR Amplification:

```
Universal PCR primer
Index PCR primer(s)
2X PCR amplification mix from iGenomX
ddH2O
```

1. Resuspend the beads in 21 µL of nuclease free water, then setup the PCR reaction by adding:
    - 2 µL **Universal PCR primer**
    - 2 µL **Index PCR primer** (1-12) (choose one barcoded primer per pool)
    - 25 µL **2X PCR Amplification Mix**
    - 50 µL total
2. Input the following program into a thermocycler

```
1 cycle:   98˚C, 2 min
11 cycles: 98˚C, 20 sec
           60˚C, 30 sec
           72˚C, 30 sec
1 cycle:   72˚C, 5 min
           4˚C, hold
```

3. Record which samples received which Index PCR primer
4. Sample can be left in the thermocycler at 4˚C overnight.
5. Briefly spin the PCR tube in a picofuge, dynamag and transfer the supernatant to new low retention eppendorf 1.5 mL tubes. Discard the PCR tubes containing the Capture Beads.

## SPRI bead cleanup and gel-isolation (~1 hour)

Consumables needed for SPRI bead cleanup and gel-isolation:

```
SPRI Beads II
1% Agarose
TAE
80% Ethanol
10 mM Tris-HCl pH 8
Gel-extraction kit
SYBR Safe
```

1. Add 70 µL of well resuspended **RT SPRI Beads II** to the samples. Mix well and incubate at RT for 10 min.
2. During this incubation, pour a 1.0% agarose gel with **SYBR Safe** for gel-extraction
3. Dynamag the beads for at least 2 min, and discard supernatant.
4. Add 200 µL of **80% ethanol** to each tube(s), wait 30 sec, then remove and discard the ethanol
     - it is unnecessary to remove from magnet for the ethanol wash
5. Repeat with another 200 µL of **80% ethanol**, carefully remove all ethanol from tube without disturbing beads.
6. Open cap and allow to air dry for 10 min on dynamag (careful not to overdry).
7. Add 25 µL of RT **10 mM Tris-HCl pH 8** to beads. Remove from dynamag and fully resuspend the beads.
8. Incubate at RT for 10min to elute, place back on dynamag and transfer supernatant to new low-retention tubes.
9. Run sample(s) on a **1.0% agarose** gel and gel extract. It should come out as a visible smear, isolate the 400-1200bp region, taking care to avoid the potential primer dimer band.
10. **At this point the protocol can be paused by putting the tube(s) at -20˚C**

# Day 2 - Library Quantification and sequencing setup (45 mins)

Consumables needed for library quantification and sequencing:
- Appropriate MiSeq Kit (see below)

Every quantification method is slightly different and each comes with its own pros and cons. These issues are exacerbated by the Riptide library product which is a heterogenous mixture of DNA fragment sizes (*see above*). The way we reduce variability in quantification is to quantify a previously run OCTOPUS library and use that as a baseline to estimate the concentration of the current libraries. In our hands, a fluorescence based assay gives sufficiently accurate quantification alongside a previously run library. If no previously run libraries are available, quantification by qPCR has been the most accurate method.
1. Select the appropriate sequencing kit: as a rule of thumb, 10,000 paired-end 150 reads/well, or 1,000,000 reads per 96-well plate (nano-kit V2 for one 96-well plate, micro-kit V2 for a 384-well plate, and a standard V2 kit for anything more than 384-well, the most we do is 3 x 384-wells due to the limited number of index primers in the riptide kit)
2. The only information necessary to generate the sample sheet is to provide the plate index sequences (see below)

**Example sample-sheet.csv:**

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

The sequencing run will take ~16-30 hours depending on the size of the kit (nano, micro, standard)

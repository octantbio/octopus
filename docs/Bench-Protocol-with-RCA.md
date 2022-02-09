Below reflects the updated protocol that replaces input DNA preparation by whole cell lysate generation on Day 2 with RCA preparation on Day 1. Advantages here are that it moves a laborious step from the library preperation on Day 2 (the longest day) to day one during colony picking. It also reduces genomic DNA contamination and increases the yield and robustness of the library prep greatly.

# Table of Contents

- [Reagents](#reagents)
- [Consumables](#consumables)
- [Equipment](#equipment)
- [Day 1 - Colony picking and RCA](#day-1---colony-picking-and-rca-1-2-hours)
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

- [Phi29 polymerase and buffer](https://lucigen.com/docs/manuals/MA114-phi29.pdf) from Lucigen
- [10 mM dNTP set](https://lucigen.com/docs/manuals/MA077-10-mM-dNTP.pdf)
- 100uM custom RCA primers: Your own custom "RCA Primers" mixed equimolarly each at 100uM, designed to bind to a common sequence on your plasmids
    - For example, we designed 8 primers that target the E1 origin of replication that is in all our plasmids.
    - We ordered these from [IDT](https://www.idtdna.com) with the last three 5' bases phosphorothiorated to make them exo resistant.
- [500 µM Exo Resistant Random Hexamers](https://www.mclab.com/Exo-Resistant-Random-Hexamer.html) from MCLAB
    - We found that having a small amount of random hexamers normalizes against amplification bias of regions of plasmid targeted by the custom RCA primers and increases yield.
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

# Day 1 - Colony Picking and RCA (~1-2 hours)

1. Aliquot 25 µl ddH2O into each well in a 96-well plate. This is your "Sample Plate".
2. For each Sample Plate, prepare a corresponding plate for glycerol stocks by filling each well of a flat-bottom 96-well plate with 110 µl of 2xYT +antibiotic media. This is your "Culture Plate".
3. Pick single colonies into the Culture Plate and stamp into the Sample Plate:
    - Pick colonies using 20 µl tips. Pipette up and down in the Culture Plate to mix, then using the same tips stamp 5 µl from the Culture Plate into the Sample Plate.
    - If sequencing minipreps is desired, deposit 1 ul of miniprep to the well.
    - Grow the Culture Plate with lid on overnight in a 37˚C shaker for at least 8 hours, but no more than 16 hours. Evaporation of edge wells may occur here; we found that taping the lid tightly to the plate reduces evaporation.
4. Seal the Sample Plate with a plastic film and heat at 95˚C for 3 minutes before bringing it down to 4˚C using a thermocycler. 
5. During heating, prepare 500 µl of the RCA enzyme master mix for one 96-well plate and keep on ice:
    - 100 µl Phi29 DNA Polymerase Buffer
    - 20 µl Phi29 Polymerase
    - 40 µl 10 mM dNTP mix
    - 4 µl 100 µM custom RCA Primers
    - 2.5 µl 500 µM Exo Resistant Random Hexamers
    - 333.5 µl ddH2O
6. Aliquot 5 µl of this RCA enzyme master mix into each well of a new 96-well plate.
7. After the heating is done for the Sample Plate. Stamp 5 µl from the Sample Plate into the RCA enzyme master mix plate, mix by pipetting.
8. Seal the plate with foil seal, spin down, incubate at 30˚C for at least 12 hours, then heat inactivate at 65˚C for 10 minutes, then hold at 4˚C in thermocycler. This RCA plate can now be used directly as the DNA input to Reaction A (See below) without further prep work.

# Day 2 - "Miniaturized" Riptide protocol (~4-6 hours)

## Post Day 1 Prep: Glycerol stocks (~1 hour)

1. To your culture plates, add 100 d ul of 30% Glycerol to each well, mix, foil, and place in -80˚C.

The following protocol is adapted from the iGenomX Riptide library prep protocol [here](https://igenomx.com/resources/)

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
4. Pipette into each well 2 µL of RCA reaction using liquidator or multichannel
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
16. Allow the beads to clear on dynamag for 2 min, then transfer elution to new low retention PCR tube(s)
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

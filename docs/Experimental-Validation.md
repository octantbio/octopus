## OCTOPUS is Reproducible and accurate

To better understand the reproducibility and accuracy of OCTOPUS, we transformed a previously validated 8.6 kilobase plasmid into E. coli and used OCTOPUS to sequence one colony in 192 different wells. With ~1,000,000 paired-end 150 reads per plate (or ~10,000 paired-end reads per well), the coverage throughout the plasmid was relatively uniform (gray ribbon is the inter-quartile range of coverage), with an average coefficient of variation of 0.485 across the wells. In this run, OCTOPUS correctly verified 157/192 (81.8%; open) wells, with at least 3x coverage across the entire plasmid. The remaining 35/192 (18.2%; blue) wells had less than 3x coverage in some parts of the plasmid, and 5 of those reported a variant (red). We thought this was unlikely given the low somatic mutation rate of E. coli, and by manually inspecting the read pileups for those wells, we found that these likely incorrectly reported variants were in regions of 0 coverage. For this reason we don’t recommend accepting plasmids with less than 3x coverage without manual inspection. 

![](https://github.com/octantbio/octopus/blob/master/img/coverage-rank.png)

## OCTOPUS Works on Diverse Sequences

We analyzed a typical 384-well OCTOPUS run (~10,000 paired-end 150 reads per well) of a pooled GPCR cloning reaction for coverage. We excluded 45/384 wells that lacked sufficient coverage to properly identify the construct and 8/384 wells that lacked an insert. Of the remaining wells, the majority (173/331; 52.3%) have at least 10x coverage and almost all (300/331; 90.6%) have at least 3x coverage across the entire plasmid. These data suggest that OCTOPUS is able to verify a diverse set of sequences with sufficient coverage.

![](https://github.com/octantbio/octopus/blob/master/img/gpcr-frac-coverage.png)

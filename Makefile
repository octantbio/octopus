SHELL := /bin/bash -euo pipefail
.DELETE_ON_ERROR:

RUNS := $(notdir $(shell find data/ -mindepth 1 -maxdepth 1 -type d))

.PHONY: denovo stats all clean
denovo: $(addprefix pipeline/, $(addsuffix /spades-contigs.fasta, $(RUNS)))
stats: $(addprefix pipeline/, $(addsuffix /aggregated-stats.tsv, $(RUNS)))
all: stats

# cleanup
clean:
	rm -rf pipeline/*

.SECONDARY:

#===============================================================================
#                        DeNovo iGenomX Pipeline
#                        ------------------------
# 0) Ensure data is linked
# 1) Demultiplex the reads using IgenomX's protocol
# 2) Preprocess the reads (trim, filter, error-correct, ...) according to JGI
# 3) De novo assemble the results
# 4) Align contigs to library to get well identity
#===============================================================================

# ---------------
# 0) Data Linking:
# ---------------
# link-data.py will traverse a standard illumina directory and clean read names
# | specifies an order only prerequisite - make will only check existence
pipeline/%/input.fasta pipeline/%/fastqs: | data/%
	@echo "Cleaning up $<"
	@python src/link-data.py $| --out-dir $(@D)

# get list of plates to process
# exclude plates that don't demultiplex
# -print0 ensures nasty filenames are handled with grace
# anything that parses plates.txt must handle null characters
pipeline/%/plates.txt: pipeline/%/fastqs
	@find -path "./$</*.fastq*" ! -regex '.*Undetermined.*' -print0 \
	    | sed --null-data -e 's/_R.*//' \
	    | sort --zero-terminated \
	    | uniq --zero-terminated \
	    | sed --null-data -e 's/.*\///'\
	    > $@

# ---------------
# 1) Demultiplex:
# ---------------
# --xapply ensures {1} {2} are processed as we expect (i.e. not a permutation)
# '{=1 ... }' enables you to run a perl expression on the first position var
# s:.*/:: strips the paths (think s/.*//); s:_.*:: strips the underscores
pipeline/%/demux: pipeline/%/fastqs
	@echo "Demultiplexing $(@D)"
	@parallel --xapply fgbio \
	    --log-level=Error \
	    DemuxFastqs \
	    --inputs {1} {2} \
	    --metadata src/igenomx-meta.txt \
	    --read-structures 8B12M134T 8M142T \
	    --output-type=Fastq \
	    --output=$(@D)/'{=1 s:.*/::;s:_.*:: =}' \
	    --max-mismatches=2 \
	    --metrics $(@D)/'{=1 s:.*/::;s:_.*:: =}'/igenomx_test_metrics.txt \
	    ::: $</*_R1.* ::: $</*_R2.* \
	    && touch $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ---------------
# 2) Preprocess:
# ---------------
# 0) see src/jgi-preproc.sh for implementation
# 1) trim adapters, filter contams (ecoli, phiX, etc)
# 2) error correct reads
# 3) get stats on well contaminants


# Grab DH5a Genome:
# -----------------
src/background.fasta:
	@echo "Grabbing e. coli genome"
	@curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&amp;id=CP017100&amp;rettype=fasta&amp;retmode=text" > $@


# Preprocessing:
# --------------
# if you want a different genome, change out here
# fgbio will demultiplex the wells into folder corresponding to whatever's in ./pipeline/fastqs
# ensure everything ends in null character to handle filenames with nasty characters
# xargs -I {} converts {} into substitution (much like parallel)
# sed can change out delimiter s|foo|bar| == s/foo/bar/
# -n2 puts two lines of the input at {}
pipeline/%/preproc: pipeline/%/plates.txt pipeline/%/demux src/background.fasta
	@echo "Preprocessing plates in $(<D)"
	@sed --null-data -e 's|^|./$(<D)/|' -e 's|$$|/*.fastq*|' $< \
	    | xargs --null -n1 -I {} find -path {} ! -regex '.*unmatched.*' -print0 \
	    | sort --zero-terminated \
	    | parallel --null -n2 src/jgi-preproc.sh {} $(lastword $^) \
	    && touch $@

# Contamination Stats:
# --------------------
# jgi's pipeline provides per-well contam stats
# recall parallel requires escaping '
pipeline/%/read-stats.tsv: pipeline/%/plates.txt pipeline/%/preproc
	@echo "Calculating well statistics for $(<D)"
	@sed --null-data -e 's|^|./$(<D)/|' $< \
	    | xargs --null -n1 -I {} find -path {} -print0 \
	    | parallel --null \
	    grep -F -B4 \'Unique 31\' {}/*.pre-proc \
	    \| sed -e \''/Total /d; /Input/d; /Unique/d; /--/d'\' \
	    -e \''s/reads/\t/g; s/bases/\t/g; s/\.pre-proc-/\t/g'\' \
	    -e \''s/\.\/pipeline\///g; s/[:()]//g'\' \
	    -e \''s/\//\t/g; s/  \+//g'\' \
	| awk '{NF = NF -2; print}' OFS="\t" \
	| mlr --tsvlite --implicit-csv-header label Run,Plate,Well,Metric,Reads,Percent \
	> $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --------------------
# 3) De novo Assembly:
# --------------------
# 0) see src/jgi-denovo.sh for implementation
# 1) merge reads
# 2) use SPAdes to assemble

pipeline/%/de-novo: pipeline/%/plates.txt pipeline/%/preproc
	@echo "De novo assembling all plates in $(<D)"
	@sed --null-data -e 's|^|./$(<D)/|' -e 's|$$|/*.ecc.fq.gz|' $< \
	    | xargs --null -n1 -I {} find -path {} -print0 \
	    | parallel --null src/jgi-denovo.sh {} \
	    && touch $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# -----------------
# 4) Align Contigs:
# -----------------
# 1) Flatten input plasmids to allow for alignments to span the junction
# 2) aggregate all assemblies into one fasta (and reorient)
# 3) align to input library
# 4) parse alignment for reference info

# Flatten input library:
# -----------------------
# Cat input references together to enable mapping through plasmid junction
pipeline/%/input-flatten.fasta: pipeline/%/input.fasta
	@echo "Flattening $<"
	@python src/flatten-fasta.py $< > $@

# Collapse Contigs:
# -----------------
# again harness parallel's string processing to output as contigs/PLATE_WELL.fasta
# {/} = basename {//} = dirname
# only take the first contig
# 1) dump all fastas into spades-contigs/Plate_Well.spades-contigs.fasta
# 2) take only the first record
pipeline/%/spades-contigs.fasta: pipeline/%/plates.txt pipeline/%/de-novo
	@echo "Aggregating all contigs in $(<D)"
	@mkdir -p $(@:.fasta=) \
	    && sed --null-data -e 's|^|./$(<D)/|' -e 's|$$|/*.spades-contig.fasta|' $< \
	    | xargs --null -n1 -I {} find -path {} -print0 \
	    | parallel --null \
	    'well="$$(basename -s .spades-contig.fasta {/})"; \
	    plate="$$(basename {//})"; \
	    cp {} $(@:.fasta=)/"$$plate"_"$$well".fasta' \
	    && python src/flatten-fasta.py --no-flat --records 1 $(@:.fasta=)/*.fasta \
	    > $@

# Align Contigs:
# --------------
# --eqx: =/X in cigar string
# -N: no secondary alignments
# -a output sam
# take only good alignments (0 for forward, 16 for rev-comp)
pipeline/%/spades-contigs.sam: pipeline/%/input-flatten.fasta pipeline/%/spades-contigs.fasta
	@echo "Mapping contigs from $(<D)"
	@minimap2 --eqx --secondary=no -x asm20 -a  -t $$(nproc --all) $^ \
	    2> $(@:.sam=.map.err) \
	    | awk '$$2 == 0 || $$2 == 16' \
	    > $@

# Parse Alignment:
# ----------------
pipeline/%/de-novo-ref-stats.tsv: pipeline/%/spades-contigs.sam
	@awk -v plate=$(<D) \
	    '$$1 !~ /@/ \
	    {len=split(plate,p,"/"); split($$1,a,"_"); print p[len],a[1],a[2],$$3,$$6,$$10}' \
	    OFS="\t" $< \
	    | mlr --tsvlite --implicit-csv-header label Run,Plate,Well,DeNovo_Ref,CIGAR,Contig \
	    > $@

#===============================================================================
#                             Variant Calling
#                             ---------------
# 0) Determine well identity with the de novo assembly contig
# 1) Prep input reference for alignment
# 2) Variant calling
# 3) Variant parsing
#===============================================================================

# -----------------------
# 1) Prep Input Reference
# -----------------------
# flatten-fasta --no-flat cleans input
pipeline/%/input-refs.fasta: pipeline/%/input.fasta
	@echo "Sanitizing input fastas"
	@python src/flatten-fasta.py --no-flat $< > $@

# split input into multiple files for downstream alignment and variant calls
pipeline/%/lib/split: pipeline/%/input-refs.fasta
	@echo "Splitting $< into individual fastas"
	@mkdir -p $(@D) \
	    && paste - - < $< \
	    | awk '{print $$1"\n"$$2 > "$(@D)/"substr($$1,2)".fasta"}' \
	    && parallel samtools faidx {} ::: $(@D)/*.fasta \
	    && touch $@

# split input into tsv for easy joining
pipeline/%/input-refs.tsv: pipeline/%/input-refs.fasta
	@echo "Splitting input library into tsv"
	@sed 's/>//g' $< \
	    | paste - - \
	    | mlr --tsvlite --implicit-csv-header label Ref,Ref_Seq \
	    > $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ------------------
# 2) Variant Calling
# ------------------
# 0) See /src/variant-calling.sh for details
# 1) get well identity
# 2) Re-map preprocessed reads to that specific reference
# 3) Call variants with Q>20, and >=1 read that makes up 50% of reads at that position

pipeline/%/guided: pipeline/%/de-novo-ref-stats.tsv pipeline/%/preproc pipeline/%/lib/split
	@echo "Calling variants for all plates in $(<D)"
	@awk 'NR > 1{print "pipeline",$$1,"lib",$$4".fasta pipeline",$$1,$$2,$$3".ecc.fq.gz"}' \
	    OFS=/ \
	    $< \
	    | parallel --col-sep ' ' src/variant-calling.sh {1} {2} \
	    && touch $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ------------------
# 3) Variant Parsing
# ------------------
# 0) See /src/vcf-parse.py for details
# 1) parse input fasta for barcodes
# 2) compare barcodes found in VCF file to expected number from input fasta
# 3) un-tidy for nice excel output

# parallel requires escaping |
pipeline/%/freebayes-tidy.tsv: pipeline/%/plates.txt pipeline/%/input-refs.fasta pipeline/%/guided
	@echo "Parsing variants from FreeBayes for $(<D)"
	@sed --null-data -e 's|^|./$(<D)/|' -e 's|$$|/*.freebayes.bcf|' $< \
	    | xargs --null -n1 -I {} find -path {} -print0 \
	    | parallel --null \
	    bcftools view {} \
	    \| python src/vcf-parse.py $(word 2, $^) - \
	    \| awk -v wells='{= s:\.freebayes\.bcf:: =}' \
	    \''{len=split(wells,a,"/"); print a[len-2],a[len-1],a[len],$$0}'\' OFS=\'\\t\' \
	    > $@

pipeline/%/freebayes.tsv: pipeline/%/freebayes-tidy.tsv
	@mlr --tsvlite --implicit-csv-header \
	    label Run,Plate,Well,Ref,idx,bc,pos,n_vars,n_barcodes,expected_bcs \
	    then reshape -i bc,pos -o foo,bar \
	    then put '$$baz = $$foo . $$idx' \
	    then cut -x -f foo,idx \
	    then reshape -s baz,bar \
	    then unsparsify --fill-with "NA" \
	    $< \
	    > $@

#===============================================================================
#                         Ancillary Analysis
#                         ------------------
# 1) Coverage
# 2) Barcode Filter
# 3) Aggregation
#===============================================================================

# -----------------
# 1) Coverage Info:
# -----------------
# 1) Calculate coverage from variant calling
# 2) Report percentages of bases with <10x, <3x coverage

# -aa -> absolutely all positions; -d0 no max depth
pipeline/%/lt-X.tsv: pipeline/%/plates.txt pipeline/%/guided
	@echo "Calculating percent bases <10x and <3x coverage for $(@D)"
	@sed --null-data -e 's|^|./$(<D)/|' -e 's|$$|/*.map.bam|' $< \
	    | xargs --null -n1 -I {} find -path {} -print0 \
	    | parallel --null \
	    samtools depth -aa -d0 {} \
	    \| awk -v path='{= s:\.map\.bam:: =}' \
	    \''{if($$3 < 3) lt3 += 1; else if($$3 < 10) lt10 += 1} END {len=split(path,a,"/"); print a[len-2],a[len-1],a[len],lt10/NR,lt3/NR}'\' \
	    OFS=\"\\t\" \
	| mlr --tsvlite --implicit-csv-header label Run,Plate,Well,LT_10,LT_3 \
	> $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# -----------------
# 2) Barcode Filter
# -----------------
# 0) see src/bc-contam.py for details
# 1) parse bcf file for barcode locations
# 2) output pileup at barcode
# 3) collapse barcodes at lev dist 1
# 4) count their frequency
# 5) <--- good ---[ ratio check ]--- contam --->
#    (0      0.04)(   2/3 > 2  )(0.1         oo)
# *) catch edge case with two barcodes but  0.04 < bc2 < 0.1 -> CONTAM
# the -e trick in mlr allows us to break up the ifs while keeping make happy
# finally clean up the pipeline/run/plate/well.map.bam into columns
pipeline/%/barcode-filter.tsv: pipeline/%/freebayes-tidy.tsv
	@echo "Checking for contaminated barcodes"
	@awk '{print "pipeline",$$1,$$2,$$3".map.bam "$$4" "$$6" "$$7}' \
	    OFS='/' $< \
	    | python src/bc-contam.py - \
	    | mlr --tsvlite --implicit-csv-header \
	    label Well,Barcode,Reads \
	    then fraction -g Well -f Reads \
	    then top -n 3 -g Well -f Reads_fraction \
	    then reshape -s top_idx,Reads_fraction_top \
	    then put \
	    -e 'if ($$2 >= 0.1){$$BC_Contam = "TRUE"}' \
	    -e 'elif ($$2 < 0.1 && $$2 > 0.04 && $$3 == ""){$$BC_Contam = "TRUE"}' \
	    -e 'elif ($$2 < 0.1 && $$2 > 0.04 && $$2 / $$3 > 2){$$BC_Contam = "TRUE"}' \
	    -e 'else{$$BC_Contam = "FALSE"}' \
	    then cut -f Well,BC_Contam \
	    then nest --explode --values --across-fields -f Well --nested-fs / \
	    then label foo,Run,Plate,Well \
	    then cut -x -f foo \
	    | sed 's/\.map\.bam//g' \
	    > $@

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ----------------------------
# 3) Aggregate All the Things:
# ----------------------------
# unsparsify may be unnecessary depending on how clean your data is
# replicates right_join
# http://johnkerl.org/miller-releases/miller-head/doc/faq.html#How_to_rectangularize_after_joins_with_unpaired?
pipeline/%/aggregated-stats.tsv: \
    pipeline/%/read-stats.tsv \
    pipeline/%/freebayes.tsv \
    pipeline/%/de-novo-ref-stats.tsv \
    pipeline/%/lt-X.tsv \
    pipeline/%/barcode-filter.tsv \
    pipeline/%/input-refs.tsv
	@echo "Aggregating everying together into $@"
	@mlr --tsvlite \
	     cut -x -f Percent \
	     then reshape -s Metric,Reads \
	     then rename Result,Leftover \
	     then join -j Run,Plate,Well --ur -f $(word 2, $^) \
	     then unsparsify --fill-with "NA" \
	     then join -j Run,Plate,Well --ur -f $(word 3, $^) \
	     then unsparsify --fill-with "NA" \
	     then join -j Run,Plate,Well --ur -f $(word 4, $^) \
	     then unsparsify --fill-with "NA" \
	     then join -j Run,Plate,Well --ur -f $(word 5, $^) \
	     then unsparsify --fill-with "NA" \
	     then join -j Ref --ur -f $(word 6, $^) \
	     then unsparsify --fill-with "NA" \
	     then cut -x -f 1,2 \
	     then put 'plate=splitnv($$Plate, "-"); well=splitnv($$Well, "-"); $$Plate_Well = plate[2]."_".well[1]' \
	     then reorder -f Run,Plate,Well,Plate_Well,Ref,DeNovo_Ref,CIGAR,LT_10,LT_3,BC_Contam,BC_Clash,n_vars \
	     then reorder -e -f Contig,Ref_Seq \
	     < $< \
	     > $@

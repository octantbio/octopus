#!/bin/bash
# bash strict mode
set -euo pipefail
IFS=$'\n\t'

#-------------------------------------------------------------------------------------------------

REF=$1
REF_DIR="$(dirname ${REF})"


READ_1=$2
WELL="$(basename -s ".ecc.fq.gz" ${READ_1})"
OUT_DIR="$(dirname ${READ_1})"

THREADS="threads=1"
MEM="-Xmx8g"

# send all temp files to the scratch directory
# trap ensures that scratch will be deleted no matter what
scratch="$(mktemp -d -t tmp.XXXXXXXXXX)"
function finish {
    mv -f ${scratch}/guided.err ${OUT_DIR}/${WELL}.guided.err
    rm -rf ${scratch}
}
trap finish EXIT

#-------------------------------------------------------------------------------------------------

bbmap.sh \
    ${MEM} \
    ${THREADS} \
    -eoom \
    in=${READ_1} \
    ref=${REF} \
    outm=${scratch}/map1.sam \
    nodisk \
    overwrite \
    ambiguous=toss \
    2> ${scratch}/guided.err

echo ">>>Map1 Complete" >> ${scratch}/guided.err

#-------------------------------------------------------------------------------------------------

awk '$1 !~ /@/{a[$3]++}END{for(ref in a) print ref, a[ref]}' OFS="\t" ${scratch}/map1.sam \
    | sort -rnk2,2 \
    > ${scratch}/refs.tsv

# if there is no alignment, bail gracefully
if [ ! -s ${scratch}/refs.tsv ]; then
    echo ">>>NO ALIGMENT!" >> ${scratch}/guided.err
    exit 0
fi

new_ref="$(awk 'NR == 1{print $1".fasta"}' ${scratch}/refs.tsv)"

echo ">>>Awk Complete" >> ${scratch}/guided.err

#-------------------------------------------------------------------------------------------------

bbmap.sh \
    ${MEM} \
    ${THREADS} \
    -eoom \
    in=${READ_1} \
    ref=\'"${REF_DIR}/lib/${new_ref}"\' \
    outm=stdout.sam \
    nodisk \
    overwrite \
    maxindel=2000 \
    2>> ${scratch}/guided.err \
    | samtools sort > ${scratch}/map2.bam

samtools index ${scratch}/map2.bam

echo ">>>Map2 Complete" >> ${scratch}/guided.err

#-------------------------------------------------------------------------------------------------

# # BCF tools caller
# bcftools mpileup -Ob \
#     -f ${REF_DIR}/lib/${new_ref} \
#     --max-depth=1000000 \
#     ${scratch}/map2.bam \
#     2>> ${scratch}/guided.err \
#     | bcftools call -Ob \
#     --ploidy=1 \
#     --multiallelic-caller \
#     --variants-only \
#     > ${scratch}/bcftools.bcf
# # bcftools index ${scratch}/bcftools.bcf
# #
# # bcftools consensus \
# #     -f ${REF_DIR}/lib/${new_ref} \
# #     ${scratch}/bcftools.bcf \
# #     > ${scratch}/cons.fa
#
# echo ">>>BCFtools Complete" >> ${scratch}/guided.err

#-------------------------------------------------------------------------------------------------

# freebayes caller
freebayes \
    -f ${REF_DIR}/lib/${new_ref} \
    --pooled-continuous \
    --ploidy 1 \
    --haplotype-length 1 \
    --min-base-quality 20 \
    --min-alternate-fraction 0.5 \
    --min-alternate-count 1 \
    ${scratch}/map2.bam \
    | bcftools convert -Ou \
    > ${scratch}/freebayes.bcf

echo ">>>freebayes Complete" >> ${scratch}/guided.err

#-------------------------------------------------------------------------------------------------

# move relevant files
if [ -e ${scratch}/refs.tsv ]; then
    mv -f ${scratch}/refs.tsv ${OUT_DIR}/${WELL}.refs.tsv
fi
if [ -e ${scratch}/map2.bam ]; then
   mv -f ${scratch}/map2.bam ${OUT_DIR}/${WELL}.map2.bam
fi
if [ -e ${scratch}/map2.bam.bai ]; then
   mv -f ${scratch}/map2.bam.bai ${OUT_DIR}/${WELL}.map2.bam.bai
fi
if [ -e ${scratch}/bcftools.bcf ]; then
   mv -f ${scratch}/bcftools.bcf ${OUT_DIR}/${WELL}.bcftools.bcf
fi
if [ -e ${scratch}/freebayes.bcf ]; then
   mv -f ${scratch}/freebayes.bcf ${OUT_DIR}/${WELL}.freebayes.bcf
fi

echo ">>>mv Complete" >> ${scratch}/guided.err

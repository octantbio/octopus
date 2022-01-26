#!/bin/bash
# bash strict mode
set -euo pipefail
IFS=$'\n\t'

#-------------------------------------------------------------------------------------------------

READ_1=$1
WELL="$(basename -s ".ecc.fq.gz" ${READ_1})"
OUT_DIR="$(dirname ${READ_1})"

THREADS=1
MEM="-Xmx8g"

# send all temp files to the scratch directory
# trap ensures that scratch will be deleted no matter what
# and the log always gets output
scratch="$(mktemp -d -t tmp.XXXXXXXXXX)"
function finish {
    cp -f ${scratch}/log ${OUT_DIR}/${WELL}.de-novo.err
    rm -rf ${scratch}
}
trap finish EXIT

#-------------------------------------------------------------------------------------------------

# attempt to merge the reads together
bbmerge-auto.sh \
    ${MEM} \
    threads=${THREADS} \
    in=${READ_1} \
    out=${scratch}/merged.fq.gz \
    outu=${scratch}/unmerged.fq.gz \
    strict \
    k=93 \
    extend2=80 \
    rem \
    ordered \
    2>> ${scratch}/log

#-------------------------------------------------------------------------------------------------

# quality trim any reads that don't merge
bbduk.sh \
    ${MEM} \
    threads=${THREADS} \
    in=${scratch}/unmerged.fq.gz \
    out=${scratch}/qtrimmed.fq.gz \
    qtrim=r \
    trimq=10 \
    minlen=70 \
    ordered \
    2>> ${scratch}/log

#-------------------------------------------------------------------------------------------------

# try a bunch of spades settings to find best assembly
spades.py \
     --merged ${scratch}/merged.fq.gz \
     --12 ${scratch}/qtrimmed.fq.gz \
     -o ${scratch}/spades \
     --threads=${THREADS} \
     -k 25,55,95,125 \
     --phred-offset 33 \
     --careful \
     >> ${scratch}/log \
     || echo "SPADES FAILED!!!" > ${OUT_DIR}/${WELL}.SPADES-FAIL

#-------------------------------------------------------------------------------------------------

# move relevant items from scratch
if [ -f ${scratch}/spades/contigs.fasta ]; then
    mv -f ${scratch}/spades/contigs.fasta ${OUT_DIR}/${WELL}.spades-contig.fasta
fi

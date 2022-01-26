#!/bin/bash
# bash strict mode
set -euo pipefail
IFS=$'\n\t'

#-------------------------------------------------------------------------------------------------

READ_1="$1"
READ_2="$2"
CONTAM_REF="$3"
WELL="$(basename -s "_R1.fastq.gz" "$1")"
OUT_DIR="$(dirname "$1")"

THREADS=1
MEM="-Xmx8g"

# send all temp files to the scratch directory
# trap ensures that scratch will be deleted no matter what
scratch="$(mktemp -d -t tmp.XXXXXXXXXX)"
function finish {
    cp -f ${scratch}/log ${OUT_DIR}/${WELL}.pre-proc
    rm -rf ${scratch}
}
trap finish EXIT

#-------------------------------------------------------------------------------------------------

# clumpify doesnt work with pipes
# no need to do optical dedupe on miseq
clumpify.sh \
    ${MEM} \
    in1=${READ_1} \
    in2=${READ_2} \
    out=${scratch}/input.fq.gz \
    threads=${THREADS} \
    dedupe \
    2> ${scratch}/log

# trim sequencing adapters, remove contaminants (phix, ecoli), and known artifacts
# then perform the first round of error correction
# NEB-5a genome is hard-coded to src
bbduk.sh \
    ${MEM} \
    in=${scratch}/input.fq.gz \
    out=stdout.fq \
    ref=adapters \
    threads=${THREADS} \
    interleaved=t \
    ktrim=r \
    k=23 \
    mink=11 \
    hdist=1 \
    tbo tpe \
    minlen=70 \
    ftm=5 \
    ordered \
    2> ${scratch}/duk-1.log \
| bbduk.sh \
    ${MEM} \
    in=stdin.fq \
    out=${scratch}/duk.fq.gz \
    threads=${THREADS} \
    interleaved=t \
    k=31 \
    ref=artifacts,phix,${CONTAM_REF} \
    stats=${scratch}/filter.stats \
    ordered \
    cardinality \
    2> ${scratch}/duk-2.log \
    && cat ${scratch}/duk-1.log ${scratch}/duk-2.log >> ${scratch}/log

#-------------------------------------------------------------------------------------------------
# error correction

bbmerge.sh \
    ${MEM} \
    in=${scratch}/duk.fq.gz \
    out=${scratch}/ecc-1.fq.gz \
    threads=${THREADS} \
    interleaved=t \
    ecco \
    mix \
    vstrict \
    ordered \
    2>> ${scratch}/log

# second stage of error correction
clumpify.sh \
    ${MEM} \
    in=${scratch}/ecc-1.fq.gz \
    out=${scratch}/ecc-2.fq.gz \
    threads=${THREADS} \
    interleaved=t \
    ecc \
    passes=4 \
    reorder \
    2>> ${scratch}/log

# thrid stage of error correction
tadpole.sh \
    ${MEM} \
    in=${scratch}/ecc-2.fq.gz \
    out=${scratch}/ecc-3.fq.gz \
    threads=${THREADS} \
    interleaved=t \
    ecc \
    k=62 \
    ordered \
    2>> ${scratch}/log

#-------------------------------------------------------------------------------------------------

mv -f ${scratch}/ecc-3.fq.gz ${OUT_DIR}/${WELL}.ecc.fq.gz

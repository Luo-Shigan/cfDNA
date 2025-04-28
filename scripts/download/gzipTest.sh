#!/bin/bash
function gzipTest(){
    dir=$1
    outfile=$2
    for f in ${1}*.gz; do
        gzip -t "$f" >/dev/null 2>&1 \
        && echo "OK:    $f" \
        || echo "ERROR: $f"
    done > ${outfile} 2>&1
}
export -f gzipTest
dir=$1
outfile=$2
gzipTest ${dir} ${outfile}
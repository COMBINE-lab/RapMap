#!/bin/bash

cmd="$@"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
bam_out=`echo $cmd | sed -n 's/.*--bamOut\s\+\(\S\+\)\s*.*/\1/p'`
bam_compress_threads=`echo $cmd | sed -n 's/.*--bamThreads\s\+\([[:digit:]]\+\)\s*.*/\1/p'`

if [ -z "$bam_out" ]
then
    #Run normally in this branch
    $DIR/rapmap ${@}
else
    new_cmd=`echo $cmd | sed 's/--bamOut\s\+\(\S\+\)\s*//'`
    # The following interleaved to split conversion is courtesy of
    # https://gist.github.com/nathanhaigh/3521724
    execmd=""
    if [ -z "$bam_compress_threads" ]
    then
        execmd="${new_cmd} -o | samtools view -Sb - > ${bam_out}"
    else
        execmd="${new_cmd} -o | samtools view -Sb -@ ${bam_compress_threads} > ${bam_out}"
    fi
    echo "Running command [$DIR/rapmap ${execmd}]"
    $DIR/rapmap ${execmd}
fi

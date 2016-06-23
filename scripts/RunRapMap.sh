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
    # from: http://stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script
    if command -v samtools >/dev/null; then 
        new_cmd=`echo $cmd | sed 's/--bamOut\s\+\(\S\+\)\s*//'`
        execmd=""
        if [ -z "$bam_compress_threads" ]
        then
            execmd="${new_cmd} -o | samtools view -Sb - > ${bam_out}"
        else
            execmd="${new_cmd} -o | samtools view -Sb -@ ${bam_compress_threads} > ${bam_out}"
        fi
        echo "Running command [$DIR/rapmap ${execmd}]"
        $DIR/rapmap ${execmd}
    else
        echo >&2 "samtools is required to convert to BAM, but it's not installed.  Aborting."; 
        exit 1; 
    fi
fi

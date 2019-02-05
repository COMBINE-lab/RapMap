[![Join the chat at https://gitter.im/COMBINE-lab/RapMap](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/COMBINE-lab/RapMap?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

# What is RapMap?

RapMap is a testing ground for ideas in quasi-mapping and selective alignment.  That means that, at this point, it is somewhat experimental.  The `develop` branch will have the latest improvements and additions, but is not guaranteed to be stable between commits.  Breaking changes to the master branch will be accompanied by a tag to the version before the breaking change.  Currently, RapMap is a stand-alone quasi-mapper that can be used with other tools.  It is also being used as part of [Salmon](https://github.com/COMBINE-lab/salmon) and [Sailfish](https://github.com/kingsfordgroup/sailfish).  Eventually, the hope is to create and stabilize an API so that it can be used as a library from other tools.

# Building RapMap

To build RapMap, you need a C++14 compliant compiler (g++ >= 4.9 and clang >= 3.4) and CMake (>= 3.9).  RapMap is built with the following steps (assuming that `path_to_rapmap` is the toplevel directory where you have cloned this repository):

```
[path_to_rapmap] > mkdir build && cd build
[path_to_rapmap/build] > cmake ..
[path_to_rapmap/build] > make
[path_to_rapmap/build] > make install
[path_to_rapmap/build] > cd ../bin
[path_to_rapmap/bin] > ./rapmap -h
```

This should output the standard help message for rapmap.

# Using RapMap

To use RapMap to map reads, you first have to index your reference transcriptome.  Once the index is created, it can be used to map many different sets of reads.  Assuming that your reference transcriptome is in the file `ref.fa`, you can produce the index as follows:

```
> rapmap quasiindex -t ref.fa -i ref_index
```

if you want to make use of a minimum perfect hash when indexing (which will lower the memory requirement during mapping), you can instead use the following command:

```
> rapmap quasiindex -t ref.fa -i ref_index -p -x 4
```

the `-p` option enables the minimum perfect hash and `-x 4` tells RapMap to use up to 4 threads when building the perfect hash (you can specify as many or as few threads as you wish).

The index itself will record whether it was built with the aid of minimum perfect hashing or not, so no extra information concerning this need be provided when mapping.  For the purposes of this example, we'll assume that we wish to map paired-end reads with the first mates in the file `r1.fq.gz` and the second mates in the file `r2.fq.gz`.  We can perform the mapping like so:

```
> rapmap quasimap -i ref_index -1 r1.fq.gz -2 r2.fq.gz -s -t 8 -o mapped_reads.sam
```

This will tell RapMap to map the paired-end reads using 8 threads, and to write the resulting `SAM` records to the file `mapped_reads.sam`.  The `-s` flag tells RapMap to use selective alignment to generate better mappings and to validate the alignment _score_ of hits. The SAM format is rather verbose, and so such output files can be rather large (and slow to write) if you're mapping many reads.  For that reason, we recommend that you use [samtools](http://www.htslib.org/) to convert the `SAM` file to a `BAM` file on-the-fly.  Assuming `samtools` is installed an in your path, that can be accomplished with the following command:

```
> rapmap quasimap -i ref_index -1 r1.fq.gz -2 r2.fq.gz -s -t 8 | samtools view -Sb -@ 4 - > mapped_reads.bam
```

This will stream the output from RapMap to standard out, and then convert it into a `BAM` file (using up to an additional 4 threads for `BAM` compression) and write the resulting output to the file `mapped_reads.bam`.  To reduce the amount that needs to be typed in the common case, and to prevent the user from having to remember invocations like the above, we inclde a simple wrapper script that simplifies this process.  After installing RapMap, there should be a script called `RunRapMap.sh` in the `bin` directory of whereever you have chosen to install RapMap.  You can issue a command equivalent to the above using this scrpt as follows:

```
> RunRapMap.sh quasimap -i ref_index -1 r1.fq.gz -2 r2.fq.gz -s -t 8 --bamOut mapped_reads.sam --bamThreads 4
```

This will run RapMap with a command equivalent to the one mentioned above.  If you leave out the `--bamThreads` argument, then a single thread will be used for compression.  The `RunRapMap.sh` script can be used even if you don't wish to write the output to `BAM` format; in that case it is simply equivalent to running whichever command you pass with the `rapmap` executable itself.

# Can I use RapMap for genomic alignment?

The index and mapping strategy employed by RapMap are highly geared toward mapping to transcriptomes.  This means that RapMap will likely use _a lot_ of memory when indexing and mapping to mammalian-sized genomes, though it's possible.  We have succesfully applied RapMap to map reads to collections of baterial and viral genomes, however.

# Caveats

RapMap is experimental, and the code, at this point, is subject to me testing out new ideas (see the description above about the master vs. develop branch). This also means that limited effort has been put into size or speed optimizaiton.  There are numerous ways that the code can be sped up and the memory footprint reduced, but that hasn't been the focus yet --- it will be eventualy.  All of this being said --- RapMap is open to the community because I'd like feedback / help / thoughts.  A contribution policy is forthcoming.  So, if you're not scared off by any of the above, please *dig in*!

# License 

Since RapMap uses Vinga's [rank implementation](https://github.com/COMBINE-lab/RapMap/blob/master/include/rank9b.h), it must be released under the GPL.  However, this is currently the only GPL dependency.  If it can be replaced, I'd like to re-license RapMap under a BSD license.  I'd be happy to accept pull-requests that replace this rank implementation with a library released under a more liberal license (BSD-compatible), but note that I will *not* accept such pull requests if they reduce the speed or increase the memory consumption over the current version.

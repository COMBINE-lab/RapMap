[![Join the chat at https://gitter.im/COMBINE-lab/RapMap](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/COMBINE-lab/RapMap?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

# What is RapMap?

RapMap is a testing ground for ideas in quasi-mapping / (lightweight / pseudo) transcriptome alignment.  That means that, at this point, it is somewhat experimental.  The `develop` branch will have the latest improvements and additions, but is not guaranteed to be stable between commits.  Breaking changes to the master branch will be accompanied by a tag to the version before the breaking change.  Currently, RapMap is a stand-alone quasi-mapper that can be used with other tools.  It is also being used as part of [Sailfish](https://github.com/kingsfordgroup/sailfish) and [Salmon](https://github.com/COMBINE-lab/salmon).  Eventually, the hope is to create and stabilize an API so that it can be used as a library from other tools.

Quasi-mapping / (lightweight / pseudo)-alignment is the term we are using here for the type of information required for certain tasks (e.g. 
transcript quantification) that is less "heavyweight" than what is provided by traditional alignment. For example, one may
only need to know the transcripts / contigs to which a read aligns and, perhaps, the position within those transcripts rather
than the optimal alignment and base-for-base `CIGAR` string that aligns the read and substring of the transcript.  For details on RapMap (quasi-mapping in particular), please check out the [associated paper](http://bioinformatics.oxfordjournals.org/content/32/12/i192.full.pdf). Note: RapMap implements both quasi-mapping and pseudo-alignment (originally introduced in [Bray et al. 2016](http://www.nature.com/nbt/journal/v34/n5/full/nbt.3519.html)), these two are not the same thing. They are distinct concepts, and RapMap simply happens to implement algorithms for computing both.

There are a number of different ways to collect such information, and the idea of RapMap (as the repository grows) will be to explore multiple different strategies in how to most rapidly determine all *feasible* / *compatible* locations for a read within the transcriptome.  In this sense, it is like an *all-mapper*; the mappings it outputs are intended to be (eventually) disambiguated (*Really, it's more like an "all-best" mapper, since it returns all hits in the top "stratum" of quasi-mapping / (lightweight / pseudo)-alignments*).  If there is a need for it, *best-mapper* functionality may be added in the future.

# Building RapMap

To build RapMap, you need a C++11 compliant compiler (g++ >= 4.7 and clang >= 3.4) and CMake.  RapMap is built with the following steps (assuming that `path_to_rapmap` is the toplevel directory where you have cloned this repository):

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
> rapmap quasimap -i ref_index -1 <(gunzip -c r1.fq.gz) -2 <(gunzip -c r2.fq.gz) -t 8 -o mapped_reads.sam
```

This will tell RapMap to map the paired-end reads using 8 threads, and to write the resulting `SAM` records to the file `mapped_reads.sam`.  The SAM format is rather verbose, and so such output files can be rather large (and slow to write) if you're mapping many reads.  For that reason, we recommend that you use [samtools](http://www.htslib.org/) to convert the `SAM` file to a `BAM` file on-the-fly.  Assuming `samtools` is installed an in your path, that can be accomplished with the following command:

```
> rapmap quasimap -i ref_index -1 <(gunzip -c r1.fq.gz) -2 <(gunzip -c r2.fq.gz) -t 8 | samtools view -Sb -@ 4 - > mapped_reads.bam
```

This will stream the output from RapMap to standard out, and then convert it into a `BAM` file (using up to an additional 4 threads for `BAM` compression) and write the resulting output to the file `mapped_reads.bam`.  To reduce the amount that needs to be typed in the common case, and to prevent the user from having to remember invocations like the above, we inclde a simple wrapper script that simplifies this process.  After installing RapMap, there should be a script called `RunRapMap.sh` in the `bin` directory of whereever you have chosen to install RapMap.  You can issue a command equivalent to the above using this scrpt as follows:

```
> RunRapMap.sh quasimap -i ref_index -1 <(gunzip -c r1.fq.gz) -2 <(gunzip -c r2.fq.gz) -t 8 --bamOut mapped_reads.sam --bamThreads 4
```

This will run RapMap with a command equivalent to the one mentioned above.  If you leave out the `--bamThreads` argument, then a single thread will be used for compression.  The `RunRapMap.sh` script can be used even if you don't wish to write the output to `BAM` format; in that case it is simply equivalent to running whichever command you pass with the `rapmap` executable itself.

# Can I use RapMap for genomic alignment?

No, at least not right now.  The index and mapping strategy employed by RapMap are highly geared toward mapping to transcriptomes.  It may be the case that some of these ideas can be successfully applied to genomic alignment, but 
this functionality is not currently suppored (and is not a high priority right now).

# How fast is RapMap?

Speed is relative, but we think it's very fast: On a synthetic test dataset comprised of 75 million 76bp paired-end reads, mapping to a human transcriptome with ~213,000 transcripts, RapMap takes ~ 10 minutes to align all of the reads *on a single core* (on an Intel Xeon E5-2690 @ 3.00 GHz) --- if you actually want to write out the alignments --- it depends on you disk speed, but for us it's ~15 minutes. Again, these mapping times are *on a single core* --- but RapMap is trivially parallelizable and can be run with multiple threads.  Additionally, there are other optimizations we are currently exploring.

# OK, that's fast, but is it accurate?

Yes; quasi-mapping seems to provide accurate mapping results. In the above mentioned synthetic dataset (generated *with* sequencing errors), the true location of origin of the read appears in the hits returned by RapMap > 97% of the time. For more details, please refer to [the paper](http://bioinformatics.oxfordjournals.org/content/32/12/i192.full.pdf).

# Caveats

RapMap is experimental, and the code, at this point, is subject to me testing out new ideas (see the description above about the master vs. develop branch). This also means that limited effort has been put into size or speed optimizaiton.  There are numerous ways that the code can be sped up and the memory footprint reduced, but that hasn't been the focus yet --- it will be eventualy.  All of this being said --- RapMap is open to the community because I'd like feedback / help / thoughts.  A contribution policy is forthcoming.  So, if you're not scared off by any of the above, please *dig in*!

# External dependencies

[tclap](http://tclap.sourceforge.net/)

[cereal](https://github.com/USCiLab/cereal)

[jellyfish](https://github.com/gmarcais/Jellyfish)

# License 

Since RapMap uses Jellyfish, it must be released under the GPL.  However, this is currently the only GPL dependency.  If it can be replaced, I'd like to re-license RapMap under the BSD license.  I'd be happy to accept pull-requests that replace the Jellyfish components with a library released under a more liberal license (BSD-compatible), but note that I will *not* accept such pull requests if they reduce the speed or increase the memory consumption over the Jellyfish-based version.

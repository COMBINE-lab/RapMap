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

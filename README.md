# What is RapMap?

RapMap is a testing ground for ideas in lightweight / pseudo / quasi transcriptome alignment.  That means that, at this point, it is **very** experimental and there are no guarantees on stability / compatibility between commits.  Eventually, I hope 
that RapMap will become a stand-alone lightweight / pseudo / quasi-aligner that can be used with other tools.

Lightweight / pseudo / quasi-alignment is the term I'm using here for the type of information required for certain tasks (e.g. 
transcript quantification) that is less "heavyweight" than what is provided by traditional alignment. For example, one may
only need to know the transcripts / contigs to which a read aligns and, perhaps, the position within those transcripts rather
than the optimal alignment and base-for-base `CIGAR` string that aligns the read and substring of the transcript.

There are a number of different ways to collect such information, and the idea of RapMap (as the repository grows) will be to explore multiple different strategies in how to most rapidly determine all *feasible* / *compatible* locations for a read within the transcriptome.  In this sense, it is like an *all-mapper*; the alignments it outputs are intended to be (eventually) disambiguated (*Really, it's more like an "all-best" mapper, since it returns all hits in the top "stratum" of lightweight/pseudo/quasi alignments*).  If there is a need for it, *best-mapper* functionality may be added in the future.

# Can I use RapMap for genomic alignment?

No, at least not right now.  The index and mapping strategy employed by RapMap are highly geared toward mapping to transcriptomes.  It may be the case that some of these ideas can be successfully applied to genomic alignment, but 
this functionality is not currently suppored (and is not a high priority right now).

# How fast is RapMap?

It's currently too early in development for a comprehensive benchmark suite, but, on a synthetic test dataset comprised of 
75 million 76bp paired-end reads, mapping to a human transcriptome with ~213,000 transcripts, RapMap takes ~ 10 minutes to 
align all of the reads *on a single core* (on an Intel Xeon E5-2690 @ 3.00 GHz) --- if you actually want to write out the alignments --- it depends on you disk speed, but for us it's ~15 minutes. Again, these mapping times are *on a single core*, and before any significant optimizations (the code is only about a week and a half old) --- but RapMap is trivially parallelizable and can already be run with multiple threads.

# OK, that's fast, but is it accurate?

Again, we're too early in development for a comprehensive benchmark or answer to this question.  However, in the above mentioned synthetic dataset (generated *with* sequencing errors), the true location of origin of the read appears in the hits returned by RapMap > 97% of the time.

# Caveats

RapMap is **highly** experimental, and the code exists primarily at this point for me to test out ideas.  This means that 
there are shortcuts and hacks abound.  It also means that I've not yet put much effort into size or speed optimizaiton.  There are numerous ways that the code can be sped up and the memory footprint reduced, but that hasn't been the focus yet --- it will be eventualy.  All of this being said --- RapMap is open to the community because I'd like feedback / help / thoughts.  So, if you're not scared off by any of the above, *dig in*!

# External dependencies

[tclap](http://tclap.sourceforge.net/)

[cereal](https://github.com/USCiLab/cereal)

# License 

Recently a hot issue in Bioinformatics :).  Since RapMap uses Jellyfish, it must be released under the GPL.  However, this is currently the only GPL dependency.  If it can be replaced, I'd like to re-license RapMap under the BSD license.  I'd be happy to accept pull-requests that replace the Jellyfish components with a library released under a more liberal license (BSD-compatible), but note that I will *not* accept such pull requests if they reduce the speed or increase the memory consumption over the Jellyfish-based version.

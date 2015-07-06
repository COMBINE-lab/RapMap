# What is RapMap?

RapMap is a testing ground for ideas in lightweight / pseudo / quasi transcriptome alignment.  That means that, at this point, it
is **incredibly** experimental and there are no guarantees on stability / compatibility between commits.  Eventually, I hope 
that RapMap *may* become a stand-alone lightweight / pseudo / quasi aligner that can be used with other tools.

Lightweight / pseudo / quasi alignment is the term I'm using here for the type of information required for certain tasks (e.g. 
transcript quantification) that is less "heavyweight" than what is provided by traditional alignment. For example, one may
only need to know the transcripts / contigs to which a read aligns and, perhaps, the position within those transcripts rather
than the optimal alignment and base-for-base `CIGAR` string that aligns the read and substring of the transcript.

RapMap will explore different ideas in how to most rapidly determine all *feasible* / *compatible* locations for a read within 
the transcriptome.  In this sense, it is like an *all-mapper*; the alignments it outputs are intended to be (eventually) 
disambiguated.  If there is a need for it, *best-mapper* functionality may be added in the future.

# How fast is RapMap?

It's currently too early in development for a comprehensive benchmark suite, but, on a synthetic test dataset comprised of 
75 million 76bp paired-end reads, mapping to a human transcriptome with ~213,000 transcripts, RapMap takes ~10 minutes to 
align all of the reads *on a single core* (on an Intel Xeon E5-2690 @ 3.00 GHz) --- if you actually want to write out a horribly
inefficient ascii file of the alignments (currently working on a better format) --- the whole process takes ~20 minutes. Again,
these mapping times are *on a single core*, and before any significant optimizations (the code is only a few days old) --- 
but RapMap is trivially parallelizable.

# Caveats

RapMap is **incredibly** experimental, and the code exists primarily at this point for me to test out ideas.  This means that 
there are shortcuts and hacks abound.  For example, there is currently no proper command line parsing (it's expected that the 
right arguments are passed in in the right order).  It currently expects only paired-end reads (this will change soon). It also means that I've not yet put much effort into size optimizaiton --- so the RapMap index on the 300M Human transcriptome with 213,121 transcripts is ~4.8G.  There are ways that this can be made **considerably** smaller, but that hasn't been the focus yet.  All of this being said --- RapMap is open to the community because I'd like feedback / help / thoughts.  So, if you're not scared off by any of the above, *dig in*!



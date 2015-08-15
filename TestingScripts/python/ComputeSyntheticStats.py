import argparse

def main(args):
    import pysam
    import sys

    alnFile = pysam.AlignmentFile(args.input, 'r')
    currQueryName = ''
    skipToNextAlignment = False

    readsWithTrueAln = 0
    readsSeen = 0

    for rec in alnFile:
        qname = rec.qname[:-2]
        # If this is a new read, remember the name
        # and increment the counter
        if qname != currQueryName:
            readsSeen += 1
            if readsSeen % 1000000 == 0:
                print("\r\rSaw {} reads --- thp = {:.2%}".format(readsSeen, \
                       float(readsWithTrueAln) / readsSeen), file=sys.stderr, end='')
            currQueryName = qname
            skipToNextAlignment = False

        # If we already found the true hit
        # for this read, don't bother processing
        # this record
        if not skipToNextAlignment:
            # Process the record to find if its the true
            # hit
            currQueryName = qname
            trueTxpName = qname.split(':')[2]
            alignedTxpName = alnFile.getrname(rec.rname)
            if (trueTxpName == alignedTxpName):
                readsWithTrueAln += 1
                skipToNextAlignment = True

    print("\nSaw {} reads, {} of them had a correct alignment: THP = {:.2%}".format(\
           readsSeen, readsWithTrueAln, float(readsWithTrueAln) / readsSeen, \
           file=sys.stderr))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute statistics from synthetic SAM file.')
    parser.add_argument('--input', type=str)
    args = parser.parse_args()
    main(args)

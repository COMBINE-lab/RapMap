from __future__ import print_function
import argparse

def main(args):
    import pysam
    import sys

    printFP = args.printFalsePositives
    alnFile = pysam.AlignmentFile(args.input, 'r')
    totalNumReads = args.totalNumReads
    currQueryName = ''
    skipToNextAlignment = False

    readsWithTrueAln = 0
    readsSeen = 0

    truePos = 0
    falsePos = 0
    falseNeg = 0
    foundTrueAlignment = False
    ## True Neg are somewhat ill-defined in our context
    prevRec = None
    useNHTag = args.useNHTag
    NHTagSum = 0
    nhData = []
    for rec in alnFile:
        if rec.is_unmapped:
            continue
        qname = rec.qname[:-2]
        # If this is a new read, remember the name
        # and increment the counter
        if qname != currQueryName:
            readsSeen += 1
            if useNHTag:
                tagVal = rec.get_tag('NH')
                NHTagSum += tagVal
                nhData.append(tagVal)
            if readsSeen % 1000000 == 0:
                print("\r\rSaw {} reads --- thp = {:.2%}".format(readsSeen, \
                       float(truePos) / readsSeen), file=sys.stderr, end='')


            currQueryName = qname
            if not foundTrueAlignment:
                falsePos += 1
                if (printFP):
                    print(prevRec)
            skipToNextAlignment = False
            foundTrueAlignment = False

        prevRec = rec
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
                truePos += 1
                skipToNextAlignment = True
                foundTrueAlignment = True

    falseNeg = totalNumReads - readsSeen
    print('\n'.join(["Total Reads = {}" ,
          "Reads Aligned  = {}" ,
          "True Pos = {}" ,
          "False Pos = {}" ,
          "False Neg = {}" ,
          "Precision = {:.2%}",
          "Recall = {:.2%}" ,
          "TPH = {:.2%}",
          "FDR = {:.2%}",
          "F1 Score = {:.2%}"]).format(totalNumReads, readsSeen,
                                  truePos, falsePos,
                                  falseNeg, truePos / float(truePos + falsePos),
                                  truePos / float(truePos + falseNeg),
                                  truePos / float(readsSeen),
                                  falsePos / float(falsePos + truePos),
                                  (2*truePos) / float(2*truePos + falsePos + falseNeg)),
          file=sys.stderr)

    if useNHTag:
        print("Average hits-per-read = {:.2}".format(
              NHTagSum / float(totalNumReads)), file=sys.stderr)
        if args.nhFreqFile:
            import pickle
            pickle.dump(nhData, args.nhFreqFile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute statistics from synthetic SAM file.')
    parser.add_argument('--input', type=str)
    parser.add_argument('--totalNumReads', type=int)
    parser.add_argument('--useNHTag', action='store_true')
    parser.add_argument('--nhFreqFile', nargs='?', type=argparse.FileType('w'))
    parser.add_argument('--printFalsePositives', action='store_true')
    args = parser.parse_args()
    main(args)

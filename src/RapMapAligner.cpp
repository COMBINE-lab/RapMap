#include <limits>
#include <iostream>
#include "RapMapAligner.hpp"
#include "RapMapUtils.hpp"

void RapMapAligner::init() {
    tm[0][0] = TraceBack::END;
    tx[0][0] = TraceBack::END;
    ty[0][0] = TraceBack::END;
    m[0][0] = 0;
    x[0][0] = 0;
    y[0][0] = 0;
    for (int i = 1; i < maxRefLen; ++i) {
        m[i][0] = std::numeric_limits<short>::min();
        x[i][0] = std::numeric_limits<short>::min();
        if (freeGapsBeforeRead) {
            y[i][0] = 0;
            ty[i][0] = TraceBack::END;
        } else {
            y[i][0] = gapStart + (i - 1) * gapExtend;
            ty[i][0] = TraceBack::GAP_IN_Y;
        }
    }
    for (int j = 1; j < maxReadLen; ++j) {
        m[0][j] = std::numeric_limits<short>::min();
        x[0][j] = gapStart + (j - 1) * gapExtend;
        tx[0][j] = TraceBack::GAP_IN_X;
        y[0][j] = std::numeric_limits<short>::min();
    }
}

int RapMapAligner::align(std::string& ref, size_t refStart, size_t refLen,
                         std::string& read, size_t readStart, size_t readLen) {
    int max, xTemp, yTemp, score;
    TraceBack trace;

    // Assumption
    assert(readLen < maxReadLen);
    assert(refLen < maxRefLen);

    for (int i = 1; i <= refLen; ++i) {
        for (int j = 1; j <= readLen; ++j) {
            // Update m
            max = m[i - 1][j - 1];
            trace = TraceBack::MATCH;
            xTemp = x[i - 1][j - 1];
            yTemp = y[i - 1][j - 1];
            if (max < xTemp) {
                max = xTemp;
                trace = TraceBack::GAP_IN_X;
            }
            if (max < yTemp) {
                max = yTemp;
                trace = TraceBack::GAP_IN_Y;
            }
            // Save match/mismatch to distinguish X vs =
            if (rapmap::utils::matches(ref[refStart + i - 1], read[readStart + j - 1])) {
                score = match;
            } else {
                score = misMatch;
                trace |= TraceBack::MISMATCH;
            }
            m[i][j] = score + max;
            tm[i][j] = trace;
            // Update x
            max = gapStart + m[i][j - 1];
            trace = TraceBack::MATCH;
            xTemp = gapExtend + x[i][j - 1];
            yTemp = gapStart + y[i][j - 1];
            if (max < xTemp) {
                max = xTemp;
                trace = TraceBack::GAP_IN_X;
            }
            if (max < yTemp) {
                max = yTemp;
                trace = TraceBack::GAP_IN_Y;
            }
            x[i][j] = max;
            tx[i][j] = trace;
            // Update y
            max = gapStart + m[i - 1][j];
            trace = TraceBack::MATCH;
            xTemp = gapStart + x[i - 1][j];
            yTemp = gapExtend + y[i - 1][j];
            if (max < xTemp) {
                max = xTemp;
                trace = TraceBack::GAP_IN_X;
            }
            if (max < yTemp) {
                max = yTemp;
                trace = TraceBack::GAP_IN_Y;
            }
            y[i][j] = max;
            ty[i][j] = trace;
        }
    }

    // Find maximum total score
    size_t maxI = refLen;
    max = m[maxI][readLen];
    trace = TraceBack::MATCH;
    xTemp = x[maxI][readLen];
    yTemp = y[maxI][readLen];
    if (max < xTemp) {
        max = xTemp;
        trace = TraceBack::GAP_IN_X;
    }
    if (max < yTemp) {
        max = yTemp;
        trace = TraceBack::GAP_IN_Y;
    }

    // Need to find highest score in last column
    if (freeGapsAfterRead) {
        for (size_t i = 0; i < refLen; ++i) {
            int mtemp = m[i][readLen];
            xTemp = x[i][readLen];
            yTemp = y[i][readLen];
            if (max < mtemp) {
                max = mtemp;
                maxI = i;
                trace = TraceBack::MATCH;
            }
            if (max < xTemp) {
                max = xTemp;
                maxI = i;
                trace = TraceBack::GAP_IN_X;
            }
            if (max < yTemp) {
                max = yTemp;
                maxI = i;
                trace = TraceBack::GAP_IN_Y;
            }
        }
    }

    // Save info for trace call
    trace_ = trace;
    maxI_ = maxI;
    maxJ_ = readLen;
    return max;
}

void RapMapAligner::trace(CigarString& cigar) {
    TraceBack trace = trace_;
    bool tracing = true;
    size_t i = maxI_, j = maxJ_;

    while (tracing and (i > 0 or j > 0)) {
        if (freeGapsBeforeRead and j == 0)
            break;
        switch (trace) {
            case TraceBack::MATCH:
                trace = tm[i][j];
                if (trace & TraceBack::MISMATCH) {
                    trace ^=  TraceBack::MISMATCH; // clear bit
                    cigar.emplace_front(CigarOp::X);
                } else {
                    cigar.emplace_front(CigarOp::EQ);
                }
                --i;
                --j;
                break;
            case TraceBack::GAP_IN_X:
                cigar.emplace_front(CigarOp::I);
                trace = tx[i][j];
                --j;
                break;
            case TraceBack::GAP_IN_Y:
                cigar.emplace_front(CigarOp::D);
                trace = ty[i][j];
                --i;
                break;
            default:
                tracing = false;
        }
    }
}


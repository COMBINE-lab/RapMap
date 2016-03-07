#ifndef __RAPMAP_RAPMAPALIGNER_HPP__
#define __RAPMAP_RAPMAPALIGNER_HPP__

#include <string>
#include <deque>

enum TraceBack : uint8_t {
    END = 0,
    MATCH = 1,
    GAP_IN_X = 2,
    GAP_IN_Y = 4,
    MISMATCH = 8,
};

inline TraceBack operator|(TraceBack a, TraceBack b) {
    return static_cast<TraceBack>(static_cast<uint8_t>(a) | static_cast<uint8_t>(b));
}

inline TraceBack operator&(TraceBack a, TraceBack b) {
    return static_cast<TraceBack>(static_cast<uint8_t>(a) & static_cast<uint8_t>(b));
}

inline TraceBack operator^(TraceBack a, TraceBack b) {
    return static_cast<TraceBack>(static_cast<uint8_t>(a) ^ static_cast<uint8_t>(b));
}

inline TraceBack& operator|=(TraceBack& a, TraceBack b) {
    a = a | b;
    return a;
}

inline TraceBack& operator&=(TraceBack& a, TraceBack b) {
    a = a & b;
    return a;
}

inline TraceBack& operator^=(TraceBack& a, TraceBack b) {
    a = a ^ b;
    return a;
}

enum CigarOp : char {
    M = 'M',  // alignment match (can be a sequence match or mismatch)
    I = 'I',  // insertion to the reference
    D = 'D',  // deletion from the reference
    N = 'N',  // skipped region from the reference
    S = 'S',  // soft clipping (clipped sequences present in SEQ)
    H = 'H',  // hard clipping (clipped sequences NOT present in SEQ)
    P = 'P',  // padding (silent deletion from padded reference)
    EQ = '=', // sequence match
    X = 'X'   // sequence mismatch
};

class CigarElement {
public:
    int count;
    CigarOp op;

    CigarElement() : count(0), op(CigarOp::M) { }

    CigarElement(int count, CigarOp op) : count(count), op(op) { }
};

class CigarString {
public:
    std::deque<CigarElement> cigar;

    CigarString() : cigar() { }

    inline void emplace_front(CigarOp op) {
        if (cigar.size() > 0 and cigar.front().op == op) {
            cigar.front().count++;
        } else {
            cigar.emplace_front(1, op);
        }
    };

    inline void emplace_front(CigarOp op, int count) {
        if (cigar.size() > 0 and cigar.front().op == op) {
            cigar.front().count += count;
        } else {
            cigar.emplace_front(count, op);
        }
    };

    inline void emplace_back(CigarOp op) {
        if (cigar.size() > 0 and cigar.back().op == op) {
            cigar.back().count++;
        } else {
            cigar.emplace_back(1, op);
        }
    };

    inline void emplace_back(CigarOp op, int count) {
        if (cigar.size() > 0 and cigar.back().op == op) {
            cigar.back().count += count;
        } else {
            cigar.emplace_back(count, op);
        }
    };

    inline void toString(std::string& out) {
        for (auto i = 0; i < cigar.size(); ++i) {
            out += std::to_string(cigar[i].count);
            out += static_cast<char>(cigar[i].op);
        }
    }
};


class RapMapAligner {
private:
    enum : int {
        maxRefLen = 192,
        maxReadLen = 64
    };
    TraceBack trace_;      // starting matrix from previous align
    size_t maxI_, maxJ_;   // Start of trace from previous align

    void init();

public:
    int m[maxRefLen][maxReadLen];
    int x[maxRefLen][maxReadLen];
    int y[maxRefLen][maxReadLen];
    TraceBack tm[maxRefLen][maxReadLen];
    TraceBack tx[maxRefLen][maxReadLen];
    TraceBack ty[maxRefLen][maxReadLen];

    int match;
    int misMatch;
    int gapExtend;
    int gapStart;
    bool freeGapsBeforeRead;
    bool freeGapsAfterRead;

    RapMapAligner() :
            match(1),
            misMatch(-3),
            gapExtend(-1),
            gapStart(-3),
            freeGapsBeforeRead(false),
            freeGapsAfterRead(false) {
        init();
    }

    RapMapAligner(bool gapBefore, bool gapAfter) :
            match(1),
            misMatch(-3),
            gapExtend(-1),
            gapStart(-3),
            freeGapsBeforeRead(gapBefore),
            freeGapsAfterRead(gapAfter) {
        init();
    }

    RapMapAligner(int matchIn, int misMatchIn, int gapExtendIn, int gapStartIn,
                  bool gapBefore, bool gapAfter) :
            match(matchIn),
            misMatch(misMatchIn),
            gapExtend(gapExtendIn),
            gapStart(gapStartIn),
            freeGapsBeforeRead(gapBefore),
            freeGapsAfterRead(gapAfter) {
        init();
    }

    int align(std::string& ref, size_t refStart, size_t refLen,
              std::string& read, size_t readStart, size_t readLen);

    inline int align(std::string& ref, size_t refStart, size_t refLen,
                     std::string& read, size_t readStart, size_t readLen,
                     std::string& cigar) {
        int score = align(ref, refStart, refLen, read, readStart, readLen);
        trace(cigar);
        return score;
    }

    inline int align(std::string& ref, size_t refStart, size_t refLen,
                     std::string& read, size_t readStart, size_t readLen,
                     CigarString& cigar) {
        int score = align(ref, refStart, refLen, read, readStart, readLen);
        trace(cigar);
        return score;
    }

    void trace(CigarString& cigar);

    inline void trace(std::string& cigarOut) {
        CigarString cigar;
        trace(cigar);
        cigar.toString(cigarOut);
    }
};


#endif //__RAPMAP_RAPMAPALIGNER_HPP__

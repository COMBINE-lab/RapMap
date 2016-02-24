#include "emphf/common.hpp"

namespace emphf {
    std::ostream& logger() {
        time_t t = std::time(nullptr);
        // XXX(ot): put_time unsupported in g++ 4.7
        // return std::cerr
        //     <<  std::put_time(std::localtime(&t), "%F %T")
        //     << ": ";
        std::locale loc;
        const std::time_put<char>& tp =
            std::use_facet<std::time_put<char>>(loc);
        const char *fmt = "%F %T";
        tp.put(std::cerr, std::cerr, ' ',
               std::localtime(&t), fmt, fmt + strlen(fmt));
        return std::cerr << ": ";
    }
    uint64_t nonzero_pairs(uint64_t x) {
        static const uint64_t ones_step_4  = 0x1111111111111111ULL;
        x = (x | (x >> 1)) & (0x5 * ones_step_4);

#if EMPHF_USE_POPCOUNT
        return (uint64_t)__builtin_popcountll(x);
#else
        static const uint64_t ones_step_8  = 0x0101010101010101ULL;
        x = (x & 3 * ones_step_4) + ((x >> 2) & 3 * ones_step_4);
        x = (x + (x >> 4)) & 0x0f * ones_step_8;
        return (x * ones_step_8) >> 56;
#endif
    }


}

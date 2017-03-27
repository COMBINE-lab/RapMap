
#include "RapMapSAIndex.hpp"
#include "RapMapUtils.hpp"
#include "SASearcher.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
#include <string>
#include "EditDistance.hpp"


typedef std::vector<int>::iterator iterv;
void expand(std::vector<int>& V_data, iterv& Vf, iterv& Vr, int& R, const int& D, const int& delta) {
    int Rp = R + (R>>1);
    V_data.resize(2 + 4*Rp);
    Vf = V_data.begin() + R;
    Vr = V_data.begin() + (3*R+1) - delta;
    iterv Vp = V_data.begin() + (3*Rp+1) - delta;
    for (int j=D+delta;  j >= -D+delta;  --j) Vp[j] = Vr[j];
    Vr = Vp;
    Vp = V_data.begin() + Rp;
    for (int j=D;  j >= -D;  --j) Vp[j] = Vf[j];
    Vf = Vp;
    R = Rp;
}

/*
 * Mayer's fast bit vector algorithm
 * for approximate string matching
 * adopted from
 */
int edit_distance(std::string& seq1, std::string& seq2, const int max_cost){
    //itr1_t seq1 = boost::begin(seq1_);
    //itr2_t seq2 = boost::begin(seq2_);
    //intype len1 = distance(seq1_);
    //intype len2 = distance(seq2_);

    int len1 = seq1.length();
    int len2 = seq2.length();
    // identify any equal suffix and/or prefix
    int eqb = 0;
    for (;  eqb < std::min(len1, len2);  ++eqb) if (seq1[eqb] != seq2[eqb]) break;
    int eqe = len1-1;
    for (int j2 = len2-1;  eqe > eqb && j2 > eqb;  --eqe,--j2) if (seq1[eqe] != seq2[j2]) break;
    eqe = len1-1-eqe;
    // sub-strings with equal suffix and/or prefix stripped
    const int S1 = eqb;
    const int L1 = len1-(eqb+eqe);
    const int S2 = eqb;
    const int L2 = len2-(eqb+eqe);
    // either or both strings are empty:
    if (L1 <= 0) return L2;
    if (L2 <= 0) return L1;
    const int delta = L1-L2;
    const bool delta_even = delta%2 == 0;
    int R = 10;
    std::vector<int> V_data(2*(1 + 2*R));
    //typedef std::vector<int>::iterator iterv;
    iterv Vf = V_data.begin()+R;
    iterv Vr = V_data.begin()+(3*R+1)-delta;

    int max_cost_check{max_cost};
    int D = 0;
    Vf[1] = 0;
    Vr[-1+delta] = L1;
    while (true) {
        // advance forward-path diagonals:
    //if(D<20)
    //std::cout<<D<<"\n";
        for (int k = -D;  k <= D;  k += 2) {
        //std::cout<<"bbb\n";
            int j1 = (k == -D  ||  (k != D  &&  Vf[k-1] < Vf[k+1]))  ?  Vf[k+1]  :  1+Vf[k-1];
            int j2 = j1-k;
            if (!delta_even  &&  (k-delta) >= -(D-1)  &&  (k-delta) <= (D-1)) {
                int r1 = Vr[k];
                int r2 = Vr[k]-k;
                if ((j1-j2) == (r1-r2)  &&  j1 >= r1){ return 2*D-1;

                }
            }
            while (j1 < L1  &&  j2 < L2  &&  (seq1[j1] == seq2[j2])) { ++j1;  ++j2; }
            Vf[k] = j1;
        //std::cout<<Vf[k] <<"a\n";
        }
        // advance the reverse-path diagonals:
        for (int k = -D+delta;  k <= D+delta;  k += 2) {
            int j1 = (k == D+delta  ||  (k != -D+delta  &&  Vr[k-1] < Vr[k+1]))  ?  Vr[k-1]  :  Vr[k+1]-1;
            int j2 = j1-k;
            if (delta_even  &&  k >= -D  &&  k <= D) {
                int f1 = Vf[k];
                int f2 = Vf[k]-k;
                if ((j1-j2) == (f1-f2)  &&  f1 >= j1) return 2*D;
            }
            while (j1 > 0  &&  j2 > 0  &&  (seq1[j1-1] == seq2[j2-1])) { --j1;  --j2; }
            Vr[k] = j1;
        }
    //if (D<10)
    //std::cout<< D <<"\n";
        if (((delta_even) ? (2*D+2) : (2*D+1))>max_cost) {
        //std::cout<< "max_cost exceeded!\n";
        //std::cout<<((delta_even))<<" " << (2*D+2)<<" " << (2*D+1) << "\n";
            return 255;
        }
        if (D >= R) expand(V_data, Vf, Vr, R, D, delta);
        ++D;
    }
    // control should not reach here
    //BOOST_ASSERT(false);*/
    return 0;
}

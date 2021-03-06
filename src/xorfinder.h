/*
 * CryptoMiniSat
 *
 * Copyright (c) 2009-2015, Mate Soos. All rights reserved.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation
 * version 2.0 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */

#ifndef _XORFINDER_H_
#define _XORFINDER_H_

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include <set>
#include "xor.h"
#include "cset.h"
#include "watcharray.h"

using std::vector;
using std::set;

namespace CMSat {

//#define VERBOSE_DEBUG_XOR_FINDER

class Solver;
class OccSimplifier;

class PossibleXor
{
    public:
        PossibleXor()
        {
        }

        void setup(
            const vector<Lit>& cl
            , const ClOffset offset
            , cl_abst_type _abst
            , vector<uint16_t>& seen
        ) {
            abst = _abst;
            size = cl.size();
            offsets.clear();
            #ifdef VERBOSE_DEBUG_XOR_FINDER
            cout << "Trying to create XOR from clause: " << cl << endl;
            #endif

            assert(cl.size() <= sizeof(origCl)/sizeof(Lit));
            for(size_t i = 0; i < size; i++) {
                origCl[i] = cl[i];
                if (i > 0)
                    assert(cl[i-1] < cl[i]);
            }
            setup_seen_rhs_foundcomb(seen);
            if (offset != std::numeric_limits<ClOffset>::max()) {
                offsets.push_back(offset);
            }
        }

        //GET-type functions
        cl_abst_type      getAbst() const;
        uint32_t          getSize() const;
        bool              getRHS() const;
        bool              foundAll() const;

        //Add
        template<class T>
        void add(const T& cl, const ClOffset offset, vector<uint32_t>& varsMissing);

        const vector<ClOffset>& get_offsets() const
        {
            return offsets;
        }

    private:
        void setup_seen_rhs_foundcomb(vector<uint16_t>& seen)
        {
            //Calculate parameters of base clause.
            //Also set 'seen' for easy check in 'findXorMatch()'
            rhs = true;
            uint32_t whichOne = 0;
            for (uint32_t i = 0; i < size; i++) {
                rhs ^= origCl[i].sign();
                whichOne += ((uint32_t)origCl[i].sign()) << i;
                seen[origCl[i].var()] = 1;
            }

            foundComb.clear();
            foundComb.resize(1UL<<size, false);
            foundComb[whichOne] = true;
        }
        uint32_t NumberOfSetBits(uint32_t i) const;
        bool     bit(const uint32_t a, const uint32_t b) const;

        //bitfield to indicate which of the following is already set
        //-1 -2 -3
        //-1  2  3
        // 1 -2  3
        // 1  2 -3
        //order the above according to sign: if sign:
        //LSB ... MSB
        // 1 1 1
        // 1 0 0
        // 0 1 0
        // 0 0 1
        vector<bool> foundComb;
        Lit origCl[7];
        cl_abst_type abst;
        uint32_t size;
        bool rhs;
        vector<ClOffset> offsets;
};

class XorFinder
{
public:
    XorFinder(OccSimplifier* occsimplifier, Solver* solver);
    void find_xors();

    struct Stats
    {
        void clear()
        {
            Stats tmp;
            *this = tmp;
        }

        Stats& operator+=(const Stats& other);
        void print_short(const Solver* solver) const;
        void print() const;

        //Time
        uint32_t numCalls = 0;
        double findTime = 0.0;
        uint32_t time_outs = 0;

        //XOR stats
        uint64_t foundXors = 0;
        uint64_t sumSizeXors = 0;
    };

    const Stats& get_stats() const;
    virtual size_t mem_used() const;
    void add_xors_to_gauss();
    void clean_up_xors();
    void recursively_xor_xors();
    bool add_new_truths_from_xors();

    vector<Xor> xors;
    vector<ClOffset> cls_of_xors;

private:
    PossibleXor poss_xor;
    void add_found_xor(const Xor& found_xor);
    void find_xors_based_on_short_clauses();
    void find_xors_based_on_long_clauses();
    void print_found_xors();
    bool xor_clause_already_inside(const Xor& xor_c);
    bool xor_has_interesting_var(const Xor& x);
    void clean_occur_from_idxs(const Lit lit, size_t idx1, size_t idx2);
    void clean_occur_from_idx(const Lit lit, size_t idx1);
    vector<uint32_t> xor_two(Xor& x1, Xor& x2, const size_t idx1, const size_t idx2, const uint32_t v);
    void clean_xors_from_empty();

    int64_t xor_find_time_limit;

    //Find XORs
    void findXor(vector<Lit>& lits, const ClOffset offset, cl_abst_type abst);

    ///Normal finding of matching clause for XOR
    void findXorMatch(watch_subarray_const occ);
    void findXorMatch(
        watch_subarray_const occ
        , const Lit lit
    );
    void findXorMatchExt(
        watch_subarray_const occ
        , const Lit lit
    );
    //TODO stamping finXorMatch with stamp
    /*void findXorMatch(
        const vector<LitExtra>& lits
        , const Lit lit
    ) const;*/

    OccSimplifier* occsimplifier;
    Solver *solver;

    //Stats
    Stats runStats;
    Stats globalStats;

    //Temporary
    vector<Lit> tmpClause;
    vector<uint32_t> varsMissing;

    //Other temporaries
    vector<uint16_t>& seen;
    vector<uint16_t>& seen2;
    vector<Lit>& toClear;
};


inline cl_abst_type PossibleXor::getAbst() const
{
    return abst;
}

inline uint32_t PossibleXor::getSize() const
{
    return size;
}

inline bool PossibleXor::getRHS() const
{
    return rhs;
}

template<class T> void PossibleXor::add(
    const T& cl
    , const ClOffset offset
    , vector<uint32_t>& varsMissing
) {
    #ifdef VERBOSE_DEBUG_XOR_FINDER
    cout << "Adding to XOR: " << cl << endl;

    cout << "FoundComb before:" << endl;
    for(size_t i = 0; i < foundComb.size(); i++) {
        cout << "foundComb[" << i << "]: " << foundComb[i] << endl;
    }
    cout << "----" << endl;
    #endif

    //It's the base clause, skip.
    if (!offsets.empty() && offset == offsets[0])
        return;

    assert(cl.size() <= size);

    //If clause covers more than one combination, this is used to calculate which ones
    varsMissing.clear();

    //Position of literal in the ORIGINAL clause.
    //This may be larger than the position in the current clause (as some literals could be missing)
    uint32_t origI = 0;

    //Position in current clause
    uint32_t i = 0;

    //Used to calculate this clause covers which combination(s)
    uint32_t whichOne = 0;

    bool thisRhs = true;

    for (typename T::const_iterator
        l = cl.begin(), end = cl.end()
        ; l != end
        ; l++, i++, origI++
    ) {
        thisRhs ^= l->sign();

        //some variables might be missing in the middle
        while(cl[i].var() != origCl[origI].var()) {
            varsMissing.push_back(origI);
            origI++;
            assert(origI < size && "cl must be sorted");
        }
        whichOne += ((uint32_t)l->sign()) << origI;
    }

    //if vars are missing from the end
    while(origI < size) {
        varsMissing.push_back(origI);
        origI++;
    }

    assert(cl.size() < size || rhs == thisRhs);

    //set to true every combination for the missing variables
    for (uint32_t i = 0; i < 1UL<<(varsMissing.size()); i++) {
        uint32_t thisWhichOne = whichOne;
        for (uint32_t i2 = 0; i2 < varsMissing.size(); i2++) {
            if (bit(i, i2)) thisWhichOne+= 1<<(varsMissing[i2]);
        }
        foundComb[thisWhichOne] = true;
    }
    if (offset != std::numeric_limits<ClOffset>::max()) {
        offsets.push_back(offset);
    }

    #ifdef VERBOSE_DEBUG_XOR_FINDER
    cout << "whichOne was:" << whichOne << endl;
    cout << "FoundComb after:" << endl;
    for(size_t i = 0; i < foundComb.size(); i++) {
        cout << "foundComb[" << i << "]: " << foundComb[i] << endl;
    }
    cout << "----" << endl;
    #endif
}

inline bool PossibleXor::foundAll() const
{
    bool OK = true;
    for (uint32_t i = 0; i < foundComb.size(); i++) {
        //Only count combinations with the correct RHS
        if ((NumberOfSetBits(i)%2) == (uint32_t)rhs) {
            continue;
        }

        //If this combination hasn't been found yet, then the XOR is not complete
        if (!foundComb[i]) {
            OK = false;
            break;
        }
    }

    #ifdef VERBOSE_DEBUG_XOR_FINDER
    if (OK) {
        cout << "Found all for this clause" << endl;
    }
    #endif

    return OK;
}

inline uint32_t PossibleXor::NumberOfSetBits(uint32_t i) const
{
    //Magic is coming! (copied from some book.... never trust code like this!)
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return (((i + (i >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
}

inline bool PossibleXor::bit(const uint32_t a, const uint32_t b) const
{
    return (((a)>>(b))&1);
}

inline const XorFinder::Stats& XorFinder::get_stats() const
{
    return globalStats;
}

} //end namespace

#endif //_XORFINDER_H_

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

#include "xorfinder.h"
#include "time_mem.h"
#include "solver.h"
#include "occsimplifier.h"
#include "clauseallocator.h"
#include "sqlstats.h"

#include <limits>

using namespace CMSat;
using std::cout;
using std::endl;

XorFinder::XorFinder(OccSimplifier* _occsimplifier, Solver* _solver) :
    occsimplifier(_occsimplifier)
    , solver(_solver)
    , seen(_solver->seen)
    , seen2(_solver->seen2)
    , toClear(_solver->toClear)
{
}

void XorFinder::find_xors_based_on_long_clauses()
{
    vector<Lit> lits;
    for (vector<ClOffset>::iterator
        it = occsimplifier->clauses.begin()
        , end = occsimplifier->clauses.end()
        ; it != end && xor_find_time_limit > 0
        ; ++it
    ) {
        ClOffset offset = *it;
        Clause* cl = solver->cl_alloc.ptr(offset);
        xor_find_time_limit -= 3;

        //Already freed
        if (cl->freed() || cl->getRemoved()) {
            continue;
        }

        //Too large -> too expensive
        if (cl->size() > solver->conf.maxXorToFind) {
            continue;
        }

        //If not tried already, find an XOR with it
        if (!cl->stats.marked_clause) {
            cl->stats.marked_clause = true;

            lits.resize(cl->size());
            std::copy(cl->begin(), cl->end(), lits.begin());
            findXor(lits, offset, cl->abst);
        }
    }
}

void XorFinder::find_xors_based_on_short_clauses()
{
    assert(solver->ok);

    vector<Lit> lits;
    for (size_t wsLit = 0, end = 2*solver->nVars()
        ; wsLit < end && xor_find_time_limit > 0
        ; wsLit++
    ) {
        const Lit lit = Lit::toLit(wsLit);
        assert(solver->watches.size() > wsLit);
        watch_subarray_const ws = solver->watches[lit];

        xor_find_time_limit -= (int64_t)ws.size()*4;

        //cannot use iterators because findXor may update the watchlist
        for (size_t i = 0, size = solver->watches[lit].size()
            ; i < size && xor_find_time_limit > 0
            ; i++
        ) {
            const Watched& w = ws[i];

            //Only care about tertiaries
            if (!w.isTri())
                continue;

            //Only bother about each tri-clause once
            if (lit > w.lit2()
                || w.lit2() > w.lit3()
            ) {
                continue;
            }

            lits.resize(3);
            lits[0] = lit;
            lits[1] = w.lit2();
            lits[2] = w.lit3();

            findXor(lits, CL_OFFSET_MAX, calcAbstraction(lits));
        }
    }
}

void XorFinder::find_xors()
{
    runStats.clear();
    runStats.numCalls = 1;

    cls_of_xors.clear();
    xors.clear();
    assert(solver->watches.get_smudged_list().empty());
    double myTime = cpuTime();
    const int64_t orig_xor_find_time_limit =
        1000LL*1000LL*solver->conf.xor_finder_time_limitM
        *solver->conf.global_timeout_multiplier;

    xor_find_time_limit = orig_xor_find_time_limit;

    #ifdef DEBUG_MARKED_CLAUSE
    assert(solver->no_marked_clauses());
    #endif

    find_xors_based_on_short_clauses();
    find_xors_based_on_long_clauses();

    //Cleanup
    for(ClOffset offset: occsimplifier->clauses) {
        Clause* cl = solver->cl_alloc.ptr(offset);
        cl->stats.marked_clause = false;
    }
    solver->clean_occur_from_idx_types_only_smudged();

    //Print stats
    const bool time_out = (xor_find_time_limit < 0);
    const double time_remain = float_div(xor_find_time_limit, orig_xor_find_time_limit);
    runStats.findTime = cpuTime() - myTime;
    runStats.time_outs += time_out;
    assert(runStats.foundXors == xors.size());
    solver->sumStats.num_xors_found_last = xors.size();
    print_found_xors();

    if (solver->conf.verbosity >= 1) {
        runStats.print_short(solver);
    }
    globalStats += runStats;

    if (solver->sqlStats) {
        solver->sqlStats->time_passed(
            solver
            , "xorfind"
            , cpuTime() - myTime
            , time_out
            , time_remain
        );
    }
}

void XorFinder::print_found_xors()
{
    if (solver->conf.verbosity >= 5) {
        cout << "c Found XORs: " << endl;
        for(vector<Xor>::const_iterator
            it = xors.begin(), end = xors.end()
            ; it != end
            ; ++it
        ) {
            cout << "c " << *it << endl;
        }
    }
}

void XorFinder::add_xors_to_gauss()
{
    solver->xorclauses = xors;
    for(ClOffset offs: cls_of_xors) {
        Clause* cl = solver->cl_alloc.ptr(offs);
        cl->set_used_in_xor(true);
    }
}

bool XorFinder::xor_clause_already_inside(const Xor& xor_c)
{
    xor_find_time_limit -= 30;
    for (const Watched ws: solver->watches.at(xor_c[0])) {
        if (ws.isIdx()
            && xors[ws.get_idx()] == xor_c
        ) {
            return true;
        }
    }
    return false;
}

void XorFinder::findXor(vector<Lit>& lits, const ClOffset offset, cl_abst_type abst)
{
    //Set this clause as the base for the XOR, fill 'seen'
    xor_find_time_limit -= 50;
    poss_xor.setup(lits, offset, abst, seen);

    bool found_all = false;
    //To be honest, -1 would do the trick fully. Or even just "i < 1", but
    //that might not find all possibilities when shorter ones are allowed for
    for (size_t i = 0; i < lits.size()-2; i++) {
        const Lit lit = lits[i];
        findXorMatch(solver->watches[lit], lit);
        findXorMatch(solver->watches[~lit], ~lit);

        //More expensive
        //findXorMatchExt(solver->watches[lit], lit);
        //findXorMatchExt(solver->watches[~lit], ~lit);

        xor_find_time_limit -= 20;
        if (poss_xor.foundAll()) {
            found_all = true;
            break;
        }
    }

    if (found_all) {
        std::sort(lits.begin(), lits.end());
        Xor found_xor(lits, poss_xor.getRHS());

        if (!xor_clause_already_inside(found_xor)) {
            add_found_xor(found_xor);
        }
        for(ClOffset off: poss_xor.get_offsets()) {
            cls_of_xors.push_back(off);
        }
    }

    //Clear 'seen'
    for (const Lit tmp_lit: lits) {
        seen[tmp_lit.var()] = 0;
    }
}

void XorFinder::add_found_xor(const Xor& found_xor)
{
    xor_find_time_limit -= 20;
    xors.push_back(found_xor);
    runStats.foundXors++;
    runStats.sumSizeXors += found_xor.size();
    uint32_t thisXorIndex = xors.size()-1;
    Lit attach_point = Lit(found_xor[0], false);
    solver->watches[attach_point].push(Watched(thisXorIndex));
    solver->watches.smudge(attach_point);
}

void XorFinder::findXorMatchExt(
    watch_subarray_const occ
    , Lit lit
) {
    //seen2 is clear

    for (watch_subarray::const_iterator
        it = occ.begin(), end = occ.end()
        ; it != end
        ; ++it
    ) {
        //TODO we should do the same as below for teritary
        if (!it->isClause()) {
            continue;
        }

        assert(it->isClause() && "This algo has not been updated to deal with TRI, sorry");

        //Deal with clause
        const ClOffset offset = it->get_offset();
        Clause& cl = *solver->cl_alloc.ptr(offset);
        if (cl.freed() || cl.getRemoved())
            continue;

        //Must not be larger than the original clauses
        if (cl.size() > poss_xor.getSize())
            continue;

        tmpClause.clear();
        //cout << "Orig clause: " << poss_xor.getOrigCl() << endl;

        bool rhs = true;
        uint32_t i = 0;
        for (const Lit *l = cl.begin(), *end = cl.end(); l != end; l++, i++) {

            //If this literal is not meant to be inside the XOR
            //then try to find a replacement for it from the cache
            if (!seen[l->var()]) {
                bool found = false;
                //TODO stamping
                const vector<LitExtra>& cache = solver->implCache[~lit].lits;
                for(vector<LitExtra>::const_iterator it2 = cache.begin(), end2 = cache.end()
                    ; it2 != end2 && !found
                    ; it2++
                ) {
                    if (seen[l->var()] && !seen2[l->var()]) {
                        found = true;
                        seen2[l->var()] = true;
                        rhs ^= l->sign();
                        tmpClause.push_back(it2->getLit());
                        //cout << "Added trans lit: " << tmpClause.back() << endl;
                    }
                }

                //Didn't find replacement
                if (!found)
                    goto end;
            }
            else
            //Fine, it's inside the orig clause, but we might have already added this lit
            {
                if (!seen2[l->var()]) {
                    seen2[l->var()] = true;
                    rhs ^= l->sign();
                    tmpClause.push_back(*l);
                    //cout << "Added lit: " << tmpClause.back() << endl;
                } else {
                    goto end; //HACK: we don't want both 'lit' and '~lit' end up in the clause
                }
            }
        }

        //either the invertedness has to match, or the size must be smaller
        if (rhs != poss_xor.getRHS() && cl.size() == poss_xor.getSize())
            goto end;

        //If the size of this clause is the same of the base clause, then
        //there is no point in using this clause as a base for another XOR
        //because exactly the same things will be found.
        if (cl.size() == poss_xor.getSize()) {
            cl.stats.marked_clause = true;
        }

        std::sort(tmpClause.begin(), tmpClause.end());
        //Note: cannot add this offset, because XOR won't cover all of the clause
        poss_xor.add(tmpClause, CL_OFFSET_MAX, varsMissing);

        end:;
        //cout << "Not OK" << endl;

        //Clear 'seen2'
        for(const Lit tmp_lit: tmpClause) {
            seen2[tmp_lit.var()] = false;
        }
    }
}

void XorFinder::findXorMatch(
    watch_subarray_const occ
    , const Lit lit
) {
    xor_find_time_limit -= (int64_t)occ.size();
    for (watch_subarray::const_iterator
        it = occ.begin(), end = occ.end()
        ; it != end
        ; ++it
    ) {
        if (it->isIdx() || it->isBin()) {
            continue;
        }

        //Deal with tertiary
        if (it->isTri()) {
            if (//Only once per tri
                lit < it->lit2() && it->lit2() < it->lit3()

                //Only for correct tri
                && seen[it->lit2().var()] && seen[it->lit3().var()]
            ) {
                bool rhs = true;
                rhs ^= lit.sign();
                rhs ^= it->lit2().sign();
                rhs ^= it->lit3().sign();

                if (rhs == poss_xor.getRHS() || poss_xor.getSize() > 3) {
                    tmpClause.clear();
                    tmpClause.push_back(lit);
                    tmpClause.push_back(it->lit2());
                    tmpClause.push_back(it->lit3());

                    xor_find_time_limit-=20;
                    poss_xor.add(tmpClause, CL_OFFSET_MAX, varsMissing);
                    if (poss_xor.foundAll())
                        break;
                }

            }
            continue;
        }

        //Deal with clause

        //Clause will be at least 4 long
        if (poss_xor.getSize() <= 3) {
            continue;
        }

        const ClOffset offset = it->get_offset();
        Clause& cl = *solver->cl_alloc.ptr(offset);
        xor_find_time_limit -= 20; //deref penalty
        if (cl.freed() || cl.getRemoved())
            continue;

        //Must not be larger than the original clause
        if (cl.size() > poss_xor.getSize())
            continue;

        //Doesn't contain variables not in the original clause
        if ((cl.abst | poss_xor.getAbst()) != poss_xor.getAbst())
            continue;

        //Check RHS, vars inside
        bool rhs = true;
        for (const Lit cl_lit :cl) {
            //early-abort, contains literals not in original clause
            if (!seen[cl_lit.var()])
                goto end;

            rhs ^= cl_lit.sign();
        }
        //either the invertedness has to match, or the size must be smaller
        if (rhs != poss_xor.getRHS() && cl.size() == poss_xor.getSize())
            continue;

        //If the size of this clause is the same of the base clause, then
        //there is no point in using this clause as a base for another XOR
        //because exactly the same things will be found.
        if (cl.size() == poss_xor.getSize()) {
            cl.stats.marked_clause = true;;
        }

        xor_find_time_limit-=50;
        poss_xor.add(cl, offset, varsMissing);
        if (poss_xor.foundAll())
            break;

        end:;
    }
}

void XorFinder::clean_up_xors()
{
    assert(toClear.empty());

    //Fill seen with vars used
    for(const Xor& x: xors) {
        for(uint32_t v: x) {
            if (seen[v] == 0) {
                toClear.push_back(Lit(v, false));
            }

            if (seen[v] < 2) {
                seen[v]++;
            }
        }
    }

    vector<Xor>::iterator it = xors.begin();
    vector<Xor>::iterator it2 = xors.begin();
    size_t num = 0;
    for(vector<Xor>::iterator end = xors.end()
        ; it != end
        ; it++
    ) {
        if (xor_has_interesting_var(*it)) {
            *it2 = *it;
            it2++;
            num++;
        }
    }
    xors.resize(num);

    for(Lit l: toClear) {
        seen[l.var()] = 0;
    }
    toClear.clear();
}

void XorFinder::recursively_xor_xors()
{
    assert(toClear.empty());
    //count how many times a var is used
    for(const Xor& x: xors) {
        for(uint32_t v: x) {
            if (seen[v] == 0) {
                toClear.push_back(Lit(v, false));
            }

            if (seen[v] < 3) {
                seen[v]++;
            }
        }
    }

    //Link in
    for(size_t i = 0; i < xors.size(); i++) {
        const Xor& x = xors[i];
        for(uint32_t v: x) {
            Lit var(v, false);
            assert(solver->watches.size() > var.toInt());
            solver->watches[var].push(Watched(i));
            solver->watches.smudge(var);
        }
    }

    //Only when a var is used exactly twice it's interesting
    vector<uint32_t> interesting;
    for(const Lit l: toClear) {
        if (seen[l.var()] == 2) {
            interesting.push_back(l.var());
        }
    }

    while(!interesting.empty()) {
        const uint32_t v = interesting.back();
        interesting.resize(interesting.size()-1);

        Xor x[2];
        size_t idxes[2];
        unsigned at = 0;
        size_t i2 = 0;
        assert(solver->watches.size() > Lit(v, false).toInt());
        watch_subarray ws = solver->watches[Lit(v, false)];
        for(size_t i = 0; i < ws.size(); i++) {
            const Watched& w = ws[i];
            if (!w.isIdx()) {
                ws[i2] = ws[i];
                i2++;
            } else if (xors[w.get_idx()] != Xor()) {
                assert(at < 2);
                x[at] = xors[w.get_idx()];
                idxes[at] = w.get_idx();
                at++;
            }
        }
        ws.resize(i2);
        if (at < 2) {
            //Has been removed thanks to some XOR-ing together, skip
            continue;
        }

        vector<uint32_t> vars = xor_two(x[0], x[1], idxes[0], idxes[1], v);
        Xor x_new(vars, x[0].rhs ^ x[1].rhs);
        xors.push_back(x_new);
        for(uint32_t v: x_new) {
            Lit var(v, false);
            solver->watches[var].push(Watched(xors.size()-1));
            solver->watches.smudge(var);
        }
        xors[idxes[0]] = Xor();
        xors[idxes[1]] = Xor();
    }

    for(const Lit l: toClear) {
        seen[l.var()] = 0;
    }
    toClear.clear();

    solver->clean_occur_from_idx_types_only_smudged();
    clean_xors_from_empty();
}

void XorFinder::clean_xors_from_empty()
{
    size_t i2 = 0;
    for(size_t i = 0;i < xors.size(); i++) {
        Xor& x = xors[i];
        if (x.size() == 0
            && x.rhs == false
        ) {
            //nothing, skip
        } else {
            xors[i2] = xors[i];
            i2++;
        }
    }
    xors.resize(i2);
}

bool XorFinder::add_new_truths_from_xors()
{
    assert(solver->ok);
    size_t i2 = 0;
    for(size_t i = 0;i < xors.size(); i++) {
        Xor& x = xors[i];
        if (x.size() > 2) {
            xors[i2] = xors[i];
            i2++;
            continue;
        }

        switch(x.size() ) {
            case 0: {
                if (x.rhs == true) {
                    solver->ok = false;
                    return false;
                }
                break;
            }

            case 1: {
                vector<Lit> lits;
                lits.push_back(Lit(x[0], !x.rhs));
                solver->add_clause_int(lits);
                if (!solver->ok) {
                    return false;
                }
                break;
            }

            case 2: {
                vector<Lit> lits{Lit(x[0], false), Lit(x[1], false)};
                solver->add_xor_clause_inter(lits, x.rhs, true);
                if (!solver->ok) {
                    return false;
                }
                break;
            }

            default: {
                assert(false && "Not possible");
            }
        }
    }
    xors.resize(i2);

    return true;
}

void XorFinder::clean_occur_from_idxs(const Lit lit, size_t idx1, size_t idx2)
{
    auto ws = solver->watches[lit];
    size_t i2 = 0;
    for(size_t i = 0; i < ws.size(); i++) {
        if (i != idx1 && i != idx2) {
            ws[i2] = ws[i];
            i2++;
        }
    }
    ws.resize(ws.size()-2);
}

void XorFinder::clean_occur_from_idx(const Lit lit, size_t idx1)
{
    auto ws = solver->watches[lit];
    size_t i2 = 0;
    for(size_t i = 0; i < ws.size(); i++) {
        if (i != idx1) {
            ws[i2] = ws[i];
            i2++;
        }
    }
    ws.resize(ws.size()-1);
}

vector<uint32_t> XorFinder::xor_two(Xor& x1, Xor& x2, const size_t idx1, const size_t idx2, const uint32_t v)
{
    x1.sort();
    x2.sort();
    vector<uint32_t> ret;
    size_t x1_at = 0;
    size_t x2_at = 0;
    while(x1_at < x1.size() || x2_at < x2.size()) {
        if (x1_at == x1.size()) {
            ret.push_back(x2[x2_at]);
            clean_occur_from_idx(Lit(x2[x2_at], false), idx2);
            x2_at++;
            continue;
        }

        if (x2_at == x2.size()) {
            ret.push_back(x1[x1_at]);
            clean_occur_from_idx(Lit(x1[x1_at], false), idx1);
            x1_at++;
            continue;
        }

        const uint32_t a = x1[x1_at];
        const uint32_t b = x2[x2_at];
        if (a == v) {
            x1_at++;
            x2_at++;
            continue;
        }

        if (a == b) {
            x1_at++;
            x2_at++;
            clean_occur_from_idxs(Lit(a, false), idx1, idx2);
            //we could/should update seen[] but in case there are too many XORs
            //we could not store the value in seen[] when counting and then
            //everything would go haywire. So this algorithm is not perfect
            //but is good enough
            continue;
        }

        if (a < b) {
            ret.push_back(a);
            x1_at++;
            continue;
        } else {
            ret.push_back(b);
            x2_at++;
            continue;
        }
    }

    return ret;
}

bool XorFinder::xor_has_interesting_var(const Xor& x)
{
    for(uint32_t v: x) {
        if (seen[v] > 1) {
            return true;
        }
    }
    return false;
}

size_t XorFinder::mem_used() const
{
    size_t mem = 0;
    mem += xors.capacity()*sizeof(Xor);

    //Temporary
    mem += tmpClause.capacity()*sizeof(Lit);
    mem += varsMissing.capacity()*sizeof(uint32_t);

    return mem;
}

void XorFinder::Stats::print_short(const Solver* solver) const
{
    cout
    << "c [occ-xor] found " << std::setw(6) << foundXors
    << " avg sz " << std::setw(4) << std::fixed << std::setprecision(1)
    << float_div(sumSizeXors, foundXors)
    << solver->conf.print_times(findTime, time_outs)
    << endl;
}

void XorFinder::Stats::print() const
{
    cout << "c --------- XOR STATS ----------" << endl;
    print_stats_line("c num XOR found on avg"
        , float_div(foundXors, numCalls)
        , "avg size"
    );

    print_stats_line("c XOR avg size"
        , float_div(sumSizeXors, foundXors)
    );

    print_stats_line("c XOR finding time"
        , findTime
        , float_div(time_outs, numCalls)*100.0
        , "time-out"
    );
    cout << "c --------- XOR STATS END ----------" << endl;
}

XorFinder::Stats& XorFinder::Stats::operator+=(const XorFinder::Stats& other)
{
    //Time
    findTime += other.findTime;

    //XOR
    foundXors += other.foundXors;
    sumSizeXors += other.sumSizeXors;

    //Usefulness
    time_outs += other.time_outs;

    return *this;
}

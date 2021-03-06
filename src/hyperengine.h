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

#include "cnf.h"
#include "propby.h"
#include "solvertypes.h"
#include <vector>
#include <set>
#include "propengine.h"
#include "mystack.h"


using std::vector;
using std::set;

namespace CMSat {
using namespace CMSat;

class HyperEngine : public PropEngine {
public:
    HyperEngine(const SolverConf *_conf, std::atomic<bool>* _must_interrupt_inter);
    size_t print_stamp_mem(size_t totalMem) const;
    size_t mem_used() const;
    size_t mem_used_stamp() const;

    bool use_depth_trick = true;
    bool perform_transitive_reduction = true;
    bool timedOutPropagateFull = false;
    Lit propagate_bfs(
        const uint64_t earlyAborTOut = std::numeric_limits<uint64_t>::max()
    );
    Lit propagate_dfs(
        StampType stampType
        , uint64_t earlyAborTOut = std::numeric_limits<uint64_t>::max()
    );
    set<BinaryClause> needToAddBinClause;       ///<We store here hyper-binary clauses to be added at the end of propagateFull()
    set<BinaryClause> uselessBin;

    ///Add hyper-binary clause given this bin clause
    void  add_hyper_bin(Lit p);

    ///Add hyper-binary clause given this tri-clause
    void  add_hyper_bin(Lit p, Lit lit1, Lit lit2);

    ///Add hyper-binary clause given this large clause
    void  add_hyper_bin(Lit p, const Clause& cl);

    void  enqueue_with_acestor_info(const Lit p, const Lit ancestor, const bool redStep);

private:
    Lit   analyzeFail(PropBy propBy);
    void  close_all_timestamps(const StampType stampType);
    Lit   remove_which_bin_due_to_trans_red(Lit conflict, Lit thisAncestor, const bool thisStepRed);
    void  remove_bin_clause(Lit lit);
    bool  is_ancestor_of(
        const Lit conflict
        , Lit thisAncestor
        , const bool thisStepRed
        , const bool onlyIrred
        , const Lit lookingForAncestor
    );

    //Find lowest common ancestor, once 'currAncestors' has been filled
    Lit deepest_common_ancestor();

    PropResult prop_bin_with_ancestor_info(
        const Lit p
        , watch_subarray::const_iterator k
        , PropBy& confl
    );
    PropResult prop_tri_clause_with_acestor_info(
        watch_subarray_const::const_iterator i
        , const Lit lit1
        , PropBy& confl
    );
    PropResult prop_normal_cl_with_ancestor_info(
        watch_subarray_const::const_iterator i
        , watch_subarray::iterator &j
        , const Lit p
        , PropBy& confl
    );
    PropResult propTriHelperComplex(
        const Lit lit1
        , const Lit lit2
        , const Lit lit3
        , const bool red
    );
    Lit prop_red_bin_dfs(
        StampType stampType
        , PropBy& confl
        , Lit& root
        , bool& restart
    );
    Lit prop_irred_bin_dfs(
       StampType stampType
        , PropBy& confl
        , const Lit root
        , bool& restart
    );
    Lit prop_larger_than_bin_cl_dfs(
        StampType stampType
        , PropBy& confl
        , Lit& root
        , bool& restart
    );
    bool need_early_abort_dfs(
        StampType stampType
        , const size_t timeout
    );

    //For proiorty propagations
    //
    MyStack<Lit> toPropNorm;
    MyStack<Lit> toPropBin;
    MyStack<Lit> toPropRedBin;

    vector<Lit> currAncestors;
};

}
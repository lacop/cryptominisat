/*
 * CryptoMiniSat
 *
 * Copyright (c) 2009-2014, Mate Soos. All rights reserved.
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


#ifndef __REDUCEDB_H__
#define __REDUCEDB_H__

#include "cleaningstats.h"
#include "clauseallocator.h"
#include "clauseusagestats.h"

namespace CMSat {

enum class CleanWhich
{
    normal, longer
};

class Solver;

class ReduceDB
{
public:
    ReduceDB(Solver* solver);
    void reduce_db_and_update_reset_stats(const CleanWhich which);
    const CleaningStats& get_cleaning_stats() const;

    uint64_t get_nbReduceDB() const
    {
        return nbReduceDB;
    }

private:
    Solver* solver;
    uint64_t nbReduceDB = 0;
    vector<ClOffset> delayed_clause_free;
    CleaningStats cleaningStats;

    vector<ClOffset> not_cleaned;
    void move_to_not_cleaned(const CleanWhich which);
    void move_from_not_cleaned();

    unsigned cl_locked;
    unsigned cl_marked;
    unsigned cl_ttl;
    unsigned cl_glue;
    unsigned cl_locked_solver;

    size_t last_reducedb_num_conflicts = 0;
    bool red_cl_too_young(const Clause* cl) const;
    void clear_clauses_stats(vector<ClOffset>& clauseset);

    bool cl_needs_removal(
        const Clause* cl
        , const ClOffset offset
        , const CleanWhich which
    ) const;
    void remove_cl_from_array_and_count_stats(
        CleaningStats& tmpStats
        , uint64_t sumConflicts
        , const CleanWhich which
    );

    CleaningStats reduceDB(const CleanWhich which);

    void sort_red_cls(ClauseCleaningTypes clean_type);
    void mark_top_N_clauses(const uint64_t keep_num);
    void print_best_red_clauses_if_required() const;
    ClauseUsageStats sumClauseData(
        const vector<ClOffset>& toprint
        , bool red
    ) const;
};

inline const CleaningStats& ReduceDB::get_cleaning_stats() const
{
    return cleaningStats;
}

}

#endif //__REDUCEDB_H__

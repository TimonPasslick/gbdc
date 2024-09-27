/*************************************************************************************************
CNFTools -- Copyright (c) 2021, Markus Iser, KIT - Karlsruhe Institute of Technology

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#ifndef ISOHASH2_H_
#define ISOHASH2_H_

#include <algorithm>
#include <list>
#include <string>
#include <vector>

#include "lib/md5/md5.h"

#include "src/util/PointerlessCNFFormula.h"


namespace CNF {
    /**
     * @brief Hashsum of formula with variables transformed to rank in ordered degree sequence of literal incidence graph
     * - literal nodes are grouped pairwise and sorted lexicographically
     * - edge weight of 1/n ==> literal node degree = occurence count
     * - sign is preserved canonically if literal node degrees are equal
     * @param filename benchmark instance
     * @return std::string isohash2
     */
    std::string isohash2(const char* filename) {
        // data structures
        struct PN {
            int p;
            int n;
            void flip() { std::swap(p, n); }
            bool operator < (PN o) const { return n != o.n ? n < o.n : p < o.p; }
            bool operator == (PN o) const { return p == o.p && n == o.n; }
        };
        enum Bool3 : unsigned {
            no = (unsigned) 0 << 30,
            yes = (unsigned) 1 << 30,
            maybe = (unsigned) 1 << 31,
        };
        struct Var {
            PN pn;
            // rank in isometrical order, Bool3 if it's flipped in first 2 bits
            unsigned rank;
        };
        std::vector<Var> vars;
        std::vector<int> order; // var numbers partially ordered by isometry
        struct TieCount {
            unsigned rank;
            unsigned count;
        };
        std::vector<TieCount> var_ties; 
        PointerlessCNFFormula cnf(filename);

        // find canonical polarities (if they exist)
        vars.resize(cnf.variables);
        for (const Lit lit : cnf.literals) {
            Var& var = vars[lit.var() - 1];
            if (lit.sign()) ++var.pn.n;
            else ++var.pn.p;
        }
        for (Var& var : vars) {
            if (var.pn.n > var.pn.p) {
                var.rank = Bool3::yes;
                var.pn.flip();
            } else if (var.pn.p > var.pn.n) var.rank = Bool3::no;
            else var.rank = Bool3::maybe;
        }

        // isohash1 order
        order.resize(vars.size());
        for (int i = 0; i < order.size(); ++i)
            order[i] = i;
        std::sort(order.begin(), order.end(), [&vars](const int a, const int b) {
            return vars[a].pn < vars[b].pn;
        });
        // rank access
        assert(order.size() <= Bool3::yes);
        PN cmp = vars[order[0]].pn;
        unsigned equal = 0u - 1u; // overflow arithmetic is defined for unsigned
        for (unsigned i = 0; i < order.size(); ++i) {
            Var& var = vars[order[i]];
            if (var.pn == cmp) ++equal;
            else {
                if (equal >= 1)
                    var_ties.push_back({i - 1u - equal, equal});
                cmp = var.pn;
                equal = 0;
            }
            var.rank |= i - equal;
        }

        // formula transformation
        // literal transformation
        for (Lit& lit : cnf.literals) {
            const bool sign = lit.sign();
            lit.x = vars[lit.var() - 1].rank;
            if (sign && lit.x & Bool3::maybe == 0) lit.x ^= Bool3::yes;
        }
        // clause sorting
        for (int i = 0; i < cnf.clauses.size(); ++i) {
            const auto slice = cnf.clauses[i];
            std::sort(
                    cnf.literals.begin() + slice.begin,
                    cnf.literals.begin() + slice.end);
        }
        // formula sorting
        std::sort(cnf.clauses.begin(), cnf.clauses.end(), [&cnf](const auto& asl, const auto& bsl) {
            const unsigned asize = asl.end - asl.begin;
            const unsigned bsize = bsl.end - bsl.begin;
            if (asize != bsize) return asize < bsize;
            for (int i = 0; i < asize; ++i) {
                const Lit alit = cnf.literals[asl.begin + i];
                const Lit blit = cnf.literals[bsl.begin + i];
                if (alit != blit) return alit < blit;
            }
            return false;
        });

        // hash
        MD5 md5;
        const auto hash = [&md5](unsigned x) {
            md5.consume(reinterpret_cast<char*>(&x), sizeof(unsigned));
        };
        for (const TieCount var_tie : var_ties) {
            hash(var_tie.rank);
            hash(var_tie.count);
        }
        for (const auto& slice : cnf.clauses) {
            hash(Bool3::maybe | Bool3::yes); // separator
            for (int i = slice.begin; i < slice.end; ++i)
                hash(cnf.literals[i].x);
        }
        return md5.produce();
    }
} // namespace CNF

#endif  // ISOHASH2_H_

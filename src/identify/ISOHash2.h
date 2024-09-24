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
        };
        enum Bool3 : unsigned {
            no = 0b00 << 29,
            yes = 0b01 << 29,
            maybe = 0b10 << 29,
        };
        struct Var {
            PN pn;
            // rank in isometrical order, Bool3 if it's flipped in first 2 bits
            unsigned rank;
        };
        std::vector<Var> vars;
        std::vector<int> order; // var numbers partially ordered by isometry

        // initialization
        PointerlessCNFFormula cnf(filename);
        vars.resize(cnf.nVars());
        for (const auto& clcb : cnf) {
            const Cl& cl = clcb.vector();
            for (int i = 0; i < cl.size(); ++i) {
                const Lit lit = cl[i];
                Var& var = vars[lit.var() - 1];
                if (lit.sign()) ++var.pn.n;
                else ++var.pn.p;
            }
        }
        order.resize(vars.size());
        for (int i = 0; i < order.size(); ++i)
            order[i] = i;

        // find canonical polarities (if they exist)
        for (Var& var : vars) {
            if (var.pn.n > var.pn.p) {
                var.rank |= Bool3::yes;
                var.pn.flip();
            } else if (var.pn.p > var.pn.n) var.rank |= Bool3::no;
            else var.rank |= Bool3::maybe;
        }
        // isohash1 order
        std::sort(order.begin(), order.end(), [&vars](const int a, const int b) {
            return vars[a].pn < vars[b].pn;
        });
        // rank access
        assert(order.size() <= Bool3::yes);
        for (unsigned i = 0; i < order.size(); ++i)
            vars[order[i]].rank |= i;

        // formula transformation
        for (auto& clcb : cnf) {
            Cl& cl = clcb.vector();
            for (Lit& lit : cl) {
                const bool sign = lit.sign();
                lit.x = vars[lit.var() - 1].rank;
                if (sign && lit.x < Bool3::maybe) lit.x ^= Bool3::yes;
            }
            std::sort(cl.begin(), cl.end());
        }
        std::sort(cnf.begin(), cnf.end(), [](const auto& acb, const auto& bcb) {
            const auto& a = acb.vector();
            const auto& b = bcb.vector();
            if (a.size() != b.size()) return a.size() < b.size();
            for (int i = 0; i < a.size(); ++i)
                if (a[i] != b[i]) return a[i] < b[i];
            return false;
        });

        // hash
        MD5 md5;
        const auto hash = [&md5](unsigned x) {
            md5.consume(reinterpret_cast<char*>(&x), sizeof(unsigned));
        };
        for (const auto& clcb : cnf) {
            for (const Lit& lit : clcb.vector())
                hash(lit.x);
            hash((unsigned) 0b11 << 29); // separator
        }
        return md5.produce();
    }
} // namespace CNF

#endif  // ISOHASH2_H_

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

#include "src/util/CNFFormula.h"


enum class Bool3 : signed char {
    yes = 1,
    no = -1,
    maybe = 0
};

namespace CNF {
    /**
     * @brief TODO Timon
     * @param filename benchmark instance
     * @return std::string isohash2
     */
    std::string isohash2(const char* filename) {
        // data structures
        struct PN {
            int p;
            int n;
            void flip() { std::swap(p, n); }
            bool operator < (PN o) { return n != o.n ? n < o.n : p < o.p; }
        };
        struct Var {
            PN pn;
            int rank; // rank in isometrical order
            Bool3 flipped;
        };
        std::vector<Var> vars;
        std::vector<int> order; // var numbers partially ordered by isometry

        // initialization
        CNFFormula cnf(filename);
        vars.resize(cnf.nVars());
        for (Cl* cl : cnf) {
            for (int i = 0; i < cl->size(); ++i) {
                const Lit lit = (*cl)[i];
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
                var.flipped = Bool3::yes;
                var.pn.flip();
            } else if (var.pn.p > var.pn.n) var.flipped = Bool3::no;
            else var.flipped = Bool3::maybe;
        }
        // isohash1 order
        std::sort(order.begin(), order.end(), [&vars](int a, int b) {
            return vars[a].pn < vars[b].pn;
        });
        // rank access
        for (int i = 0; i < order.size(); ++i)
            vars[order[i]].rank = i;

        // TODO Timon: formula transformation

        // hash
        MD5 md5;
        char buffer[64];
        // TODO Timon
        return md5.produce();
    }
} // namespace CNF

namespace WCNF {
    std::string isohash2(const char* filename) {
        // TODO Timon
        // hash
        MD5 md5;
        char buffer[64];
        // TODO Timon
        return md5.produce();
    }
} // namespace CNF

#endif  // ISOHASH2_H_

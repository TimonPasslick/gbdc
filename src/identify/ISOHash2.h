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
        struct Var {
            int n;
            int p;
            void flip() { std::swap(p, n); }
            bool operator < (Var o) const { return n != o.n ? n < o.n : p < o.p; }
            bool operator != (Var o) const { return n != o.n || p != o.p; }
        };
        std::vector<Var> vars;
        std::vector<std::vector<Var>*> normal_form;
        CNFFormula cnf(filename);

        vars.resize(cnf.nVars());
        for (const Cl* cl : cnf) {
            for (const Lit lit : *cl) {
                Var& var = vars[lit.var() - 1];
                if (lit.sign()) ++var.n;
                else ++var.p;
            }
        }

        for (Var& var : vars)
            if (var.n > var.p)
                var.flip();

        // formula transformation
        // literal transformation
        for (const Cl* cl : cnf) {
            normal_form.push_back(new std::vector<Var>);
            for (const Lit lit : *cl) {
                auto& normal_cl = *normal_form.back();
                normal_cl.push_back(vars[lit.var() - 1]);
                if (lit.sign())
                    normal_cl.back().flip(); // flip
            }
        }
        // clause sorting
        for (auto* cl : normal_form)
            std::sort(cl->begin(), cl->end());
        // formula sorting
        std::sort(normal_form.begin(), normal_form.end(), [](const auto* a, const auto* b) {
            if (a->size() != b->size()) return a->size() < b->size();
            for (int i = 0; i < a->size(); ++i)
                if ((*a)[i] != (*b)[i]) return (*a)[i] < (*b)[i];
            return false;
        });

        // hash
        MD5 md5;
        const auto hash = [&md5](int x) {
            md5.consume(reinterpret_cast<char*>(&x), sizeof(int));
        };
        for (const auto* cl : normal_form) {
            for (const Var var : *cl) {
                hash(var.n);
                hash(var.p);
            }
            hash(-1); // separator
        }

        // cleanup
        for (const auto* cl : normal_form)
            delete cl;

        return md5.produce();
    }
} // namespace CNF

#endif  // ISOHASH2_H_

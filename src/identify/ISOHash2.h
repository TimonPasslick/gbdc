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

struct Iv { // interval
    int start;
    int length;
};

namespace CNF {
    /**
     * @brief TODO Timon
     * @param filename benchmark instance
     * @return std::string isohash2
     */
    std::string isohash2(const char* filename) {
        // data structures
        struct Var {
            std::vector<Cl*> pos;
            std::vector<Cl*> neg;
            int rank; // rank in isometrical order
            bool in_interval; // whether the rank is shared
            Bool3 flipped;
            bool normalized() {
                return rank != 0 && flipped != Bool3::maybe;
            }
        };
        std::vector<Var> vars;
        std::vector<int> order; // var numbers that get partially ordered by isometry
        std::list<Iv> unordered;
        std::list<Var*> maybe_flipped;

        // initialization
        CNFFormula cnf(filename);
        vars.resize(cnf.nVars());
        for (Cl* cl : cnf) {
            for (int i = 0; i < cl->size(); ++i) {
                const Lit lit = (*cl)[i];
                Var& var = vars[lit.var() - 1];
                auto& cls = lit.sign() ? var.neg : var.pos;
                cls.push_back(cl);
            }
        }
        order.resize(vars.size());
        for (int i = 0; i < order.size(); ++i)
            order[i] = i;

        { // isohash1 calculation to be at least as good
            for (Var& var : vars) {
                int pos = var.pos.size();
                int neg = var.neg.size();
                if (neg > pos) var.flipped = Bool3::yes;
                else if (pos > neg) var.flipped = Bool3::no;
                else maybe_flipped.push_back(&var);
            }
            struct PN {
                int p;
                int n;
                void flip() { std::swap(p, n); }
                bool operator == (PN o) { return p == o.p && n == o.n; }
                bool operator < (PN o) { return n != o.n ? n < o.n : p < o.p; }
            };
            const auto pn = [&vars](int i) {
                Var& var = vars[i];
                PN pn {(int) var.pos.size(), (int) var.neg.size()};
                if (var.flipped == Bool3::yes) pn.flip();
                return pn;
            };
            std::sort(order.begin(), order.end(), [&pn](int a, int b) { return pn(a) < pn(b); });
            // check order
            Iv iv {0, 0};
            PN start = pn(0);
            for (int i = 0; i < order.size(); ++i) {
                int j = order[i];
                PN here = pn(j);
                if (here == start) ++iv.length;
                else {
                    if (iv.length > 1) {
                        unordered.push_back(iv);
                        for (int k = 0; k < iv.length; ++i)
                            vars[order[i - k]].in_interval = true;
                    }
                    iv = {i, 1};
                    start = here;
                }
                vars[j].rank = iv.start;
            }
            if (iv.length > 1) unordered.push_back(iv);
        }

        { // TODO Timon: fixed depth recursion to improve order and find polarities
            const auto size_count = [](const std::vector<Cl*>& cls) {
                std::vector<int> size_count;
                for (Cl* cl : cls) {
                    int size = cl->size();
                    if (size >= size_count.size())
                        size_count.reserve(size * 2); // logarithmic reallocation
                        size_count.resize(size);
                    ++size_count[size];
                }
                return size_count;
            };
            const auto flagship = [](const std::vector<int>& a, const std::vector<int>& b) {
                int result = a.size() - b.size();
                for (int i = 0; i < a.size(); ++i) {
                    if (result != 0) return result;
                    result = a[i] - b[i];
                }
                return result;
            };
            auto fi = maybe_flipped.begin();
            while (fi != maybe_flipped.end()) {
                Var& var = **fi;
                int cmp = flagship(size_count(var.pos), size_count(var.neg));
                auto prev = fi++;
                if (cmp == 0) continue;
                else if (cmp < 0) var.flipped = Bool3::yes;
                else var.flipped = Bool3::no;
                maybe_flipped.erase(prev);
            }
            auto ii = unordered.begin();
            while (ii != unordered.end()) {
                // TODO Timon
            }
        }

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

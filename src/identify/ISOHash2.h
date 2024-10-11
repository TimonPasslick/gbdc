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
    template <typename Var>
    struct IHData {
        CNFFormula cnf;
        std::vector<Var> vars;
        std::vector<std::vector<Var>*> normal_form;

        IHData(const char* filename) : cnf(filename) {
            // variable normalization
            vars.resize(cnf.nVars());
            for (const Cl* cl : cnf) {
                for (const Lit lit : *cl) {
                    Var& var = vars[lit.var() - 1];
                    if (lit.sign()) ++var.n;
                    else ++var.p;
                }
            }
            // variable flipping
            for (Var& var : vars)
                if (var.n > var.p)
                    var.flip();
            // literal transformation
            for (const Cl* cl : cnf) {
                normal_form.push_back(new std::vector<Var>);
                for (const Lit lit : *cl) {
                    auto& normal_cl = *normal_form.back();
                    normal_cl.push_back(vars[lit.var() - 1]);
                    if (lit.sign())
                        normal_cl.back().flip();
                }
            }
            // clause sorting
            for (auto* cl : normal_form)
                std::sort(cl->begin(), cl->end());
        }

        ~IHData() {
            for (const auto* cl : normal_form)
                delete cl;
        }

        std::string final_hash() {
            // formula sorting
            std::sort(normal_form.begin(), normal_form.end(), [](const auto* a, const auto* b) {
                if (a->size() != b->size()) return a->size() < b->size();
                for (int i = 0; i < a->size(); ++i)
                    if ((*a)[i] != (*b)[i]) return (*a)[i] < (*b)[i];
                return false;
            });
            // normal form hashing
            MD5 md5;
            for (const auto* cl : normal_form) {
                for (const Var var : *cl)
                    md5.consume(reinterpret_cast<const char*>(&var), sizeof(Var));
                constexpr int separator = -1;
                md5.consume(reinterpret_cast<const char*>(&separator), sizeof(int));
            }
            return md5.produce();
        }
    };

    /**
     * @brief Hashsum of formula with variables transformed to rank in ordered degree sequence of literal incidence graph
     * - literal nodes are grouped pairwise and sorted lexicographically
     * - edge weight of 1/n ==> literal node degree = occurence count
     * - sign is preserved canonically if literal node degrees are equal
     * @param filename benchmark instance
     * @return std::string isohash2
     */
    std::string isohash2(const char* filename) {
        struct Var {
            int n;
            int p;
            void flip() { std::swap(p, n); }
            bool operator < (Var o) const { return n != o.n ? n < o.n : p < o.p; }
            bool operator != (Var o) const { return n != o.n || p != o.p; }
        };
        return IHData<Var>(filename).final_hash();
    }

    std::string isohash3plus(const unsigned n, const char* filename) {
        struct SigVar {
            MD5::Signature n;
            MD5::Signature p;
            void flip() { std::swap(p, n); }
            bool operator < (SigVar o) const { return n != o.n ? n < o.n : p < o.p; }
            bool operator != (SigVar o) const { return n != o.n || p != o.p; }
        };
        IHData<SigVar> data(filename);
        auto& cnf = data.cnf;
        auto& vars = data.vars;
        auto& normal_form = data.normal_form;

        for (int i = 0; i < n + 1; ++i) {
            // variable normalization
            for (int i = 0; i < normal_form.size(); ++i) {
                // clause hashing
                MD5 md5;
                for (const SigVar var : *normal_form[i])
                    md5.consume(reinterpret_cast<const char*>(&var), sizeof(SigVar));
                const auto sig = md5.finish();

                // variable hashing
                for (const Lit lit : *cnf[i]) {
                    SigVar& var = vars[lit.var() - 1];
                    if (lit.sign()) var.n += sig;
                    else var.p += sig;
                }
            }
            // variable flipping
            for (SigVar& var : vars)
                if (var.n > var.p)
                    var.flip();
            // literal transformation
            for (int i = 0; i < cnf.nClauses(); ++i) {
                for (int j = 0; j < cnf[i]->size(); ++i) {
                    const Lit lit = (*cnf[i])[j];
                    SigVar& var = (*normal_form[i])[j];
                    var = vars[lit.var() - 1];
                    if (lit.sign())
                        var.flip();
                }
            }
            // clause sorting
            for (auto* cl : normal_form)
                std::sort(cl->begin(), cl->end());
        }
        return data.final_hash();
    }
} // namespace CNF

#endif  // ISOHASH2_H_

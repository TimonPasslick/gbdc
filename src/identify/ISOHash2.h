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
#include <map>
#include <set>
#include <string>
#include <vector>

#include "lib/md5/md5.h"

#include "src/util/CNFFormula.h"


namespace CNF {
    template <typename T>
    struct Var {
        T n;
        T p;
        void flip() { std::swap(p, n); }
        bool operator < (Var o) const { return n != o.n ? n < o.n : p < o.p; }
        bool operator != (Var o) const { return n != o.n || p != o.p; }
    };

    template <typename T> // needs to be flat, no pointers or heap data
    MD5::Signature sort_and_hash(std::vector<T>* v) {
        std::sort(v->begin(), v->end());
        MD5 md5;
        for (const T t : *v)
            md5.consume_binary(t);
        return md5.finish();
    }

    template <typename T>
    struct IHData {
        const CNFFormula cnf;
        std::vector<Var<T>> vars;
        std::vector<std::vector<Var<T>>*> normal_form;

        IHData(const char* filename) : cnf(filename) {
            vars.resize(cnf.nVars());
            normal_form.resize(cnf.nClauses());
            for (int i = 0; i < cnf.nClauses(); ++i)
                normal_form[i]->resize(cnf[i]->size());

            iteration_step([](const int) { return 1; });
        }

        void iteration_step(const std::function<T(const int i)>& summand) {
            // variable normalization
            for (int i = 0; i < cnf.nClauses(); ++i) {
                const auto s = summand(i);
                for (const Lit lit : *cnf[i]) {
                    Var<T>& var = vars[lit.var() - 1];
                    if (lit.sign()) var.n += s;
                    else var.p += s;
                }
            }
            for (Var<T>& var : vars)
                if (var.n > var.p)
                    var.flip();
            // literal transformation
            for (int i = 0; i < cnf.nClauses(); ++i) {
                for (int j = 0; j < cnf[i]->size(); ++i) {
                    const Lit lit = (*cnf[i])[j];
                    auto& var = (*normal_form[i])[j];
                    var = vars[lit.var() - 1];
                    if (lit.sign())
                        var.flip();
                }
            }
        }

        ~IHData() {
            for (const auto* cl : normal_form)
                delete cl;
        }
        std::string final_hash() {
            MD5::Signature result;
            for (auto* cl : normal_form)
                result += sort_and_hash(cl);

            char str[MD5_STRING_SIZE];
            md5::sig_to_string(result.data, str, sizeof(str));
            return std::string(str);
        }
    };

    /**
     * @brief Hashsum of literal interaction graph with variables transformed to degrees of literal incidence graph
     * - literal nodes are grouped pairwise and the pair is sorted lexicographically
     * - edge weight of 1/n ==> literal node degree = occurence count
     * - sign is preserved canonically if literal node degrees are unequal
     * @param filename benchmark instance
     * @return std::string isohash3
     */
    std::string isohash2(const char* filename) {
        // TODO
        return std::string();
    }

    /**
     * @brief Hashsum of formula with variables transformed to degrees of literal incidence graph
     * - literal nodes are grouped pairwise and the pair is sorted lexicographically
     * - edge weight of 1/n ==> literal node degree = occurence count
     * - sign is preserved canonically if literal node degrees are unequal
     * @param filename benchmark instance
     * @return std::string isohash3
     */
    std::string isohash3(const char* filename) {
        return IHData<int>(filename).final_hash();
    }

    /**
     * @brief Hashsum of formula with variables transformed to degrees of literal incidence graph, then occurence clause hashes n+1 times
     * First normal form:
     * - literal nodes are grouped pairwise and the pair is sorted lexicographically
     * - edge weight of 1/n ==> literal node degree = occurence count
     * - sign is preserved canonically if literal node degrees are unequal
     * Variable hashes aren't cleared between normal form iterations, but they are always canonical, so it's not necessary.
     * @param filename benchmark instance
     * @return std::string isohash3
     */
    std::string isohash4plus(const unsigned n, const char* filename) {
        IHData<MD5::Signature> data(filename);
        auto& normal_form = data.normal_form;
        for (int i = 0; i < n + 1; ++i)
            data.iteration_step([&normal_form](const int clause_index) {
                return sort_and_hash(normal_form[clause_index]);
            });
        return data.final_hash();
    }
} // namespace CNF

#endif  // ISOHASH2_H_

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
    MD5::Signature hash(const T t) {
        MD5 md5;
        md5.consume_binary(t);
        return md5.finish();
    }
    template <typename T> // needs to be flat, no pointers or heap data
    MD5::Signature sort_and_hash(std::vector<T>* v) {
        std::sort(v->begin(), v->end());
        MD5 md5;
        for (const T t : *v)
            md5.consume_binary(t);
        return md5.finish();
    }

    template <typename T>
    struct WLData {
        const CNFFormula cnf;
        std::vector<Var<T>> vars;
        std::vector<std::vector<Var<T>>*> normal_form;

        WLData(const char* filename) : cnf(filename) {
            vars.resize(cnf.nVars());
            normal_form.resize(cnf.nClauses());
            for (unsigned i = 0; i < cnf.nClauses(); ++i)
                normal_form[i]->resize(cnf[i]->size());

            normalize_variables([]() {}, [](const unsigned) { return 1; });
        }
        ~WLData() {
            for (const auto* cl : normal_form)
                delete cl;
        }
        void iteration_step(const std::function<void()>& preprocess_vars, const std::function<T(const unsigned i)>& summand) {
            transform_literals();
            normalize_variables(preprocess_vars, summand);
        }
        void transform_literals() {
            for (unsigned i = 0; i < cnf.nClauses(); ++i) {
                for (unsigned j = 0; j < cnf[i]->size(); ++i) {
                    const Lit lit = (*cnf[i])[j];
                    auto& var = (*normal_form[i])[j];
                    var = vars[lit.var() - 1];
                    if (lit.sign())
                        var.flip();
                }
            }
        }
        void normalize_variables(const std::function<void()>& preprocess_vars, const std::function<T(const unsigned i)>& summand) {
            preprocess_vars();
            for (unsigned i = 0; i < cnf.nClauses(); ++i) {
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
        }
        std::string final_hash() {
            std::sort(vars.begin(), vars.end());
            MD5 md5;
            for (const Var<T> var : vars)
                md5.consume_binary(var);
            return md5.produce();
        }
        std::string half_iteration_final_hash() {
            transform_literals();
            MD5::Signature result;
            for (auto* cl : normal_form)
                result += sort_and_hash(cl);

            char str[MD5_STRING_SIZE];
            md5::sig_to_string(result.data, str, sizeof(str));
            return std::string(str);
        }
    };

    /**
     * @brief Comparing Weisfeiler-Leman hashes with depth 2*h is approximately as strong as
     * running the Weisfeiler-Leman algorithm on the literal hypergraph and stopping after h iterations.
     * Runtime O(depth*nlogn), space O(n).
     * Note that depth 0 is not exactly isohash because this algorithm hashes binary data directly and not strings.
     * @param filename benchmark instance
     * @return std::string Weisfeiler-Leman hash
     */
    std::string weisfeiler_leman_hash(const unsigned depth, const char* filename) {
        if (depth == 0) return WLData<int>(filename).final_hash();
        if (depth == 1) return WLData<int>(filename).half_iteration_final_hash();

        WLData<MD5::Signature> data(filename);

        auto& normal_form = data.normal_form;
        const auto summand = [&normal_form](const unsigned clause_index) {
            return sort_and_hash(normal_form[clause_index]);
        };

        auto& vars = data.vars;
        const auto preprocess_vars = [&vars](){
            for (Var<MD5::Signature>& var : vars) {
                // double hash to avoid ambiguity with unit clauses
                const MD5::Signature previous_color_hash = hash(hash(var));
                var = {previous_color_hash, previous_color_hash};
            }
        };

        for (unsigned i = 0; i < depth / 2; ++i)
            data.iteration_step(preprocess_vars, summand);

        if (depth % 2 == 0) return data.final_hash();
        return data.half_iteration_final_hash();
    }

    /**
     * @brief Like weisfeiler_leman_hash, but with the literal interaction graph.
     * @param filename benchmark instance
     * @return std::string isohash2
     */
    std::string isohash2(const char* filename) {
        // TODO
        // unclear what I should do because first iteration of Weisfeiler-Leman on LIG is not isohash1, so it is not necessarily stronger
        return std::string();
    }
} // namespace CNF

#endif  // ISOHASH2_H_

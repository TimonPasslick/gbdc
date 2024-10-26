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
#include <functional>
#include <vector>

#include "src/external/md5/md5.h"

#include "src/util/PointerlessCNFFormula.h"


namespace CNF {
    struct WeisfeilerLemanHasher {
        const PointerlessCNFFormula cnf;
        using Hash = MD5::Signature;
        using Clause = PointerlessCNFFormula::Clause;
        struct LitColors {
            Hash p;
            Hash n;
            void flip() { std::swap(n, p); }
            void cross_reference() {
                const Hash ncr = hash(*this);
                flip();
                const Hash pcr = hash(*this);
                n = ncr;
                p = pcr;
            }
            Hash variable_hash() {
                if (n > p) flip();
                return hash(*this);
            }
        };
        struct ColorFunction {
            std::vector<LitColors> colors;
            explicit ColorFunction(const std::size_t n) : colors(n, LitColors {1, 1}) {}
            Hash& operator () (const Lit lit) { return reinterpret_cast<Hash*>(&colors[0])[lit]; }
        };
        // old and new color function, swapping in each iteration
        ColorFunction color_functions[2];
        unsigned iteration = 0;

        template <typename T> // needs to be flat, no pointers or heap data
        static Hash hash(const T t) {
            MD5 md5;
            md5.consume_binary(t);
            return md5.finish();
        }
        template <typename T, typename C>
        static Hash hash_sum(const C& c, const std::function<Hash(const T&)> f) {
            Hash h;
            for (const T& t : c)
                h += f(t);
            return h;
        }
        ColorFunction& old_color() { return color_functions[iteration % 2]; }
        ColorFunction& new_color() { return color_functions[(iteration + 1) % 2]; }

        WeisfeilerLemanHasher(const char* filename)
                : cnf(filename)
                , color_functions {ColorFunction(cnf.nVars()), ColorFunction(cnf.nVars())}
        {}
        Hash clause_hash(const Clause cl) {
            // hash again to preserve clause structure and avoid collisions of unit clauses with old color
            return hash(hash_sum<Lit>(cl, old_color()));
        }
        void iteration_step() {
            if (iteration > 0) {
                for (unsigned i = 0; i < cnf.nVars(); ++i) {
                    LitColors& olc = old_color().colors[i];
                    olc.cross_reference();
                    LitColors& nlc = new_color().colors[i];
                    nlc = olc;
                }
            }
            for (const Clause cl : cnf.clauses()) {
                const Hash clh = clause_hash(cl);
                for (const Lit lit : cl)
                    new_color()(lit) += clh;
            }
            ++iteration;
        }
        Hash variable_hash() {
            return hash_sum<LitColors>(old_color().colors, [](LitColors lc) { return lc.variable_hash(); });
        }
        Hash cnf_hash() {
            for (LitColors& lc : old_color().colors)
                lc.cross_reference();
            return hash_sum<Clause>(cnf.clauses(), [this](const Clause cl) { return clause_hash(cl); });
        }
        std::string operator () (const unsigned depth) {
            while (iteration < depth / 2)
                iteration_step();
            const Hash h = depth % 2 == 0 ? variable_hash() : cnf_hash();
            return std::to_string(h);
        }
    };

    /**
     * @brief Comparing Weisfeiler-Leman hashes with depth 2h is approximately as strong as
     * running the Weisfeiler-Leman algorithm on the literal hypergraph and stopping after h iterations.
     * Runtime O(h*n), space O(n).
     * Note that zero iterations do not conform to isohash directly because this algorithm uses a different hashing technique to achieve linear runtime.
     * @param filename benchmark instance
     * @return std::string Weisfeiler-Leman hash
     */
    std::string weisfeiler_leman_hash(const unsigned depth, const char* filename) {
        return WeisfeilerLemanHasher(filename)(depth);
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

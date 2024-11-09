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
#include <optional>
#include <unordered_set>
#include <vector>

#define XXH_INLINE_ALL
#include "xxhash.h"

#include "src/util/PointerlessCNFFormula.h"

namespace CNF {
    struct WeisfeilerLemanHasher {
        using Hash = XXH64_hash_t;
        constexpr static bool debug_output = false;
        const std::string file; // just for debugging
        const PointerlessCNFFormula cnf;
        using Clause = PointerlessCNFFormula::Clause;
        struct LitColors {
            Hash p = 0;
            Hash n = 0;
            inline void flip() { std::swap(n, p); }
            inline void cross_reference() {
                const Hash pcr = hash(*this);
                flip();
                const Hash ncr = hash(*this);
                p = pcr;
                n = ncr;
            }
            inline Hash variable_hash() const {
                LitColors copy = *this;
                if (n > p) copy.flip();
                return hash(copy);
            }
        };
        struct ColorFunction {
            std::vector<LitColors> colors;
            explicit ColorFunction(const std::size_t n) : colors(n) {}
            inline Hash& operator () (const Lit lit) { return reinterpret_cast<Hash*>(&colors[0])[lit]; }
        };
        // old and new color function, swapping in each iteration
        ColorFunction color_functions[2];
        unsigned iteration = 0;
        std::unordered_set<Hash> unique_hashes;
        unsigned previous_unique_hashes = 1;

        inline ColorFunction& old_color() { return color_functions[iteration % 2]; }
        inline ColorFunction& new_color() { return color_functions[(iteration + 1) % 2]; }

        template <typename T> // needs to be flat, no pointers or heap data
        static inline Hash hash(const T t) {
            return {XXH3_64bits(&t, sizeof(T))};
        }
        static inline void combine(Hash* acc, const Hash in) {
            *acc += in + (*acc > std::numeric_limits<Hash>::max() - in);
        }
        template <typename T, typename C>
        static inline Hash hash_sum(const C& c, const std::function<Hash(const T&)>& f) {
            Hash h = 0;
            for (const T& t : c)
                combine(&h, f(t));
            return h;
        }

        WeisfeilerLemanHasher(const char* filename)
                : cnf(filename)
                , color_functions {ColorFunction(cnf.nVars()), ColorFunction(cnf.nVars())}
                , file(filename)
        {
            // optimized first iteration (theoretically all colors are 1 initially)
            for (const Clause cl : cnf.clauses()) {
                const Hash clh = hash((unsigned) cl.size());
                for (const Lit lit : cl)
                    combine(&new_color()(lit), clh);
            }
            ++iteration;
        }
        Hash clause_hash(const Clause cl) {
            // hash again to preserve clause structure and avoid collisions of unit clauses with old color
            return hash(hash_sum<const Lit>(cl, [this](const Lit lit) { return old_color()(lit); }));
        }
        void iteration_step() {
            for (unsigned i = 0; i < cnf.nVars(); ++i) {
                LitColors& olc = old_color().colors[i];
                olc.cross_reference();
                LitColors& nlc = new_color().colors[i];
                nlc = olc;
            }
            for (const Clause cl : cnf.clauses()) {
                const Hash clh = clause_hash(cl);
                for (const Lit lit : cl)
                    combine(&new_color()(lit), clh);
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
        std::optional<std::string> check_progress() {
            // few hits at the start
            if (iteration <= 2) return std::nullopt;

            unique_hashes.reserve(previous_unique_hashes);
            const Hash vh = hash_sum<LitColors>(old_color().colors, [this](LitColors lc) {
                const Hash vh = lc.variable_hash();
                unique_hashes.insert(vh);
                return vh;
            });
            if (unique_hashes.size() <= previous_unique_hashes) {
                if constexpr (debug_output) std::cout << iteration << " iterations for " << file << std::endl;
                return std::to_string(vh);
            }
            previous_unique_hashes = unique_hashes.size();
            unique_hashes.clear();
            return std::nullopt;
        }
        std::string operator () (const unsigned depth = std::numeric_limits<unsigned>::max()) {
            while (iteration <= depth / 2) {
                if (const auto result = check_progress())
                    return *result;
                iteration_step();
            }
            if constexpr (debug_output) std::cout << "iteration limit (" << ((double) depth) / 2 + 1 << ") reached for " << file << std::endl;
            const Hash h = depth % 2 == 0 ? variable_hash() : cnf_hash();
            return std::to_string(h);
        }
    };

    /**
     * @brief Comparing Weisfeiler-Leman hashes with depth 2h is approximately as strong as
     * running the Weisfeiler-Leman algorithm on the literal hypergraph and stopping after h+1 iterations.
     * Runtime O(h*n), space O(n).
     * Note that zero iterations do not conform to isohash directly because this algorithm uses a different hashing technique to achieve linear runtime.
     * @param filename benchmark instance
     * @return std::string Weisfeiler-Leman hash
     */
    std::string weisfeiler_leman_hash(const unsigned depth, const char* filename) {
        return WeisfeilerLemanHasher(filename)(depth);
    }
} // namespace CNF

#endif  // ISOHASH2_H_

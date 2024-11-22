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
#include <chrono>
#include <functional>
#include <optional>
#include <unordered_set>
#include <vector>

#define XXH_INLINE_ALL
#include "xxhash.h"

#include "src/util/PointerlessCNFFormula.h"

namespace CNF {
    struct WLHRuntimeConfig {
        const unsigned depth;
        const bool cross_reference_literals;
        const unsigned first_progress_check_iteration;

        const bool return_iteration_count;
        const bool return_time;
    };
    template < // compile time config
        bool use_xxh3 = true, // MD5 otherwise
        unsigned hash_size = 64,
        bool carrying_addition = true,
        unsigned ring_prime_offset = 0, // 0 means no modulo calculations

        unsigned formula_optimization_level = 2, // 0 means classic, 1 bounds array, 2 size concatenation

        bool optimize_first_iteration = true
    >
    struct WeisfeilerLemanHasher {
        const WLHRuntimeConfig cfg;
        using Hash = XXH64_hash_t;
        const std::string file; // just for debugging
        using Clock = std::chrono::high_resolution_clock;
        Clock::time_point start_time;
        Clock::time_point parsing_start_time;
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
            explicit ColorFunction(const std::size_t n) : colors(n, {1, 1}) {}
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

        WeisfeilerLemanHasher(const char* filename, const WLHRuntimeConfig cfg)
                : file(filename)
                , cfg(cfg)
                , parsing_start_time(Clock::now())
                , cnf(filename)
                , start_time(Clock::now())
                , color_functions {ColorFunction(cnf.nVars()), ColorFunction(cnf.nVars())}
        {
        }
        void cross_reference() {
            for (LitColors& lc : old_color().colors)
                lc.cross_reference();
        }
        Hash clause_hash(const Clause cl) {
            // hash again to preserve clause structure
            return hash(hash_sum<const Lit>(cl, [this](const Lit lit) { return old_color()(lit); }));
        };
        void iteration_step() {
            if (iteration != 0) cross_reference();
            for (const Clause cl : cnf.clauses()) {
                const Hash clh = (iteration != 0) ?
                    clause_hash(cl)
                    : hash((unsigned) cl.size());
                for (const Lit lit : cl)
                    combine(&new_color()(lit), clh);
            }
            ++iteration;
        }
        Hash variable_hash() {
            return hash_sum<LitColors>(old_color().colors, [](LitColors lc) { return lc.variable_hash(); });
        }
        Hash cnf_hash() {
            cross_reference();
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
            if (unique_hashes.size() <= previous_unique_hashes)
                return std::to_string(vh);
            previous_unique_hashes = unique_hashes.size();
            unique_hashes.clear();
            return std::nullopt;
        }
        std::string operator () () {
            while (iteration < cfg.depth / 2) {
                if (const auto result = check_progress())
                    return *result;
                iteration_step();
            }
            const Hash h = cfg.depth % 2 == 0 ? variable_hash() : cnf_hash();
            return std::to_string(h);
        }
    };

    /**
     * @brief Comparing Weisfeiler-Leman hashes is approximately as strong as
     * running the Weisfeiler-Leman algorithm on the literal hypergraph.
     * Runtime O(h*n), space O(n).
     * @param filename benchmark instance
     * @param depth maximum iterations / 2, half iterations hash clause labels
     * @param cross_reference_literals whether the information which literals
     * belong to the same variable should be used in the calculation
     * @param first_progress_check_iteration the first iteration in which the
     * progress check runs
     * @param return_iteration_count whether the amount of iterations that were
     * calculated (possibly half) should be returned
     * @param return_time whether the calculation time and parsing time should
     * be returned
     * @return comma separated list, std::string Weisfeiler-Leman hash,
     * possibly iteration count, calculation time and parsing time
     */
    std::string weisfeiler_leman_hash(
        const char* filename,

        const unsigned depth = 13,
        const bool cross_reference_literals = true,
        const unsigned first_progress_check_iteration = 3,

        const bool return_iteration_count = true,
        const bool return_time = true
    ) {
        WeisfeilerLemanHasher hasher(filename, WLHRuntimeConfig {
            depth,
            cross_reference_literals,
            first_progress_check_iteration,
            return_iteration_count,
            return_time
        });
        std::string result = hasher();
        const auto elapsed = WeisfeilerLemanHasher<>::Clock::now() - hasher.start_time;
        const auto parsing_time = hasher.start_time - hasher.parsing_start_time;
        result += "," + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count())
                + "," + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(parsing_time).count());
        return result;
    }
} // namespace CNF

#endif  // ISOHASH2_H_

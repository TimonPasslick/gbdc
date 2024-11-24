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
#include <type_traits>
#include <unordered_set>
#include <vector>

#include "src/external/md5/md5.h"
#define XXH_INLINE_ALL
#include "xxhash.h"

#include "src/util/CNFFormula.h"
#include "src/util/IntervalCNFFormula.h"
#include "src/util/SizeGroupedCNFFormula.h"

namespace CNF {
    struct WLHRuntimeConfig {
        unsigned depth;
        bool cross_reference_literals;
        bool rehash_clauses;
        bool optimize_first_iteration;
        unsigned first_progress_check_iteration;
        bool return_measurements;
    };
    template < // compile time config
        typename CNF = SizeGroupedCNFFormula,
        bool use_xxh3 = true, // MD5 otherwise
        unsigned hash_size = 64, // can be either 32 or 64
        unsigned ring_prime_offset = 0 // 0 means no modulo calculations
    >
    struct WeisfeilerLemanHasher {
        const WLHRuntimeConfig cfg;
        using Hash = std::conditional_t<hash_size == 32, std::uint32_t, std::uint64_t>;
        using Clock = std::chrono::high_resolution_clock;
        Clock::time_point start_time;
        Clock::time_point parsing_start_time;
        const CNF cnf;
        using Clause = typename CNF::Clause;
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
            if (use_xxh3)
                return XXH3_64bits(&t, sizeof(T));

            MD5 md5;
            md5.consume_binary<T>(t);
            return md5.finish();
        }
        static inline void combine(Hash* acc, Hash in) {
            if constexpr (ring_prime_offset > 0) {
                constexpr Hash ring_size = ((Hash) 0) - ring_prime_offset;
                in %= ring_size;
                const Hash first_overflow_acc = ring_size - in;
                if (*acc >= first_overflow_acc) {
                    *acc -= first_overflow_acc;
                    return;
                }
            }
            *acc += in;
        }
        template <typename T, typename C>
        static inline Hash hash_sum(const C& c, const std::function<Hash(const T&)>& f) {
            Hash h = 0;
            for (const T& t : c)
                combine(&h, f(t));
            return h;
        }

        WeisfeilerLemanHasher(const char* filename, const WLHRuntimeConfig cfg)
                : cfg(cfg)
                , parsing_start_time(Clock::now())
                , cnf(filename)
                , start_time(Clock::now())
                , color_functions {ColorFunction(cnf.nVars()), ColorFunction(cnf.nVars())}
        {
        }
        inline bool in_optimized_iteration() {
            return iteration == 0 && cfg.optimize_first_iteration;
        }
        void cross_reference() {
            if (!cfg.cross_reference_literals || in_optimized_iteration())
                return;
            for (LitColors& lc : old_color().colors)
                lc.cross_reference();
        }
        Hash clause_hash(const Clause cl) {
            Hash h = hash_sum<const Lit>(cl, [this](const Lit lit) { return old_color()(lit); });
            // hash again to preserve clause structure
            if (cfg.rehash_clauses) h = hash(h);
            return h;
        };
        void iteration_step() {
            cross_reference();
            for (const Clause cl : cnf.clauses()) {
                const Hash clh = (!in_optimized_iteration()) ?
                    clause_hash(cl)
                    : hash((unsigned) cl.size());
                for (const Lit lit : cl)
                    combine(&new_color()(lit), clh);
            }
            ++iteration;
        }
        Hash variable_hash() {
            if (cfg.cross_reference_literals)
                return hash_sum<LitColors>(old_color().colors, [](LitColors lc) { return lc.variable_hash(); });

            Hash h = 0;
            for (Lit lit {}; lit != cnf.nVars() * 2; ++lit)
                combine(&h, old_color()(lit));
            return h;
        }
        Hash cnf_hash() {
            cross_reference();
            return hash_sum<Clause>(cnf.clauses(), [this](const Clause cl) { return clause_hash(cl); });
        }
        std::optional<Hash> check_progress() {
            // few hits at the start
            if (iteration < cfg.first_progress_check_iteration) return std::nullopt;

            unique_hashes.reserve(previous_unique_hashes);
            const Hash vh = hash_sum<LitColors>(old_color().colors, [this](LitColors lc) {
                const Hash vh = lc.variable_hash();
                unique_hashes.insert(vh);
                return vh;
            });
            if (unique_hashes.size() <= previous_unique_hashes)
                return vh;
            previous_unique_hashes = unique_hashes.size();
            unique_hashes.clear();
            return std::nullopt;
        }
        Hash run() {
            while (iteration < cfg.depth / 2) {
                if (const auto result = check_progress())
                    return *result;
                iteration_step();
            }
            return cfg.depth % 2 == 0 ? variable_hash() : cnf_hash();
        }
        std::string operator () () {
            std::string result = std::to_string(run());
            if (cfg.return_measurements) {
                const auto elapsed = Clock::now() - start_time;
                const auto parsing_time = start_time - parsing_start_time;
                result += "," + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count())
                        + "," + std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(parsing_time).count());
                double iteration_count = iteration;
                if (iteration * 2 >= cfg.depth)
                    iteration_count = cfg.depth / 2.;
                result += "," + std::to_string(iteration_count);
            }
            return result;
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
     * @param optimize_first_iteration whether the first iteration should be
     * optimized
     * @param first_progress_check_iteration the first iteration in which the
     * progress check runs
     * @param return_measurements whether the calculation time and parsing time
     * and the amount of iterations that were calculated (possibly half) should
     * be returned
     * @return comma separated list, std::string Weisfeiler-Leman hash,
     * possibly iteration count, calculation time and parsing time
     */
    std::string weisfeiler_leman_hash(
        const char* filename,

        const unsigned depth = 13,
        const bool cross_reference_literals = true,
        const bool rehash_clauses = true,
        const bool optimize_first_iteration = true,
        const unsigned first_progress_check_iteration = 3,
        const bool return_measurements = true
    ) {
        return WeisfeilerLemanHasher(filename, WLHRuntimeConfig {
            depth,
            cross_reference_literals,
            rehash_clauses,
            optimize_first_iteration,
            first_progress_check_iteration,
            return_measurements
        })();
    }
} // namespace CNF

#endif  // ISOHASH2_H_

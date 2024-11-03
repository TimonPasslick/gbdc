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

#ifndef SRC_UTIL_POINTERLESSCNFFORMULA_H_
#define SRC_UTIL_POINTERLESSCNFFORMULA_H_

#include <vector>
#include <algorithm>
#include <memory>
#include <string>

#include "src/util/StreamBuffer.h"
#include "src/util/SolverTypes.h"

class PointerlessCNFFormula {
    std::vector<std::vector<Lit>> clause_length_literals;
    unsigned variables = 0;

 public:
    explicit PointerlessCNFFormula(const char* filename) {
        readDimacsFromFile(filename);
    }

    inline size_t nVars() const {
        return variables;
    }

    struct Clause {
        using It = std::vector<Lit>::const_iterator;
        const It begin_, end_;
        inline It begin() const {
            return begin_;
        }
        inline It end() const {
            return end_;
        }
        std::size_t size() const {
            return end_ - begin_;
        }
    };
    struct ClauseIt {
        const std::vector<std::vector<Lit>>& clause_length_literals;
        unsigned length;
        unsigned clause_start;
        inline ClauseIt& operator ++ () {
            clause_start += length;
            if (clause_start >= clause_length_literals[length].size()) {
                while (++length < clause_length_literals.size() && clause_length_literals[length].size() == 0)
                    ;
                clause_start = 0;
            }
            return *this;
        }
        Clause operator * () {
            const Clause::It begin = clause_length_literals[length].begin() + clause_start;
            return Clause {begin, begin + length};
        }
        bool operator != (ClauseIt o) {
            return length != o.length || clause_start != o.clause_start;
        }
    };
    struct Clauses {
        const ClauseIt begin_, end_;
        inline ClauseIt begin() const {
            return begin_;
        }
        inline ClauseIt end() const {
            return end_;
        }
    };
    Clauses clauses() const {
        return Clauses {
            ++ClauseIt{clause_length_literals, 0, 0},
            {clause_length_literals, (unsigned) clause_length_literals.size(), 0}
        };
    }

private:
    // create gapless representation of variables
    void normalizeVariableNames() {
        std::vector<unsigned> name;
        constexpr unsigned empty = ~0U;
        name.resize(variables+1, empty);
        unsigned max = 0;
        for (std::vector<Lit>& clause_length : clause_length_literals) {
            for (Lit& lit : clause_length) {
                if (name[lit.var()] == empty) name[lit.var()] = max++;
                lit = Lit(name[lit.var()], lit.sign());
            }
        }
        variables = max;
    }

    // https://stackoverflow.com/a/1322548/27720282
    static unsigned next_power_of_2(unsigned n) {
        n--;
        n |= n >> 0b00001;
        n |= n >> 0b00010;
        n |= n >> 0b00100;
        n |= n >> 0b01000;
        n |= n >> 0b10000;
        n++;
        return n;
    }
    void readDimacsFromFile(const char* filename) {
        StreamBuffer in(filename);
        while (in.skipWhitespace()) {
            if (*in == 'p' || *in == 'c') {
                if (!in.skipLine()) break;
            } else {
                std::vector<Lit> clause;
                int plit;
                while (in.readInteger(&plit)) {
                    if (plit == 0) break;
                    const unsigned var = abs(plit);
                    clause.push_back(Lit(var, plit < 0));
                    if (var > variables) variables = var;
                }
                if (clause.size() >= clause_length_literals.size()) {
                    const unsigned new_size = clause.size() + 1;
                    clause_length_literals.reserve(next_power_of_2(new_size));
                    clause_length_literals.resize(new_size);
                }
                std::vector<Lit>& insert_here = clause_length_literals[clause.size()];
                insert_here.reserve(next_power_of_2(insert_here.size() + clause.size()));
                insert_here.insert(insert_here.end(), clause.begin(), clause.end());
            }
        }
        normalizeVariableNames();
    }
};

#endif  // SRC_UTIL_POINTERLESSCNFFORMULA_H_


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
    std::vector<unsigned> clause_bounds;
    std::vector<Lit> literals;
    unsigned variables;

 public:
    PointerlessCNFFormula() : clause_bounds(), literals(), variables(0) {
        clause_bounds.push_back(0);
    }

    explicit PointerlessCNFFormula(const char* filename) : PointerlessCNFFormula() {
        readDimacsFromFile(filename);
    }

    inline size_t nVars() const {
        return variables;
    }

    inline const std::vector<Lit>& literal() {
        return literals;
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
    };
    struct ClauseIt {
        std::vector<unsigned>::const_iterator lower_clause_bound;
        const std::vector<Lit>::const_iterator literals_begin;
        inline ClauseIt& operator ++ () {
            ++lower_clause_bound;
            return *this;
        }
        Clause operator * () {
            return Clause {
                literals_begin + *lower_clause_bound,
                literals_begin + *(lower_clause_bound + 1)
            };
        }
        bool operator != (ClauseIt o) {
            return lower_clause_bound != o.lower_clause_bound;
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
            {clause_bounds.begin(), literals.begin()},
            {clause_bounds.end() - 1, literals.begin()}
        };
    }

    // create gapless representation of variables
    void normalizeVariableNames() {
        std::vector<unsigned> name;
        name.resize(variables+1, 0);
        unsigned int max = 0;
        for (Lit& lit : literals) {
            if (name[lit.var()] == 0) name[lit.var()] = max++;
            lit = Lit(name[lit.var()], lit.sign());
        }
        variables = max;
    }

    void readDimacsFromFile(const char* filename) {
        StreamBuffer in(filename);
        while (in.skipWhitespace()) {
            if (*in == 'p' || *in == 'c') {
                if (!in.skipLine()) break;
            } else {
                int plit;
                while (in.readInteger(&plit)) {
                    if (plit == 0) break;
                    literals.push_back(Lit(abs(plit), plit < 0));
                }
                clause_bounds.push_back(literals.size());
            }
        }
        normalizeVariableNames();
    }
};

#endif  // SRC_UTIL_POINTERLESSCNFFORMULA_H_


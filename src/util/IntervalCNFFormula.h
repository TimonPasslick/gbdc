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

#ifndef SRC_UTIL_INTERVALCNFFORMULA
#define SRC_UTIL_INTERVALCNFFORMULA

#include <vector>
#include <algorithm>
#include <memory>
#include <string>

#include "src/util/StreamBuffer.h"
#include "src/util/SolverTypes.h"

class IntervalCNFFormula {
    std::vector<Lit> literals;
    std::vector<unsigned short> clause_lengths;
    unsigned variables = 0;

 public:
    explicit inline IntervalCNFFormula(const char* filename) {
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
        inline unsigned short size() const {
            return end_ - begin_;
        }
    };
    struct ClauseIt {
        std::vector<Lit>::const_iterator clause_begin;
        std::vector<unsigned short>::const_iterator length_it;
        inline ClauseIt& operator ++ () {
            clause_begin += *(length_it++);
            return *this;
        }
        inline Clause operator * () {
            return Clause {clause_begin, clause_begin + *length_it};
        }
        inline bool operator != (ClauseIt o) {
            return clause_begin != o.clause_begin;
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
            {literals.begin(), clause_lengths.begin()},
            {literals.end(), clause_lengths.end()}
        };
    }

 private:
    // create gapless representation of variables
    void normalizeVariableNames() {
        std::vector<unsigned> name;
        constexpr unsigned empty = ~0U;
        name.resize(variables+1, empty);
        unsigned max = 0;
        for (Lit& lit : literals) {
            if (name[lit.var()] == empty) name[lit.var()] = max++;
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
                unsigned short length = 0;
                int plit;
                while (in.readInteger(&plit)) {
                    if (plit == 0) break;
                    const unsigned var = abs(plit);
                    literals.push_back(Lit(var, plit < 0));
                    ++length;
                    if (var > variables) variables = var;
                }
                clause_lengths.push_back(length);
            }
        }
        normalizeVariableNames();
    }
};

#endif  // SRC_UTIL_INTERVALCNFFORMULA


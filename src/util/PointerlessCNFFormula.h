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

struct PointerlessCNFFormula {
    struct Slice {
        unsigned begin;
        unsigned end;
    };
    std::vector<Lit> literals;
    std::vector<Slice> clauses;
    unsigned variables;

    PointerlessCNFFormula() : literals(), clauses(), variables(0) { }

    explicit PointerlessCNFFormula(const char* filename) : PointerlessCNFFormula() {
        readDimacsFromFile(filename);
    }

    inline size_t nClauses() const {
        return clauses.size();
    }

    inline int newVar() {
        return ++variables;
    }

    inline void clear() {
        literals.clear();
        clauses.clear();
        variables = 0;
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
        Cl clause;
        // TODO Timon: literals.reserve, clauses.reserve
        while (in.skipWhitespace()) {
            if (*in == 'p' || *in == 'c') {
                if (!in.skipLine()) break;
            } else {
                int plit;
                while (in.readInteger(&plit)) {
                    if (plit == 0) break;
                    clause.push_back(Lit(abs(plit), plit < 0));
                }
                // TODO Timon: avoid double allocation
                readClause(clause.begin(), clause.end());
                clause.clear();
            }
        }
    }

    void readClause(std::initializer_list<Lit> list) {
        readClause(list.begin(), list.end());
    }

    void readClauses(const For& formula) {
        for (Cl* clause : formula) {
            readClause(clause->begin(), clause->end());
        }
    }

    template <typename Iterator>
    void readClause(Iterator begin, Iterator end) {
        const unsigned size = end - begin;
        if (size > 0) {
            literals.insert(literals.end(), begin, end);
            // remove redundant literals
            std::sort(literals.end() - size, literals.end());
            unsigned dup = 0;
            const auto begin = literals.end() - size
            for (auto it = begin, jt = begin + 1; jt != literals.end(); ++jt) {
                if (*it != *jt) {  // unique
                    if (it->var() == jt->var()) {
                        literals.resize(literals.size() - size);
                        return;  // no tautologies
                    }
                    ++it;
                    *it = *jt;
                } else {
                    ++dup;
                }
            }
            literals.resize(literals.size() - dup);
            clauses.push_back({
                    (unsigned) literals.size() - size,
                    (unsigned) literals.size()});
            variables = std::max(variables, (unsigned int)literals.back().var());
        }
    }
};

#endif  // SRC_UTIL_POINTERLESSCNFFORMULA_H_


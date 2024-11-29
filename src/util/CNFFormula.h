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

#ifndef SRC_UTIL_CNFFORMULA_H_
#define SRC_UTIL_CNFFORMULA_H_

#include <vector>
#include <algorithm>
#include <memory>
#include <string>

#include "src/util/StreamBuffer.h"
#include "src/util/SolverTypes.h"

template <bool normalized = false>
class TemplateCNFFormula {
    For formula;
    unsigned variables;

 public:
    TemplateCNFFormula() : formula(), variables(0) { }

    explicit TemplateCNFFormula(const char* filename) : TemplateCNFFormula() {
        readDimacsFromFile(filename);
    }

    ~TemplateCNFFormula() {
        for (Cl* clause : formula) {
            delete clause;
        }
    }

    typedef For::const_iterator const_iterator;

    inline const_iterator begin() const {
        return formula.begin();
    }

    inline const_iterator end() const {
        return formula.end();
    }

    inline const Cl* operator[] (int i) const {
        return formula[i];
    }

    inline size_t nVars() const {
        return variables;
    }

    inline size_t nClauses() const {
        return formula.size();
    }

    inline int newVar() {
        return ++variables;
    }

    inline void clear() {
        formula.clear();
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
        For::const_iterator it;
        inline ClauseIt& operator ++ () {
            ++it;
            return *this;
        }
        inline Clause operator * () {
            return Clause {(*it)->begin(), (*it)->end()};
        }
        inline bool operator != (ClauseIt o) {
            return it != it;
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
            {formula.begin()},
            {formula.end()}
        };
    }

    // create gapless representation of variables
    void normalizeVariableNames() {
        std::vector<unsigned> name;
        constexpr unsigned empty = ~0U;
        name.resize(variables+1, empty);
        unsigned int max = 0;
        for (Cl* clause : formula) {
            for (Lit& lit : *clause) {
                if (name[lit.var()] == empty) name[lit.var()] = max++;
                lit = Lit(name[lit.var()], lit.sign());
            }
        }
        variables = max;
    }

    void readDimacsFromFile(const char* filename) {
        StreamBuffer in(filename);
        Cl clause;
        while (in.skipWhitespace()) {
            if (*in == 'p' || *in == 'c') {
                if (!in.skipLine()) break;
            } else {
                int plit;
                while (in.readInteger(&plit)) {
                    if (plit == 0) break;
                    clause.push_back(Lit(abs(plit), plit < 0));
                }
                readClause(clause.begin(), clause.end());
                clause.clear();
            }
        }
        if constexpr (normalized)
            normalizeVariableNames();
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
        Cl* clause = new Cl { begin, end };
        if (clause->size() > 0) {
            // remove redundant literals
            std::sort(clause->begin(), clause->end());
            unsigned dup = 0;
            for (auto it = clause->begin(), jt = clause->begin()+1; jt != clause->end(); ++jt) {
                if (*it != *jt) {  // unique
                    if (it->var() == jt->var()) {
                        delete clause;
                        return;  // no tautologies
                    }
                    ++it;
                    *it = *jt;
                } else {
                    ++dup;
                }
            }
            clause->resize(clause->size() - dup);
            clause->shrink_to_fit();
            variables = std::max(variables, (unsigned int)clause->back().var());
        }
        formula.push_back(clause);
    }
};

using CNFFormula = TemplateCNFFormula<>;
using NormalizedCNFFormula = TemplateCNFFormula<true>;

#endif  // SRC_UTIL_CNFFORMULA_H_


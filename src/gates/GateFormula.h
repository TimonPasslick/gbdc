/*************************************************************************************************
CNFTools -- Copyright (c) 2015, Markus Iser, KIT - Karlsruhe Institute of Technology

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

#ifndef SRC_GATES_GATEFORMULA_H_
#define SRC_GATES_GATEFORMULA_H_

#include <cstdlib>
#include <algorithm>
#include <vector>
#include <set>

#include "src/util/CNFFormula.h"
#include "src/util/Stamp.h"


struct Gate {
    Lit out = lit_Undef;
    For fwd, bwd;
    bool notMono = false;
    std::vector<Lit> inp;

    inline bool isDefined() const { return out != lit_Undef; }
    inline bool hasNonMonotonicParent() const { return notMono; }
};


class GateFormula {
 public:
    std::vector<Cl*> roots;  // top-level clauses
    std::vector<char> inputs;  // mark literals which are used as input to a gate (used in detection of monotonicity)
    std::vector<Gate> gates;  // stores gate-struct for every output
    For remainder;  // stores clauses remaining outside of recognized gate-structure
    Cl* artificialRoot;  // top-level unit-clause that can be generated by normalizeRoots()

    explicit GateFormula(unsigned nVars) :
            roots(), gates(), artificialRoot() {
        inputs.resize(2 + 2*nVars, false);
        gates.resize(2 + nVars);
        artificialRoot = new Cl();
    }

    ~GateFormula() {
        // delete artificialRoot;
    }

    inline void setUsedAsInput(Lit lit) {
        inputs[lit] = true;
    }

    inline bool isUsedAsInput(Lit lit) {
        return inputs[lit];
    }

    inline bool isNestedMonotonic(Lit lit) {
        return !(isUsedAsInput(lit) && isUsedAsInput(~lit));
    }

    void addGate(Lit o, const For& fwd, const For& bwd) {
        Gate& gate = gates[o.var()];
        gate.out = o;
        gate.fwd.insert(gate.fwd.end(), fwd.begin(), fwd.end());
        gate.bwd.insert(gate.bwd.end(), bwd.begin(), bwd.end());
        gate.notMono = !isNestedMonotonic(o);
        for (Cl* c : fwd) for (Lit lit : *c) {
            if (lit != ~o) gate.inp.push_back(lit);
        }
        sort(gate.inp.begin(), gate.inp.end());
        gate.inp.erase(std::unique(gate.inp.begin(), gate.inp.end()), gate.inp.end());
        for (Lit lit : gate.inp) {
            setUsedAsInput(lit);
            if (gate.notMono) setUsedAsInput(~lit);
        }
    }

    Gate& getGate(Lit output) {
        return gates[output.var()];
    }

    inline bool isGateOutput(Lit output) const {
        return gates[output.var()].isDefined();
    }

    typedef std::vector<Gate>::const_iterator const_iterator;

    inline const_iterator begin() const {
        return gates.begin();
    }

    inline const_iterator end() const {
        return gates.end();
    }

    inline const Gate& operator[] (Var var) const {
        return gates[var];
    }

    inline unsigned nGates() const {
        return std::count_if(gates.begin(), gates.end(), [] (const Gate& gate) {
            return gate.isDefined();
        });
    }

    inline unsigned nMonotonicGates() const {
        return std::count_if(gates.begin(), gates.end(), [] (const Gate& gate) {
            return gate.isDefined() && !gate.hasNonMonotonicParent();
        });
    }

    inline unsigned nRoots() const {
        return roots.size();
    }

    const std::vector<Cl*> getRoots() const {
        return roots;
    }

    // create unique list of root literals base on root clauses
    std::vector<Lit> getRootLiterals() {
        std::vector<Lit> literals;

        for (const Cl* c : getRoots()) {
            literals.insert(literals.end(), c->begin(), c->end());
        }
        std::sort(literals.begin(), literals.end());
        auto last = std::unique(literals.begin(), literals.end());
        literals.erase(last, literals.end());

        return literals;
    }

    /**
     * @brief GateAnalyzer::getPrunedProblem
     * @param model
     * @return clauses of all satisfied branches
     */
    For getPrunedProblem(const std::vector<uint8_t>& model) {
        For result(roots.begin(), roots.end());

        std::vector<Lit> literals = getRootLiterals();
        Stamp<uint8_t> visited { gates.size() };

        while (literals.size() > 0) {
            Lit o = literals.back();
            literals.pop_back();
            Gate gate = gates[o.var()];

            if (!gate.isDefined()) continue;

            // std::cout << "output " << gate.out << ", not-mono " << gate.notMono << std::endl;
            // std::cout << "satisfied " << model.satisfies(o) << std::endl;
            // std::cout << "inputs " << gate.inp << std::endl;
            // std::cout << "visited " << visited[o.var()] << std::endl;

            if (!visited[o.var()] && (gate.hasNonMonotonicParent() || model[o])) {  // Skip "don't cares"
                std::copy(gate.fwd.begin(), gate.fwd.end(), result.end());
                if (gate.hasNonMonotonicParent()) {  // BCE
                    std::copy(gate.bwd.begin(), gate.bwd.end(), result.end());
                }
                literals.insert(literals.end(), gate.inp.begin(), gate.inp.end());
                visited.set(o.var());
            }
        }

        result.insert(result.end(), remainder.begin(), remainder.end());

        return result;
    }

    bool hasArtificialRoot() const {
        return artificialRoot->size() > 0;
    }

    Cl* getArtificialRoot() const {
        return artificialRoot;
    }

    void printGates() {
        std::cerr << "Found " << nGates() << " gates of which " << nMonotonicGates() << " are monotonic" << std::endl;
        std::cerr << "Number of root clauses is " << nRoots() << std::endl;
        // std::vector<Lit> outputs;
        // std::vector<bool> done(gates.size());
        // for (Cl* root : roots) {
        //     outputs.insert(outputs.end(), root->begin(), root->end());
        // }
        // for (size_t i = 0; i < outputs.size(); ++i) {
        //     Gate& gate = getGate(outputs[i]);

        //     if (gate.isDefined() && !done[outputs[i].var()]) {
        //         done[outputs[i].var()] = true;
        //         std::cout << gate.out << " <-> f(" << gate.inp << ")" << std::endl;
        //         outputs.insert(outputs.end(), gate.inp.begin(), gate.inp.end());
        //     }
        // }
    }

    Lit getRoot() const {
        assert(roots.size() == 1 && roots.front()->size() == 1);
        return roots.front()->front();
    }

    /**
     * Execute after analysis in order to
     * tronsform many roots to one big and gate with one output
     * Side-effect: introduces a fresh variable
     */
    void normalizeRoots() {
        Var root = Var(gates.size());
        gates.resize(gates.size() + 1);
        // inputs.resize(2 * problem.nVars() + 2, false);
        // create gate
        gates[root].out = Lit(root, false);
        gates[root].notMono = false;
        std::set<Lit> inp;
        roots.insert(roots.end(), remainder.begin(), remainder.end());
        remainder.clear();
        for (Cl* c : roots) {
            inp.insert(c->begin(), c->end());
            c->push_back(Lit(root, true));
            gates[root].fwd.push_back(new Cl(*c));
        }
        gates[root].inp.insert(gates[root].inp.end(), inp.begin(), inp.end());
        this->roots.clear();
        artificialRoot->push_back(gates[root].out);
        this->roots.push_back(artificialRoot);
    }
};

#endif  // SRC_GATES_GATEFORMULA_H_
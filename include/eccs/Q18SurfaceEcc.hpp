/*
 * This file is part of JKQ QFR library which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum/ for more information.
 */

#ifndef QFR_Q18SurfaceEcc_HPP
#define QFR_Q18SurfaceEcc_HPP

#include "Ecc.hpp"
#include "QuantumComputation.hpp"

class Q18SurfaceEcc: public Ecc {
public:
    Q18SurfaceEcc(qc::QuantumComputation& qc, int measureFq, bool decomposeMC, bool cliffOnly);

    static const std::string getName() {
        return "Q18Surface";
    }

protected:
    void initMappedCircuit() override;

    void writeEncoding() override;

    void measureAndCorrect() override;

    void writeDecoding() override;

    void mapGate(const std::unique_ptr<qc::Operation>& gate, qc::QuantumComputation& qc) override;
};

#endif //QFR_Q18SurfaceEcc_HPP

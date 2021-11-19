/*
 * This file is part of JKQ QFR library which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum/ for more information.
 */

#include "eccs/Q7SteaneEcc.hpp"

//7 data qubits, 6 for measuring -> 13 qubits per physical qubit (6 classical for measuring at end)
Q7SteaneEcc::Q7SteaneEcc(qc::QuantumComputation& qc, int measureFq, bool decomposeMC): Ecc({ID::Q7Steane, 7, 6, Q7SteaneEcc::getName()}, qc, measureFq, decomposeMC) {}

void Q7SteaneEcc::initMappedCircuit() {
//method is overridden because we need 2 kinds of classical measurement output registers
    qc.stripIdleQubits(true, false);
    statistics.nInputQubits = qc.getNqubits();
    statistics.nInputClassicalBits = (int)qc.getNcbits();
	statistics.nOutputQubits = qc.getNqubits()*ecc.nRedundantQubits+ecc.nCorrectingBits;
	statistics.nOutputClassicalBits = statistics.nInputClassicalBits+ecc.nCorrectingBits;
	qcMapped.addQubitRegister(statistics.nOutputQubits);
	qcMapped.addClassicalRegister(statistics.nInputClassicalBits);
	qcMapped.addClassicalRegister(3, "qeccX");
	qcMapped.addClassicalRegister(3, "qeccZ");
}

void Q7SteaneEcc::writeEncoding() {
    measureAndCorrect();
    decodingDone = false;
}

void Q7SteaneEcc::measureAndCorrect() {
    const int nQubits = qc.getNqubits();
    const int ancStart = nQubits * ecc.nRedundantQubits;
    const int clAncStart = qc.getNcbits();

    for(int i=0;i<nQubits;i++) {
        qcMapped.reset(ancStart);
        qcMapped.reset(ancStart+1);
        qcMapped.reset(ancStart+2);
        qcMapped.reset(ancStart+3);
        qcMapped.reset(ancStart+4);
        qcMapped.reset(ancStart+5);

        qcMapped.h(ancStart);
        qcMapped.h(ancStart+1);
        qcMapped.h(ancStart+2);
        qcMapped.h(ancStart+3);
        qcMapped.h(ancStart+4);
        qcMapped.h(ancStart+5);

        auto c0 = dd::Control{dd::Qubit(ancStart), dd::Control::Type::pos};
        auto c1 = dd::Control{dd::Qubit(ancStart+1), dd::Control::Type::pos};
        auto c2 = dd::Control{dd::Qubit(ancStart+2), dd::Control::Type::pos};
        auto c3 = dd::Control{dd::Qubit(ancStart+3), dd::Control::Type::pos};
        auto c4 = dd::Control{dd::Qubit(ancStart+4), dd::Control::Type::pos};
        auto c5 = dd::Control{dd::Qubit(ancStart+5), dd::Control::Type::pos};

        //K1: XIXIXIX
        qcMapped.x(i+nQubits*0, c0);
        qcMapped.x(i+nQubits*2, c0);
        qcMapped.x(i+nQubits*4, c0);
        qcMapped.x(i+nQubits*6, c0);

        //K2: IXXIIXX
        qcMapped.x(i+nQubits*1, c1);
        qcMapped.x(i+nQubits*2, c1);
        qcMapped.x(i+nQubits*5, c1);
        qcMapped.x(i+nQubits*6, c1);

        //K3: IIIXXXX
        qcMapped.x(i+nQubits*3, c2);
        qcMapped.x(i+nQubits*4, c2);
        qcMapped.x(i+nQubits*5, c2);
        qcMapped.x(i+nQubits*6, c2);

        //K2: ZIZIZIZ
        qcMapped.z(i+nQubits*0, c3);
        qcMapped.z(i+nQubits*2, c3);
        qcMapped.z(i+nQubits*4, c3);
        qcMapped.z(i+nQubits*6, c3);

        //K3: IZZIIZZ
        qcMapped.z(i+nQubits*1, c4);
        qcMapped.z(i+nQubits*2, c4);
        qcMapped.z(i+nQubits*5, c4);
        qcMapped.z(i+nQubits*6, c4);

        //K1: IIIZZZZ
        qcMapped.z(i+nQubits*3, c5);
        qcMapped.z(i+nQubits*4, c5);
        qcMapped.z(i+nQubits*5, c5);
        qcMapped.z(i+nQubits*6, c5);


        qcMapped.h(ancStart);
        qcMapped.h(ancStart+1);
        qcMapped.h(ancStart+2);
        qcMapped.h(ancStart+3);
        qcMapped.h(ancStart+4);
        qcMapped.h(ancStart+5);

        qcMapped.measure(ancStart, clAncStart);
        qcMapped.measure(ancStart+1, clAncStart+1);
        qcMapped.measure(ancStart+2, clAncStart+2);
        qcMapped.measure(ancStart+3, clAncStart+3);
        qcMapped.measure(ancStart+4, clAncStart+4);
        qcMapped.measure(ancStart+5, clAncStart+5);

        //correct Z_i for i+1 = c0*1+c1*2+c2*4
        //correct X_i for i+1 = c3*1+c4*2+c5*4
        for(unsigned int j=0;j<7;j++) {
            std::unique_ptr<qc::Operation> opZ = std::make_unique<qc::StandardOperation>(qcMapped.getNqubits(), i+j*nQubits, qc::Z);
            const auto pairZ = std::make_pair(dd::Qubit(clAncStart), dd::QubitCount(3));
            qcMapped.emplace_back<qc::ClassicControlledOperation>(opZ, pairZ, j+1U);

            std::unique_ptr<qc::Operation> opX = std::make_unique<qc::StandardOperation>(qcMapped.getNqubits(), i+j*nQubits, qc::X);
            const auto pairX = std::make_pair(dd::Qubit(clAncStart+3), dd::QubitCount(3));
            qcMapped.emplace_back<qc::ClassicControlledOperation>(opX, pairX, j+1U);
        }

    }
}

void Q7SteaneEcc::writeDecoding() {
    const int nQubits = qc.getNqubits();
    const int ancStart = nQubits * ecc.nRedundantQubits;

    for(int i=0;i<nQubits;i++) {
        qcMapped.reset(ancStart);
        qcMapped.h(ancStart);

        auto c = dd::Control{dd::Qubit(ancStart), dd::Control::Type::pos};
        for(int j=0;j<7;j++) {
            qcMapped.z(i+j*nQubits, c);
        }
        qcMapped.h(ancStart);

        qcMapped.measure(ancStart, i);
        qcMapped.swap(ancStart, i);
    }
    decodingDone = true;
}

void Q7SteaneEcc::mapGate(const std::unique_ptr<qc::Operation> &gate) {
    if(decodingDone && gate.get()->getType()!=qc::Measure) {
        writeEncoding();
    }
    const int nQubits = qc.getNqubits();
    qc::NonUnitaryOperation *measureGate=nullptr;
    switch(gate.get()->getType()) {
    case qc::I: break;
    case qc::X:
    case qc::H:
    case qc::Y:
    case qc::Z:
        for(std::size_t t=0;t<gate.get()->getNtargets();t++) {
            int i = gate.get()->getTargets()[t];
            if(gate.get()->getNcontrols()==2 && decomposeMultiControlledGates) {

                auto& ctrls = gate.get()->getControls();
                int idx=0;
                int ctrl2[2] = {-1, -1};
                bool ctrl2T[2] = {true, true};
                for(const auto &ct: ctrls) {
                    ctrl2[idx] = ct.qubit;
                    ctrl2T[idx] = ct.type == dd::Control::Type::pos;
                    idx++;
                }
                if(gate.get()->getType()==qc::X) {
                    for(int j=0;j<7;j++) {
                        writeToffoli(i+j*nQubits, ctrl2[0]+j*nQubits, ctrl2T[0], ctrl2[1]+j*nQubits, ctrl2T[1]);
                    }
                } else if(gate.get()->getType()==qc::Z) {
                    for(int j=0;j<7;j++) {
                        qcMapped.h(i+j*nQubits);
                        writeToffoli(i+j*nQubits, ctrl2[0]+j*nQubits, ctrl2T[0], ctrl2[1]+j*nQubits, ctrl2T[1]);
                        qcMapped.h(i+j*nQubits);
                    }
                } else if(gate.get()->getType()==qc::Y) {
                    for(int j=0;j<7;j++) {
                        writeToffoli(i+j*nQubits, ctrl2[0]+j*nQubits, ctrl2T[0], ctrl2[1]+j*nQubits, ctrl2T[1]);
                        qcMapped.h(i+j*nQubits);
                        writeToffoli(i+j*nQubits, ctrl2[0]+j*nQubits, ctrl2T[0], ctrl2[1]+j*nQubits, ctrl2T[1]);
                        qcMapped.h(i+j*nQubits);
                    }
                } else {
                    gateNotAvailableError(gate);
                }
            } else if(gate.get()->getNcontrols()>2 && decomposeMultiControlledGates) {
                gateNotAvailableError(gate);
            } else if(gate.get()->getNcontrols()) {
                auto& ctrls = gate.get()->getControls();
                for(int j=0;j<7;j++) {
                    dd::Controls ctrls2;
                    for(const auto &ct: ctrls) {
                        ctrls2.insert(dd::Control{dd::Qubit(ct.qubit+j*nQubits), ct.type});
                    }
                    qcMapped.emplace_back<qc::StandardOperation>(nQubits*ecc.nRedundantQubits, ctrls2, i+j*nQubits, gate.get()->getType());
                }
            } else {
                for(int j=0;j<7;j++) {
                    qcMapped.emplace_back<qc::StandardOperation>(nQubits*ecc.nRedundantQubits, i+j*nQubits, gate.get()->getType());
                }
            }
        }
        break;
    //locigal S = 3 physical S's
    case qc::S:
    case qc::Sdag:
        for(std::size_t t=0;t<gate.get()->getNtargets();t++) {
            int i = gate.get()->getTargets()[t];
            if(gate.get()->getNcontrols()>1 && decomposeMultiControlledGates) {
                gateNotAvailableError(gate);
            } else if(gate.get()->getNcontrols()) {
                auto& ctrls = gate.get()->getControls();
                for(int j=0;j<7;j++) {
                    dd::Controls ctrls2;
                    for(const auto &ct: ctrls) {
                        ctrls2.insert(dd::Control{dd::Qubit(ct.qubit+j*nQubits), ct.type});
                    }
                    qcMapped.emplace_back<qc::StandardOperation>(nQubits*ecc.nRedundantQubits, ctrls2, i+j*nQubits, gate.get()->getType());
                    qcMapped.emplace_back<qc::StandardOperation>(nQubits*ecc.nRedundantQubits, ctrls2, i+j*nQubits, gate.get()->getType());
                    qcMapped.emplace_back<qc::StandardOperation>(nQubits*ecc.nRedundantQubits, ctrls2, i+j*nQubits, gate.get()->getType());
                }
            } else {
                for(int j=0;j<7;j++) {
                    qcMapped.emplace_back<qc::StandardOperation>(nQubits*ecc.nRedundantQubits, i+j*nQubits, gate.get()->getType());
                    qcMapped.emplace_back<qc::StandardOperation>(nQubits*ecc.nRedundantQubits, i+j*nQubits, gate.get()->getType());
                    qcMapped.emplace_back<qc::StandardOperation>(nQubits*ecc.nRedundantQubits, i+j*nQubits, gate.get()->getType());
                }
            }
        }
        break;
    case qc::Measure:
        if(!decodingDone) {
            measureAndCorrect();
            writeDecoding();
        }
        measureGate = (qc::NonUnitaryOperation*)gate.get();
        for(std::size_t j=0;j<measureGate->getNclassics();j++) {
            qcMapped.measure(measureGate->getClassics()[j], measureGate->getTargets()[j]);
        }
        break;
    case qc::T:
    case qc::Tdag:
    case qc::V:
    case qc::Vdag:
    case qc::U3:
    case qc::U2:
    case qc::Phase:
    case qc::SX:
    case qc::SXdag:
    case qc::RX:
    case qc::RY:
    case qc::RZ:
    case qc::SWAP:
    case qc::iSWAP:
    case qc::Peres:
    case qc::Peresdag:
    case qc::Compound:
    case qc::ClassicControlled:
    default:
        gateNotAvailableError(gate);
    }
}


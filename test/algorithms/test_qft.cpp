/*
 * This file is part of JKQ QFR library which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum/ for more information.
 */

#include "algorithms/QFT.hpp"

#include "gtest/gtest.h"
#include <cmath>
#include <iostream>

class QFT: public testing::TestWithParam<dd::QubitCount> {
protected:
    void TearDown() override {
        if (!sim.isTerminal())
            dd->decRef(sim);
        if (!func.isTerminal())
            dd->decRef(func);
        dd->garbageCollect(true);

        // number of complex table entries after clean-up should equal initial number of entries
        EXPECT_EQ(dd->cn.magnitudeTable.getCount(), initialMagCount);
        EXPECT_EQ(dd->cn.phaseTable.getCount(), initialPhaseCount);
        // number of available cache entries after clean-up should equal initial number of entries
        EXPECT_EQ(dd->cn.complexCache.getCount(), initialCacheCount);
    }

    void SetUp() override {
        nqubits           = GetParam();
        dd                = std::make_unique<dd::Package>(nqubits);
        initialCacheCount = dd->cn.complexCache.getCount();
        initialMagCount   = dd->cn.magnitudeTable.getCount();
        initialPhaseCount = dd->cn.phaseTable.getCount();
    }

    dd::QubitCount               nqubits = 0;
    std::unique_ptr<dd::Package> dd;
    std::unique_ptr<qc::QFT>     qc;
    std::size_t                  initialCacheCount = 0;
    std::size_t                  initialMagCount   = 0;
    std::size_t                  initialPhaseCount = 0;
    qc::VectorDD                 sim{};
    qc::MatrixDD                 func{};
};

/// Findings from the QFT Benchmarks:
/// The DDpackage has to be able to represent all 2^n different amplitudes in order to produce correct results
/// The smallest entry seems to be closely related to '1-cos(pi/2^(n-1))'
/// The following CN::TOLERANCE values suffice up until a certain number of qubits:
/// 	10e-10	..	18 qubits
///		10e-11	..	20 qubits
///		10e-12	..	22 qubits
///		10e-13	..	23 qubits
/// The accuracy of double floating points allows for a minimal CN::TOLERANCE value of 10e-15
///	Utilizing more qubits requires the use of fp=long double
constexpr dd::QubitCount QFT_MAX_QUBITS = 20;

INSTANTIATE_TEST_SUITE_P(QFT, QFT,
                         testing::Range(static_cast<dd::QubitCount>(0), static_cast<dd::QubitCount>(QFT_MAX_QUBITS + 1), 3),
                         [](const testing::TestParamInfo<QFT::ParamType>& info) {
			dd::QubitCount nqubits = info.param;
			std::stringstream ss{};
			ss << static_cast<std::size_t>(nqubits);
			if (nqubits == 1) {
				ss << "_qubit";
			} else {
				ss << "_qubits";
			}
			return ss.str(); });

TEST_P(QFT, Functionality) {
    // there should be no error constructing the circuit
    ASSERT_NO_THROW({ qc = std::make_unique<qc::QFT>(nqubits); });

    // there should be no error building the functionality
    ASSERT_NO_THROW({ func = qc->buildFunctionality(dd); });

    qc->printStatistics(std::cout);
    // QFT DD should consist of 2^n nodes
    ASSERT_EQ(dd->size(func), std::pow(2, nqubits));

    // Force garbage collection of compute table and complex table
    dd->garbageCollect(true);

    // the final DD should store all 2^n different amplitudes
    // since only positive real values are stored in the complex table
    // this number has to be divided by 4
    EXPECT_NEAR(dd->cn.magnitudeTable.getCount(), 2, 1);
    EXPECT_EQ(dd->cn.phaseTable.getCount(), std::max(2U, static_cast<unsigned int>(std::pow(2, nqubits - 2))) - 1);

    // top edge weight should equal sqrt(0.5)^n
    EXPECT_NEAR(dd::MagnitudeTable<>::Entry::val(func.w.mag), static_cast<dd::fp>(std::pow(1.L / std::sqrt(2.L), nqubits)), dd->cn.magnitudeTable.tolerance());

    // first row and first column should consist only of (1/sqrt(2))**nqubits
    for (unsigned long long i = 0; i < std::pow(static_cast<long double>(2), nqubits); ++i) {
        auto c = dd->getValueByPath(func, 0, i);
        EXPECT_NEAR(c.mag, static_cast<dd::fp>(std::pow(1.L / std::sqrt(2.L), nqubits)), dd->cn.magnitudeTable.tolerance());
        EXPECT_NEAR(c.phase, 0, dd->cn.phaseTable.tolerance());
        c = dd->getValueByPath(func, i, 0);
        EXPECT_NEAR(c.mag, static_cast<dd::fp>(std::pow(1.L / std::sqrt(2.L), nqubits)), dd->cn.magnitudeTable.tolerance());
        EXPECT_NEAR(c.phase, 0, dd->cn.phaseTable.tolerance());
    }
}

TEST_P(QFT, FunctionalityRecursive) {
    // there should be no error constructing the circuit
    ASSERT_NO_THROW({ qc = std::make_unique<qc::QFT>(nqubits); });

    // there should be no error building the functionality
    ASSERT_NO_THROW({ func = qc->buildFunctionalityRecursive(dd); });

    qc->printStatistics(std::cout);
    // QFT DD should consist of 2^n nodes
    ASSERT_EQ(dd->size(func), std::pow(2, nqubits));

    // Force garbage collection of compute table and complex table
    dd->garbageCollect(true);

    // the final DD should store all 2^n different amplitudes
    // since only positive real values are stored in the complex table
    // this number has to be divided by 4
    EXPECT_NEAR(dd->cn.magnitudeTable.getCount(), 2, 1);
    EXPECT_EQ(dd->cn.phaseTable.getCount(), std::max(2U, static_cast<unsigned int>(std::pow(2, nqubits - 2))) - 1);

    // top edge weight should equal sqrt(0.5)^n
    EXPECT_NEAR(dd::MagnitudeTable<>::Entry::val(func.w.mag), static_cast<dd::fp>(std::pow(1.L / std::sqrt(2.L), nqubits)), dd->cn.magnitudeTable.tolerance());

    // first row and first column should consist only of (1/sqrt(2))**nqubits
    for (unsigned long long i = 0; i < std::pow(static_cast<long double>(2), nqubits); ++i) {
        auto c = dd->getValueByPath(func, 0, i);
        EXPECT_NEAR(c.mag, static_cast<dd::fp>(std::pow(1.L / std::sqrt(2.L), nqubits)), dd->cn.magnitudeTable.tolerance());
        EXPECT_NEAR(c.phase, 0, dd->cn.phaseTable.tolerance());
        c = dd->getValueByPath(func, i, 0);
        EXPECT_NEAR(c.mag, static_cast<dd::fp>(std::pow(1.L / std::sqrt(2.L), nqubits)), dd->cn.magnitudeTable.tolerance());
        EXPECT_NEAR(c.phase, 0, dd->cn.phaseTable.tolerance());
    }
}

TEST_P(QFT, Simulation) {
    // there should be no error constructing the circuit
    ASSERT_NO_THROW({ qc = std::make_unique<qc::QFT>(nqubits); });

    // there should be no error building the functionality
    ASSERT_NO_THROW({
        auto in = dd->makeZeroState(nqubits);
        sim     = qc->simulate(in, dd);
    });
    qc->printStatistics(std::cout);

    // QFT DD |0...0> sim should consist of n nodes
    ASSERT_EQ(dd->size(sim), nqubits + 1);

    // Force garbage collection of compute table and complex table
    dd->garbageCollect(true);

    // top edge weight should equal 1
    EXPECT_NEAR(dd::MagnitudeTable<>::Entry::val(sim.w.mag), 1, dd->cn.magnitudeTable.tolerance());
    EXPECT_NEAR(dd::PhaseTable<>::Entry::val(sim.w.phase), 0, dd->cn.phaseTable.tolerance());

    // first column should consist only of sqrt(0.5)^n's
    for (unsigned long long i = 0; i < std::pow(static_cast<long double>(2), nqubits); ++i) {
        auto c = dd->getValueByPath(sim, i);
        EXPECT_NEAR(c.mag, static_cast<dd::fp>(std::pow(1.L / std::sqrt(2.L), nqubits)), dd->cn.magnitudeTable.tolerance());
        EXPECT_NEAR(c.phase, 0, dd->cn.phaseTable.tolerance());
    }
}

TEST_P(QFT, FunctionalityRecursiveEquality) {
    // there should be no error constructing the circuit
    ASSERT_NO_THROW({ qc = std::make_unique<qc::QFT>(nqubits); });
    std::cout << *qc << std::endl;

    // there should be no error building the functionality recursively
    ASSERT_NO_THROW({ func = qc->buildFunctionalityRecursive(dd); });

    // there should be no error building the functionality regularly
    qc::MatrixDD funcRec{};
    ASSERT_NO_THROW({ funcRec = qc->buildFunctionality(dd); });

    ASSERT_EQ(func, funcRec);
    dd->decRef(funcRec);
}

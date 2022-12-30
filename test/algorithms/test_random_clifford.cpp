/*
 * This file is part of MQT QFR library which is released under the MIT license.
 * See file README.md or go to https://www.cda.cit.tum.de/research/quantum/ for more information.
 */

#include "algorithms/RandomCliffordCircuit.hpp"
#include "dd/FunctionalityConstruction.hpp"
#include "dd/Simulation.hpp"

#include "gtest/gtest.h"
#include <string>

class RandomClifford: public testing::TestWithParam<dd::QubitCount> {
protected:
    void TearDown() override {}
    void SetUp() override {}
};

INSTANTIATE_TEST_SUITE_P(RandomClifford, RandomClifford,
                         testing::Range(static_cast<dd::QubitCount>(1), static_cast<dd::QubitCount>(9)),
                         [](const testing::TestParamInfo<RandomClifford::ParamType>& inf) {
                             // Generate names for test cases
                             dd::QubitCount    nqubits = inf.param;
                             std::stringstream ss{};
                             ss << static_cast<std::size_t>(nqubits) << "_qubits";
                             return ss.str();
                         });

TEST_P(RandomClifford, simulate) {
    const dd::QubitCount nq = GetParam();

    auto dd = std::make_unique<dd::Package<>>(nq);
    auto qc = qc::RandomCliffordCircuit(nq, nq * nq, 12345);
    auto in = dd->makeZeroState(nq);

    std::cout << qc << std::endl;
    ASSERT_NO_THROW({ simulate(&qc, in, dd); });
    qc.printStatistics(std::cout);
}

TEST_P(RandomClifford, buildFunctionality) {
    const dd::QubitCount nq = GetParam();

    auto dd = std::make_unique<dd::Package<>>(nq);
    auto qc = qc::RandomCliffordCircuit(nq, nq * nq, 12345);
    std::cout << qc << std::endl;
    ASSERT_NO_THROW({ buildFunctionality(&qc, dd); });
    qc.printStatistics(std::cout);
}

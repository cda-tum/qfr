#include "CircuitOptimizer.hpp"
#include "algorithms/BernsteinVazirani.hpp"
#include "algorithms/Entanglement.hpp"
#include "algorithms/Grover.hpp"
#include "algorithms/QFT.hpp"
#include "algorithms/QPE.hpp"
#include "algorithms/RandomCliffordCircuit.hpp"
#include "algorithms/WState.hpp"
#include "dd/Benchmark.hpp"
#include "dd/FunctionalityConstruction.hpp"
#include "dd/statistics/PackageStatistics.hpp"

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include <string>
#include <utility>

namespace dd {

std::string runCLI(const char* cmd) {
  std::array<char, 128> buffer{};
  std::string result;
  const std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
  if (!pipe) {
    throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
  }
  std::string realResult = result.erase(result.size() - 1);
  return realResult;
}

[[maybe_unused]] static const std::string CURRENT_BRANCH =
    runCLI("git symbolic-ref --short -q HEAD");
[[maybe_unused]] static const std::string CURRENT_COMMIT =
    runCLI("git rev-parse --short HEAD");
static const std::string FILENAME = "results.json";
static const std::string FILENAME_REDUCED = "results_reduced.json";

static constexpr std::size_t SEED = 42U;

void transposeAndReduceJson(const std::string& fileName,
                            const std::string& outFilename) {
  std::ifstream ifs(fileName);
  nlohmann::json j;
  ifs >> j;
  ifs.close();

  nlohmann::json k;

  for (const auto& [algorithm, resultsA] : j.items()) {
    for (const auto& [type, resultsT] : resultsA.items()) {
      for (const auto& [nqubits, resultsN] : resultsT.items()) {
        for (const auto& [branch, resultsB] : resultsN.items()) {
          const auto& runtime = resultsB["runtime"];
          k[algorithm][type][nqubits]["runtime"][branch] = runtime;

          const auto& gateCount = resultsB["gate_count"];
          k[algorithm][type][nqubits]["gate_count"][branch] = gateCount;

          const auto& dd = resultsB["dd"];
          const auto& activeMemoryMiB = dd["active_memory_mib"];
          k[algorithm][type][nqubits]["dd"]["active_memory_mib"][branch] =
              activeMemoryMiB;
          const auto& peakMemoryMiB = dd["peak_memory_mib"];
          k[algorithm][type][nqubits]["dd"]["peak_memory_mib"][branch] =
              peakMemoryMiB;
          for (const auto& stat : {"matrix", "vector", "density_matrix",
                                   "real_numbers", "compute_tables"}) {
            for (const auto& [key, value] : dd[stat].items()) {
              if (value == "unused") {
                k[algorithm][type][nqubits]["dd"][stat][key][branch] = value;
                continue;
              }

              for (const auto& [key2, value2] : value.items()) {
                if ((std::strcmp(stat, "matrix") != 0 ||
                     std::strcmp(stat, "vector") != 0 ||
                     std::strcmp(stat, "density_matrix") != 0) &&
                    key == "unique_table") {
                  for (const auto& [key3, value3] : value2.items()) {
                    if (!key3.empty()) {
                      k[algorithm][type][nqubits]["dd"][stat][key][key2][key3]
                       [branch] = value3;
                    } else {
                      k[algorithm][type][nqubits]["dd"][stat][key][key2]
                       [branch] = value3;
                    }
                  }
                  continue;
                }
                k[algorithm][type][nqubits]["dd"][stat][key][key2][branch] =
                    value2;
              }
            }
          }
        }
      }
    }
  }

  for (const auto& [algorithm, resultsA] : k.items()) {
    for (const auto& [type, resultsT] : resultsA.items()) {
      for (const auto& [nqubits, resultsN] : resultsT.items()) {
        auto& dd = resultsN["dd"];
        dd.erase("density_matrix");

        auto& computeTables = dd["compute_tables"];
        computeTables.erase("density_matrix_add");
        computeTables.erase("density_density_mult");
        computeTables.erase("density_noise_operations");
        computeTables.erase("stochastic_noise_operations");
        computeTables.erase("matrix_kronecker");
        computeTables.erase("vector_kronecker");
        computeTables.erase("vector_inner_product");
        computeTables.erase("matrix_conjugate_transpose");

        if (type == "Functionality") {
          dd.erase("vector");
          computeTables.erase("vector_add");
          computeTables.erase("matrix_vector_mult");
        }
      }
    }
  }

  std::ofstream ofs(outFilename);
  ofs << k.dump(2U);
  ofs.close();
}

// a function that parses a nlohmann::json from a file "results.json", populates
// it with the results of the current run and writes it back to the file

void verifyAndSave(const std::string& name, const std::string& type,
                   qc::QuantumComputation& qc, const Experiment& exp) {
  EXPECT_TRUE(exp.success());

  nlohmann::json j;
  std::fstream file(FILENAME, std::ios::in | std::ios::out | std::ios::ate);
  if (!file.is_open()) {
    std::ofstream outputFile(FILENAME);
    outputFile << nlohmann::json();
  } else if (file.tellg() == 0) {
    file << nlohmann::json();
  }
  file.close();

  std::ifstream ifs(FILENAME);
  ifs >> j;
  ifs.close();

  auto& entry = j[name][type][std::to_string(qc.getNqubits())][CURRENT_BRANCH];
  // Change this line to CURRENT_COMMIT when comparing commits, or anything else
  // to distinguish between runs

  entry["gate_count"] = qc.getNindividualOps();
  entry["runtime"] = exp.runtime.count();

  // collect statistics from DD package
  entry["dd"] = exp.stats;

  std::ofstream ofs(FILENAME);
  ofs << j.dump(2U);
  ofs.close();

  transposeAndReduceJson(FILENAME, FILENAME_REDUCED);
}

class GHZEval : public testing::TestWithParam<std::size_t> {
protected:
  void TearDown() override {}
  void SetUp() override {
    nqubits = GetParam();
    qc = std::make_unique<qc::Entanglement>(nqubits);
  }

  std::size_t nqubits = 0;
  std::unique_ptr<qc::Entanglement> qc;
};

INSTANTIATE_TEST_SUITE_P(GHZ, GHZEval,
                         testing::Values(256U, 512U, 1024U, 2048U, 4096U));

TEST_P(GHZEval, GHZSimulation) {
  const auto out = benchmarkSimulate(*qc);
  verifyAndSave("GHZ", "Simulation", *qc, *out);
}

TEST_P(GHZEval, GHZFunctionality) {
  const auto out = benchmarkFunctionalityConstruction(*qc);
  verifyAndSave("GHZ", "Functionality", *qc, *out);
}

class WStateEval : public testing::TestWithParam<std::size_t> {
protected:
  void TearDown() override {}
  void SetUp() override {
    nqubits = GetParam();
    qc = std::make_unique<qc::WState>(nqubits);
  }

  std::size_t nqubits = 0;
  std::unique_ptr<qc::WState> qc;
};

INSTANTIATE_TEST_SUITE_P(WState, WStateEval,
                         testing::Values(256U, 512U, 1024U, 2048U, 4096U));

TEST_P(WStateEval, WStateSimulation) {
  const auto out = benchmarkSimulate(*qc);
  verifyAndSave("WState", "Simulation", *qc, *out);
}

TEST_P(WStateEval, WStateFunctionality) {
  const auto out = benchmarkFunctionalityConstruction(*qc);
  verifyAndSave("WState", "Functionality", *qc, *out);
}

// add dynamic
class BVEval : public testing::TestWithParam<std::size_t> {
protected:
  void TearDown() override {}
  void SetUp() override {
    nqubits = GetParam();
    qc = std::make_unique<qc::BernsteinVazirani>(nqubits);
    qc::CircuitOptimizer::removeFinalMeasurements(*qc);
  }

  std::size_t nqubits = 0;
  std::unique_ptr<qc::BernsteinVazirani> qc;
};

INSTANTIATE_TEST_SUITE_P(BV, BVEval,
                         testing::Values(255U, 511U, 1023U, 2047U, 4095U));

TEST_P(BVEval, BVSimulation) {
  const auto out = benchmarkSimulate(*qc);
  verifyAndSave("BV", "Simulation", *qc, *out);
}

TEST_P(BVEval, BVFunctionality) {
  const auto out = benchmarkFunctionalityConstruction(*qc);
  verifyAndSave("BV", "Functionality", *qc, *out);
}

class QFTEval : public testing::TestWithParam<std::size_t> {
protected:
  void TearDown() override {}
  void SetUp() override {
    nqubits = GetParam();
    qc = std::make_unique<qc::QFT>(nqubits, false);
  }

  std::size_t nqubits = 0;
  std::unique_ptr<qc::QFT> qc;
};

INSTANTIATE_TEST_SUITE_P(QFT, QFTEval,
                         testing::Values(256U, 512U, 1024U, 2048U, 4096U));

TEST_P(QFTEval, QFTSimulation) {
  const auto out = benchmarkSimulate(*qc);
  verifyAndSave("QFT", "Simulation", *qc, *out);
}

class QFTEvalFunctionality : public testing::TestWithParam<std::size_t> {
protected:
  void TearDown() override {}
  void SetUp() override {
    nqubits = GetParam();
    qc = std::make_unique<qc::QFT>(nqubits, false);
  }

  std::size_t nqubits = 0;
  std::unique_ptr<qc::QFT> qc;
};

INSTANTIATE_TEST_SUITE_P(QFT, QFTEvalFunctionality,
                         testing::Values(18U, 19U, 20U, 21U, 22U));

TEST_P(QFTEvalFunctionality, QFTFunctionality) {
  const auto out = benchmarkFunctionalityConstruction(*qc);
  verifyAndSave("QFT", "Functionality", *qc, *out);
}

class GroverEval : public testing::TestWithParam<qc::Qubit> {
protected:
  void TearDown() override {}
  void SetUp() override {
    nqubits = GetParam();
    qc = std::make_unique<qc::Grover>(nqubits, 12345U);
    dd = std::make_unique<dd::Package<>>(qc->getNqubits());
  }

  qc::Qubit nqubits = 0;
  std::unique_ptr<dd::Package<>> dd;
  std::unique_ptr<qc::Grover> qc;
};

INSTANTIATE_TEST_SUITE_P(Grover, GroverEval,
                         testing::Values(27U, 31U, 35U, 39U, 41U));

TEST_P(GroverEval, GroverSimulator) {
  const auto start = std::chrono::high_resolution_clock::now();

  // apply state preparation setup
  qc::QuantumComputation statePrep(qc->getNqubits());
  qc->setup(statePrep);
  auto s = buildFunctionality(&statePrep, dd);
  auto e = dd->multiply(s, dd->makeZeroState(qc->getNqubits()));
  dd->incRef(e);
  dd->decRef(s);

  qc::QuantumComputation groverIteration(qc->getNqubits());
  qc->oracle(groverIteration);
  qc->diffusion(groverIteration);

  auto iter = buildFunctionalityRecursive(&groverIteration, dd);
  std::bitset<128U> iterBits(qc->iterations);
  auto msb = static_cast<std::size_t>(std::floor(std::log2(qc->iterations)));
  auto f = iter;
  dd->incRef(f);
  for (std::size_t j = 0U; j <= msb; ++j) {
    if (iterBits[j]) {
      auto g = dd->multiply(f, e);
      dd->incRef(g);
      dd->decRef(e);
      e = g;
      dd->garbageCollect();
    }
    if (j < msb) {
      auto tmp = dd->multiply(f, f);
      dd->incRef(tmp);
      dd->decRef(f);
      f = tmp;
    }
  }
  dd->decRef(f);
  const auto end = std::chrono::high_resolution_clock::now();
  const auto runtime =
      std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
  std::unique_ptr<SimulationExperiment> exp =
      std::make_unique<SimulationExperiment>();
  exp->dd = std::move(dd);
  exp->sim = e;
  exp->runtime = runtime;
  exp->stats = dd::getStatistics(exp->dd.get());

  verifyAndSave("Grover", "Simulation", *qc, *exp);
}

TEST_P(GroverEval, GroverFunctionality) {
  const auto out = benchmarkFunctionalityConstruction(*qc, true);
  verifyAndSave("Grover", "Functionality", *qc, *out);
}

class QPEEval : public testing::TestWithParam<std::size_t> {
protected:
  void TearDown() override {}
  void SetUp() override {
    nqubits = GetParam();
    qc = std::make_unique<qc::QPE>(nqubits, false);
    qc::CircuitOptimizer::removeFinalMeasurements(*qc);
  }

  std::size_t nqubits = 0;
  std::unique_ptr<qc::QPE> qc;
};

INSTANTIATE_TEST_SUITE_P(QPE, QPEEval,
                         testing::Values(14U, 15U, 16U, 17U, 18U));

TEST_P(QPEEval, QPESimulation) {
  const auto out = benchmarkSimulate(*qc);
  verifyAndSave("QPE", "Simulation", *qc, *out);
}

class QPEEvalFunctionality : public testing::TestWithParam<std::size_t> {
protected:
  void TearDown() override {}
  void SetUp() override {
    nqubits = GetParam();
    qc = std::make_unique<qc::QPE>(nqubits, false);
    qc::CircuitOptimizer::removeFinalMeasurements(*qc);
  }

  std::size_t nqubits = 0;
  std::unique_ptr<qc::QPE> qc;
};

INSTANTIATE_TEST_SUITE_P(QPE, QPEEvalFunctionality,
                         testing::Values(7U, 8U, 9U, 10U, 11U));

TEST_P(QPEEvalFunctionality, QPEFunctionality) {
  const auto out = benchmarkFunctionalityConstruction(*qc);
  verifyAndSave("QPE", "Functionality", *qc, *out);
}

class RandomCliffordEval : public testing::TestWithParam<std::size_t> {
protected:
  void TearDown() override {}
  void SetUp() override {
    nqubits = GetParam();
    qc = std::make_unique<qc::RandomCliffordCircuit>(nqubits, nqubits * nqubits,
                                                     SEED);
  }

  std::size_t nqubits = 0;
  std::unique_ptr<qc::RandomCliffordCircuit> qc;
};

INSTANTIATE_TEST_SUITE_P(RandomCliffordCircuit, RandomCliffordEval,
                         testing::Values(14U, 15U, 16U, 17U, 18U));

TEST_P(RandomCliffordEval, RandomCliffordSimulation) {
  const auto out = benchmarkSimulate(*qc);
  verifyAndSave("RandomClifford", "Simulation", *qc, *out);
}

class RandomCliffordEvalFunctionality
    : public testing::TestWithParam<std::size_t> {
protected:
  void TearDown() override {}
  void SetUp() override {
    nqubits = GetParam();
    qc = std::make_unique<qc::RandomCliffordCircuit>(nqubits, nqubits * nqubits,
                                                     SEED);
  }

  std::size_t nqubits = 0;
  std::unique_ptr<qc::RandomCliffordCircuit> qc;
};

INSTANTIATE_TEST_SUITE_P(RandomCliffordCircuit, RandomCliffordEvalFunctionality,
                         testing::Values(7U, 8U, 9U, 10U, 11U));

TEST_P(RandomCliffordEvalFunctionality, RandomCliffordFunctionality) {
  const auto out = benchmarkFunctionalityConstruction(*qc);
  verifyAndSave("RandomClifford", "Functionality", *qc, *out);
}

} // namespace dd

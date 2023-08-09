#include "algorithms/WState.hpp"

namespace qc {
void fGate(QuantumComputation& qc, const std::size_t i, const std::size_t j, const std::size_t n,
           const std::size_t k) {
  const auto theta = std::acos(std::sqrt(1.0 / static_cast<double>(n - k + 1)));
  qc.ry(static_cast<Qubit>(j), -theta);
  qc.z(static_cast<Qubit>(j), qc::Control{static_cast<Qubit>(i)});
  qc.ry(static_cast<Qubit>(j), theta);
}

WState::WState(const std::size_t nq) : QuantumComputation(nq) {
  name = "wstate_" + std::to_string(nq);
  const auto top = static_cast<Qubit>(nq - 1);

  x(top);

  for (std::size_t m = 1; m < nq; m++) {
    fGate(*this, nq - m, nq - m - 1, nq, m);
  }

  for (std::size_t k = nq - 1; k > 0; k--) {
    x(static_cast<Qubit>(k), qc::Control{static_cast<Qubit>(k - 1)});
  }
}
} // namespace qc
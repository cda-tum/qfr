//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#include "Layer.hpp"

#include "UndirectedGraph.hpp"

#include <memory>
#include <set>
#include <unordered_set>

namespace qc {

auto Layer::constructDAG(const QuantumComputation& qc) -> void {
  const auto nQubits = qc.getNqubits();
  // For a pair of self-canceling operations like two consecutive X operations
  // or RY rotations with opposite angles the first operations is a
  // destructive operation that disables operations until the consecutive
  // constructive operation enables them again
  // ---
  // those that add a (+) edge to the current group members
  std::vector<std::vector<std::shared_ptr<DAGVertex>>> constructive(nQubits);
  // those that add a (-) edge to the current group members
  std::vector<std::vector<std::shared_ptr<DAGVertex>>> destructive(nQubits);
  // those that are already in the current group where all gates commute on
  // this qubit
  std::vector<std::vector<std::shared_ptr<DAGVertex>>> currentGroup(nQubits);
  // lookahead of 1 serves as a buffer for the next operation on each qubit
  std::vector<std::shared_ptr<DAGVertex>> lookahead(nQubits);
  // the predecessor of the current group members
  std::vector<std::vector<std::shared_ptr<DAGVertex>>> predecessorGroup(
      nQubits);
  // all operations acting on a qubit (processed so far) excluding
  // constructive and destructive operations
  std::vector<std::vector<std::shared_ptr<DAGVertex>>> qubitOperations(nQubits);
  // iterate over all operations in the quantum circuit
  for (const auto& op : qc) {
    // create a vertex for the current operation
    std::shared_ptr<DAGVertex> const vertex =
        DAGVertex::create(&op, &executableSet);
    // iterate over all qubits the operation acts on
    for (const auto& qubit : op->getUsedQubits()) {
      // check whether the lookahead is empty
      if (lookahead[qubit] == nullptr) {
        // here: the lookahead is empty
        // add the current operation to the lookahead
        lookahead[qubit] = vertex;
      } else {
        // here: the lookahead is not empty
        // get the current vertex from the lookahead and store the new
        // vertex in the lookahead
        // Note: this might seem odd and one might think that the gates can be
        // eliminated, however, since we also allow global gates, those might
        // cancel out each other on a particular qubit but not on all qubits
        // and, hence, cannot be eliminated.
        std::shared_ptr<DAGVertex> const current = lookahead[qubit];
        lookahead[qubit] = vertex;
        // check whether the current operation is the inverse of the
        // lookahead
        if ((*current->getOperation())
                ->isInverseOf(**lookahead[qubit]->getOperation())) {
          // here: the current operation is the inverse of the lookahead
          // add an enabling edge from the lookahead to all operations on this
          // qubit including the destructive ones
          for (const auto& qubitOperation : qubitOperations[qubit]) {
            lookahead[qubit]->addEnabledSuccessor(qubitOperation);
          }
          for (const auto& qubitOperation : destructive[qubit]) {
            lookahead[qubit]->addEnabledSuccessor(qubitOperation);
          }
          // add the lookahead to the constructive group
          constructive[qubit].emplace_back(lookahead[qubit]);
          // add a disabling edge to all operations on this qubit including
          // the destructive ones
          for (const auto& qubitOperation : qubitOperations[qubit]) {
            current->addDisabledSuccessor(qubitOperation);
          }
          for (const auto& qubitOperation : destructive[qubit]) {
            current->addDisabledSuccessor(qubitOperation);
          }
          // add an enabling edge to the lookahead
          current->addEnabledSuccessor(lookahead[qubit]);
          // add the current vertex to the destructive group
          destructive[qubit].emplace_back(current);
          // clear the lookahead
          lookahead[qubit] = nullptr;
        } else {
          // add an enabling edge from each constructive operation
          for (const auto& constructiveOp : constructive[qubit]) {
            constructiveOp->addEnabledSuccessor(current);
          }
          // add a disabling edge from each destructive operation
          for (const auto& destructiveOp : destructive[qubit]) {
            destructiveOp->addDisabledSuccessor(current);
          }
          // check whether the current operation commutes with the current
          // group members
          if (!currentGroup[qubit].empty() and
              !(*currentGroup[qubit][0]->getOperation())
                   ->commutesAtQubit(**current->getOperation(), qubit)) {
            // here: the current operation does not commute with the current
            // group members and is not the inverse of the lookahead
            // --> start a new group
            predecessorGroup[qubit].clear();
            predecessorGroup[qubit] = currentGroup[qubit];
            currentGroup[qubit].clear();
          }
          // add an enabling edge from each predecessor
          for (const auto& predecessor : predecessorGroup[qubit]) {
            predecessor->addEnabledSuccessor(current);
          }
          // add the current vertex to the current group
          currentGroup[qubit].emplace_back(current);
          qubitOperations[qubit].emplace_back(current);
        }
      }
    }
  }
  // process the remaining lookahead for every qubit
  for (Qubit qubit = 0; qubit < nQubits; ++qubit) {
    if (lookahead[qubit] != nullptr) {
      auto const current = lookahead[qubit];
      lookahead[qubit] = nullptr;
      // add an enabling edge from each constructive operation
      for (const auto& constructiveOp : constructive[qubit]) {
        constructiveOp->addEnabledSuccessor(current);
      }
      // add a disabling edge from each destructive operation
      for (const auto& destructiveOp : destructive[qubit]) {
        destructiveOp->addDisabledSuccessor(current);
      }
      // check whether the current operation commutes with the current
      // group members
      if (!currentGroup[qubit].empty() and
          !(*currentGroup[qubit][0]->getOperation())
               ->commutesAtQubit(**current->getOperation(), qubit)) {
        // here: the current operation does not commute with the current
        // group members and is not the inverse of the lookahead
        // --> start a new group
        predecessorGroup[qubit].clear();
        predecessorGroup[qubit] = currentGroup[qubit];
        currentGroup[qubit].clear();
      }
      // add an enabling edge from each predecessor
      for (const auto& predecessor : predecessorGroup[qubit]) {
        predecessor->addEnabledSuccessor(current);
      }
      // add the current vertex to the current group
      currentGroup[qubit].emplace_back(current);
      qubitOperations[qubit].emplace_back(current);
    }
  }
}
auto Layer::constructInteractionGraph(OpType opType, std::size_t nctrl) const
    -> UndirectedGraph<Qubit, std::shared_ptr<DAGVertex>> {
  switch (opType) {
  case X:
  case Y:
  case Z:
  case RX:
  case RY:
  case RZ:
    if (nctrl == 1) {
      break;
    }
    [[fallthrough]];
  default:
    std::stringstream ss;
    ss << "The operation type ";
    for (std::size_t i = 0; i < nctrl; ++i) {
      ss << "c";
    }
    ss << opType << " is not supported for constructing an interaction graph.";
    throw std::invalid_argument(ss.str());
  }
  UndirectedGraph<Qubit, std::shared_ptr<DAGVertex>> graph;
  for (const auto& vertex : *executableSet) {
    const auto& gate = *vertex->getOperation();
    if (gate->getType() == opType and gate->getNcontrols() == nctrl) {
      const auto& usedQubits = gate->getUsedQubits();
      if (usedQubits.size() != 2) {
        throw std::invalid_argument(
            "The interaction graph can only be constructed for two-qubit "
            "gates.");
      }
      std::vector<Qubit> q(usedQubits.cbegin(), usedQubits.cend());
      graph.addEdge(q[0], q[1], vertex);
    }
  }
  return graph;
}
auto Layer::getExecutablesOfType(OpType opType, std::size_t nctrl) const
    -> std::vector<std::shared_ptr<DAGVertex>> {
  std::vector<std::shared_ptr<DAGVertex>> executables;
  for (const auto& vertex : *executableSet) {
    if ((*vertex->getOperation())->getType() == opType and
        (*vertex->getOperation())->getNcontrols() == nctrl) {
      executables.emplace_back(vertex);
    }
  }
  return executables;
}
} // namespace qc
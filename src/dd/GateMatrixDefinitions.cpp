/*
 * Copyright (c) 2024 Chair for Design Automation, TUM
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#include "dd/GateMatrixDefinitions.hpp"

#include "dd/DDDefinitions.hpp"
#include "ir/operations/OpType.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <vector>

namespace {
using namespace dd;

inline GateMatrix uMat(const fp lambda, const fp phi, const fp theta) {
  return GateMatrix{{{std::cos(theta / 2.), 0.},
                     {-std::cos(lambda) * std::sin(theta / 2.),
                      -std::sin(lambda) * std::sin(theta / 2.)},
                     {std::cos(phi) * std::sin(theta / 2.),
                      std::sin(phi) * std::sin(theta / 2.)},
                     {std::cos(lambda + phi) * std::cos(theta / 2.),
                      std::sin(lambda + phi) * std::cos(theta / 2.)}}};
}

inline GateMatrix u2Mat(const fp lambda, const fp phi) {
  return GateMatrix{
      SQRT2_2,
      {-std::cos(lambda) * SQRT2_2, -std::sin(lambda) * SQRT2_2},
      {std::cos(phi) * SQRT2_2, std::sin(phi) * SQRT2_2},
      {std::cos(lambda + phi) * SQRT2_2, std::sin(lambda + phi) * SQRT2_2}};
}

inline GateMatrix pMat(const fp lambda) {
  return GateMatrix{1, 0, 0, {std::cos(lambda), std::sin(lambda)}};
}

inline GateMatrix rxMat(const fp lambda) {
  return GateMatrix{{{std::cos(lambda / 2.), 0.},
                     {0., -std::sin(lambda / 2.)},
                     {0., -std::sin(lambda / 2.)},
                     {std::cos(lambda / 2.), 0.}}};
}

inline GateMatrix ryMat(const fp lambda) {
  return GateMatrix{{{std::cos(lambda / 2.), 0.},
                     {-std::sin(lambda / 2.), 0.},
                     {std::sin(lambda / 2.), 0.},
                     {std::cos(lambda / 2.), 0.}}};
}

inline GateMatrix rzMat(const fp lambda) {
  return GateMatrix{{{std::cos(lambda / 2.), -std::sin(lambda / 2.)},
                     0,
                     0,
                     {std::cos(lambda / 2.), std::sin(lambda / 2.)}}};
}

inline TwoQubitGateMatrix rxxMat(const fp theta) {
  const auto cosTheta = std::cos(theta / 2.);
  const auto sinTheta = std::sin(theta / 2.);

  return TwoQubitGateMatrix{{{cosTheta, 0, 0, {0., -sinTheta}},
                             {0, cosTheta, {0., -sinTheta}, 0},
                             {0, {0., -sinTheta}, cosTheta, 0},
                             {std::complex{0., -sinTheta}, 0, 0, cosTheta}}};
}

inline TwoQubitGateMatrix ryyMat(const fp theta) {
  const auto cosTheta = std::cos(theta / 2.);
  const auto sinTheta = std::sin(theta / 2.);

  return TwoQubitGateMatrix{{{cosTheta, 0, 0, {0., sinTheta}},
                             {0, cosTheta, {0., -sinTheta}, 0},
                             {0, {0., -sinTheta}, cosTheta, 0},
                             {std::complex{0., sinTheta}, 0, 0, cosTheta}}};
}

inline TwoQubitGateMatrix rzzMat(const fp theta) {
  const auto cosTheta = std::cos(theta / 2.);
  const auto sinTheta = std::sin(theta / 2.);

  return TwoQubitGateMatrix{{{std::complex{cosTheta, -sinTheta}, 0, 0, 0},
                             {0, {cosTheta, sinTheta}, 0, 0},
                             {0, 0, {cosTheta, sinTheta}, 0},
                             {0, 0, 0, {cosTheta, -sinTheta}}}};
}

inline TwoQubitGateMatrix rzxMat(const fp theta) {
  const auto cosTheta = std::cos(theta / 2.);
  const auto sinTheta = std::sin(theta / 2.);

  return TwoQubitGateMatrix{{{cosTheta, {0., -sinTheta}, 0, 0},
                             {std::complex{0., -sinTheta}, cosTheta, 0, 0},
                             {0, 0, cosTheta, {0., sinTheta}},
                             {0, 0, {0., sinTheta}, cosTheta}}};
}

inline TwoQubitGateMatrix xxMinusYYMat(const fp theta, const fp beta = 0.) {
  const auto cosTheta = std::cos(theta / 2.);
  const auto sinTheta = std::sin(theta / 2.);
  const auto cosBeta = std::cos(beta);
  const auto sinBeta = std::sin(beta);

  return TwoQubitGateMatrix{
      {{cosTheta, 0, 0, {-sinBeta * sinTheta, -cosBeta * sinTheta}},
       {0, 1, 0, 0},
       {0, 0, 1, 0},
       {std::complex{sinBeta * sinTheta, -cosBeta * sinTheta}, 0, 0,
        cosTheta}}};
}

inline TwoQubitGateMatrix xxPlusYYMat(const fp theta, const fp beta = 0.) {
  const auto cosTheta = std::cos(theta / 2.);
  const auto sinTheta = std::sin(theta / 2.);
  const auto cosBeta = std::cos(beta);
  const auto sinBeta = std::sin(beta);

  return TwoQubitGateMatrix{
      {{1, 0, 0, 0},
       {0, cosTheta, {sinBeta * sinTheta, -cosBeta * sinTheta}, 0},
       {0, {-sinBeta * sinTheta, -cosBeta * sinTheta}, cosTheta, 0},
       {0, 0, 0, 1}}};
}
} // namespace

namespace dd {

GateMatrix opToSingleQubitGateMatrix(const qc::OpType t,
                                     const std::vector<fp>& params) {
  switch (t) {
  case qc::I:
    return {1, 0, 0, 1};
  case qc::H:
    return {SQRT2_2, SQRT2_2, SQRT2_2, -SQRT2_2};
  case qc::X:
    return {0, 1, 1, 0};
  case qc::Y:
    return {0, {0, -1}, {0, 1}, 0};
  case qc::Z:
    return {1, 0, 0, -1};
  case qc::S:
    return {1, 0, 0, {0, 1}};
  case qc::Sdg:
    return {1, 0, 0, {0, -1}};
  case qc::T:
    return {1, 0, 0, {SQRT2_2, SQRT2_2}};
  case qc::Tdg:
    return {1, 0, 0, {SQRT2_2, -SQRT2_2}};
  case qc::SX:
    return {std::complex{0.5, 0.5}, std::complex{0.5, -0.5},
            std::complex{0.5, -0.5}, std::complex{0.5, 0.5}};
  case qc::SXdg:
    return {std::complex{0.5, -0.5}, std::complex{0.5, 0.5},
            std::complex{0.5, 0.5}, std::complex{0.5, -0.5}};
  case qc::V:
    return {SQRT2_2, {0., -SQRT2_2}, {0., -SQRT2_2}, SQRT2_2};
  case qc::Vdg:
    return {SQRT2_2, {0., SQRT2_2}, {0., SQRT2_2}, SQRT2_2};
  case qc::U:
    return uMat(params[0], params[1], params[2]);
  case qc::U2:
    return u2Mat(params[0], params[1]);
  case qc::P:
    return pMat(params[0]);
  case qc::RX:
    return rxMat(params[0]);
  case qc::RY:
    return ryMat(params[0]);
  case qc::RZ:
    return rzMat(params[0]);
  default:
    assert(false && "Invalid single-qubit gate type");
    return {}; // unreachable
  }
}

TwoQubitGateMatrix opToTwoQubitGateMatrix(const qc::OpType t,
                                          const std::vector<fp>& params) {
  switch (t) {
  case qc::SWAP:
    return {{{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}}};
  case qc::iSWAP:
    return {{{1, 0, 0, 0}, {0, 0, {0, 1}, 0}, {0, {0, 1}, 0, 0}, {0, 0, 0, 1}}};
  case qc::iSWAPdg:
    return {
        {{1, 0, 0, 0}, {0, 0, {0, -1}, 0}, {0, {0, -1}, 0, 0}, {0, 0, 0, 1}}};
  case qc::ECR:
    return {{{0, 0, SQRT2_2, {0, SQRT2_2}},
             {0, 0, {0, SQRT2_2}, SQRT2_2},
             {SQRT2_2, {0, -SQRT2_2}, 0, 0},
             {std::complex{0., -SQRT2_2}, SQRT2_2, 0, 0}}};
  case qc::DCX:
    return {{{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {0, 1, 0, 0}}};
  case qc::Peres:
    return {{{0, 0, 0, 1}, {0, 0, 1, 0}, {1, 0, 0, 0}, {0, 1, 0, 0}}};
  case qc::Peresdg:
    return {{{0, 0, 1, 0}, {0, 0, 0, 1}, {0, 1, 0, 0}, {1, 0, 0, 0}}};
  case qc::RXX:
    return rxxMat(params[0]);
  case qc::RYY:
    return ryyMat(params[0]);
  case qc::RZZ:
    return rzzMat(params[0]);
  case qc::RZX:
    return rzxMat(params[0]);
  case qc::XXminusYY:
    return xxMinusYYMat(params[0], params[1]);
  case qc::XXplusYY:
    return xxPlusYYMat(params[0], params[1]);
  default:
    assert(false && "Invalid two-qubit gate type");
    return {{{}, {}, {}, {}}}; // unreachable
  }
}
} // namespace dd

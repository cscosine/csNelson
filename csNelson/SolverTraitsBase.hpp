#pragma once

#include <csBlockMatrix/Global.hpp>

#include <csNelson/OrderingTraits.hpp>

#include <csNelson/SolverCholeskyDense.hpp>
#include <csNelson/SolverCholeskySparse.hpp>

namespace csNelson {

constexpr int solverCholeskyDense = 1;
constexpr int solverCholeskySparse = 2;

template <int solverType> struct SolverTraits {};

template <> struct SolverTraits<solverCholeskyDense> {
  template <class HessianTraits, int choleskyOrderingIgnored = -1>
  using SolverBase = SolverCholeskyDenseBase<HessianTraits::matType, typename HessianTraits::Type, HessianTraits::B,
                                             HessianTraits::NB>;

  template <class HessianTraits, int choleskyOrderingIgnored = -1>
  using Solver =
      SolverCholeskyDense<HessianTraits::matType, typename HessianTraits::Type, HessianTraits::B, HessianTraits::NB>;

  template <class EigenMat, int choleskyOrderingIgnored = -1> using SolverEigen = SolverCholeskyEigenDense<EigenMat>;
};

template <> struct SolverTraits<solverCholeskySparse> {
  template <class HessianTraits, int choleskyOrdering>
  using SolverBase =
      SolverCholeskySparseBase<HessianTraits::matType, typename HessianTraits::Type, HessianTraits::B,
                               HessianTraits::NB, typename CholeskyOrderingTraits<choleskyOrdering>::Ordering>;

  template <class HessianTraits, int choleskyOrdering>
  using Solver = SolverCholeskySparse<HessianTraits::matType, typename HessianTraits::Type, HessianTraits::B,
                                      HessianTraits::NB, typename CholeskyOrderingTraits<choleskyOrdering>::Ordering>;

  template <class EigenMat, int choleskyOrdering>
  using SolverEigen = SolverCholeskyEigenSparse<EigenMat, typename CholeskyOrderingTraits<choleskyOrdering>::Ordering>;
};

} // namespace csNelson

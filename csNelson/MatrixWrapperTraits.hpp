#pragma once
#include <csNelson/Global.hpp>

#include <csBlockMatrix/Global.hpp>

#include <csNelson/MatrixDenseWrapper.hpp>
#include <csNelson/MatrixSparseWrapper.hpp>

namespace csNelson {

constexpr int matrixWrapperDense = 1;
constexpr int matrixWrapperSparse = 2;

template <int wrapperType>
struct MatrixWrapperTraits {};

template <>
struct MatrixWrapperTraits<matrixWrapperDense> {
  template <int matType, class T, int Ordering, int BR, int BC, int NBR = csBlockMatrix::Dynamic, int NBC = csBlockMatrix::Dynamic>
  using Wrapper = DenseWrapper<matType, T, Ordering, BR, BC, NBR, NBC>;
};

template <>
struct MatrixWrapperTraits<matrixWrapperSparse> {
  template <int matType, class T, int Ordering, int BR, int BC, int NBR = csBlockMatrix::Dynamic, int NBC = csBlockMatrix::Dynamic>
  using Wrapper = SparseWrapper<matType, T, Ordering, BR, BC, NBR, NBC>;
};

} // namespace csNelson

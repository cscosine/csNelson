#pragma once
#include "SingleSectionHessian.h"

#include "csBlockMatrix/DenseMatrixBlock.hpp"
#include "csBlockMatrix/DiagonalMatrixBlock.hpp"
#include "csBlockMatrix/SparseCoeffDiagonalMatrixBlock.hpp"
#include "csBlockMatrix/SparseCoeffMatrixBlock.hpp"
#include "csBlockMatrix/SparseMatrixBlock.hpp"

#include "csBlockMatrix/VectorBlock.hpp"

namespace csNelson {

template <int matTypeV, class T, int B, int NB>
SingleSectionHessian<matTypeV, T, B, NB>::SingleSectionHessian()
    : _chi2(0) {}

template <int matTypeV, class T, int B, int NB>
SingleSectionHessian<matTypeV, T, B, NB>::~SingleSectionHessian() {}

template <int matTypeV, class T, int B, int NB>
void SingleSectionHessian<matTypeV, T, B, NB>::resize(BlockSizeTypePar blockSizes, int nBlocks,
                                                      const csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>::CSPtr& sp) {
  auto blockDescriptor = MatType::BlockDescriptor::squareMatrix(blockSizes, nBlocks);
  _H.resize(blockDescriptor, sp);
  _b.resize(blockDescriptor.rowDescriptionCSPtr());
}

template <int matTypeV, class T, int B, int NB>
void SingleSectionHessian<matTypeV, T, B, NB>::clearAll() {
  _chi2 = 0;
  _H.setZero();
  _b.setZero();
}

template <int matTypeV, class T, int B, int NB>
T SingleSectionHessian<matTypeV, T, B, NB>::maxAbsValBVect() const {
  return _b.mat().cwiseAbs().maxCoeff();
}

} // namespace csNelson

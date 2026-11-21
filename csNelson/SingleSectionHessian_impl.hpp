#pragma once
#include <csNelson/SingleSectionHessian.hpp>

#include <csBlockMatrix/DenseMatrixBlock_impl.hpp>
#include <csBlockMatrix/DiagonalMatrixBlock_impl.hpp>
#include <csBlockMatrix/SparseCoeffDiagonalMatrixBlock_impl.hpp>
#include <csBlockMatrix/SparseCoeffMatrixBlock_impl.hpp>
#include <csBlockMatrix/SparseMatrixBlock_impl.hpp>

#include <csBlockMatrix/VectorBlock_impl.hpp>

namespace csNelson {

template <int matTypeV, class T, int B, int NB>
SingleSectionHessian<matTypeV, T, B, NB>::SingleSectionHessian() : _chi2(0) {}

template <int matTypeV, class T, int B, int NB> SingleSectionHessian<matTypeV, T, B, NB>::~SingleSectionHessian() {}

template <int matTypeV, class T, int B, int NB>
void SingleSectionHessian<matTypeV, T, B, NB>::resize(
    BlockSizeTypePar blockSizes, int nBlocks,
    const csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>::CSPtr &sp) {
  auto blockDescriptor = MatType::BlockDescriptor::squareMatrix(blockSizes, nBlocks);
  _H.resize(blockDescriptor, sp);
  _b.resize(blockDescriptor.rowDescriptionCSPtr());
}

template <int matTypeV, class T, int B, int NB> void SingleSectionHessian<matTypeV, T, B, NB>::clearAll() {
  _chi2 = 0;
  _H.setZero();
  _b.setZero();
}

template <int matTypeV, class T, int B, int NB> T SingleSectionHessian<matTypeV, T, B, NB>::maxAbsValBVect() const {
  return _b.mat().cwiseAbs().maxCoeff();
}

} // namespace csNelson

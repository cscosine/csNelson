#pragma once
#include <csNelson/DoubleSectionHessian.hpp>

#include <csBlockMatrix/DenseMatrixBlock_impl.hpp>
#include <csBlockMatrix/DiagonalMatrixBlock_impl.hpp>
#include <csBlockMatrix/SparseCoeffDiagonalMatrixBlock_impl.hpp>
#include <csBlockMatrix/SparseCoeffMatrixBlock_impl.hpp>
#include <csBlockMatrix/SparseMatrixBlock_impl.hpp>

#include <csBlockMatrix/VectorBlock_impl.hpp>

#include <csNelson/DoubleSectionHessianMatrices_impl.hpp>

namespace csNelson {

template <int matTypeUv, int matTypeWv, int matTypeVv, class Tv, int BUv, int BVv, int NBUv, int NBVv>
DoubleSectionHessian<matTypeUv, matTypeWv, matTypeVv, Tv, BUv, BVv, NBUv, NBVv>::DoubleSectionHessian() : _chi2(0) {}

template <int matTypeUv, int matTypeWv, int matTypeVv, class Tv, int BUv, int BVv, int NBUv, int NBVv>
DoubleSectionHessian<matTypeUv, matTypeWv, matTypeVv, Tv, BUv, BVv, NBUv, NBVv>::~DoubleSectionHessian() {}

template <int matTypeUv, int matTypeWv, int matTypeVv, class Tv, int BUv, int BVv, int NBUv, int NBVv>
void DoubleSectionHessian<matTypeUv, matTypeWv, matTypeVv, Tv, BUv, BVv, NBUv, NBVv>::resize(
    BlockSizeTypeParU blockSizesU, int nBlocksU, BlockSizeTypeParV blockSizesV, int nBlocksV,
    const csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>::CSPtr &spU,
    const csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>::CSPtr &spV,
    const csBlockMatrix::SparsityPattern<csBlockMatrix::RowMajor>::CSPtr &spW) {
  auto blockDescriptorU = MatTypeU::BlockDescriptor::squareMatrix(blockSizesU, nBlocksU);
  _H.U().resize(blockDescriptorU, spU);
  _b.bU().resize(blockDescriptorU.rowDescriptionCSPtr());

  auto blockDescriptorV = MatTypeV::BlockDescriptor::squareMatrix(blockSizesV, nBlocksV);
  _H.V().resize(blockDescriptorV, spV);
  _b.bV().resize(blockDescriptorV.rowDescriptionCSPtr());

  auto blockDescriptorW = typename MatTypeW::BlockDescriptor(blockDescriptorU.rowDescriptionCSPtr(),
                                                             blockDescriptorV.colDescriptionCSPtr());
  _H.W().resize(blockDescriptorW, spW);
}

template <int matTypeUv, int matTypeWv, int matTypeVv, class Tv, int BUv, int BVv, int NBUv, int NBVv>
void DoubleSectionHessian<matTypeUv, matTypeWv, matTypeVv, Tv, BUv, BVv, NBUv, NBVv>::clearAll() {
  _chi2 = 0;
  _H.U().setZero();
  _H.V().setZero();
  _H.W().setZero();
  _b.bU().setZero();
  _b.bV().setZero();
}

template <int matTypeUv, int matTypeWv, int matTypeVv, class Tv, int BUv, int BVv, int NBUv, int NBVv>
Tv DoubleSectionHessian<matTypeUv, matTypeWv, matTypeVv, Tv, BUv, BVv, NBUv, NBVv>::maxAbsValBVect() const {
  return std::max(_b.bU().mat().cwiseAbs().maxCoeff(), _b.bV().mat().cwiseAbs().maxCoeff());
}

} // namespace csNelson

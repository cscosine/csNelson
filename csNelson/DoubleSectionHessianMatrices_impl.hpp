#pragma once
#include <csNelson/DoubleSectionHessianMatrices.hpp>

#include <csBlockMatrix/DenseMatrixBlock_impl.hpp>
#include <csBlockMatrix/DiagonalMatrixBlock_impl.hpp>
#include <csBlockMatrix/SparseCoeffDiagonalMatrixBlock_impl.hpp>
#include <csBlockMatrix/SparseCoeffMatrixBlock_impl.hpp>
#include <csBlockMatrix/SparseMatrixBlock_impl.hpp>
#include <csBlockMatrix/VectorBlock_impl.hpp>

namespace csNelson {

template <int matTypeUv, int matTypeVv, int matTypeWv, class Tv, int BUv, int BVv, int NBUv, int NBVv>
DoubleSectionHessianMatrices<matTypeUv, matTypeVv, matTypeWv, Tv, BUv, BVv, NBUv,
                             NBVv>::DoubleSectionHessianMatrices() {}

template <int matTypeUv, int matTypeVv, int matTypeWv, class Tv, int BUv, int BVv, int NBUv, int NBVv>
DoubleSectionHessianMatrices<matTypeUv, matTypeVv, matTypeWv, Tv, BUv, BVv, NBUv,
                             NBVv>::~DoubleSectionHessianMatrices() {}

//----------------------------------------------------------------------
template <class Tv, int BUv, int BVv, int NBUv, int NBVv>
DoubleSectionHessianVectors<Tv, BUv, BVv, NBUv, NBVv>::DoubleSectionHessianVectors() {}

template <class Tv, int BUv, int BVv, int NBUv, int NBVv>
DoubleSectionHessianVectors<Tv, BUv, BVv, NBUv, NBVv>::~DoubleSectionHessianVectors() {}

} // namespace csNelson

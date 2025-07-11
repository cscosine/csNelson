#pragma once
#include "DoubleSectionHessianMatrices.h"

#include "csBlockMatrix/DenseMatrixBlock.hpp"
#include "csBlockMatrix/DiagonalMatrixBlock.hpp"
#include "csBlockMatrix/SparseCoeffDiagonalMatrixBlock.hpp"
#include "csBlockMatrix/SparseCoeffMatrixBlock.hpp"
#include "csBlockMatrix/SparseMatrixBlock.hpp"
#include "csBlockMatrix/VectorBlock.hpp"

namespace csNelson {

template <int matTypeUv, int matTypeVv, int matTypeWv, class Tv, int BUv, int BVv, int NBUv, int NBVv>
DoubleSectionHessianMatrices<matTypeUv, matTypeVv, matTypeWv, Tv, BUv, BVv, NBUv, NBVv>::DoubleSectionHessianMatrices() {}

template <int matTypeUv, int matTypeVv, int matTypeWv, class Tv, int BUv, int BVv, int NBUv, int NBVv>
DoubleSectionHessianMatrices<matTypeUv, matTypeVv, matTypeWv, Tv, BUv, BVv, NBUv, NBVv>::~DoubleSectionHessianMatrices() {}

//----------------------------------------------------------------------
template <class Tv, int BUv, int BVv, int NBUv, int NBVv>
DoubleSectionHessianVectors<Tv, BUv, BVv, NBUv, NBVv>::DoubleSectionHessianVectors() {}

template <class Tv, int BUv, int BVv, int NBUv, int NBVv>
DoubleSectionHessianVectors<Tv, BUv, BVv, NBUv, NBVv>::~DoubleSectionHessianVectors() {}

} // namespace csNelson

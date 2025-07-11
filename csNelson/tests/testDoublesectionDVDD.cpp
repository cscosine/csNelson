#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "testDoubleSectionCommon.hpp"

#include "csNelson/GaussNewton.hpp"
#include "csNelson/LevenbergMarquardt.hpp"

#include <array>
#include <iostream>

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3d
    : public Points2d3dBase<Points2d3d<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv, csBlockMatrix::Dynamic,
                            csBlockMatrix::Variable, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic> {
  std::vector<int> v_sizes;

public:
  Points2d3d()
      : v_sizes(numPoints3d, Point3d::blockSize) {}

  int parameterUSize(void) const override {
    return Point2d::blockSize;
  }

  const std::vector<int>& parameterVSize(void) const override {
    return v_sizes;
  }

  int numParametersU() const override {
    return numPoints2d;
  }
};

using PointsSection_DenseDenseDense = Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, csBlockMatrix::BlockDense>;
using PointsSection_DenseDenseDiago = Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal>;
using PointsSection_DenseDenseSpars = Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse>;
using PointsSection_DenseDenseSpaco =
    Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse>;
using PointsSection_DenseDiagoDense = Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSection_DenseDiagoDiago =
    Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSection_DenseDiagoSpars = Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSection_DenseDiagoSpaco =
    Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSection_DenseSparsDense = Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense>;
using PointsSection_DenseSparsDiago = Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal>;
using PointsSection_DenseSparsSpars = Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse>;
using PointsSection_DenseSparsSpaco =
    Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse>;
using PointsSection_DenseSpacoDense =
    Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense>;
using PointsSection_DenseSpacoDiago =
    Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal>;
using PointsSection_DenseSpacoSpars =
    Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse>;
using PointsSection_DenseSpacoSpaco =
    Points2d3d<csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse>;

using PointsSection_DiagoDenseDense = Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, csBlockMatrix::BlockDense>;
using PointsSection_DiagoDenseDiago =
    Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal>;
using PointsSection_DiagoDenseSpars = Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse>;
using PointsSection_DiagoDenseSpaco =
    Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse>;
using PointsSection_DiagoDiagoDense =
    Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSection_DiagoDiagoDiago =
    Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSection_DiagoDiagoSpars =
    Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSection_DiagoDiagoSpaco =
    Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSection_DiagoSparsDense = Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense>;
using PointsSection_DiagoSparsDiago =
    Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal>;
using PointsSection_DiagoSparsSpars = Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse>;
using PointsSection_DiagoSparsSpaco =
    Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse>;
using PointsSection_DiagoSpacoDense =
    Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense>;
using PointsSection_DiagoSpacoDiago =
    Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal>;
using PointsSection_DiagoSpacoSpars =
    Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse>;
using PointsSection_DiagoSpacoSpaco =
    Points2d3d<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse>;

using PointsSection_SparsDenseDense = Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, csBlockMatrix::BlockDense>;
using PointsSection_SparsDenseDiago = Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal>;
using PointsSection_SparsDenseSpars = Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse>;
using PointsSection_SparsDenseSpaco =
    Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse>;
using PointsSection_SparsDiagoDense = Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSection_SparsDiagoDiago =
    Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSection_SparsDiagoSpars = Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSection_SparsDiagoSpaco =
    Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSection_SparsSparsDense = Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense>;
using PointsSection_SparsSparsDiago = Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal>;
using PointsSection_SparsSparsSpars = Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse>;
using PointsSection_SparsSparsSpaco =
    Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse>;
using PointsSection_SparsSpacoDense =
    Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense>;
using PointsSection_SparsSpacoDiago =
    Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal>;
using PointsSection_SparsSpacoSpars =
    Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse>;
using PointsSection_SparsSpacoSpaco =
    Points2d3d<csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse>;

using PointsSection_SpacoDenseDense =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, csBlockMatrix::BlockDense>;
using PointsSection_SpacoDenseDiago =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal>;
using PointsSection_SpacoDenseSpars =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse>;
using PointsSection_SpacoDenseSpaco =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse>;
using PointsSection_SpacoDiagoDense =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSection_SpacoDiagoDiago =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSection_SpacoDiagoSpars =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSection_SpacoDiagoSpaco =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSection_SpacoSparsDense =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense>;
using PointsSection_SpacoSparsDiago =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal>;
using PointsSection_SpacoSparsSpars =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse>;
using PointsSection_SpacoSparsSpaco =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse>;
using PointsSection_SpacoSpacoDense =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense>;
using PointsSection_SpacoSpacoDiago =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal>;
using PointsSection_SpacoSpacoSpars =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse>;
using PointsSection_SpacoSpacoSpaco =
    Points2d3d<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE(
    "DoubleSection-", "[DoubleSection-]", PointsSection_DenseDenseDense, PointsSection_DenseDenseDiago, PointsSection_DenseDenseSpars,
    PointsSection_DenseDenseSpaco, PointsSection_DenseDiagoDense, PointsSection_DenseDiagoDiago, PointsSection_DenseDiagoSpars,
    PointsSection_DenseDiagoSpaco, PointsSection_DenseSparsDense, PointsSection_DenseSparsDiago, PointsSection_DenseSparsSpars,
    PointsSection_DenseSparsSpaco, PointsSection_DenseSpacoDense, PointsSection_DenseSpacoDiago, PointsSection_DenseSpacoSpars,
    PointsSection_DenseSpacoSpaco,
    //------------------------------------
    PointsSection_DiagoDenseDense, PointsSection_DiagoDenseDiago, PointsSection_DiagoDenseSpars, PointsSection_DiagoDenseSpaco,
    PointsSection_DiagoDiagoDense, PointsSection_DiagoDiagoDiago, PointsSection_DiagoDiagoSpars, PointsSection_DiagoDiagoSpaco,
    PointsSection_DiagoSparsDense, PointsSection_DiagoSparsDiago, PointsSection_DiagoSparsSpars, PointsSection_DiagoSparsSpaco,
    PointsSection_DiagoSpacoDense, PointsSection_DiagoSpacoDiago, PointsSection_DiagoSpacoSpars, PointsSection_DiagoSpacoSpaco,
    //------------------------------------
    PointsSection_SparsDenseDense, PointsSection_SparsDenseDiago, PointsSection_SparsDenseSpars, PointsSection_SparsDenseSpaco,
    PointsSection_SparsDiagoDense, PointsSection_SparsDiagoDiago, PointsSection_SparsDiagoSpars, PointsSection_SparsDiagoSpaco,
    PointsSection_SparsSparsDense, PointsSection_SparsSparsDiago, PointsSection_SparsSparsSpars, PointsSection_SparsSparsSpaco,
    PointsSection_SparsSpacoDense, PointsSection_SparsSpacoDiago, PointsSection_SparsSpacoSpars, PointsSection_SparsSpacoSpaco,
    //------------------------------------
    PointsSection_SpacoDenseDense, PointsSection_SpacoDenseDiago, PointsSection_SpacoDenseSpars, PointsSection_SpacoDenseSpaco,
    PointsSection_SpacoDiagoDense, PointsSection_SpacoDiagoDiago, PointsSection_SpacoDiagoSpars, PointsSection_SpacoDiagoSpaco,
    PointsSection_SpacoSparsDense, PointsSection_SpacoSparsDiago, PointsSection_SpacoSparsSpars, PointsSection_SpacoSparsSpaco,
    PointsSection_SpacoSpacoDense, PointsSection_SpacoSpacoDiago, PointsSection_SpacoSpacoSpars, PointsSection_SpacoSpacoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

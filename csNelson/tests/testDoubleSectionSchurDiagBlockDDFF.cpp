#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "testDoubleSectionSchurCommon.hpp"

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dDDFF : public Points2d3dBase<Points2d3dDDFF<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                                             csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, numPoints2d, numPoints3d> {

public:
  int parameterUSize(void) const override {
    return Point2d::blockSize;
  }

  int parameterVSize(void) const override {
    return Point3d::blockSize;
  }
};

using PointsSectionDDFF_DenseDiagoDense =
    Points2d3dDDFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDDFF_DenseDiagoDiago =
    Points2d3dDDFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDDFF_DenseDiagoSpars =
    Points2d3dDDFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDDFF_DenseDiagoSpaco =
    Points2d3dDDFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDDFF_DiagoDiagoDense =
    Points2d3dDDFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDDFF_DiagoDiagoDiago =
    Points2d3dDDFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDDFF_DiagoDiagoSpars =
    Points2d3dDDFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDDFF_DiagoDiagoSpaco =
    Points2d3dDDFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDDFF_SparsDiagoDense =
    Points2d3dDDFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDDFF_SparsDiagoDiago =
    Points2d3dDDFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDDFF_SparsDiagoSpars =
    Points2d3dDDFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDDFF_SparsDiagoSpaco =
    Points2d3dDDFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDDFF_SpacoDiagoDense =
    Points2d3dDDFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDDFF_SpacoDiagoDiago =
    Points2d3dDDFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDDFF_SpacoDiagoSpars =
    Points2d3dDDFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDDFF_SpacoDiagoSpaco =
    Points2d3dDDFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-DDFF", "[DoubleSection-DDFF]", PointsSectionDDFF_DenseDiagoDense, PointsSectionDDFF_DenseDiagoDiago,
                   PointsSectionDDFF_DenseDiagoSpars, PointsSectionDDFF_DenseDiagoSpaco, PointsSectionDDFF_DiagoDiagoDense,
                   PointsSectionDDFF_DiagoDiagoDiago, PointsSectionDDFF_DiagoDiagoSpars, PointsSectionDDFF_DiagoDiagoSpaco,
                   PointsSectionDDFF_SparsDiagoDense, PointsSectionDDFF_SparsDiagoDiago, PointsSectionDDFF_SparsDiagoSpars,
                   PointsSectionDDFF_SparsDiagoSpaco, PointsSectionDDFF_SpacoDiagoDense, PointsSectionDDFF_SpacoDiagoDiago,
                   PointsSectionDDFF_SpacoDiagoSpars, PointsSectionDDFF_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

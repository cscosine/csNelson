#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "testDoubleSectionSchurCommon.hpp"

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dFDFF : public Points2d3dBase<Points2d3dFDFF<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                                             Point2d::blockSize, csBlockMatrix::Dynamic, numPoints2d, numPoints3d> {

public:
  int parameterVSize(void) const override {
    return Point3d::blockSize;
  }
};

using PointsSectionFDFF_DenseDiagoDense =
    Points2d3dFDFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFDFF_DenseDiagoDiago =
    Points2d3dFDFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFDFF_DenseDiagoSpars =
    Points2d3dFDFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFDFF_DenseDiagoSpaco =
    Points2d3dFDFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFDFF_DiagoDiagoDense =
    Points2d3dFDFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFDFF_DiagoDiagoDiago =
    Points2d3dFDFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFDFF_DiagoDiagoSpars =
    Points2d3dFDFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFDFF_DiagoDiagoSpaco =
    Points2d3dFDFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFDFF_SparsDiagoDense =
    Points2d3dFDFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFDFF_SparsDiagoDiago =
    Points2d3dFDFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFDFF_SparsDiagoSpars =
    Points2d3dFDFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFDFF_SparsDiagoSpaco =
    Points2d3dFDFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFDFF_SpacoDiagoDense =
    Points2d3dFDFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFDFF_SpacoDiagoDiago =
    Points2d3dFDFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFDFF_SpacoDiagoSpars =
    Points2d3dFDFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFDFF_SpacoDiagoSpaco =
    Points2d3dFDFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-FDFF", "[DoubleSection-FDFF]", PointsSectionFDFF_DenseDiagoDense, PointsSectionFDFF_DenseDiagoDiago,
                   PointsSectionFDFF_DenseDiagoSpars, PointsSectionFDFF_DenseDiagoSpaco, PointsSectionFDFF_DiagoDiagoDense,
                   PointsSectionFDFF_DiagoDiagoDiago, PointsSectionFDFF_DiagoDiagoSpars, PointsSectionFDFF_DiagoDiagoSpaco,
                   PointsSectionFDFF_SparsDiagoDense, PointsSectionFDFF_SparsDiagoDiago, PointsSectionFDFF_SparsDiagoSpars,
                   PointsSectionFDFF_SparsDiagoSpaco, PointsSectionFDFF_SpacoDiagoDense, PointsSectionFDFF_SpacoDiagoDiago,
                   PointsSectionFDFF_SpacoDiagoSpars, PointsSectionFDFF_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

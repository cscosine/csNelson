#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "testDoubleSectionSchurCommon.hpp"

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dVFFD : public Points2d3dBase<Points2d3dVFFD<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                                             csBlockMatrix::Variable, Point3d::blockSize, numPoints2d, csBlockMatrix::Dynamic> {
  std::vector<int> u_sizes;

public:
  Points2d3dVFFD()
      : u_sizes(numPoints2d, Point2d::blockSize) {}
  const std::vector<int>& parameterUSize(void) const override {
    return u_sizes;
  }
  int numParametersV() const override {
    return numPoints3d;
  }
};

using PointsSectionVFFD_DenseDiagoDense =
    Points2d3dVFFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVFFD_DenseDiagoDiago =
    Points2d3dVFFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVFFD_DenseDiagoSpars =
    Points2d3dVFFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVFFD_DenseDiagoSpaco =
    Points2d3dVFFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVFFD_DiagoDiagoDense =
    Points2d3dVFFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVFFD_DiagoDiagoDiago =
    Points2d3dVFFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVFFD_DiagoDiagoSpars =
    Points2d3dVFFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVFFD_DiagoDiagoSpaco =
    Points2d3dVFFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVFFD_SparsDiagoDense =
    Points2d3dVFFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVFFD_SparsDiagoDiago =
    Points2d3dVFFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVFFD_SparsDiagoSpars =
    Points2d3dVFFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVFFD_SparsDiagoSpaco =
    Points2d3dVFFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVFFD_SpacoDiagoDense =
    Points2d3dVFFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVFFD_SpacoDiagoDiago =
    Points2d3dVFFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVFFD_SpacoDiagoSpars =
    Points2d3dVFFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVFFD_SpacoDiagoSpaco =
    Points2d3dVFFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-VFFD", "[DoubleSection-VFFD]", PointsSectionVFFD_DenseDiagoDense, PointsSectionVFFD_DenseDiagoDiago,
                   PointsSectionVFFD_DenseDiagoSpars, PointsSectionVFFD_DenseDiagoSpaco, PointsSectionVFFD_DiagoDiagoDense,
                   PointsSectionVFFD_DiagoDiagoDiago, PointsSectionVFFD_DiagoDiagoSpars, PointsSectionVFFD_DiagoDiagoSpaco,
                   PointsSectionVFFD_SparsDiagoDense, PointsSectionVFFD_SparsDiagoDiago, PointsSectionVFFD_SparsDiagoSpars,
                   PointsSectionVFFD_SparsDiagoSpaco, PointsSectionVFFD_SpacoDiagoDense, PointsSectionVFFD_SpacoDiagoDiago,
                   PointsSectionVFFD_SpacoDiagoSpars, PointsSectionVFFD_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

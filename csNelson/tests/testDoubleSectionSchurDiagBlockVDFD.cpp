#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "testDoubleSectionSchurCommon.hpp"

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dVDFD : public Points2d3dBase<Points2d3dVDFD<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                                             csBlockMatrix::Variable, csBlockMatrix::Dynamic, numPoints2d, csBlockMatrix::Dynamic> {
  std::vector<int> u_sizes;

public:
  Points2d3dVDFD()
      : u_sizes(numPoints2d, Point2d::blockSize) {}
  const std::vector<int>& parameterUSize(void) const override {
    return u_sizes;
  }

  int parameterVSize(void) const override {
    return Point3d::blockSize;
  }

  int numParametersV() const override {
    return numPoints3d;
  }
};

using PointsSectionVDFD_DenseDiagoDense =
    Points2d3dVDFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDFD_DenseDiagoDiago =
    Points2d3dVDFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDFD_DenseDiagoSpars =
    Points2d3dVDFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDFD_DenseDiagoSpaco =
    Points2d3dVDFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVDFD_DiagoDiagoDense =
    Points2d3dVDFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDFD_DiagoDiagoDiago =
    Points2d3dVDFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDFD_DiagoDiagoSpars =
    Points2d3dVDFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDFD_DiagoDiagoSpaco =
    Points2d3dVDFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVDFD_SparsDiagoDense =
    Points2d3dVDFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDFD_SparsDiagoDiago =
    Points2d3dVDFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDFD_SparsDiagoSpars =
    Points2d3dVDFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDFD_SparsDiagoSpaco =
    Points2d3dVDFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVDFD_SpacoDiagoDense =
    Points2d3dVDFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDFD_SpacoDiagoDiago =
    Points2d3dVDFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDFD_SpacoDiagoSpars =
    Points2d3dVDFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDFD_SpacoDiagoSpaco =
    Points2d3dVDFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-VDFD", "[DoubleSection-VDFD]", PointsSectionVDFD_DenseDiagoDense, PointsSectionVDFD_DenseDiagoDiago,
                   PointsSectionVDFD_DenseDiagoSpars, PointsSectionVDFD_DenseDiagoSpaco, PointsSectionVDFD_DiagoDiagoDense,
                   PointsSectionVDFD_DiagoDiagoDiago, PointsSectionVDFD_DiagoDiagoSpars, PointsSectionVDFD_DiagoDiagoSpaco,
                   PointsSectionVDFD_SparsDiagoDense, PointsSectionVDFD_SparsDiagoDiago, PointsSectionVDFD_SparsDiagoSpars,
                   PointsSectionVDFD_SparsDiagoSpaco, PointsSectionVDFD_SpacoDiagoDense, PointsSectionVDFD_SpacoDiagoDiago,
                   PointsSectionVDFD_SpacoDiagoSpars, PointsSectionVDFD_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

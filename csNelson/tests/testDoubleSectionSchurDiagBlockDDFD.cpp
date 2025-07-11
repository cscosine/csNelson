#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "testDoubleSectionSchurCommon.hpp"

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dDDFD : public Points2d3dBase<Points2d3dDDFD<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                                             csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, numPoints2d, csBlockMatrix::Dynamic> {

public:
  int parameterUSize(void) const override {
    return Point2d::blockSize;
  }

  int parameterVSize(void) const override {
    return Point3d::blockSize;
  }

  int numParametersV() const override {
    return numPoints3d;
  }
};

using PointsSectionDDFD_DenseDiagoDense =
    Points2d3dDDFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDDFD_DenseDiagoDiago =
    Points2d3dDDFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDDFD_DenseDiagoSpars =
    Points2d3dDDFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDDFD_DenseDiagoSpaco =
    Points2d3dDDFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDDFD_DiagoDiagoDense =
    Points2d3dDDFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDDFD_DiagoDiagoDiago =
    Points2d3dDDFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDDFD_DiagoDiagoSpars =
    Points2d3dDDFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDDFD_DiagoDiagoSpaco =
    Points2d3dDDFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDDFD_SparsDiagoDense =
    Points2d3dDDFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDDFD_SparsDiagoDiago =
    Points2d3dDDFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDDFD_SparsDiagoSpars =
    Points2d3dDDFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDDFD_SparsDiagoSpaco =
    Points2d3dDDFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDDFD_SpacoDiagoDense =
    Points2d3dDDFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDDFD_SpacoDiagoDiago =
    Points2d3dDDFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDDFD_SpacoDiagoSpars =
    Points2d3dDDFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDDFD_SpacoDiagoSpaco =
    Points2d3dDDFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-DDFD", "[DoubleSection-DDFD]", PointsSectionDDFD_DenseDiagoDense, PointsSectionDDFD_DenseDiagoDiago,
                   PointsSectionDDFD_DenseDiagoSpars, PointsSectionDDFD_DenseDiagoSpaco, PointsSectionDDFD_DiagoDiagoDense,
                   PointsSectionDDFD_DiagoDiagoDiago, PointsSectionDDFD_DiagoDiagoSpars, PointsSectionDDFD_DiagoDiagoSpaco,
                   PointsSectionDDFD_SparsDiagoDense, PointsSectionDDFD_SparsDiagoDiago, PointsSectionDDFD_SparsDiagoSpars,
                   PointsSectionDDFD_SparsDiagoSpaco, PointsSectionDDFD_SpacoDiagoDense, PointsSectionDDFD_SpacoDiagoDiago,
                   PointsSectionDDFD_SpacoDiagoSpars, PointsSectionDDFD_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

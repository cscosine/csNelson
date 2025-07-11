#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "testDoubleSectionSchurCommon.hpp"

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dFDFD : public Points2d3dBase<Points2d3dFDFD<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                                             Point2d::blockSize, csBlockMatrix::Dynamic, numPoints2d, csBlockMatrix::Dynamic> {

public:
  int parameterVSize(void) const override {
    return Point3d::blockSize;
  }

  int numParametersV() const override {
    return numPoints3d;
  }
};

using PointsSectionFDFD_DenseDiagoDense =
    Points2d3dFDFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFDFD_DenseDiagoDiago =
    Points2d3dFDFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFDFD_DenseDiagoSpars =
    Points2d3dFDFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFDFD_DenseDiagoSpaco =
    Points2d3dFDFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFDFD_DiagoDiagoDense =
    Points2d3dFDFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFDFD_DiagoDiagoDiago =
    Points2d3dFDFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFDFD_DiagoDiagoSpars =
    Points2d3dFDFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFDFD_DiagoDiagoSpaco =
    Points2d3dFDFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFDFD_SparsDiagoDense =
    Points2d3dFDFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFDFD_SparsDiagoDiago =
    Points2d3dFDFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFDFD_SparsDiagoSpars =
    Points2d3dFDFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFDFD_SparsDiagoSpaco =
    Points2d3dFDFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFDFD_SpacoDiagoDense =
    Points2d3dFDFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFDFD_SpacoDiagoDiago =
    Points2d3dFDFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFDFD_SpacoDiagoSpars =
    Points2d3dFDFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFDFD_SpacoDiagoSpaco =
    Points2d3dFDFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-FDFD", "[DoubleSection-FDFD]", PointsSectionFDFD_DenseDiagoDense, PointsSectionFDFD_DenseDiagoDiago,
                   PointsSectionFDFD_DenseDiagoSpars, PointsSectionFDFD_DenseDiagoSpaco, PointsSectionFDFD_DiagoDiagoDense,
                   PointsSectionFDFD_DiagoDiagoDiago, PointsSectionFDFD_DiagoDiagoSpars, PointsSectionFDFD_DiagoDiagoSpaco,
                   PointsSectionFDFD_SparsDiagoDense, PointsSectionFDFD_SparsDiagoDiago, PointsSectionFDFD_SparsDiagoSpars,
                   PointsSectionFDFD_SparsDiagoSpaco, PointsSectionFDFD_SpacoDiagoDense, PointsSectionFDFD_SpacoDiagoDiago,
                   PointsSectionFDFD_SpacoDiagoSpars, PointsSectionFDFD_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

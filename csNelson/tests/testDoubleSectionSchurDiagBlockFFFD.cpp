#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/testDoubleSectionSchurCommon_impl.hpp>

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dFFFD : public Points2d3dBase<Points2d3dFFFD<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                                             Point2d::blockSize, Point3d::blockSize, numPoints2d, csBlockMatrix::Dynamic> {

public:
  int numParametersV() const override {
    return numPoints3d;
  }
};

using PointsSectionFFFD_DenseDiagoDense =
    Points2d3dFFFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFFD_DenseDiagoDiago =
    Points2d3dFFFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFFD_DenseDiagoSpars =
    Points2d3dFFFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFFD_DenseDiagoSpaco =
    Points2d3dFFFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFFFD_DiagoDiagoDense =
    Points2d3dFFFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFFD_DiagoDiagoDiago =
    Points2d3dFFFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFFD_DiagoDiagoSpars =
    Points2d3dFFFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFFD_DiagoDiagoSpaco =
    Points2d3dFFFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFFFD_SparsDiagoDense =
    Points2d3dFFFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFFD_SparsDiagoDiago =
    Points2d3dFFFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFFD_SparsDiagoSpars =
    Points2d3dFFFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFFD_SparsDiagoSpaco =
    Points2d3dFFFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFFFD_SpacoDiagoDense =
    Points2d3dFFFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFFD_SpacoDiagoDiago =
    Points2d3dFFFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFFD_SpacoDiagoSpars =
    Points2d3dFFFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFFD_SpacoDiagoSpaco =
    Points2d3dFFFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-FFFD", "[DoubleSection-FFFD]", PointsSectionFFFD_DenseDiagoDense, PointsSectionFFFD_DenseDiagoDiago,
                   PointsSectionFFFD_DenseDiagoSpars, PointsSectionFFFD_DenseDiagoSpaco, PointsSectionFFFD_DiagoDiagoDense,
                   PointsSectionFFFD_DiagoDiagoDiago, PointsSectionFFFD_DiagoDiagoSpars, PointsSectionFFFD_DiagoDiagoSpaco,
                   PointsSectionFFFD_SparsDiagoDense, PointsSectionFFFD_SparsDiagoDiago, PointsSectionFFFD_SparsDiagoSpars,
                   PointsSectionFFFD_SparsDiagoSpaco, PointsSectionFFFD_SpacoDiagoDense, PointsSectionFFFD_SpacoDiagoDiago,
                   PointsSectionFFFD_SpacoDiagoSpars, PointsSectionFFFD_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

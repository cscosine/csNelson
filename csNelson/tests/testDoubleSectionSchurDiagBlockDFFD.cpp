#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/testDoubleSectionSchurCommon_impl.hpp>

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dDFFD
    : public Points2d3dBase<Points2d3dDFFD<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                            csBlockMatrix::Dynamic, Point3d::blockSize, numPoints2d, csBlockMatrix::Dynamic> {

public:
  int parameterUSize(void) const override { return Point2d::blockSize; }

  int numParametersV() const override { return numPoints3d; }
};

using PointsSectionDFFD_DenseDiagoDense =
    Points2d3dDFFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDFFD_DenseDiagoDiago =
    Points2d3dDFFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDFFD_DenseDiagoSpars =
    Points2d3dDFFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDFFD_DenseDiagoSpaco =
    Points2d3dDFFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDFFD_DiagoDiagoDense =
    Points2d3dDFFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDFFD_DiagoDiagoDiago =
    Points2d3dDFFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDFFD_DiagoDiagoSpars =
    Points2d3dDFFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDFFD_DiagoDiagoSpaco =
    Points2d3dDFFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDFFD_SparsDiagoDense =
    Points2d3dDFFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDFFD_SparsDiagoDiago =
    Points2d3dDFFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDFFD_SparsDiagoSpars =
    Points2d3dDFFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDFFD_SparsDiagoSpaco =
    Points2d3dDFFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDFFD_SpacoDiagoDense =
    Points2d3dDFFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDFFD_SpacoDiagoDiago =
    Points2d3dDFFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDFFD_SpacoDiagoSpars =
    Points2d3dDFFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDFFD_SpacoDiagoSpaco =
    Points2d3dDFFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-DFFD", "[DoubleSection-DFFD]", PointsSectionDFFD_DenseDiagoDense,
                   PointsSectionDFFD_DenseDiagoDiago, PointsSectionDFFD_DenseDiagoSpars,
                   PointsSectionDFFD_DenseDiagoSpaco, PointsSectionDFFD_DiagoDiagoDense,
                   PointsSectionDFFD_DiagoDiagoDiago, PointsSectionDFFD_DiagoDiagoSpars,
                   PointsSectionDFFD_DiagoDiagoSpaco, PointsSectionDFFD_SparsDiagoDense,
                   PointsSectionDFFD_SparsDiagoDiago, PointsSectionDFFD_SparsDiagoSpars,
                   PointsSectionDFFD_SparsDiagoSpaco, PointsSectionDFFD_SpacoDiagoDense,
                   PointsSectionDFFD_SpacoDiagoDiago, PointsSectionDFFD_SpacoDiagoSpars,
                   PointsSectionDFFD_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

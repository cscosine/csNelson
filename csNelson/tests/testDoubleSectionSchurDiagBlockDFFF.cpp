#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/testDoubleSectionSchurCommon_impl.hpp>

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dDFFF : public Points2d3dBase<Points2d3dDFFF<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                                             csBlockMatrix::Dynamic, Point3d::blockSize, numPoints2d, numPoints3d> {

public:
  int parameterUSize(void) const override {
    return Point2d::blockSize;
  }
};

using PointsSectionDFFF_DenseDiagoDense =
    Points2d3dDFFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDFFF_DenseDiagoDiago =
    Points2d3dDFFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDFFF_DenseDiagoSpars =
    Points2d3dDFFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDFFF_DenseDiagoSpaco =
    Points2d3dDFFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDFFF_DiagoDiagoDense =
    Points2d3dDFFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDFFF_DiagoDiagoDiago =
    Points2d3dDFFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDFFF_DiagoDiagoSpars =
    Points2d3dDFFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDFFF_DiagoDiagoSpaco =
    Points2d3dDFFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDFFF_SparsDiagoDense =
    Points2d3dDFFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDFFF_SparsDiagoDiago =
    Points2d3dDFFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDFFF_SparsDiagoSpars =
    Points2d3dDFFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDFFF_SparsDiagoSpaco =
    Points2d3dDFFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDFFF_SpacoDiagoDense =
    Points2d3dDFFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDFFF_SpacoDiagoDiago =
    Points2d3dDFFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDFFF_SpacoDiagoSpars =
    Points2d3dDFFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDFFF_SpacoDiagoSpaco =
    Points2d3dDFFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-DFFF", "[DoubleSection-DFFF]", PointsSectionDFFF_DenseDiagoDense, PointsSectionDFFF_DenseDiagoDiago,
                   PointsSectionDFFF_DenseDiagoSpars, PointsSectionDFFF_DenseDiagoSpaco, PointsSectionDFFF_DiagoDiagoDense,
                   PointsSectionDFFF_DiagoDiagoDiago, PointsSectionDFFF_DiagoDiagoSpars, PointsSectionDFFF_DiagoDiagoSpaco,
                   PointsSectionDFFF_SparsDiagoDense, PointsSectionDFFF_SparsDiagoDiago, PointsSectionDFFF_SparsDiagoSpars,
                   PointsSectionDFFF_SparsDiagoSpaco, PointsSectionDFFF_SpacoDiagoDense, PointsSectionDFFF_SpacoDiagoDiago,
                   PointsSectionDFFF_SpacoDiagoSpars, PointsSectionDFFF_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

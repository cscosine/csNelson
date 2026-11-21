#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/testDoubleSectionSchurCommon_impl.hpp>

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dVDDF
    : public Points2d3dBase<Points2d3dVDDF<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                            csBlockMatrix::Variable, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, numPoints3d> {

  std::vector<int> u_sizes;

public:
  Points2d3dVDDF() : u_sizes(numPoints2d, Point2d::blockSize) {}
  const std::vector<int> &parameterUSize(void) const override { return u_sizes; }
  int parameterVSize(void) const override { return Point3d::blockSize; }
};

using PointsSectionVDDF_DenseDiagoDense =
    Points2d3dVDDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDDF_DenseDiagoDiago =
    Points2d3dVDDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDDF_DenseDiagoSpars =
    Points2d3dVDDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDDF_DenseDiagoSpaco =
    Points2d3dVDDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVDDF_DiagoDiagoDense =
    Points2d3dVDDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDDF_DiagoDiagoDiago =
    Points2d3dVDDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDDF_DiagoDiagoSpars =
    Points2d3dVDDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDDF_DiagoDiagoSpaco =
    Points2d3dVDDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVDDF_SparsDiagoDense =
    Points2d3dVDDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDDF_SparsDiagoDiago =
    Points2d3dVDDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDDF_SparsDiagoSpars =
    Points2d3dVDDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDDF_SparsDiagoSpaco =
    Points2d3dVDDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVDDF_SpacoDiagoDense =
    Points2d3dVDDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDDF_SpacoDiagoDiago =
    Points2d3dVDDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDDF_SpacoDiagoSpars =
    Points2d3dVDDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDDF_SpacoDiagoSpaco =
    Points2d3dVDDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-VDDF", "[DoubleSection-VDDF]", PointsSectionVDDF_DenseDiagoDense,
                   PointsSectionVDDF_DenseDiagoDiago, PointsSectionVDDF_DenseDiagoSpars,
                   PointsSectionVDDF_DenseDiagoSpaco, PointsSectionVDDF_DiagoDiagoDense,
                   PointsSectionVDDF_DiagoDiagoDiago, PointsSectionVDDF_DiagoDiagoSpars,
                   PointsSectionVDDF_DiagoDiagoSpaco, PointsSectionVDDF_SparsDiagoDense,
                   PointsSectionVDDF_SparsDiagoDiago, PointsSectionVDDF_SparsDiagoSpars,
                   PointsSectionVDDF_SparsDiagoSpaco, PointsSectionVDDF_SpacoDiagoDense,
                   PointsSectionVDDF_SpacoDiagoDiago, PointsSectionVDDF_SpacoDiagoSpars,
                   PointsSectionVDDF_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

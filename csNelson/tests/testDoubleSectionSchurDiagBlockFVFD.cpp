#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/testDoubleSectionSchurCommon_impl.hpp>

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dFVFD
    : public Points2d3dBase<Points2d3dFVFD<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                            Point2d::blockSize, csBlockMatrix::Variable, numPoints2d, csBlockMatrix::Dynamic> {
  std::vector<int> v_sizes;

public:
  Points2d3dFVFD() : v_sizes(numPoints3d, Point3d::blockSize) {}

  const std::vector<int> &parameterVSize(void) const override { return v_sizes; }
};

using PointsSectionFVFD_DenseDiagoDense =
    Points2d3dFVFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFVFD_DenseDiagoDiago =
    Points2d3dFVFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFVFD_DenseDiagoSpars =
    Points2d3dFVFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFVFD_DenseDiagoSpaco =
    Points2d3dFVFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFVFD_DiagoDiagoDense =
    Points2d3dFVFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFVFD_DiagoDiagoDiago =
    Points2d3dFVFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFVFD_DiagoDiagoSpars =
    Points2d3dFVFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFVFD_DiagoDiagoSpaco =
    Points2d3dFVFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFVFD_SparsDiagoDense =
    Points2d3dFVFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFVFD_SparsDiagoDiago =
    Points2d3dFVFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFVFD_SparsDiagoSpars =
    Points2d3dFVFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFVFD_SparsDiagoSpaco =
    Points2d3dFVFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFVFD_SpacoDiagoDense =
    Points2d3dFVFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFVFD_SpacoDiagoDiago =
    Points2d3dFVFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFVFD_SpacoDiagoSpars =
    Points2d3dFVFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFVFD_SpacoDiagoSpaco =
    Points2d3dFVFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-FVFD", "[DoubleSection-FVFD]", PointsSectionFVFD_DenseDiagoDense,
                   PointsSectionFVFD_DenseDiagoDiago, PointsSectionFVFD_DenseDiagoSpars,
                   PointsSectionFVFD_DenseDiagoSpaco, PointsSectionFVFD_DiagoDiagoDense,
                   PointsSectionFVFD_DiagoDiagoDiago, PointsSectionFVFD_DiagoDiagoSpars,
                   PointsSectionFVFD_DiagoDiagoSpaco, PointsSectionFVFD_SparsDiagoDense,
                   PointsSectionFVFD_SparsDiagoDiago, PointsSectionFVFD_SparsDiagoSpars,
                   PointsSectionFVFD_SparsDiagoSpaco, PointsSectionFVFD_SpacoDiagoDense,
                   PointsSectionFVFD_SpacoDiagoDiago, PointsSectionFVFD_SpacoDiagoSpars,
                   PointsSectionFVFD_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/testDoubleSectionSchurCommon_impl.hpp>

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dVVFD
    : public Points2d3dBase<Points2d3dVVFD<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                            csBlockMatrix::Variable, csBlockMatrix::Variable, numPoints2d, csBlockMatrix::Dynamic> {
  std::vector<int> v_sizes;
  std::vector<int> u_sizes;

public:
  Points2d3dVVFD() : v_sizes(numPoints3d, Point3d::blockSize), u_sizes(numPoints2d, Point2d::blockSize) {}

  const std::vector<int> &parameterUSize(void) const override { return u_sizes; }

  const std::vector<int> &parameterVSize(void) const override { return v_sizes; }
};

using PointsSectionVVFD_DenseDiagoDense =
    Points2d3dVVFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVVFD_DenseDiagoDiago =
    Points2d3dVVFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVVFD_DenseDiagoSpars =
    Points2d3dVVFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVVFD_DenseDiagoSpaco =
    Points2d3dVVFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVVFD_DiagoDiagoDense =
    Points2d3dVVFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVVFD_DiagoDiagoDiago =
    Points2d3dVVFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVVFD_DiagoDiagoSpars =
    Points2d3dVVFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVVFD_DiagoDiagoSpaco =
    Points2d3dVVFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVVFD_SparsDiagoDense =
    Points2d3dVVFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVVFD_SparsDiagoDiago =
    Points2d3dVVFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVVFD_SparsDiagoSpars =
    Points2d3dVVFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVVFD_SparsDiagoSpaco =
    Points2d3dVVFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVVFD_SpacoDiagoDense =
    Points2d3dVVFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVVFD_SpacoDiagoDiago =
    Points2d3dVVFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVVFD_SpacoDiagoSpars =
    Points2d3dVVFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVVFD_SpacoDiagoSpaco =
    Points2d3dVVFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-VVFD", "[DoubleSection-VVFD]", PointsSectionVVFD_DenseDiagoDense,
                   PointsSectionVVFD_DenseDiagoDiago, PointsSectionVVFD_DenseDiagoSpars,
                   PointsSectionVVFD_DenseDiagoSpaco, PointsSectionVVFD_DiagoDiagoDense,
                   PointsSectionVVFD_DiagoDiagoDiago, PointsSectionVVFD_DiagoDiagoSpars,
                   PointsSectionVVFD_DiagoDiagoSpaco, PointsSectionVVFD_SparsDiagoDense,
                   PointsSectionVVFD_SparsDiagoDiago, PointsSectionVVFD_SparsDiagoSpars,
                   PointsSectionVVFD_SparsDiagoSpaco, PointsSectionVVFD_SpacoDiagoDense,
                   PointsSectionVVFD_SpacoDiagoDiago, PointsSectionVVFD_SpacoDiagoSpars,
                   PointsSectionVVFD_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

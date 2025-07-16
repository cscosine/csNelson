#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/testDoubleSectionSchurCommon_impl.hpp>

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dDVFD : public Points2d3dBase<Points2d3dDVFD<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                                             csBlockMatrix::Dynamic, csBlockMatrix::Variable, numPoints2d, csBlockMatrix::Dynamic> {
  std::vector<int> v_sizes;

public:
  Points2d3dDVFD()
      : v_sizes(numPoints3d, Point3d::blockSize) {}

  int parameterUSize(void) const override {
    return Point2d::blockSize;
  }

  const std::vector<int>& parameterVSize(void) const override {
    return v_sizes;
  }
};

using PointsSectionDVFD_DenseDiagoDense =
    Points2d3dDVFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDVFD_DenseDiagoDiago =
    Points2d3dDVFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDVFD_DenseDiagoSpars =
    Points2d3dDVFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDVFD_DenseDiagoSpaco =
    Points2d3dDVFD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDVFD_DiagoDiagoDense =
    Points2d3dDVFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDVFD_DiagoDiagoDiago =
    Points2d3dDVFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDVFD_DiagoDiagoSpars =
    Points2d3dDVFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDVFD_DiagoDiagoSpaco =
    Points2d3dDVFD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDVFD_SparsDiagoDense =
    Points2d3dDVFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDVFD_SparsDiagoDiago =
    Points2d3dDVFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDVFD_SparsDiagoSpars =
    Points2d3dDVFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDVFD_SparsDiagoSpaco =
    Points2d3dDVFD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDVFD_SpacoDiagoDense =
    Points2d3dDVFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDVFD_SpacoDiagoDiago =
    Points2d3dDVFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDVFD_SpacoDiagoSpars =
    Points2d3dDVFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDVFD_SpacoDiagoSpaco =
    Points2d3dDVFD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-DVFD", "[DoubleSection-DVFD]", PointsSectionDVFD_DenseDiagoDense, PointsSectionDVFD_DenseDiagoDiago,
                   PointsSectionDVFD_DenseDiagoSpars, PointsSectionDVFD_DenseDiagoSpaco, PointsSectionDVFD_DiagoDiagoDense,
                   PointsSectionDVFD_DiagoDiagoDiago, PointsSectionDVFD_DiagoDiagoSpars, PointsSectionDVFD_DiagoDiagoSpaco,
                   PointsSectionDVFD_SparsDiagoDense, PointsSectionDVFD_SparsDiagoDiago, PointsSectionDVFD_SparsDiagoSpars,
                   PointsSectionDVFD_SparsDiagoSpaco, PointsSectionDVFD_SpacoDiagoDense, PointsSectionDVFD_SpacoDiagoDiago,
                   PointsSectionDVFD_SpacoDiagoSpars, PointsSectionDVFD_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

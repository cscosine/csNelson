#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/testDoubleSectionSchurCommon_impl.hpp>

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dVDFF
    : public Points2d3dBase<Points2d3dVDFF<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                            csBlockMatrix::Variable, csBlockMatrix::Dynamic, numPoints2d, numPoints3d> {

  std::vector<int> u_sizes;

public:
  Points2d3dVDFF() : u_sizes(numPoints2d, Point2d::blockSize) {}
  const std::vector<int> &parameterUSize(void) const override { return u_sizes; }

  int parameterVSize(void) const override { return Point3d::blockSize; }
};

using PointsSectionVDFF_DenseDiagoDense =
    Points2d3dVDFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDFF_DenseDiagoDiago =
    Points2d3dVDFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDFF_DenseDiagoSpars =
    Points2d3dVDFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDFF_DenseDiagoSpaco =
    Points2d3dVDFF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVDFF_DiagoDiagoDense =
    Points2d3dVDFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDFF_DiagoDiagoDiago =
    Points2d3dVDFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDFF_DiagoDiagoSpars =
    Points2d3dVDFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDFF_DiagoDiagoSpaco =
    Points2d3dVDFF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVDFF_SparsDiagoDense =
    Points2d3dVDFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDFF_SparsDiagoDiago =
    Points2d3dVDFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDFF_SparsDiagoSpars =
    Points2d3dVDFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDFF_SparsDiagoSpaco =
    Points2d3dVDFF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVDFF_SpacoDiagoDense =
    Points2d3dVDFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDFF_SpacoDiagoDiago =
    Points2d3dVDFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDFF_SpacoDiagoSpars =
    Points2d3dVDFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDFF_SpacoDiagoSpaco =
    Points2d3dVDFF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-VDFF", "[DoubleSection-VDFF]", PointsSectionVDFF_DenseDiagoDense,
                   PointsSectionVDFF_DenseDiagoDiago, PointsSectionVDFF_DenseDiagoSpars,
                   PointsSectionVDFF_DenseDiagoSpaco, PointsSectionVDFF_DiagoDiagoDense,
                   PointsSectionVDFF_DiagoDiagoDiago, PointsSectionVDFF_DiagoDiagoSpars,
                   PointsSectionVDFF_DiagoDiagoSpaco, PointsSectionVDFF_SparsDiagoDense,
                   PointsSectionVDFF_SparsDiagoDiago, PointsSectionVDFF_SparsDiagoSpars,
                   PointsSectionVDFF_SparsDiagoSpaco, PointsSectionVDFF_SpacoDiagoDense,
                   PointsSectionVDFF_SpacoDiagoDiago, PointsSectionVDFF_SpacoDiagoSpars,
                   PointsSectionVDFF_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

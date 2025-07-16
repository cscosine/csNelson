#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/testDoubleSectionSchurCommon_impl.hpp>

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dVDDD
    : public Points2d3dBase<Points2d3dVDDD<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv, csBlockMatrix::Variable,
                            csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic> {

  std::vector<int> u_sizes;

public:
  Points2d3dVDDD()
      : u_sizes(numPoints2d, Point2d::blockSize) {}
  const std::vector<int>& parameterUSize(void) const override {
    return u_sizes;
  }
  int parameterVSize(void) const override {
    return Point3d::blockSize;
  }

  int numParametersV() const override {
    return numPoints3d;
  }
};

using PointsSectionVDDD_DenseDiagoDense =
    Points2d3dVDDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDDD_DenseDiagoDiago =
    Points2d3dVDDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDDD_DenseDiagoSpars =
    Points2d3dVDDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDDD_DenseDiagoSpaco =
    Points2d3dVDDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVDDD_DiagoDiagoDense =
    Points2d3dVDDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDDD_DiagoDiagoDiago =
    Points2d3dVDDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDDD_DiagoDiagoSpars =
    Points2d3dVDDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDDD_DiagoDiagoSpaco =
    Points2d3dVDDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVDDD_SparsDiagoDense =
    Points2d3dVDDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDDD_SparsDiagoDiago =
    Points2d3dVDDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDDD_SparsDiagoSpars =
    Points2d3dVDDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDDD_SparsDiagoSpaco =
    Points2d3dVDDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionVDDD_SpacoDiagoDense =
    Points2d3dVDDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionVDDD_SpacoDiagoDiago =
    Points2d3dVDDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionVDDD_SpacoDiagoSpars =
    Points2d3dVDDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionVDDD_SpacoDiagoSpaco =
    Points2d3dVDDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-VDDD", "[DoubleSection-VDDD]", PointsSectionVDDD_DenseDiagoDense, PointsSectionVDDD_DenseDiagoDiago,
                   PointsSectionVDDD_DenseDiagoSpars, PointsSectionVDDD_DenseDiagoSpaco, PointsSectionVDDD_DiagoDiagoDense,
                   PointsSectionVDDD_DiagoDiagoDiago, PointsSectionVDDD_DiagoDiagoSpars, PointsSectionVDDD_DiagoDiagoSpaco,
                   PointsSectionVDDD_SparsDiagoDense, PointsSectionVDDD_SparsDiagoDiago, PointsSectionVDDD_SparsDiagoSpars,
                   PointsSectionVDDD_SparsDiagoSpaco, PointsSectionVDDD_SpacoDiagoDense, PointsSectionVDDD_SpacoDiagoDiago,
                   PointsSectionVDDD_SpacoDiagoSpars, PointsSectionVDDD_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

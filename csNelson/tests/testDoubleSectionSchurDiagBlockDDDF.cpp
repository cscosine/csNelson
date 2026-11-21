#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/testDoubleSectionSchurCommon_impl.hpp>

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dDDDF
    : public Points2d3dBase<Points2d3dDDDF<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                            csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, numPoints3d> {

public:
  int parameterVSize(void) const override { return Point3d::blockSize; }

  int parameterUSize(void) const override { return Point2d::blockSize; }

  int numParametersU() const override { return numPoints2d; }
};

using PointsSectionDDDF_DenseDiagoDense =
    Points2d3dDDDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDDDF_DenseDiagoDiago =
    Points2d3dDDDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDDDF_DenseDiagoSpars =
    Points2d3dDDDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDDDF_DenseDiagoSpaco =
    Points2d3dDDDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDDDF_DiagoDiagoDense =
    Points2d3dDDDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDDDF_DiagoDiagoDiago =
    Points2d3dDDDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDDDF_DiagoDiagoSpars =
    Points2d3dDDDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDDDF_DiagoDiagoSpaco =
    Points2d3dDDDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDDDF_SparsDiagoDense =
    Points2d3dDDDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDDDF_SparsDiagoDiago =
    Points2d3dDDDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDDDF_SparsDiagoSpars =
    Points2d3dDDDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDDDF_SparsDiagoSpaco =
    Points2d3dDDDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDDDF_SpacoDiagoDense =
    Points2d3dDDDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDDDF_SpacoDiagoDiago =
    Points2d3dDDDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDDDF_SpacoDiagoSpars =
    Points2d3dDDDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDDDF_SpacoDiagoSpaco =
    Points2d3dDDDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-DDDF", "[DoubleSection-DDDF]", PointsSectionDDDF_DenseDiagoDense,
                   PointsSectionDDDF_DenseDiagoDiago, PointsSectionDDDF_DenseDiagoSpars,
                   PointsSectionDDDF_DenseDiagoSpaco, PointsSectionDDDF_DiagoDiagoDense,
                   PointsSectionDDDF_DiagoDiagoDiago, PointsSectionDDDF_DiagoDiagoSpars,
                   PointsSectionDDDF_DiagoDiagoSpaco, PointsSectionDDDF_SparsDiagoDense,
                   PointsSectionDDDF_SparsDiagoDiago, PointsSectionDDDF_SparsDiagoSpars,
                   PointsSectionDDDF_SparsDiagoSpaco, PointsSectionDDDF_SpacoDiagoDense,
                   PointsSectionDDDF_SpacoDiagoDiago, PointsSectionDDDF_SpacoDiagoSpars,
                   PointsSectionDDDF_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

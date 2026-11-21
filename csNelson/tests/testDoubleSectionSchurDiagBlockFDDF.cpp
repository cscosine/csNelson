#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/testDoubleSectionSchurCommon_impl.hpp>

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dFDDF
    : public Points2d3dBase<Points2d3dFDDF<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                            Point2d::blockSize, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, numPoints3d> {

public:
  int parameterVSize(void) const override { return Point3d::blockSize; }

  int numParametersU() const override { return numPoints2d; }
};

using PointsSectionFDDF_DenseDiagoDense =
    Points2d3dFDDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFDDF_DenseDiagoDiago =
    Points2d3dFDDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFDDF_DenseDiagoSpars =
    Points2d3dFDDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFDDF_DenseDiagoSpaco =
    Points2d3dFDDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFDDF_DiagoDiagoDense =
    Points2d3dFDDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFDDF_DiagoDiagoDiago =
    Points2d3dFDDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFDDF_DiagoDiagoSpars =
    Points2d3dFDDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFDDF_DiagoDiagoSpaco =
    Points2d3dFDDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFDDF_SparsDiagoDense =
    Points2d3dFDDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFDDF_SparsDiagoDiago =
    Points2d3dFDDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFDDF_SparsDiagoSpars =
    Points2d3dFDDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFDDF_SparsDiagoSpaco =
    Points2d3dFDDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFDDF_SpacoDiagoDense =
    Points2d3dFDDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFDDF_SpacoDiagoDiago =
    Points2d3dFDDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFDDF_SpacoDiagoSpars =
    Points2d3dFDDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFDDF_SpacoDiagoSpaco =
    Points2d3dFDDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-FDDF", "[DoubleSection-FDDF]", PointsSectionFDDF_DenseDiagoDense,
                   PointsSectionFDDF_DenseDiagoDiago, PointsSectionFDDF_DenseDiagoSpars,
                   PointsSectionFDDF_DenseDiagoSpaco, PointsSectionFDDF_DiagoDiagoDense,
                   PointsSectionFDDF_DiagoDiagoDiago, PointsSectionFDDF_DiagoDiagoSpars,
                   PointsSectionFDDF_DiagoDiagoSpaco, PointsSectionFDDF_SparsDiagoDense,
                   PointsSectionFDDF_SparsDiagoDiago, PointsSectionFDDF_SparsDiagoSpars,
                   PointsSectionFDDF_SparsDiagoSpaco, PointsSectionFDDF_SpacoDiagoDense,
                   PointsSectionFDDF_SpacoDiagoDiago, PointsSectionFDDF_SpacoDiagoSpars,
                   PointsSectionFDDF_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/testDoubleSectionSchurCommon_impl.hpp>

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dFFDF
    : public Points2d3dBase<Points2d3dFFDF<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                            Point2d::blockSize, Point3d::blockSize, csBlockMatrix::Dynamic, numPoints3d> {

public:
  int numParametersU() const override { return numPoints2d; }
};

using PointsSectionFFDF_DenseDiagoDense =
    Points2d3dFFDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFDF_DenseDiagoDiago =
    Points2d3dFFDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFDF_DenseDiagoSpars =
    Points2d3dFFDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFDF_DenseDiagoSpaco =
    Points2d3dFFDF<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFFDF_DiagoDiagoDense =
    Points2d3dFFDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFDF_DiagoDiagoDiago =
    Points2d3dFFDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFDF_DiagoDiagoSpars =
    Points2d3dFFDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFDF_DiagoDiagoSpaco =
    Points2d3dFFDF<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFFDF_SparsDiagoDense =
    Points2d3dFFDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFDF_SparsDiagoDiago =
    Points2d3dFFDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFDF_SparsDiagoSpars =
    Points2d3dFFDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFDF_SparsDiagoSpaco =
    Points2d3dFFDF<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFFDF_SpacoDiagoDense =
    Points2d3dFFDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFDF_SpacoDiagoDiago =
    Points2d3dFFDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFDF_SpacoDiagoSpars =
    Points2d3dFFDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFDF_SpacoDiagoSpaco =
    Points2d3dFFDF<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-FFDF", "[DoubleSection-FFDF]", PointsSectionFFDF_DenseDiagoDense,
                   PointsSectionFFDF_DenseDiagoDiago, PointsSectionFFDF_DenseDiagoSpars,
                   PointsSectionFFDF_DenseDiagoSpaco, PointsSectionFFDF_DiagoDiagoDense,
                   PointsSectionFFDF_DiagoDiagoDiago, PointsSectionFFDF_DiagoDiagoSpars,
                   PointsSectionFFDF_DiagoDiagoSpaco, PointsSectionFFDF_SparsDiagoDense,
                   PointsSectionFFDF_SparsDiagoDiago, PointsSectionFFDF_SparsDiagoSpars,
                   PointsSectionFFDF_SparsDiagoSpaco, PointsSectionFFDF_SpacoDiagoDense,
                   PointsSectionFFDF_SpacoDiagoDiago, PointsSectionFFDF_SpacoDiagoSpars,
                   PointsSectionFFDF_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/testDoubleSectionSchurCommon_impl.hpp>

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dFFDD : public Points2d3dBase<Points2d3dFFDD<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv,
                                             Point2d::blockSize, Point3d::blockSize, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic> {

public:
  int numParametersV() const override {
    return numPoints3d;
  }

  int numParametersU() const override {
    return numPoints2d;
  }
};

using PointsSectionFFDD_DenseDiagoDense =
    Points2d3dFFDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFDD_DenseDiagoDiago =
    Points2d3dFFDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFDD_DenseDiagoSpars =
    Points2d3dFFDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFDD_DenseDiagoSpaco =
    Points2d3dFFDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFFDD_DiagoDiagoDense =
    Points2d3dFFDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFDD_DiagoDiagoDiago =
    Points2d3dFFDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFDD_DiagoDiagoSpars =
    Points2d3dFFDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFDD_DiagoDiagoSpaco =
    Points2d3dFFDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFFDD_SparsDiagoDense =
    Points2d3dFFDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFDD_SparsDiagoDiago =
    Points2d3dFFDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFDD_SparsDiagoSpars =
    Points2d3dFFDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFDD_SparsDiagoSpaco =
    Points2d3dFFDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFFDD_SpacoDiagoDense =
    Points2d3dFFDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFDD_SpacoDiagoDiago =
    Points2d3dFFDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFDD_SpacoDiagoSpars =
    Points2d3dFFDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFDD_SpacoDiagoSpaco =
    Points2d3dFFDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

using PointsSectionFFDD_DenseDiagoDense =
    Points2d3dFFDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFDD_DenseDiagoDiago =
    Points2d3dFFDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFDD_DenseDiagoSpars =
    Points2d3dFFDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFDD_DenseDiagoSpaco =
    Points2d3dFFDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFFDD_DiagoDiagoDense =
    Points2d3dFFDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFDD_DiagoDiagoDiago =
    Points2d3dFFDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFDD_DiagoDiagoSpars =
    Points2d3dFFDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFDD_DiagoDiagoSpaco =
    Points2d3dFFDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFFDD_SparsDiagoDense =
    Points2d3dFFDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFDD_SparsDiagoDiago =
    Points2d3dFFDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFDD_SparsDiagoSpars =
    Points2d3dFFDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFDD_SparsDiagoSpaco =
    Points2d3dFFDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionFFDD_SpacoDiagoDense =
    Points2d3dFFDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionFFDD_SpacoDiagoDiago =
    Points2d3dFFDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionFFDD_SpacoDiagoSpars =
    Points2d3dFFDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionFFDD_SpacoDiagoSpaco =
    Points2d3dFFDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-FFDD", "[DoubleSection-FFDD]", PointsSectionFFDD_DenseDiagoDense, PointsSectionFFDD_DenseDiagoDiago,
                   PointsSectionFFDD_DenseDiagoSpars, PointsSectionFFDD_DenseDiagoSpaco, PointsSectionFFDD_DiagoDiagoDense,
                   PointsSectionFFDD_DiagoDiagoDiago, PointsSectionFFDD_DiagoDiagoSpars, PointsSectionFFDD_DiagoDiagoSpaco,
                   PointsSectionFFDD_SparsDiagoDense, PointsSectionFFDD_SparsDiagoDiago, PointsSectionFFDD_SparsDiagoSpars,
                   PointsSectionFFDD_SparsDiagoSpaco, PointsSectionFFDD_SpacoDiagoDense, PointsSectionFFDD_SpacoDiagoDiago,
                   PointsSectionFFDD_SpacoDiagoSpars, PointsSectionFFDD_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

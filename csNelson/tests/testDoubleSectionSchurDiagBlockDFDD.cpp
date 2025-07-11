#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "testDoubleSectionSchurCommon.hpp"

template <int matTypeUv, int matTypeVv, int matTypeWv>
class Points2d3dDFDD
    : public Points2d3dBase<Points2d3dDFDD<matTypeUv, matTypeVv, matTypeWv>, matTypeUv, matTypeVv, matTypeWv, csBlockMatrix::Dynamic,
                            Point3d::blockSize, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic> {

public:
  int parameterUSize(void) const override {
    return Point2d::blockSize;
  }

  int numParametersV() const override {
    return numPoints3d;
  }

  int numParametersU() const override {
    return numPoints2d;
  }
};

using PointsSectionDFDD_DenseDiagoDense =
    Points2d3dDFDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDFDD_DenseDiagoDiago =
    Points2d3dDFDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDFDD_DenseDiagoSpars =
    Points2d3dDFDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDFDD_DenseDiagoSpaco =
    Points2d3dDFDD<csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDFDD_DiagoDiagoDense =
    Points2d3dDFDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDFDD_DiagoDiagoDiago =
    Points2d3dDFDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDFDD_DiagoDiagoSpars =
    Points2d3dDFDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDFDD_DiagoDiagoSpaco =
    Points2d3dDFDD<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDFDD_SparsDiagoDense =
    Points2d3dDFDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDFDD_SparsDiagoDiago =
    Points2d3dDFDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDFDD_SparsDiagoSpars =
    Points2d3dDFDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDFDD_SparsDiagoSpaco =
    Points2d3dDFDD<csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;
using PointsSectionDFDD_SpacoDiagoDense =
    Points2d3dDFDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense>;
using PointsSectionDFDD_SpacoDiagoDiago =
    Points2d3dDFDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal>;
using PointsSectionDFDD_SpacoDiagoSpars =
    Points2d3dDFDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>;
using PointsSectionDFDD_SpacoDiagoSpaco =
    Points2d3dDFDD<csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>;

TEMPLATE_TEST_CASE("DoubleSection-DFDD", "[DoubleSection-DFDD]", PointsSectionDFDD_DenseDiagoDense, PointsSectionDFDD_DenseDiagoDiago,
                   PointsSectionDFDD_DenseDiagoSpars, PointsSectionDFDD_DenseDiagoSpaco, PointsSectionDFDD_DiagoDiagoDense,
                   PointsSectionDFDD_DiagoDiagoDiago, PointsSectionDFDD_DiagoDiagoSpars, PointsSectionDFDD_DiagoDiagoSpaco,
                   PointsSectionDFDD_SparsDiagoDense, PointsSectionDFDD_SparsDiagoDiago, PointsSectionDFDD_SparsDiagoSpars,
                   PointsSectionDFDD_SparsDiagoSpaco, PointsSectionDFDD_SpacoDiagoDense, PointsSectionDFDD_SpacoDiagoDiago,
                   PointsSectionDFDD_SpacoDiagoSpars, PointsSectionDFDD_SpacoDiagoSpaco) {

  testFunction<TestType, csNelson::GaussNewton>();
  testFunction<TestType, csNelson::LevenbergMarquardt>();
}

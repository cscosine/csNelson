#pragma once
#include "testDoubleSectionCommon.h"

#include "csNelson/GaussNewton.hpp"
#include "csNelson/LevenbergMarquardt.hpp"

template <class TestType, template <typename> class Solver>
void testFunction() {
  std::cout << "-------------------------------------------------------" << std::endl;
  TestType pss;
  pss.parametersReady();
  REQUIRE(pss.numParametersU() == numPoints2d);
  REQUIRE(pss.numParametersV() == numPoints3d);

  // unary edge first section
  for (int i = 0; i < numPoints2d; i++) {
    pss.addEdge(i, new EdgeUnaryPoint2d<TestType>(i, pss.parameterU(i).p2d));
  }
  // unary edge second section
  for (int i = 0; i < numPoints3d; i++) {
    pss.addEdge(i, new EdgeUnaryPoint3d<TestType>(i, pss.parameterV(i).p3d));
  }
  // binary edge first section
  if (pss.matTypeU() != csBlockMatrix::BlockDiagonal && pss.matTypeU() != csBlockMatrix::SparseCoeffBlockDiagonal) {
    for (int i = 0; i < numPoints2d; i++) {
      for (int j = i + 1; j < numPoints2d; j++) {
        pss.addEdge(i, j, new EdgeBinaryPoint2d<TestType>(i, j, pss.parameterU(i).p2d - pss.parameterU(j).p2d));
      }
    }
  }
  // binary edge section section
  if (pss.matTypeV() != csBlockMatrix::BlockDiagonal && pss.matTypeV() != csBlockMatrix::SparseCoeffBlockDiagonal) {
    for (int i = 0; i < numPoints3d; i++) {
      for (int j = i + 1; j < numPoints3d; j++) {
        pss.addEdge(i, j, new EdgeBinaryPoint3d<TestType>(i, j, pss.parameterV(i).p3d - pss.parameterV(j).p3d));
      }
    }
  }
  // binary edge first section to second section
  if (pss.matTypeW() != csBlockMatrix::BlockDiagonal && pss.matTypeW() != csBlockMatrix::SparseCoeffBlockDiagonal) {
    for (int i = 0; i < numPoints2d; i++) {
      for (int j = 0; j < numPoints3d; j++) {
        pss.addEdge(i, j, new EdgeBinaryPoint2d3d<TestType>(i, j, pss.parameterU(i).p2d - pss.parameterV(j).p3d.template head<2>()));
      }
    }
  } else {
    for (int i = 0; i < std::min(numPoints2d, numPoints3d); i++) {
      pss.addEdge(i, i, new EdgeBinaryPoint2d3d<TestType>(i, i, pss.parameterU(i).p2d - pss.parameterV(i).p3d.template head<2>()));
    }
  }

  pss.structureReady();

  pss.update(true);

  {
    using SolverAlgorithm = Solver<typename csNelson::SolverTraits<csNelson::solverCholeskySchur>::Solver<
        typename TestType::Hessian::Traits, csNelson::solverCholeskyDense, csNelson::matrixWrapperDense, csNelson::matrixWrapperDense,
        csNelson::choleskyAMDOrdering, csNelson::choleskyAMDOrdering>>;
    SolverAlgorithm gn;
    pss.addNoise(0.5);
    auto tc = gn.solve(pss);
    std::cout << SolverAlgorithm::Utils::toString(tc) << std::endl;
    std::cout << "stats " << gn.stats().toString() << std::endl;
    REQUIRE(pss.hessian().chi2() < Eigen::NumTraits<double>::dummy_precision());
  }

  {
    using SolverAlgorithm = Solver<typename csNelson::SolverTraits<csNelson::solverCholeskySchur>::Solver<
        typename TestType::Hessian::Traits, csNelson::solverCholeskyDense, csNelson::matrixWrapperDense, csNelson::matrixWrapperSparse,
        csNelson::choleskyAMDOrdering, csNelson::choleskyAMDOrdering>>;
    SolverAlgorithm gn;
    pss.addNoise(0.5);
    auto tc = gn.solve(pss);
    std::cout << SolverAlgorithm::Utils::toString(tc) << std::endl;
    std::cout << "stats " << gn.stats().toString() << std::endl;
    REQUIRE(pss.hessian().chi2() < Eigen::NumTraits<double>::dummy_precision());
  }
  {
    using SolverAlgorithm = Solver<typename csNelson::SolverTraits<csNelson::solverCholeskySchur>::Solver<
        typename TestType::Hessian::Traits, csNelson::solverCholeskyDense, csNelson::matrixWrapperSparse, csNelson::matrixWrapperDense,
        csNelson::choleskyAMDOrdering, csNelson::choleskyAMDOrdering>>;
    SolverAlgorithm gn;
    pss.addNoise(0.5);
    auto tc = gn.solve(pss);
    std::cout << SolverAlgorithm::Utils::toString(tc) << std::endl;
    std::cout << "stats " << gn.stats().toString() << std::endl;
    REQUIRE(pss.hessian().chi2() < Eigen::NumTraits<double>::dummy_precision());
  }
  {
    using SolverAlgorithm = Solver<typename csNelson::SolverTraits<csNelson::solverCholeskySchur>::Solver<
        typename TestType::Hessian::Traits, csNelson::solverCholeskyDense, csNelson::matrixWrapperSparse,
        csNelson::matrixWrapperSparse, csNelson::choleskyAMDOrdering, csNelson::choleskyAMDOrdering>>;
    SolverAlgorithm gn;
    pss.addNoise(0.5);
    auto tc = gn.solve(pss);
    std::cout << SolverAlgorithm::Utils::toString(tc) << std::endl;
    std::cout << "stats " << gn.stats().toString() << std::endl;
    REQUIRE(pss.hessian().chi2() < Eigen::NumTraits<double>::dummy_precision());
  }
  {
    using SolverAlgorithm = Solver<typename csNelson::SolverTraits<csNelson::solverCholeskySchur>::Solver<
        typename TestType::Hessian::Traits, csNelson::solverCholeskySparse, csNelson::matrixWrapperDense, csNelson::matrixWrapperDense,
        csNelson::choleskyAMDOrdering, csNelson::choleskyAMDOrdering>>;
    SolverAlgorithm gn;
    pss.addNoise(0.5);
    auto tc = gn.solve(pss);
    std::cout << SolverAlgorithm::Utils::toString(tc) << std::endl;
    std::cout << "stats " << gn.stats().toString() << std::endl;
    REQUIRE(pss.hessian().chi2() < Eigen::NumTraits<double>::dummy_precision());
  }
  {
    using SolverAlgorithm = Solver<typename csNelson::SolverTraits<csNelson::solverCholeskySchur>::Solver<
        typename TestType::Hessian::Traits, csNelson::solverCholeskySparse, csNelson::matrixWrapperDense,
        csNelson::matrixWrapperSparse, csNelson::choleskyAMDOrdering, csNelson::choleskyAMDOrdering>>;
    SolverAlgorithm gn;
    pss.addNoise(0.5);
    auto tc = gn.solve(pss);
    std::cout << SolverAlgorithm::Utils::toString(tc) << std::endl;
    std::cout << "stats " << gn.stats().toString() << std::endl;
    REQUIRE(pss.hessian().chi2() < Eigen::NumTraits<double>::dummy_precision());
  }
  {
    using SolverAlgorithm = Solver<typename csNelson::SolverTraits<csNelson::solverCholeskySchur>::Solver<
        typename TestType::Hessian::Traits, csNelson::solverCholeskySparse, csNelson::matrixWrapperSparse,
        csNelson::matrixWrapperDense, csNelson::choleskyAMDOrdering, csNelson::choleskyAMDOrdering>>;
    SolverAlgorithm gn;
    pss.addNoise(0.5);
    auto tc = gn.solve(pss);
    std::cout << SolverAlgorithm::Utils::toString(tc) << std::endl;
    std::cout << "stats " << gn.stats().toString() << std::endl;
    REQUIRE(pss.hessian().chi2() < Eigen::NumTraits<double>::dummy_precision());
  }
  {
    using SolverAlgorithm = Solver<typename csNelson::SolverTraits<csNelson::solverCholeskySchur>::Solver<
        typename TestType::Hessian::Traits, csNelson::solverCholeskySparse, csNelson::matrixWrapperSparse,
        csNelson::matrixWrapperSparse, csNelson::choleskyAMDOrdering, csNelson::choleskyAMDOrdering>>;
    SolverAlgorithm gn;
    pss.addNoise(0.5);
    auto tc = gn.solve(pss);
    std::cout << SolverAlgorithm::Utils::toString(tc) << std::endl;
    std::cout << "stats " << gn.stats().toString() << std::endl;
    REQUIRE(pss.hessian().chi2() < Eigen::NumTraits<double>::dummy_precision());
  }
}

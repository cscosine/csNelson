#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "csBlockMatrix/DenseMatrixBlock.hpp"
#include "csBlockMatrix/SparsityPattern.h"

#include "csNelson/SingleSectionHessian.hpp"

static constexpr int secSizeFix = 3;

template <int secType>
class SecType {
public:
  static constexpr int secSize = secSizeFix;
};

template <>
class SecType<csBlockMatrix::Variable> {
public:
  static const std::vector<int> secSize;
};
const std::vector<int> SecType<csBlockMatrix::Variable>::secSize = {1, 2, 3, 4, 5}; // sum = 15

static constexpr int numBlocks = 5;

TEMPLATE_TEST_CASE_SIG(
    "SingleSectionHessian", "[SingleSectionHessian]", ((int matType, int B, int NB), matType, B, NB),
    (csBlockMatrix::BlockDense, secSizeFix, numBlocks), (csBlockMatrix::BlockDense, secSizeFix, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, numBlocks),
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDense, csBlockMatrix::Variable, numBlocks),
    (csBlockMatrix::BlockDense, csBlockMatrix::Variable, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, secSizeFix, numBlocks), (csBlockMatrix::BlockDiagonal, secSizeFix, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, numBlocks),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Variable, numBlocks),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Variable, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, secSizeFix, numBlocks), (csBlockMatrix::BlockSparse, secSizeFix, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, numBlocks),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Variable, numBlocks),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Variable, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, secSizeFix, numBlocks), (csBlockMatrix::BlockCoeffSparse, secSizeFix, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, numBlocks),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Variable, numBlocks),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Variable, csBlockMatrix::Dynamic)) {
  auto sp = std::make_shared<csBlockMatrix::SparsityPatternColMajor>(numBlocks, numBlocks);
  sp->setDiagonal();
  if (matType != csBlockMatrix::BlockDiagonal) {
    for (int i = 0; i < 10; i++) {
      int r = rand() % numBlocks;
      int c = rand() % numBlocks;
      sp->add(r, c);
    }
  }

  csNelson::SingleSectionHessian<matType, double, B, NB> sec;
  sec.resize(SecType<B>::secSize, numBlocks, sp);
  sec.clearAll();
  REQUIRE(sec.chi2() == 0);
}

TEST_CASE("temp", "temp") {
  // just compilation test

  csNelson::SingleSectionHessian<csBlockMatrix::BlockDense, double, secSizeFix, numBlocks> sec_01;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockDense, double, secSizeFix, csBlockMatrix::Dynamic> sec_02;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockDense, double, csBlockMatrix::Dynamic, secSizeFix> sec_03;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockDense, double, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic> sec_04;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockDense, double, csBlockMatrix::Variable, secSizeFix> sec_05;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockDense, double, csBlockMatrix::Variable, csBlockMatrix::Dynamic> sec_06;

  csNelson::SingleSectionHessian<csBlockMatrix::BlockDiagonal, double, secSizeFix, numBlocks> sec_11;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockDiagonal, double, secSizeFix, csBlockMatrix::Dynamic> sec_12;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockDiagonal, double, csBlockMatrix::Dynamic, secSizeFix> sec_13;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockDiagonal, double, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic> sec_14;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockDiagonal, double, csBlockMatrix::Variable, secSizeFix> sec_15;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockDiagonal, double, csBlockMatrix::Variable, csBlockMatrix::Dynamic> sec_16;

  csNelson::SingleSectionHessian<csBlockMatrix::BlockSparse, double, secSizeFix, numBlocks> sec_21;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockSparse, double, secSizeFix, csBlockMatrix::Dynamic> sec_22;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockSparse, double, csBlockMatrix::Dynamic, secSizeFix> sec_23;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockSparse, double, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic> sec_24;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockSparse, double, csBlockMatrix::Variable, secSizeFix> sec_25;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockSparse, double, csBlockMatrix::Variable, csBlockMatrix::Dynamic> sec_26;

  csNelson::SingleSectionHessian<csBlockMatrix::BlockCoeffSparse, double, secSizeFix, numBlocks> sec_31;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockCoeffSparse, double, secSizeFix, csBlockMatrix::Dynamic> sec_32;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockCoeffSparse, double, csBlockMatrix::Dynamic, secSizeFix> sec_33;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockCoeffSparse, double, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic> sec_34;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockCoeffSparse, double, csBlockMatrix::Variable, secSizeFix> sec_35;
  csNelson::SingleSectionHessian<csBlockMatrix::BlockCoeffSparse, double, csBlockMatrix::Variable, csBlockMatrix::Dynamic> sec_36;
}

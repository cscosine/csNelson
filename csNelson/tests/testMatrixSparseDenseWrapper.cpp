#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csBlockMatrix/DenseMatrixBlock_impl.hpp>
#include <csBlockMatrix/SparsityPattern.hpp>

#include <csNelson/MatrixDenseWrapper_impl.hpp>
#include <csNelson/MatrixSparseWrapper_impl.hpp>
#include <csNelson/SingleSectionHessian_impl.hpp>

static constexpr int secSizeFix = 3;

template <int secType> class SecType {
public:
  static constexpr int secSize = secSizeFix;
};

template <> class SecType<csBlockMatrix::Variable> {
public:
  static const std::vector<int> secSize;
};
const std::vector<int> SecType<csBlockMatrix::Variable>::secSize = {1, 2, 3, 4, 5}; // sum = 15

static constexpr int numBlocks = 5;

TEMPLATE_TEST_CASE_SIG("DenseSparseWrapper", "[MatrixSparseDenseWrapper]",
                       ((int matType, int B, int NB), matType, B, NB),
                       (csBlockMatrix::BlockDense, secSizeFix, numBlocks),
                       (csBlockMatrix::BlockDense, secSizeFix, csBlockMatrix::Dynamic),
                       (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, numBlocks),
                       (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
                       (csBlockMatrix::BlockDense, csBlockMatrix::Variable, numBlocks),
                       (csBlockMatrix::BlockDense, csBlockMatrix::Variable, csBlockMatrix::Dynamic),
                       (csBlockMatrix::BlockDiagonal, secSizeFix, numBlocks),
                       (csBlockMatrix::BlockDiagonal, secSizeFix, csBlockMatrix::Dynamic),
                       (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, numBlocks),
                       (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
                       (csBlockMatrix::BlockDiagonal, csBlockMatrix::Variable, numBlocks),
                       (csBlockMatrix::BlockDiagonal, csBlockMatrix::Variable, csBlockMatrix::Dynamic),
                       (csBlockMatrix::SparseCoeffBlockDiagonal, secSizeFix, numBlocks),
                       (csBlockMatrix::SparseCoeffBlockDiagonal, secSizeFix, csBlockMatrix::Dynamic),
                       (csBlockMatrix::SparseCoeffBlockDiagonal, csBlockMatrix::Dynamic, numBlocks),
                       (csBlockMatrix::SparseCoeffBlockDiagonal, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
                       (csBlockMatrix::SparseCoeffBlockDiagonal, csBlockMatrix::Variable, numBlocks),
                       (csBlockMatrix::SparseCoeffBlockDiagonal, csBlockMatrix::Variable, csBlockMatrix::Dynamic),
                       (csBlockMatrix::BlockSparse, secSizeFix, numBlocks),
                       (csBlockMatrix::BlockSparse, secSizeFix, csBlockMatrix::Dynamic),
                       (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, numBlocks),
                       (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
                       (csBlockMatrix::BlockSparse, csBlockMatrix::Variable, numBlocks),
                       (csBlockMatrix::BlockSparse, csBlockMatrix::Variable, csBlockMatrix::Dynamic),
                       (csBlockMatrix::BlockCoeffSparse, secSizeFix, numBlocks),
                       (csBlockMatrix::BlockCoeffSparse, secSizeFix, csBlockMatrix::Dynamic),
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

  {
    csNelson::DenseWrapper<matType, double, csBlockMatrix::ColMajor, B, B, NB, NB> wrapperDense;
    wrapperDense.set(&sec.H());

    {
      const auto &wr = wrapperDense;
      auto mat2 = wr.mat() * 2;
    }

    csNelson::SparseWrapper<matType, double, csBlockMatrix::ColMajor, B, B, NB, NB> wrapperSparse;
    wrapperSparse.set(&sec.H());

    {
      const auto &wr = wrapperSparse;
      auto mat2 = wr.mat() * 2;
    }
  }

  // dense square
  {
    csNelson::DenseSquareWrapper<matType, double, csBlockMatrix::ColMajor, B, NB> wrapperSquareDense;
    wrapperSquareDense.set(&sec.H());

    typename csNelson::DenseSquareWrapper<matType, double, csBlockMatrix::ColMajor, B, NB>::DiagType diagonal =
        wrapperSquareDense.diagonal();
    diagonal.setConstant(3);
    wrapperSquareDense.setDiagonal(diagonal);
    wrapperSquareDense.diagonalCopy(diagonal);
    REQUIRE((diagonal.array() == wrapperSquareDense.mat().diagonal().array()).all());

    wrapperSquareDense.mat().diagonal().setRandom();
    wrapperSquareDense.diagonalCopy(diagonal);
    REQUIRE((diagonal.array() == wrapperSquareDense.mat().diagonal().array()).all());
  }
  // sparse square
  {
    csNelson::SparseSquareWrapper<matType, double, csBlockMatrix::ColMajor, B, NB> wrapperSquareSparse;
    wrapperSquareSparse.set(&sec.H());

    typename csNelson::SparseSquareWrapper<matType, double, csBlockMatrix::ColMajor, B, NB>::DiagType diagonal =
        wrapperSquareSparse.diagonal();
    diagonal.setConstant(3);
    wrapperSquareSparse.setDiagonal(diagonal);
    wrapperSquareSparse.diagonalCopy(diagonal);
    REQUIRE((diagonal.array() == wrapperSquareSparse.mat().diagonal().array()).all());

    wrapperSquareSparse.mat().diagonal().setRandom();
    wrapperSquareSparse.diagonalCopy(diagonal);
    REQUIRE((diagonal.array() == wrapperSquareSparse.mat().diagonal().array()).all());
  }
}

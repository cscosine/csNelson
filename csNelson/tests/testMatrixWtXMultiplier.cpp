#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/MatrixDenseWrapper_impl.hpp>
#include <csNelson/MatrixWtXMultiplier_impl.hpp>

#include <iostream>

constexpr int BRv = 6;
constexpr int NBRv = 3;

constexpr int BCv = 3;
constexpr int NBCv = 5;

#define DEBUGOUT if (false)

template <int matUType, int matType, int BR, int BC, int NBR, int NBC, int matOutType> void templateTestFunc() {
  typename csBlockMatrix::MatrixBlockIterableTypeTraits<matType, double, csBlockMatrix::RowMajor, BR, BC, NBR,
                                                        NBC>::MatrixType W;

  csBlockMatrix::SparsityPattern<csBlockMatrix::RowMajor>::SPtr sp(
      new csBlockMatrix::SparsityPattern<csBlockMatrix::RowMajor>(NBRv, NBCv));
  if (matType != csBlockMatrix::BlockDiagonal) {
    // 0 1 2 3 4
    // X - - X X
    sp->add(0, 0);
    sp->add(0, 3);
    sp->add(0, 4);
    // X X - X -
    sp->add(1, 0);
    sp->add(1, 1);
    sp->add(1, 3);
    // - X X X X
    sp->add(2, 1);
    sp->add(2, 2); // sp->add(2, 3); sp->add(2, 4);

    //   0 1 2 3 4       0 1 2
    // 0 X - - X X     0 X X -
    // 1 X X - X -  *  1 - X X
    // 2 - X X - -     2 - - X
    //                 3 X X -
    //                 4 X - -

    // [0,0] = 0,0 * (0,0)' + 0,3 * (0,3)' + 0,4 * (0,4)'
    // [0,1] = 0,0 * (1,0)' + 0,3 * (1,3)'
    // [0,2] = 0,3 * (2,3)' + 0,4 * (2,4)'
    // etc

    // result patter (tri up)
    // X X -
    //   X X
    //     X
  } else {
    sp->setDiagonal();
  }

  W.resize(csBlockMatrix::MatrixBlockDescriptor<BR, BC, NBR, NBC>(BRv, BCv, NBRv, NBCv), sp);
  W.setZero();
  for (int bi = 0; bi < W.nonZeroBlocks(); bi++) {
    W.blockByUID(bi).setRandom();
  }

  csBlockMatrix::VectorBlock<double, BR, NBR> xU;
  xU.resize(W.blockDescriptor().rowDescriptionCSPtr());
  xU.mat().setRandom();

  csBlockMatrix::VectorBlock<double, BC, NBC> bV, bVOrig;
  bV.resize(W.blockDescriptor().colDescriptionCSPtr());
  bV.mat().setRandom();
  bVOrig.mat() = bV.mat();

  Eigen::SparseMatrix<int, Eigen::RowMajor> spMat = sp->toSparseMatrix();

  csNelson::MatrixWtXMultiplier<matType, double, BR, BC, NBR, NBC> mWtX;

  SECTION("SINGLE THREAD") { mWtX.settings().setSingleThread(); }
  SECTION("MULTI THREAD") { mWtX.settings().setNumThreadsMax(); }

  mWtX.prepare(W, &spMat);
  mWtX.rightMultVectorSub(W, xU, bV);

  csNelson::DenseWrapper<matType, double, csBlockMatrix::RowMajor, BR, BC, NBR, NBC> W_Wrap;
  W_Wrap.set(&W);
  DEBUGOUT std::cout << "W_Wrap" << std::endl << W_Wrap.mat() << std::endl << std::endl;

  Eigen::MatrixXd diff = bV.mat() - (bVOrig.mat() - W_Wrap.mat().transpose() * xU.mat());
  DEBUGOUT std::cout << "diff" << std::endl << diff << std::endl << std::endl;

  bool ok = (diff.array().abs() < Eigen::NumTraits<double>::dummy_precision()).all();
  REQUIRE(ok);
}

TEMPLATE_TEST_CASE_SIG(
    "MatrixWWtMultiplier-Dense", "[MatrixWWtMultiplier-Base]",
    ((int matUType, int matType, int BR, int BC, int NBR, int NBC, int matOutType), matUType, matType, BR, BC, NBR, NBC,
     matOutType),
    //***************************************************************************
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDense, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse)) {
  templateTestFunc<matUType, matType, BR, BC, NBR, NBC, matOutType>();
}

TEMPLATE_TEST_CASE_SIG(
    "MatrixWWtMultiplier-Sparse", "[MatrixWWtMultiplier-Base]",
    ((int matUType, int matType, int BR, int BC, int NBR, int NBC, int matOutType), matUType, matType, BR, BC, NBR, NBC,
     matOutType),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse)) {
  templateTestFunc<matUType, matType, BR, BC, NBR, NBC, matOutType>();
}

TEMPLATE_TEST_CASE_SIG(
    "MatrixWWtMultiplier-SparseCoeff", "[MatrixWWtMultiplier-Base]",
    ((int matUType, int matType, int BR, int BC, int NBR, int NBC, int matOutType), matUType, matType, BR, BC, NBR, NBC,
     matOutType),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockCoeffSparse)) {
  templateTestFunc<matUType, matType, BR, BC, NBR, NBC, matOutType>();
}

TEMPLATE_TEST_CASE_SIG(
    "MatrixWWtMultiplier-Diagonal", "[MatrixWWtMultiplier-Base]",
    ((int matUType, int matType, int BR, int BC, int NBR, int NBC, int matOutType), matUType, matType, BR, BC, NBR, NBC,
     matOutType),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, NBCv, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, NBCv,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockCoeffSparse),
    //---------------------------------------------------------------------------------
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockDense),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockSparse),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic, csBlockMatrix::BlockCoeffSparse)) {
  templateTestFunc<matUType, matType, BR, BC, NBR, NBC, matOutType>();
}

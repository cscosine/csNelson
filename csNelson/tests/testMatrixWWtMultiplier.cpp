#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/MatrixDenseWrapper_impl.hpp>
#include <csNelson/MatrixWWtMultiplier_impl.hpp>

#include <iostream>

constexpr int BRv = 6;
constexpr int NBRv = 3;

constexpr int BCv = 3;
constexpr int NBCv = 5;

#define DEBUGOUT if (false)

template <int matUType, int matType, int BR, int BC, int NBR, int NBC, int matOutType> void templateTestFunc() {
  typename csBlockMatrix::MatrixBlockIterableTypeTraits<matUType, double, csBlockMatrix::ColMajor, BR, BR, NBR,
                                                        NBR>::MatrixType U;
  csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>::SPtr spU(
      new csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>(NBRv, NBRv));
  spU->setDiagonal();
  if (matUType != csBlockMatrix::BlockDiagonal) {
    spU->add(0, 2);
  }

  U.resize(csBlockMatrix::MatrixBlockDescriptor<BR, BR, NBR, NBR>(BRv, NBRv, BRv, NBRv), spU);
  U.setZero();
  for (int bi = 0; bi < U.nonZeroBlocks(); bi++) {
    U.blockByUID(bi).setRandom();
  }

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

  csNelson::MatrixWWtMultiplier<matOutType, double, BR, NBR> mwwt;

  SECTION("SINGLE THREAD") { mwwt.settings().setSingleThread(); }
  SECTION("MULTI THREAD") { mwwt.settings().setNumThreadsMax(); }

  mwwt.prepare(U, W);
  mwwt.multiply(U, W, W);

  csNelson::DenseWrapper<matType, double, csBlockMatrix::RowMajor, BR, BC, NBR, NBC> W_Wrap;
  W_Wrap.set(&W);
  DEBUGOUT std::cout << "W_Wrap" << std::endl << W_Wrap.mat() << std::endl << std::endl;
  csNelson::DenseWrapper<matUType, double, csBlockMatrix::ColMajor, BR, BR, NBR, NBR> U_Wrap;
  U_Wrap.set(&U);
  DEBUGOUT std::cout << "U_Wrap" << std::endl << U_Wrap.mat() << std::endl << std::endl;

  Eigen::MatrixXd WWt = U_Wrap.mat() - W_Wrap.mat() * W_Wrap.mat().transpose();
  DEBUGOUT std::cout << "WWt" << std::endl << WWt << std::endl << std::endl;

  csNelson::DenseWrapper<matOutType, double, csBlockMatrix::ColMajor, BR, BR, NBR, NBR> WWt_Wrap;
  WWt_Wrap.set(&mwwt.result());

  DEBUGOUT std::cout << "WWt_Wrap" << std::endl << WWt_Wrap.mat() << std::endl << std::endl;

  Eigen::MatrixXd diff = (WWt_Wrap.mat() - WWt).template triangularView<Eigen::Upper>();
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

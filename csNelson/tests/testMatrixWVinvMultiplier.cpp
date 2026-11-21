#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/MatrixDenseWrapper_impl.hpp>
#include <csNelson/MatrixWVinvMultiplier_impl.hpp>

#include <iostream>

constexpr int BRv = 6;
constexpr int NBRv = 3;

constexpr int BCv = 3;
constexpr int NBCv = 5;

#define DEBUGOUT if (false)

template <int matWType, int matVType, int BR, int BC, int NBR, int NBC> void templateTestFunc() {
  typename csBlockMatrix::MatrixBlockIterableTypeTraits<matVType, double, csBlockMatrix::ColMajor, BC, BC, NBC,
                                                        NBC>::MatrixType Vinv;
  csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>::SPtr spVinv(
      new csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>(NBCv, NBCv));
  spVinv->setDiagonal();

  auto blockDescriptorV = csBlockMatrix::MatrixBlockDescriptor<BC, BC, NBC, NBC>(BCv, NBCv, BCv, NBCv);
  auto blockDescriptorW = csBlockMatrix::MatrixBlockDescriptor<BR, BC, NBR, NBC>(BRv, BCv, NBRv, NBCv);

  Vinv.resize(blockDescriptorV, spVinv);
  Vinv.setZero();
  for (int bi = 0; bi < Vinv.nonZeroBlocks(); bi++) {
    Vinv.blockByUID(bi).setRandom();
  }

  csBlockMatrix::VectorBlock<double, BC, NBC> b;
  b.resize(blockDescriptorV.rowDescriptionCSPtr());
  b.mat().setRandom();

  csBlockMatrix::VectorBlock<double, BR, NBR> bv, bvOrig;
  bv.resize(blockDescriptorW.rowDescriptionCSPtr());
  bvOrig.resize(blockDescriptorW.rowDescriptionCSPtr());
  bvOrig.mat().setRandom();
  bv.mat() = bvOrig.mat();

  typename csBlockMatrix::MatrixBlockIterableTypeTraits<matWType, double, csBlockMatrix::RowMajor, BR, BC, NBR,
                                                        NBC>::MatrixType W;

  csBlockMatrix::SparsityPattern<csBlockMatrix::RowMajor>::SPtr sp(
      new csBlockMatrix::SparsityPattern<csBlockMatrix::RowMajor>(NBRv, NBCv));
  if (matWType != csBlockMatrix::BlockDiagonal) {
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

    //   0 1 2 3 4
    // 0 X - - X X
    // 1 X X - X -
    // 2 - X X - -
    //
    //

  } else {
    sp->setDiagonal();
  }

  W.resize(blockDescriptorW, sp);
  W.setZero();
  for (int bi = 0; bi < W.nonZeroBlocks(); bi++) {
    W.blockByUID(bi).setRandom();
  }

  csNelson::MatrixWVinvMultiplier<matWType, matVType, double, BR, BC, NBR, NBC> mwvinv;

  SECTION("SINGLE THREAD") {
    mwvinv.settings().multiplication.setSingleThread();
    mwvinv.settings().rightVectorMult.setSingleThread();
  }
  SECTION("MULTI THREAD") {
    mwvinv.settings().multiplication.setNumThreadsMax();
    mwvinv.settings().rightVectorMult.setSingleThread();
  }

  mwvinv.prepare(W);
  mwvinv.multiply(W, Vinv);
  mwvinv.rightMultVectorSub(b, bv);

  csNelson::DenseWrapper<matWType, double, csBlockMatrix::RowMajor, BR, BC, NBR, NBC> W_Wrap;

  W_Wrap.set(&W);
  DEBUGOUT std::cout << "W_Wrap" << std::endl << W_Wrap.mat() << std::endl << std::endl;
  csNelson::DenseWrapper<matVType, double, csBlockMatrix::ColMajor, BC, BC, NBC, NBC> Vinv_Wrap;
  Vinv_Wrap.set(&Vinv);
  DEBUGOUT std::cout << "Vinv_Wrap" << std::endl << Vinv_Wrap.mat() << std::endl << std::endl;

  Eigen::MatrixXd WVinv = W_Wrap.mat() * Vinv_Wrap.mat();
  DEBUGOUT std::cout << "WVinv" << std::endl << WVinv << std::endl << std::endl;

  csNelson::DenseWrapper<matWType, double, csBlockMatrix::RowMajor, BR, BC, NBR, NBC> WVinv_Wrap;
  WVinv_Wrap.set(&mwvinv.result());

  DEBUGOUT std::cout << "WVinv_Wrap" << std::endl << WVinv_Wrap.mat() << std::endl << std::endl;

  Eigen::MatrixXd diff = (WVinv_Wrap.mat() - WVinv).template triangularView<Eigen::Upper>();
  DEBUGOUT std::cout << "diff" << std::endl << diff << std::endl << std::endl;

  bool ok = (diff.array().abs() < Eigen::NumTraits<double>::dummy_precision()).all();
  REQUIRE(ok);

  Eigen::VectorXd diff2 = (bvOrig.mat() - WVinv_Wrap.mat() * b.mat()) - bv.mat();
  DEBUGOUT std::cout << "diff" << std::endl << diff.transpose() << std::endl << std::endl;

  bool ok2 = (diff2.array().abs() < Eigen::NumTraits<double>::dummy_precision()).all();
  REQUIRE(ok2);
}

TEMPLATE_TEST_CASE_SIG(
    "MatrixWVinvMultiplier-VDiagBlock", "[MatrixWVinvMultiplier-VDiagBlock]",
    ((int matWType, int BR, int BC, int NBR, int NBC), matWType, BR, BC, NBR, NBC),
    (csBlockMatrix::BlockDense, BRv, BCv, NBRv, NBCv), (csBlockMatrix::BlockSparse, BRv, BCv, NBRv, NBCv),
    (csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, NBCv), (csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, BRv, BCv, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, csBlockMatrix::Dynamic),
    //----------------------------------------
    (csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    //***************************************************************************
    (csBlockMatrix::BlockDense, BRv, csBlockMatrix::Dynamic, NBRv, NBCv),
    (csBlockMatrix::BlockSparse, BRv, csBlockMatrix::Dynamic, NBRv, NBCv),
    (csBlockMatrix::BlockCoeffSparse, BRv, csBlockMatrix::Dynamic, NBRv, NBCv),
    (csBlockMatrix::BlockDiagonal, BRv, csBlockMatrix::Dynamic, NBRv, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, BRv, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, BRv, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, BRv, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, BRv, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    //----------------------------------------
    (csBlockMatrix::BlockDense, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockSparse, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockCoeffSparse, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockDiagonal, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    //***************************************************************************
    //***************************************************************************
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, BCv, NBRv, NBCv),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, BCv, NBRv, NBCv),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, BCv, NBRv, NBCv),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, BCv, NBRv, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, BCv, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, BCv, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, BCv, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, BCv, NBRv, csBlockMatrix::Dynamic),
    //----------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    //***************************************************************************
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, NBCv),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, NBCv),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, NBCv),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    //----------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic)) {
  templateTestFunc<matWType, csBlockMatrix::BlockDiagonal, BR, BC, NBR, NBC>();
}

TEMPLATE_TEST_CASE_SIG(
    "MatrixWVinvMultiplier-VDiagCoeff", "[MatrixWVinvMultiplier-VDiagCoeff]",
    ((int matWType, int BR, int BC, int NBR, int NBC), matWType, BR, BC, NBR, NBC),
    (csBlockMatrix::BlockDense, BRv, BCv, NBRv, NBCv), (csBlockMatrix::BlockSparse, BRv, BCv, NBRv, NBCv),
    (csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, NBCv), (csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, BRv, BCv, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, BRv, BCv, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, BRv, BCv, NBRv, csBlockMatrix::Dynamic),
    //----------------------------------------
    (csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, BRv, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    //***************************************************************************
    (csBlockMatrix::BlockDense, BRv, csBlockMatrix::Dynamic, NBRv, NBCv),
    (csBlockMatrix::BlockSparse, BRv, csBlockMatrix::Dynamic, NBRv, NBCv),
    (csBlockMatrix::BlockCoeffSparse, BRv, csBlockMatrix::Dynamic, NBRv, NBCv),
    (csBlockMatrix::BlockDiagonal, BRv, csBlockMatrix::Dynamic, NBRv, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, BRv, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, BRv, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, BRv, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, BRv, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    //----------------------------------------
    (csBlockMatrix::BlockDense, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockSparse, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockCoeffSparse, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockDiagonal, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, BRv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    //***************************************************************************
    //***************************************************************************
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, BCv, NBRv, NBCv),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, BCv, NBRv, NBCv),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, BCv, NBRv, NBCv),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, BCv, NBRv, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, BCv, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, BCv, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, BCv, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, BCv, NBRv, csBlockMatrix::Dynamic),
    //----------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, BCv, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
    //***************************************************************************
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, NBCv),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, NBCv),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, NBCv),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBRv, csBlockMatrix::Dynamic),
    //----------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, NBCv),
    //----------------------------------------
    (csBlockMatrix::BlockDense, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockCoeffSparse, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic),
    (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic,
     csBlockMatrix::Dynamic)) {
  templateTestFunc<matWType, csBlockMatrix::SparseCoeffBlockDiagonal, BR, BC, NBR, NBC>();
}

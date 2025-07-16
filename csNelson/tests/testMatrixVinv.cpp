#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <csNelson/MatrixDenseWrapper_impl.hpp>
#include <csNelson/MatrixDiagInv_impl.hpp>

#include <csBlockMatrix/DiagonalMatrixBlock_impl.hpp>
#include <csBlockMatrix/MatrixBlockDescriptor_impl.hpp>
#include <csBlockMatrix/SparseCoeffDiagonalMatrixBlock_impl.hpp>
#include <csBlockMatrix/VectorBlock_impl.hpp>

#include <iostream>

constexpr int BCv = 3;
constexpr int NBCv = 5;

#define DEBUGOUT if (false)

template <int matType, int BC, int NBC>
void templateTestFunc() {
  typename csBlockMatrix::MatrixBlockIterableTypeTraits<matType, double, csBlockMatrix::ColMajor, BC, BC, NBC, NBC>::MatrixType V;
  csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>::SPtr spV(
      new csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>(NBCv, NBCv));
  spV->setDiagonal();

  auto blockDescriptorV = csBlockMatrix::MatrixBlockDescriptor<BC, BC, NBC, NBC>(BCv, NBCv, BCv, NBCv);
  V.resize(blockDescriptorV, spV);
  V.setZero();
  for (int bi = 0; bi < V.nonZeroBlocks(); bi++) {
    V.blockByUID(bi).setRandom();
  }

  csBlockMatrix::VectorBlock<double, BC, NBC> b, bv;
  b.resize(blockDescriptorV.rowDescriptionCSPtr());
  b.mat().setRandom();
  bv.resize(blockDescriptorV.rowDescriptionCSPtr());

  csNelson::MatrixDiagInv<matType, double, BC, NBC, csBlockMatrix::BlockDiagonal> mvinv;

  SECTION("SINGLE THREAD") {
    mvinv.settings().blockInversion.setSingleThread();
    mvinv.settings().rightVectorMult.setSingleThread();
  }
  SECTION("MULTI THREAD") {
    mvinv.settings().blockInversion.setNumThreadsMax();
    mvinv.settings().rightVectorMult.setNumThreadsMax();
  }

  mvinv.init(V);
  mvinv.compute(V, 0, 0);
  mvinv.rightMultVector(b, bv);

  csNelson::DenseWrapper<csBlockMatrix::BlockDiagonal, double, csBlockMatrix::ColMajor, BC, BC, NBC, NBC> Vinv_Wrap;
  Vinv_Wrap.set(&mvinv.Vinv());
  DEBUGOUT std::cout << "Vinv_Wrap" << std::endl << Vinv_Wrap.mat() << std::endl << std::endl;

  Eigen::VectorXd diff = Vinv_Wrap.mat() * b.mat() - bv.mat();
  DEBUGOUT std::cout << "diff" << std::endl << diff.transpose() << std::endl << std::endl;

  bool ok = (diff.array().abs() < Eigen::NumTraits<double>::dummy_precision()).all();
  REQUIRE(ok);
}

TEMPLATE_TEST_CASE_SIG("MatrixVinvMultiplier", "[MatrixVinvMultiplier]", ((int matType, int BC, int NBC), matType, BC, NBC),
                       (csBlockMatrix::BlockDiagonal, BCv, NBCv), (csBlockMatrix::BlockDiagonal, BCv, csBlockMatrix::Dynamic),
                       (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, NBCv),
                       (csBlockMatrix::BlockDiagonal, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic),
                       (csBlockMatrix::SparseCoeffBlockDiagonal, BCv, NBCv),
                       (csBlockMatrix::SparseCoeffBlockDiagonal, BCv, csBlockMatrix::Dynamic),
                       (csBlockMatrix::SparseCoeffBlockDiagonal, csBlockMatrix::Dynamic, NBCv),
                       (csBlockMatrix::SparseCoeffBlockDiagonal, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic)) {
  templateTestFunc<matType, BC, NBC>();
}

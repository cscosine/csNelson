#include <iostream>
#include <vector>

#include <samples/Problems.hpp>

#include <chrono>
#include <fstream>

#include <csNelson/MatrixDenseWrapper_impl.hpp>
#include <csNelson/MatrixSparseWrapper_impl.hpp>
#include <csNelson/MatrixWVinvMultiplier_impl.hpp>

constexpr int BR = 6;
constexpr int BC = 3;

template <int matVType, int matWType> void testFunction(const Problem &problem) {
  typename csBlockMatrix::MatrixBlockIterableTypeTraits<matVType, double, csBlockMatrix::ColMajor, BC, BC,
                                                        csBlockMatrix::Dynamic, csBlockMatrix::Dynamic>::MatrixType V;
  csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>::SPtr spV(
      new csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>(problem.nPoints, problem.nPoints));
  spV->setDiagonal();

  V.resize(csBlockMatrix::MatrixBlockDescriptor<BC, BC, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic>(
               BC, problem.nPoints, BC, problem.nPoints),
           spV);
  V.setZero();
  for (int bi = 0; bi < V.nonZeroBlocks(); bi++) {
    V.blockByUID(bi).setRandom();
  }

  typename csBlockMatrix::MatrixBlockIterableTypeTraits<matWType, double, csBlockMatrix::RowMajor, BR, BC,
                                                        csBlockMatrix::Dynamic, csBlockMatrix::Dynamic>::MatrixType W;
  csBlockMatrix::SparsityPattern<csBlockMatrix::RowMajor>::SPtr sp(
      new csBlockMatrix::SparsityPattern<csBlockMatrix::RowMajor>(problem.nCameras, problem.nPoints));
  for (int i = 0; i < problem.edges.size(); i++) {
    const auto &e = problem.edges[i];
    sp->add(e.first, e.second);
  }

  W.resize(csBlockMatrix::MatrixBlockDescriptor<BR, BC, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic>(
               BR, problem.nCameras, BC, problem.nPoints),
           sp);
  W.setZero();
  for (int bi = 0; bi < W.nonZeroBlocks(); bi++) {
    W.blockByUID(bi).setRandom();
    // W.blockByUID(bi).setZero();
  }

  csBlockMatrix::VectorBlock<double, BC, csBlockMatrix::Dynamic> b;
  b.resize(V.blockDescriptor().rowDescriptionCSPtr());
  b.mat().setRandom();

  csBlockMatrix::VectorBlock<double, BR, csBlockMatrix::Dynamic> bv;
  bv.resize(W.blockDescriptor().rowDescriptionCSPtr());
  bv.mat().setRandom();

  csNelson::MatrixWVinvMultiplier<matWType, matVType, double, BR, BC, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic>
      mwvinv;
  const int half_maxthreads = csNelson::ParallelExecSettings::maxSupportedThreads() / 2;
  const int maxthreads = csNelson::ParallelExecSettings::maxSupportedThreads();

  // with multiplier
  {
    auto t0 = std::chrono::steady_clock::now();
    mwvinv.prepare(W);
    auto t1 = std::chrono::steady_clock::now();
    mwvinv.settings().multiplication.setSingleThread();
    mwvinv.multiply(W, V);
    auto t2 = std::chrono::steady_clock::now();
    mwvinv.settings().multiplication.setNumThreads(half_maxthreads);
    mwvinv.multiply(W, V);
    auto t3 = std::chrono::steady_clock::now();
    mwvinv.settings().multiplication.setNumThreadsMax();
    mwvinv.multiply(W, V);
    auto t4 = std::chrono::steady_clock::now();

    double t_prep = std::chrono::duration<double>(t1 - t0).count();
    double t_mult = std::chrono::duration<double>(t2 - t1).count();
    double t_tot = std::chrono::duration<double>(t2 - t0).count();
    double t_mult2 = std::chrono::duration<double>(t3 - t2).count();
    double t_mult3 = std::chrono::duration<double>(t4 - t3).count();

    std::cout << "t_prep " << t_prep << std::endl;
    std::cout << "t_mult " << t_mult << std::endl;
    std::cout << "t_tot  " << t_tot << std::endl;
    std::cout << "t_mult(" << half_maxthreads << ") " << t_mult2 << std::endl;
    std::cout << "t_mult(" << maxthreads << ") " << t_mult3 << std::endl;
  }

  // standard
  csNelson::SparseWrapper<matVType, double, csBlockMatrix::ColMajor, BC, BC, csBlockMatrix::Dynamic,
                          csBlockMatrix::Dynamic>
      V_Wrap;
  V_Wrap.set(&V);
  csNelson::SparseWrapper<matWType, double, csBlockMatrix::RowMajor, BR, BC, csBlockMatrix::Dynamic,
                          csBlockMatrix::Dynamic>
      W_Wrap;
  W_Wrap.set(&W);
  Eigen::SparseMatrix<double, csBlockMatrix::RowMajor> WVinv;
  {
    auto t0 = std::chrono::steady_clock::now();
    WVinv = W_Wrap.mat() * V_Wrap.mat();
    auto t1 = std::chrono::steady_clock::now();
    double t_tot = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "t_tot  " << t_tot << std::endl;

    auto t20 = std::chrono::steady_clock::now();
    WVinv = W_Wrap.mat() * V_Wrap.mat();
    auto t21 = std::chrono::steady_clock::now();
    double t_tot2 = std::chrono::duration<double>(t21 - t20).count();
    std::cout << "t_tot(2)  " << t_tot2 << std::endl;
    WVinv = W_Wrap.mat() * V_Wrap.mat();
    auto t22 = std::chrono::steady_clock::now();
    double t_tot3 = std::chrono::duration<double>(t22 - t21).count();
    std::cout << "t_tot(3)  " << t_tot3 << std::endl;
  }

  // checks
  csNelson::SparseWrapper<matWType, double, csBlockMatrix::RowMajor, BR, BC, csBlockMatrix::Dynamic,
                          csBlockMatrix::Dynamic>
      WVinv_wrap;
  WVinv_wrap.set(&mwvinv.result());

  Eigen::SparseMatrix<double, csBlockMatrix::RowMajor> diff =
      (WVinv_wrap.mat() - WVinv).template triangularView<Eigen::Upper>();

  std::cout << "Eigen::NumTraits<double>::dummy_precision() " << Eigen::NumTraits<double>::dummy_precision()
            << std::endl;

  double maxErr = diff.coeffs().cwiseAbs().maxCoeff();
  bool ok = (diff.coeffs().cwiseAbs() < 16 * Eigen::NumTraits<double>::dummy_precision()).all();
  if (false) {
    std::ofstream f;
    f.open("loadWWt.m");
    f << "WWt = [" << std::endl;
    f << WVinv_wrap.mat() << std::endl;
    f << "];" << std::endl;
  }
  if (!ok) {
    std::cerr << "ERROR!! max is " << maxErr << std::endl;
    {
      std::ofstream f;
      f.open("loadDiffMat.m");
      f << "D = [" << std::endl;
      f << diff << std::endl;
      f << "];" << std::endl;
    }
    std::exit(-1);
  } else {
    std::cerr << "OK" << std::endl;
  }

  // with multiplier, vector
  {
    auto t1 = std::chrono::steady_clock::now();
    mwvinv.settings().rightVectorMult.setSingleThread();
    mwvinv.rightMultVectorSub(b, bv);
    auto t2 = std::chrono::steady_clock::now();
    mwvinv.settings().rightVectorMult.setNumThreads(half_maxthreads);
    mwvinv.rightMultVectorSub(b, bv);
    auto t3 = std::chrono::steady_clock::now();
    mwvinv.settings().rightVectorMult.setNumThreadsMax();
    mwvinv.rightMultVectorSub(b, bv);
    auto t4 = std::chrono::steady_clock::now();

    double t_mult = std::chrono::duration<double>(t2 - t1).count();
    double t_mult2 = std::chrono::duration<double>(t3 - t2).count();
    double t_mult3 = std::chrono::duration<double>(t4 - t3).count();

    std::cout << "**** VECTOR ****" << std::endl;
    std::cout << "t_mult " << t_mult << std::endl;
    std::cout << "t_mult(" << half_maxthreads << ") " << t_mult2 << std::endl;
    std::cout << "t_mult(" << maxthreads << ") " << t_mult3 << std::endl;
  }
}

int main(int argc, char *argv[]) {
  const Problem *problemPtr = &ProblemCollections::ladybug_49_7776;
  std::unique_ptr<Problem> newProblem;
  if (argc == 2) {
    newProblem.reset(new Problem());
    bool ok = newProblem->load(argv[1]);
    if (!ok) {
      std::cerr << "error reading " << argv[1] << std::endl;
      std::exit(-1);
    } else {
      problemPtr = newProblem.get();
    }
  }
  const Problem &problem = *problemPtr;

  std::cout << "--------- V BlockDiagonal, W BlockCoeffSparse" << std::endl;
  testFunction<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockCoeffSparse>(problem);
  std::cout << std::endl;

  std::cout << "--------- V SparseCoeffBlockDiagonal, W BlockCoeffSparse" << std::endl;
  testFunction<csBlockMatrix::SparseCoeffBlockDiagonal, csBlockMatrix::BlockCoeffSparse>(problem);
  std::cout << std::endl;

  std::cout << "--------- V BlockDiagonal, W BlockSparse" << std::endl;
  testFunction<csBlockMatrix::BlockDiagonal, csBlockMatrix::BlockSparse>(problem);
  std::cout << std::endl;

  std::cout << "--------- V SparseCoeffBlockDiagonal, W BlockSparse" << std::endl;
  testFunction<csBlockMatrix::SparseCoeffBlockDiagonal, csBlockMatrix::BlockSparse>(problem);
  std::cout << std::endl;

  return 0;
}

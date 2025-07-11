#include <iostream>
#include <vector>

#include "Problems.h"

#include <chrono>
#include <fstream>

#include "csNelson/MatrixDenseWrapper.hpp"
#include "csNelson/MatrixDiagInv.hpp"
#include "csNelson/MatrixSparseWrapper.hpp"

constexpr int BC = 3;

template <int matVType>
void testFunction(const Problem& problem) {
  typename csBlockMatrix::MatrixBlockIterableTypeTraits<matVType, double, csBlockMatrix::ColMajor, BC, BC, csBlockMatrix::Dynamic,
                                                        csBlockMatrix::Dynamic>::MatrixType V;
  csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>::SPtr spV(
      new csBlockMatrix::SparsityPattern<csBlockMatrix::ColMajor>(problem.nPoints, problem.nPoints));
  spV->setDiagonal();
  csBlockMatrix::MatrixBlockDescriptor<BC, BC, csBlockMatrix::Dynamic, csBlockMatrix::Dynamic> blockDescriptorV(BC, problem.nPoints,
                                                                                                                BC, problem.nPoints);
  V.resize(blockDescriptorV, spV);
  V.setZero();
  for (int bi = 0; bi < V.nonZeroBlocks(); bi++) {
    V.blockByUID(bi).setRandom();
  }

  csBlockMatrix::VectorBlock<double, BC, csBlockMatrix::Dynamic> b, bv;
  b.resize(blockDescriptorV.rowDescriptionCSPtr());
  b.mat().setRandom();
  bv.resize(blockDescriptorV.rowDescriptionCSPtr());

  csNelson::MatrixDiagInv<matVType, double, BC, csBlockMatrix::Dynamic, csBlockMatrix::BlockDiagonal> mwvinv;
  const int half_maxthreads = csNelson::ParallelExecSettings::maxSupportedThreads() / 2;
  const int maxthreads = csNelson::ParallelExecSettings::maxSupportedThreads();

  {
    auto t0 = std::chrono::steady_clock::now();
    mwvinv.init(V);
    auto t1 = std::chrono::steady_clock::now();
    mwvinv.settings().blockInversion.setSingleThread();
    mwvinv.compute(V, 0, 0);
    auto t2 = std::chrono::steady_clock::now();
    mwvinv.settings().blockInversion.setNumThreads(half_maxthreads);
    mwvinv.compute(V, 0, 0);
    auto t3 = std::chrono::steady_clock::now();
    mwvinv.settings().blockInversion.setNumThreadsMax();
    mwvinv.compute(V, 0, 0);
    auto t4 = std::chrono::steady_clock::now();

    double t_prep = std::chrono::duration<double>(t1 - t0).count();
    double t_mult = std::chrono::duration<double>(t2 - t1).count();
    double t_tot = std::chrono::duration<double>(t2 - t0).count();
    double t_mult2 = std::chrono::duration<double>(t3 - t2).count();
    double t_mult3 = std::chrono::duration<double>(t4 - t3).count();

    std::cout << "*** inversion time **" << std::endl;
    std::cout << "t_prep " << t_prep << std::endl;
    std::cout << "t_mult " << t_mult << std::endl;
    std::cout << "t_tot  " << t_tot << std::endl;
    std::cout << "t_mult(" << half_maxthreads << ") " << t_mult2 << std::endl;
    std::cout << "t_mult(" << maxthreads << ") " << t_mult3 << std::endl;
  }

  {
    auto t1 = std::chrono::steady_clock::now();
    mwvinv.settings().rightVectorMult.setSingleThread();
    mwvinv.rightMultVector(b, bv);
    auto t2 = std::chrono::steady_clock::now();
    mwvinv.settings().rightVectorMult.setNumThreads(half_maxthreads);
    mwvinv.rightMultVector(b, bv);
    auto t3 = std::chrono::steady_clock::now();
    mwvinv.settings().rightVectorMult.setNumThreadsMax();
    mwvinv.rightMultVector(b, bv);
    auto t4 = std::chrono::steady_clock::now();

    double t_mult = std::chrono::duration<double>(t2 - t1).count();
    double t_mult2 = std::chrono::duration<double>(t3 - t2).count();
    double t_mult3 = std::chrono::duration<double>(t4 - t3).count();

    std::cout << "*** vector mult time **" << std::endl;
    std::cout << "t_mult " << t_mult << std::endl;
    std::cout << "t_mult(" << half_maxthreads << ") " << t_mult2 << std::endl;
    std::cout << "t_mult(" << maxthreads << ") " << t_mult3 << std::endl;
  }
}

int main(int argc, char* argv[]) {
  const Problem* problemPtr = &ProblemCollections::ladybug_49_7776;
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
  const Problem& problem = *problemPtr;
  std::cout << problem.toString() << std::endl;

  std::cout << "******* csBlockMatrix::BlockDiagonal" << std::endl;
  testFunction<csBlockMatrix::BlockDiagonal>(problem);

  std::cout << "******* csBlockMatrix::SparseCoeffBlockDiagonal" << std::endl;
  testFunction<csBlockMatrix::SparseCoeffBlockDiagonal>(problem);

  return 0;
}

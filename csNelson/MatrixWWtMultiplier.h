#pragma once
#include "Global.h"

#include "csBlockMatrix/MatrixTypeTraits.h"

#include "ParallelExecHelper.h"

#include <Eigen/Sparse>

namespace csNelson {

namespace _private {
struct UIDPair {
  int uid_1, uid_2;
};
} // namespace _private

using MatrixWWtMultiplierSettings = ParallelExecSettings;

template <int matOutputType, class T, int BR, int NBR>
class MatrixWWtMultiplier {
public:
  using MatOuputType =
      typename csBlockMatrix::MatrixBlockIterableTypeTraits<matOutputType, T, csBlockMatrix::ColMajor, BR, BR, NBR, NBR>::MatrixType;

  using Settings = MatrixWWtMultiplierSettings;

private:
  using UIDPair = _private::UIDPair;
  std::vector<std::vector<UIDPair>> _blockPairs;
  std::vector<int> _blockPairsU;
  MatOuputType _matOutput;

  Settings _settings;

public:
  Settings& settings() {
    return _settings;
  }
  const Settings& settings() const {
    return _settings;
  }

  template <class MatrixBlockU, class MatrixBlockW>
  void prepare(const MatrixBlockU& U, const MatrixBlockW& Wpattern, const Eigen::SparseMatrix<int, Eigen::RowMajor>* spWmat = nullptr);

  template <class MatrixBlockU, class MatrixBlockW>
  void multiply(const MatrixBlockU& U, const MatrixBlockW& A, const MatrixBlockW& B);

  const MatOuputType& result() const {
    return _matOutput;
  }
  MatOuputType& result() {
    return _matOutput;
  }
};

} // namespace csNelson

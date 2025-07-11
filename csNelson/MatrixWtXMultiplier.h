#pragma once
#include "Global.h"

#include "csBlockMatrix/MatrixTypeTraits.h"
#include "csBlockMatrix/VectorBlock.h"

#include "ParallelExecHelper.h"

#include <Eigen/Sparse>
#include <vector>

namespace csNelson {

namespace _private {
struct UIDColPair {
  int uid;
  int xid;
};
} // namespace _private

using MatrixWtXMultiplierSettings = ParallelExecSettings;

template <int matWType, class T, int BR, int BC, int NBR, int NBC>
class MatrixWtXMultiplier {
public:
  using MatType =
      typename csBlockMatrix::MatrixBlockIterableTypeTraits<matWType, T, csBlockMatrix::RowMajor, BR, BC, NBR, NBC>::MatrixType;

  using Settings = MatrixWtXMultiplierSettings;

private:
  Settings _settings;

  using UIDColPair = _private::UIDColPair;
  std::vector<std::vector<UIDColPair>> _blocks;

public:
  Settings& settings() {
    return _settings;
  }
  const Settings& settings() const {
    return _settings;
  }

  void prepare(const MatType& W, const Eigen::SparseMatrix<int, Eigen::RowMajor>* spWmat = nullptr);

  void rightMultVectorSub(const MatType& W, const csBlockMatrix::VectorBlock<T, BR, NBR>& v,
                          csBlockMatrix::VectorBlock<T, BC, NBC>& res) const;
};

} // namespace csNelson

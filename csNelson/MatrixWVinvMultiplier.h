#pragma once
#include "Global.h"

#include "csBlockMatrix/MatrixTypeTraits.h"
#include "csBlockMatrix/VectorBlock.h"

#include "ParallelExecHelper.h"

namespace csNelson {

struct MatrixWVinvMultiplierSettings {
  ParallelExecSettings multiplication;
  ParallelExecSettings rightVectorMult;
};

template <int matWType, int matVType, class T, int BR, int BC, int NBR, int NBC>
class MatrixWVinvMultiplier {
public:
  using MatType =
      typename csBlockMatrix::MatrixBlockIterableTypeTraits<matWType, T, csBlockMatrix::RowMajor, BR, BC, NBR, NBC>::MatrixType;
  using MatTypeV =
      typename csBlockMatrix::MatrixBlockIterableTypeTraits<matVType, T, csBlockMatrix::ColMajor, BC, BC, NBC, NBC>::MatrixType;

  using Settings = MatrixWVinvMultiplierSettings;

private:
  MatType _matOutput;

  Settings _settings;

public:
  Settings& settings() {
    return _settings;
  }
  const Settings& settings() const {
    return _settings;
  }

  // only one input, no need for the other matrix since it is the diagonal block
  void prepare(const MatType& W);

  void multiply(const MatType& W, const MatTypeV& Vinv);

  const MatType& result() const {
    return _matOutput;
  }
  MatType& result() {
    return _matOutput;
  }

  void rightMultVectorSub(const csBlockMatrix::VectorBlock<T, BC, NBC>& v, csBlockMatrix::VectorBlock<T, BR, NBR>& res) const;
};

} // namespace csNelson
